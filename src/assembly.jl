# Compile the symbolic stencils into fast Julia functions and assemble the
# global sparse SSA operator.

"""
    CompiledStencil

Holds the in-place coefficient function, scalar right-hand-side function, and
the offset metadata needed to map a local stencil into the global state vector.
"""
struct CompiledStencil
    coeffs!::Function       # (out::Vector{T}, params::Vector{T}) -> nothing
    rhs::Function           # (params::Vector{T}) -> T
    nunknowns::Int
    unknown_offsets::Vector{Tuple{Symbol,Int,Int}}
    H_offsets::Vector{Tuple{Int,Int}}
    η_offsets::Vector{Tuple{Int,Int}}
    β_offsets::Vector{Tuple{Int,Int}}
    s_offsets::Vector{Tuple{Int,Int}}
    nparams::Int
    p_H::UnitRange{Int}
    p_η::UnitRange{Int}
    p_β::UnitRange{Int}
    p_s::UnitRange{Int}
    p_scalar::UnitRange{Int}   # dx, dy, ρ, g
end

"""
    compile_stencil(spec) -> CompiledStencil

Use Symbolics.jl to extract the linear stencil coefficients (Jacobian w.r.t.
the velocity unknowns) and the constant RHS term, then `build_function` both
into compiled Julia closures.
"""
function compile_stencil(spec::StencilSpec)
    nH = length(spec.H_vars); nη = length(spec.η_vars)
    nβ = length(spec.β_vars); ns = length(spec.s_vars)
    p_H = 1:nH
    p_η = (nH+1):(nH+nη)
    p_β = (nH+nη+1):(nH+nη+nβ)
    p_s = (nH+nη+nβ+1):(nH+nη+nβ+ns)
    p_scalar = (nH+nη+nβ+ns+1):(nH+nη+nβ+ns+4)

    params = vcat(spec.H_vars, spec.η_vars, spec.β_vars, spec.s_vars,
                  Num[spec.dx, spec.dy, spec.ρ, spec.g])

    J = Symbolics.jacobian([spec.residual], spec.unknowns)
    coeffs = collect(J[1, :])

    # RHS = -(residual with all unknowns set to zero)
    sub = Dict{Num,Num}(uk => Num(0) for uk in spec.unknowns)
    rhs_expr = -Symbolics.substitute(spec.residual, sub)

    _, coeffs_ip = build_function(coeffs, params; expression=Val{false})
    rhs_oop = build_function(rhs_expr, params; expression=Val{false})

    return CompiledStencil(coeffs_ip, rhs_oop,
                           length(spec.unknowns),
                           spec.unknown_offsets,
                           spec.H_offsets, spec.η_offsets,
                           spec.β_offsets, spec.s_offsets,
                           length(params),
                           p_H, p_η, p_β, p_s, p_scalar)
end

"""
    SSAStencils

Cached compiled stencils for both momentum equations.
"""
struct SSAStencils
    x::CompiledStencil
    y::CompiledStencil
end

const _CACHED_STENCILS = Ref{Union{Nothing,SSAStencils}}(nothing)

"""
    ssa_stencils() -> SSAStencils

Lazily build (and cache) the compiled SSA stencils. The first call pays the
Symbolics + `build_function` cost; subsequent calls are free.
"""
function ssa_stencils()
    if _CACHED_STENCILS[] === nothing
        _CACHED_STENCILS[] = SSAStencils(
            compile_stencil(derive_x_residual()),
            compile_stencil(derive_y_residual()),
        )
    end
    return _CACHED_STENCILS[]
end

"""
    SSAFields{T}

Cell-centred input fields on a grid of size `(Nx, Ny)`:
- `H`  ice thickness
- `η`  effective viscosity (held fixed; outer Picard iteration lives elsewhere)
- `β`  basal-friction coefficient
- `s`  surface elevation
"""
struct SSAFields{T}
    H::Matrix{T}
    η::Matrix{T}
    β::Matrix{T}
    s::Matrix{T}
    function SSAFields(H::Matrix{T}, η::Matrix{T}, β::Matrix{T}, s::Matrix{T}) where {T}
        size(H) == size(η) == size(β) == size(s) ||
            throw(ArgumentError("H, η, β, s must share the same size"))
        new{T}(H, η, β, s)
    end
end

Base.size(f::SSAFields) = size(f.H)

@inline _u_index(i, j, Nx, Ny) = i + (j - 1) * Nx
@inline _v_index(i, j, Nx, Ny) = Nx * Ny + i + (j - 1) * Nx

@inline function _global_col(kind::Symbol, i::Int, j::Int, di::Int, dj::Int, Nx::Int, Ny::Int)
    ip = mod1(i + di, Nx)
    jp = mod1(j + dj, Ny)
    return kind === :u ? _u_index(ip, jp, Nx, Ny) : _v_index(ip, jp, Nx, Ny)
end

# Fill `pbuf` with the cell-centred parameter values pulled from the field
# arrays at offsets `offs` relative to (i, j), with periodic wrapping.
@inline function _gather!(pbuf, range, F, i, j, offs, Nx, Ny)
    @inbounds for (k, off) in enumerate(offs)
        pbuf[range[k]] = F[mod1(i + off[1], Nx), mod1(j + off[2], Ny)]
    end
end

"""
    assemble_ssa(fields::SSAFields; dx, dy, ρ=910.0, g=9.81,
                 stencils::SSAStencils=ssa_stencils())
        -> (A::SparseMatrixCSC, b::Vector)

Assemble the linear system `A * [u; v] = b` for the SSA momentum balance on a
periodic grid, using the C-grid stencil derived symbolically by Symbolics.jl.

The unknown vector is `[u_flat; v_flat]` with column-major flattening of each
`Nx × Ny` array.
"""
function assemble_ssa(fields::SSAFields{T};
                      dx::Real, dy::Real,
                      ρ::Real = 910.0, g::Real = 9.81,
                      stencils::SSAStencils = ssa_stencils()) where {T}
    Nx, Ny = size(fields)
    N = Nx * Ny
    Tp = float(T)

    # Per-stencil scratch buffers.
    pbuf_x = Vector{Tp}(undef, stencils.x.nparams)
    pbuf_y = Vector{Tp}(undef, stencils.y.nparams)
    cbuf_x = Vector{Tp}(undef, stencils.x.nunknowns)
    cbuf_y = Vector{Tp}(undef, stencils.y.nunknowns)

    # Set the dx, dy, ρ, g entries once (same across all rows).
    for (buf, sten) in ((pbuf_x, stencils.x), (pbuf_y, stencils.y))
        sc = sten.p_scalar
        buf[sc[1]] = dx; buf[sc[2]] = dy; buf[sc[3]] = ρ; buf[sc[4]] = g
    end

    # Triplet storage. Each row contributes `nunknowns` nonzeros.
    nnz_x = N * stencils.x.nunknowns
    nnz_y = N * stencils.y.nunknowns
    Is = Vector{Int}(undef, nnz_x + nnz_y)
    Js = Vector{Int}(undef, nnz_x + nnz_y)
    Vs = Vector{Tp}(undef, nnz_x + nnz_y)
    b  = Vector{Tp}(undef, 2 * N)

    idx = 0
    for j in 1:Ny, i in 1:Nx
        # ---- x-momentum row ----
        row = _u_index(i, j, Nx, Ny)
        _gather!(pbuf_x, stencils.x.p_H, fields.H, i, j, stencils.x.H_offsets, Nx, Ny)
        _gather!(pbuf_x, stencils.x.p_η, fields.η, i, j, stencils.x.η_offsets, Nx, Ny)
        _gather!(pbuf_x, stencils.x.p_β, fields.β, i, j, stencils.x.β_offsets, Nx, Ny)
        _gather!(pbuf_x, stencils.x.p_s, fields.s, i, j, stencils.x.s_offsets, Nx, Ny)
        stencils.x.coeffs!(cbuf_x, pbuf_x)
        b[row] = stencils.x.rhs(pbuf_x)
        @inbounds for k in 1:stencils.x.nunknowns
            kind, di, dj = stencils.x.unknown_offsets[k]
            idx += 1
            Is[idx] = row
            Js[idx] = _global_col(kind, i, j, di, dj, Nx, Ny)
            Vs[idx] = cbuf_x[k]
        end

        # ---- y-momentum row ----
        row = _v_index(i, j, Nx, Ny)
        _gather!(pbuf_y, stencils.y.p_H, fields.H, i, j, stencils.y.H_offsets, Nx, Ny)
        _gather!(pbuf_y, stencils.y.p_η, fields.η, i, j, stencils.y.η_offsets, Nx, Ny)
        _gather!(pbuf_y, stencils.y.p_β, fields.β, i, j, stencils.y.β_offsets, Nx, Ny)
        _gather!(pbuf_y, stencils.y.p_s, fields.s, i, j, stencils.y.s_offsets, Nx, Ny)
        stencils.y.coeffs!(cbuf_y, pbuf_y)
        b[row] = stencils.y.rhs(pbuf_y)
        @inbounds for k in 1:stencils.y.nunknowns
            kind, di, dj = stencils.y.unknown_offsets[k]
            idx += 1
            Is[idx] = row
            Js[idx] = _global_col(kind, i, j, di, dj, Nx, Ny)
            Vs[idx] = cbuf_y[k]
        end
    end

    A = sparse(Is, Js, Vs, 2 * N, 2 * N)
    return A, b
end
