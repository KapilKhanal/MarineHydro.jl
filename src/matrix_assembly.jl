# Below are several variants of the same function.
# TODO: refactor with multiple dispatch?

function assemble_matrices_comprehension(green_functions, mesh, wavenumber; direct=true, all_normals=nothing, include_identity=true)
    
    free_surface = 0.0

    if eltype(mesh.vertices)<: ForwardDiff.Dual
        T = eltype(mesh.vertices)
    else
        T = eltype(wavenumber)
    end

    S = @inbounds [-1/2τ̅ * Complex{T}(integral(green_functions, element(mesh, i), element(mesh, j), wavenumber)) for i in 1:mesh.nfaces, j in 1:mesh.nfaces]

    D = @inbounds [begin
            element_i = element(mesh, i)
            element_j = element(mesh, j)

            norm_vec = direct ? normal(element_j) : normal(element_i)
        
            n = isnothing(all_normals) ? norm_vec : all_normals                      

            c = i == j ? Complex{T}(1.0, 0.0) : Complex{T}(0.0, 0.0) # if diagonal

            constant = abs(mesh.centers[i,3]-free_surface) < 1e-8 ? c : c/2 # if panel on surface

            # (n' * norm_vec)=1 when n==panel normal vector. When all_normals is provided (for forward speed problems), this is not always the case.
            jump = include_identity ? (constant * (n' * norm_vec)) : zero(Complex{T})
            jump - 1/2τ̅ * Complex{T}(n' * integral_gradient(green_functions, element_i, element_j, wavenumber; with_respect_to_first_variable=!direct))
        end for i in 1:mesh.nfaces, j in 1:mesh.nfaces]

    return S, D
end

function assemble_matrices_explicit_both(green_functions, mesh, wavenumber; direct=true)
    # Variant of the above using `both_integral_and_integral_gradient`
    # Probably not differentiable with Zygote
    S = Array{ComplexF64, 2}(undef, (mesh.nfaces, mesh.nfaces))
    D = Array{ComplexF64, 2}(undef, (mesh.nfaces, mesh.nfaces))
    for i in 1:mesh.nfaces
        for j in 1:mesh.nfaces
            element_i = element(mesh, i)
            element_j = element(mesh, j)
            Sij, Dij = both_integral_and_integral_gradient(green_functions, element_i, element_j, wavenumber; with_respect_to_first_variable=!direct)
            S[i, j] = -1/4π * Complex(Sij)
            n = direct ? normal(element_j) : normal(element_i)
            D[i, j] = -1/4π * Complex(n' * Dij)
        end
        D[i, i] = D[i, i] + 0.5
    end
    return S, D
end

function assemble_matrices_broadcasting(greens_functions, mesh, wavenumber; direct=true, arrtype=Array, include_identity=true)
    # Broadcasting variant: `arrtype=CuArray` (or another GPU array type) uploads
    # isbits StaticElements and compiles the Green's kernels as a GPU broadcast.
    # Mesh geometry is constant for frequency AD; dropping it from the tape avoids
    # a Zygote crash in `element(::StaticArraysMesh)` vertex indexing.
    elements = ChainRulesCore.ignore_derivatives() do
        arrtype([element(mesh, i) for i in 1:mesh.nfaces])
    end
    return assemble_matrices_broadcasting(greens_functions, elements, wavenumber; direct, include_identity)
end

# Autodispatch when the caller already holds a GPU (or other) vector of elements.
# Delhommeau `both_integral` S differs from `integral()` in the near field, so
# S and D stay on those kernels. Rankine CPU assemble uses `assemble_birk_rankine`.
function assemble_matrices_broadcasting(greens_functions, elements::AbstractVector, wavenumber; direct=true, include_identity=true)
    co_elements = reshape(elements, (1, length(elements)))
    S_kernel(e1, e2) = integral(greens_functions, e1, e2, wavenumber)
    S_matrix = (-1 / 4π) .* S_kernel.(elements, co_elements)

    if direct
        D_kernel(e1, e2) = normal(e2)' * integral_gradient(greens_functions, e1, e2, wavenumber; with_respect_to_first_variable=false)
        D_matrix = (-1 / 4π) .* D_kernel.(elements, co_elements)
    else
        K_kernel(e1, e2) = normal(e1)' * integral_gradient(greens_functions, e1, e2, wavenumber; with_respect_to_first_variable=true)
        D_matrix = (-1 / 4π) .* K_kernel.(elements, co_elements)
    end
    T = eltype(D_matrix)
    include_identity || return S_matrix, D_matrix
    return S_matrix, D_matrix + T(0.5) * I
end

# Vendor-agnostic panel-panel assembly. `get_backend(elements)` selects CPU,
# CUDA, Metal, ROC, or the JLArrays dummy GPU from the array type.
@kernel function _assemble_S_D!(S, D, @Const(elements), wavenumber, gfs, direct, include_identity)
    i, j = @index(Global, NTuple)
    @inbounds begin
        e1 = elements[i]
        e2 = elements[j]
        S[i, j] = -(1 / 4π) * integral(gfs, e1, e2, wavenumber)
        if direct
            n = normal(e2)
            grad = integral_gradient(gfs, e1, e2, wavenumber; with_respect_to_first_variable=false)
        else
            n = normal(e1)
            grad = integral_gradient(gfs, e1, e2, wavenumber; with_respect_to_first_variable=true)
        end
        D[i, j] = -(1 / 4π) * (n' * grad)
        if include_identity && i == j
            D[i, j] += convert(eltype(D), 0.5)
        end
    end
end

function assemble_matrices_ka(greens_functions, mesh, wavenumber; direct=true, arrtype=Array, include_identity=true)
    elements = arrtype([element(mesh, i) for i in 1:mesh.nfaces])
    return assemble_matrices_ka(greens_functions, elements, wavenumber; direct, include_identity)
end

function assemble_matrices_ka(greens_functions, elements::AbstractVector, wavenumber; direct=true, include_identity=true)
    backend = KernelAbstractions.get_backend(elements)
    n = length(elements)
    RT = typeof(float(real(wavenumber)))
    T = Complex{RT}
    S = KernelAbstractions.zeros(backend, T, n, n)
    D = KernelAbstractions.zeros(backend, T, n, n)
    kernel = _assemble_S_D!(backend)
    kernel(S, D, elements, wavenumber, greens_functions, direct, include_identity; ndrange=(n, n))
    _ka_synchronize(backend)
    return S, D
end

function _ka_synchronize(backend)
    try
        KernelAbstractions.synchronize(backend)
    catch e
        (e isa MethodError) || rethrow()
    end
end


"""
    assemble_matrices(green_functions, mesh, wavenumber; direct=true)

Assembles the influence matrices based on the tuple of provided Green's functions, mesh, and wavenumber.

# Arguments
- `green_functions`: Iterable of `GreensFunction` objects.
- `mesh`: Floating body mesh with panel information such as vertices, faces, normals, areas etc.
- `wavenumber`: Incoming ocean wavenumber
- `direct=true`: A flag to specify whether to use direct BEM vs Indirect BEM.

# Returns
- A tuple of assembled matrices. S and (D or K) depending on the flag.
"""
# Default: comprehension is Zygote-friendly on dense `Mesh`.
# `StaticArraysMesh` uses broadcasting so `arrtype=CuArray` (or a `CuArray` of
# elements) autodispatches to the GPU path.
assemble_matrices(green_functions, mesh, wavenumber; kwargs...) =
    assemble_matrices_comprehension(green_functions, mesh, wavenumber; kwargs...)

function assemble_matrices(green_functions, mesh::StaticArraysMesh, wavenumber; direct=true, arrtype=Array, all_normals=nothing, include_identity=true)
    if !isnothing(all_normals)
        error("all_normals is not supported on the broadcasting / GPU assembly path; use a dense Mesh.")
    end
    arrtype === Array || return assemble_matrices_broadcasting(green_functions, mesh, wavenumber; direct, arrtype, include_identity)
    return assemble_vectorized(green_functions, mesh, wavenumber; direct, include_identity)
end

# Split k-independent Rankine from the wave term, precompute source-panel Birk
# frames once per column, and evaluate S and D from the same pair kernel.
# Geometry is Float64 on this path; Zygote d/dω only needs the wave broadcast.
function assemble_vectorized(green_functions::Tuple{Rankine, RankineReflected, GFWu}, mesh::StaticArraysMesh, wavenumber; direct=true, include_identity=true)
    SD = ChainRulesCore.ignore_derivatives() do
        assemble_birk_rankine((Rankine(), RankineReflected()), mesh; direct, include_identity)
    end
    Sw, Dw = assemble_wu_centers(mesh, wavenumber, Val(direct); include_identity=false)
    return SD[1] .+ Sw, SD[2] .+ Dw
end

function assemble_vectorized(green_functions::Tuple{Rankine, RankineReflected}, mesh::StaticArraysMesh, wavenumber; direct=true, include_identity=true)
    return ChainRulesCore.ignore_derivatives() do
        assemble_birk_rankine(green_functions, mesh; direct, include_identity)
    end
end

function assemble_vectorized(green_functions, mesh::StaticArraysMesh, wavenumber; direct=true, include_identity=true)
    gfs = green_functions isa Tuple ? green_functions : (green_functions...,)
    static_gfs, wave_gfs = partition_greens_functions(gfs)
    SD_static = ChainRulesCore.ignore_derivatives() do
        isempty(static_gfs) ? nothing : assemble_static_vectorized(static_gfs, mesh; direct, include_identity)
    end
    if isempty(wave_gfs)
        return SD_static
    end
    wave_identity = isnothing(SD_static) && include_identity
    Sw, Dw = assemble_wave_vectorized(wave_gfs, mesh, wavenumber; direct, include_identity=wave_identity)
    isnothing(SD_static) && return Sw, Dw
    return SD_static[1] .+ Sw, SD_static[2] .+ Dw
end

function assemble_static_vectorized(gfs::Tuple, mesh::StaticArraysMesh; direct=true, include_identity=true)
    if _birk_rankine_tuple(gfs)
        return assemble_birk_rankine(gfs, mesh; direct, include_identity)
    end
    return assemble_matrices_broadcasting(gfs, mesh, 0.0; direct, include_identity)
end

_birk_rankine_tuple(gfs::Tuple) = all(_is_birk_rankine, gfs)
_is_birk_rankine(::Rankine) = true
_is_birk_rankine(::RankineReflected) = true
_is_birk_rankine(::RankineReflectedNegative) = true
_is_birk_rankine(::GreensFunction) = false

function assemble_wave_vectorized(gfs::Tuple{GFWu}, mesh::StaticArraysMesh, wavenumber; direct=true, include_identity=true)
    return assemble_wu_centers(mesh, wavenumber, Val(direct); include_identity)
end

function assemble_wave_vectorized(gfs::Tuple, mesh::StaticArraysMesh, wavenumber; direct=true, include_identity=true)
    length(gfs) == 1 && gfs[1] isa GFWu && return assemble_wu_centers(mesh, wavenumber, Val(direct); include_identity)
    return assemble_matrices_broadcasting(gfs, mesh, wavenumber; direct, include_identity)
end

@inline function _face_smatrix(mesh::StaticArraysMesh, j)
    f = mesh.faces[j]
    v1 = mesh.vertices[f[1]]
    v2 = mesh.vertices[f[2]]
    v3 = mesh.vertices[f[3]]
    v4 = mesh.vertices[f[4]]
    return SMatrix{4,3,Float64,12}(
        v1[1], v2[1], v3[1], v4[1],
        v1[2], v2[2], v3[2], v4[2],
        v1[3], v2[3], v3[3], v4[3],
    )
end

@inline _reflect_z(c::StaticVector{3}) = typeof(c)(c[1], c[2], -c[3])

@inline function _source_birk_geom(mesh::StaticArraysMesh, j)
    Tmat, qgc, local_corners, _ = birk_panel_geometry(_face_smatrix(mesh, j))
    return (T=Tmat, qgc=qgc, local_corners=local_corners)
end

function assemble_birk_rankine(gfs::Tuple, mesh::StaticArraysMesh; direct=true, include_identity=true)
    n = mesh.nfaces
    centers = mesh.centers
    normals = mesh.normals
    areas = mesh.areas
    radii = mesh.radii
    need_R = any(gf -> gf isa Rankine, gfs)
    need_I = any(gf -> gf isa RankineReflected, gfs)
    need_N = any(gf -> gf isa RankineReflectedNegative, gfs)
    image = need_I | need_N
    image_sign = (need_I ? 1.0 : 0.0) + (need_N ? -1.0 : 0.0)

    geoms = [_source_birk_geom(mesh, j) for j in 1:n]
    scale = -1 / 4π
    S = Matrix{ComplexF64}(undef, n, n)
    D = Matrix{ComplexF64}(undef, n, n)
    @inbounds for j in 1:n
        c2 = centers[j]
        r2 = radii[j]
        a2 = areas[j]
        n2 = normals[j]
        geom = geoms[j]
        Tmat, qgc, lc = geom.T, geom.qgc, geom.local_corners
        for i in 1:n
            c1 = centers[i]
            nvec = direct ? n2 : normals[i]
            s = 0.0
            d = 0.0
            if need_R
                φ, gsrc = _vrankine_from_geom(c1, c2, r2, a2, Tmat, qgc, lc)
                g = direct ? gsrc : -gsrc
                s += φ
                d += nvec' * _as_center_type(c1, g)
            end
            if image
                c1r = _reflect_z(c1)
                φr, gsrcr = _vrankine_from_geom(c1r, c2, r2, a2, Tmat, qgc, lc)
                if direct
                    gr = _as_center_type(c1r, gsrcr)
                else
                    gr = vertical_reflection(_as_center_type(c1r, -gsrcr))
                end
                s += image_sign * φr
                d += image_sign * (nvec' * gr)
            end
            S[i, j] = scale * s
            D[i, j] = scale * d
        end
    end
    if include_identity
        @inbounds for i in 1:n
            D[i, i] += 0.5
        end
    end
    return S, D
end

function assemble_wu_centers(mesh::StaticArraysMesh, wavenumber; direct=true, include_identity=true)
    assemble_wu_centers(mesh, wavenumber, Val(direct); include_identity)
end

function assemble_wu_centers(mesh::StaticArraysMesh, wavenumber, ::Val{direct}; include_identity=true) where {direct}
    _assemble_wu_centers_fused(mesh, wavenumber, Val(direct), include_identity)
end

# Fused S/D primal. ForwardDiff seeds Dual into `k` and runs this loop.
# Zygote uses the rrule below (mutation is not on the reverse tape).
function _assemble_wu_centers_fused(mesh, k, direct::Val, include_identity::Bool)
    _assemble_wu_centers_loop(mesh, k, direct, include_identity)
end

function _assemble_wu_centers_loop(mesh, k, ::Val{direct}, include_identity::Bool) where {direct}
    n = mesh.nfaces
    centers = mesh.centers
    areas = mesh.areas
    normals = mesh.normals
    RT = typeof(float(real(k)))
    T = Complex{RT}
    scale = convert(T, -1 / 4π)
    S = Matrix{T}(undef, n, n)
    D = Matrix{T}(undef, n, n)
    @inbounds for j in 1:n
        c2 = centers[j]
        a2 = areas[j]
        n2 = normals[j]
        for i in 1:n
            nvec = direct ? n2 : normals[i]
            sij, dij = _wu_sd_centers(centers[i], c2, nvec, a2, k, Val(direct))
            S[i, j] = scale * sij
            D[i, j] = scale * dij
        end
    end
    if include_identity
        half = convert(T, 0.5)
        @inbounds for i in 1:n
            D[i, i] += half
        end
    end
    return S, D
end

@inline function _real_inner(Ā, A)
    Ā isa AbstractZero && return zero(real(A[1])) * false
    s = zero(real(A[1])) * zero(real(first(Ā)))
    @inbounds for i in eachindex(A)
        gi = Ā[i]
        ai = A[i]
        s += real(gi) * real(ai) + imag(gi) * imag(ai)
    end
    return s
end

function ChainRulesCore.rrule(::typeof(_assemble_wu_centers_fused), mesh, k, direct::Val, include_identity::Bool)
    y = _assemble_wu_centers_loop(mesh, k, direct, include_identity)
    function assemble_wu_centers_pullback(ȳ)
        ȳu = unthunk(ȳ)
        if ȳu isa AbstractZero
            return (NoTangent(), NoTangent(), zero(k), NoTangent(), NoTangent())
        end
        S̄ = unthunk(ȳu[1])
        D̄ = unthunk(ȳu[2])
        ∂k = ForwardDiff.derivative(k) do κ
            S, D = _assemble_wu_centers_loop(mesh, κ, direct, include_identity)
            _real_inner(S̄, S) + _real_inner(D̄, D)
        end
        return (NoTangent(), NoTangent(), ∂k, NoTangent(), NoTangent())
    end
    return y, assemble_wu_centers_pullback
end


function assemble_matrix_wu(mesh, wavenumber; direct=true, all_normals=nothing)
    return assemble_matrices((Rankine(), RankineReflected(), GFWu()), mesh, wavenumber; direct, all_normals)
end


# ImplicitAD is CPU / Zygote reverse-mode. Dual numbers use native `\`
# (ForwardDiff). GPU arrays use `gpu_linsolve` (cuSOLVER + explicit rrule).
# Do not replace `implicit_linear` with `\` for Array — Zygote cannot
# reverse dense LU without that implicit function theorem wrapper.
function linsolve(A, b)
    if eltype(real(A)) <: ForwardDiff.Dual
        return A \ b
    elseif A isa Array
        return implicit_linear(A, b)
    else
        return gpu_linsolve(A, b)
    end
end

gpu_linsolve(A, b) = _backend_ldiv(A, b)

# Prefer native `\` (cuSOLVER). Dummy GPU arrays (JLArray) and some Metal
# types hit generic LAPACK via scalar indexing, so factorize on the host.
function _backend_ldiv(A, b)
    if A isa Array || _uses_device_ldiv(A)
        return A \ b
    end
    x = Array(A) \ Array(b)
    out = similar(b, eltype(x), size(x))
    copyto!(out, x)
    return out
end

function _uses_device_ldiv(A)
    return string(nameof(typeof(A))) == "CuArray"
end

function ChainRulesCore.rrule(::typeof(gpu_linsolve), A, b)
    x = _backend_ldiv(A, b)
    function gpu_linsolve_pullback(ȳ)
        ȳu = unthunk(ȳ)
        if ȳu isa AbstractZero
            return (NoTangent(), ZeroTangent(), ZeroTangent())
        end
        λ = _backend_ldiv(A', ȳu)
        return (NoTangent(), -λ * x', λ)
    end
    return x, gpu_linsolve_pullback
end


function solve(D, S, bc; direct::Bool=true)
    if direct
        ϕ = linsolve(D,S*bc)
        sources = nothing
    else
        K = D
        sources = linsolve(K,bc)
        ϕ = S * sources
    end
    return ϕ, sources
end

