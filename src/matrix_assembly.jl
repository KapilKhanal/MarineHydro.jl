# Below are several variants of the same function.
# TODO: refactor with multiple dispatch?

function assemble_matrices_comprehension(green_functions, mesh, wavenumber; direct=true, all_normals=nothing)
    
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
            (constant * (n' * norm_vec)) - 1/2τ̅ * Complex{T}(n' * integral_gradient(green_functions, element_i, element_j, wavenumber; with_respect_to_first_variable=!direct))
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

function assemble_matrices_broadcasting(greens_functions, mesh, wavenumber; direct=true, arrtype=Array)
    # Broadcasting variant: `arrtype=CuArray` (or another GPU array type) uploads
    # isbits StaticElements and compiles the Green's kernels as a GPU broadcast.
    elements = arrtype([element(mesh, i) for i in 1:mesh.nfaces])
    return assemble_matrices_broadcasting(greens_functions, elements, wavenumber; direct)
end

# Autodispatch when the caller already holds a GPU (or other) vector of elements.
function assemble_matrices_broadcasting(greens_functions, elements::AbstractVector, wavenumber; direct=true)
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
    return S_matrix, D_matrix + T(0.5) * one(D_matrix)
end

# Vendor-agnostic panel-panel assembly. `get_backend(elements)` selects CPU,
# CUDA, Metal, ROC, or the JLArrays dummy GPU from the array type.
@kernel function _assemble_S_D!(S, D, @Const(elements), wavenumber, gfs, direct)
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
        if i == j
            D[i, j] += convert(eltype(D), 0.5)
        end
    end
end

function assemble_matrices_ka(greens_functions, mesh, wavenumber; direct=true, arrtype=Array)
    elements = arrtype([element(mesh, i) for i in 1:mesh.nfaces])
    return assemble_matrices_ka(greens_functions, elements, wavenumber; direct)
end

function assemble_matrices_ka(greens_functions, elements::AbstractVector, wavenumber; direct=true)
    backend = KernelAbstractions.get_backend(elements)
    n = length(elements)
    RT = typeof(float(real(wavenumber)))
    T = Complex{RT}
    S = KernelAbstractions.zeros(backend, T, n, n)
    D = KernelAbstractions.zeros(backend, T, n, n)
    kernel = _assemble_S_D!(backend)
    kernel(S, D, elements, wavenumber, greens_functions, direct; ndrange=(n, n))
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

function assemble_matrices(green_functions, mesh::StaticArraysMesh, wavenumber; direct=true, arrtype=Array, all_normals=nothing)
    if !isnothing(all_normals)
        error("all_normals is not supported on the broadcasting / GPU assembly path; use a dense Mesh.")
    end
    return assemble_matrices_broadcasting(green_functions, mesh, wavenumber; direct, arrtype)
end


function assemble_matrix_wu(mesh, wavenumber; direct=true, all_normals=nothing)
    return assemble_matrices((Rankine(), RankineReflected(), GFWu()), mesh, wavenumber; direct, all_normals)
end


# ImplicitAD is CPU / Zygote oriented. Dual numbers and GPU arrays (CuArray, …)
# use native `\` (cuSOLVER on CUDA). GPU reverse-mode goes through `gpu_linsolve`.
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

