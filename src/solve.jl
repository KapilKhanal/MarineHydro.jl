
# Place hull BCs onto hull+lid without setindex! (Zygote) and without BLAS
# `ger` (Enzyme reverse cannot load `zger_64_` on some platforms).
function _pad_hull_bc(bc_on_hull, hull_mask)
    v = vec(bc_on_hull)
    n = length(hull_mask)
    length(v) == n && return v
    T = eltype(v)
    idx = findall(hull_mask)
    return T[sum(i == idx[k] ? v[k] : zero(T) for k in eachindex(idx)) for i in 1:n]
end

function _problem_omega_k(problem)
    if problem.forward_speed == 0
        return problem.omega, problem.wavenumber
    else
        return problem.encountered_omega, problem.encountered_wavenumber
    end
end

function _mesh_including_lid(problem)
    fb = problem.floatingbody
    if isnothing(fb.lid_mesh)
        return fb.mesh, trues(fb.mesh.nfaces)
    end
    mesh = fb.mesh + fb.lid_mesh
    hull_mask = eachrow(mesh.faces) .∈ Ref(eachrow(fb.mesh.faces))
    return mesh, hull_mask
end

# Use the mesh already on the body. `FloatingBody{<:Mesh}` stays on the dense
# comprehension path; `FloatingBody{<:StaticArraysMesh}` uses the vectorized
# kernels. Convert once when you build the body — not inside the solve.
_assembly_mesh(mesh) = mesh

# Add a precomputed Rankine (k-independent) pair to the wave-only assembly.
# Identity / free-surface jump lives on the static matrices so it is not applied twice.
function _add_wave_matrices(static_SD, wave_gfs, mesh, wavenumber; direct, all_normals=nothing)
    if isempty(wave_gfs)
        return static_SD
    end
    include_identity = isnothing(static_SD)
    amesh = isnothing(all_normals) ? _assembly_mesh(mesh) : mesh
    Sw, Dw = assemble_matrices(wave_gfs, amesh, wavenumber; direct, all_normals, include_identity)
    isnothing(static_SD) && return Sw, Dw
    return static_SD[1] .+ Sw, static_SD[2] .+ Dw
end

# `solve(prob)` / `solve(prob, DirectBEM())` — same name as the linear-system
# `solve(D, S, bc)`, dispatched on the first argument.
function solve(problem::LinearPotentialFlowProblem, alg::BEMAlgorithm=DirectBEM();
        static_SD=nothing, static_SD_normals=nothing, wave_gfs=nothing)
    omega, wavenumber = _problem_omega_k(problem)
    return _solve_problem(problem, omega, wavenumber;
        direct=is_direct(alg), greens_functions=alg.greens,
        static_SD, static_SD_normals, wave_gfs)
end

function solve_problem(problem::LinearPotentialFlowProblem;
        alg=nothing, direct::Bool=true, gf::String="Wu", greens_functions=nothing,
        static_SD=nothing, static_SD_normals=nothing, wave_gfs=nothing)
    return solve(problem, bem_algorithm(; alg, direct, gf, greens_functions);
        static_SD, static_SD_normals, wave_gfs)
end

@noinline function _solve_problem(problem, omega::T, wavenumber::T;
        direct, greens_functions, static_SD, static_SD_normals, wave_gfs) where T
    bc_on_hull = vec(compute_bc(problem))
    gfs = greens_functions
    mesh_including_lid, hull_mask = _mesh_including_lid(problem)
    bc = _pad_hull_bc(bc_on_hull, hull_mask)

    if isnothing(static_SD)
        S, D = assemble_matrices(gfs, _assembly_mesh(mesh_including_lid), wavenumber; direct=direct)
    else
        S, D = _add_wave_matrices(static_SD, wave_gfs, mesh_including_lid, wavenumber; direct)
    end

    potential, sources = solve(D, S, bc; direct=direct)

    pressure = 1im * SETTINGS.rho * omega * potential

    pressure_on_hull = pressure[hull_mask]

    if problem.forward_speed!=0
        if isnothing(static_SD_normals)
            S, D = assemble_matrices(gfs, mesh_including_lid, wavenumber; direct=direct, all_normals=[1,0,0])
        else
            S, D = _add_wave_matrices(static_SD_normals, wave_gfs, mesh_including_lid, wavenumber; direct, all_normals=[1,0,0])
        end
        if direct
            error("MarineHydro.jl has yet to be developed for nonzero forward speeds 
            with the direct method. Try changing direct to false.")
        else
            K = D
            partial_phi_partial_x = K * sources
        end
        additional_pressure = SETTINGS.rho * problem.forward_speed * partial_phi_partial_x
        pressure_on_hull = pressure_on_hull + additional_pressure[hull_mask]
    end

    forces = integrate_pressure(problem.floatingbody, problem.influenced_dofs, pressure_on_hull)

    return make_result(problem, forces)
end

function _same_geometry(problems)
    fb0 = problems[1].floatingbody
    return all(p -> p.floatingbody.mesh === fb0.mesh && p.floatingbody.lid_mesh === fb0.lid_mesh, problems)
end

# Rankine is assembled once per mesh and reused across wavenumbers.
function solve(problems::AbstractVector{<:LinearPotentialFlowProblem}, alg::BEMAlgorithm=DirectBEM())
    isempty(problems) && return LinearPotentialFlowResult[]
    gfs = alg.greens
    direct = is_direct(alg)
    static_gfs, wave_gfs = partition_greens_functions(gfs)

    if isempty(static_gfs) || !_same_geometry(problems)
        return [solve(problem, alg) for problem in problems]
    end

    mesh, _ = _mesh_including_lid(problems[1])
    static_SD = assemble_matrices(static_gfs, _assembly_mesh(mesh), 1.0; direct)
    static_SD_normals = any(p -> p.forward_speed != 0, problems) ?
        assemble_matrices(static_gfs, mesh, 1.0; direct, all_normals=[1,0,0]) : nothing

    return [solve(problem, alg; static_SD, static_SD_normals, wave_gfs) for problem in problems]
end

function solve_all_problems(problems::AbstractVector{<:LinearPotentialFlowProblem};
        alg=nothing, direct::Bool=true, gf::String="Wu", greens_functions=nothing)
    return solve(problems, bem_algorithm(; alg, direct, gf, greens_functions))
end
