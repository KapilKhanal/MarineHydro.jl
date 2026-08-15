

# Place hull BCs onto hull+lid without setindex! (Zygote cannot mutate arrays).
function _pad_hull_bc(bc_on_hull, hull_mask)
    v = vec(bc_on_hull)
    n = length(hull_mask)
    length(v) == n && return v
    idx = findall(hull_mask)
    return Matrix{eltype(v)}(I, n, n)[:, idx] * v
end

function default_greens_functions(gf::String)
    if gf == "Wu"
        wave_gf = GFWu()
    elseif gf == "ExactGuevelDelhommeau"
        wave_gf = ExactGuevelDelhommeau()
    else
        error("Unknown Green's function \"$gf\". Use \"Wu\" or \"ExactGuevelDelhommeau\".")
    end
    return (Rankine(), RankineReflected(), wave_gf)
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

# Primal Float64 solves use StaticArraysMesh broadcasting (Capytaine-class panel loops).
# Dual / Number meshes stay dense so Zygote can see `mesh.vertices`.
# `all_normals` (forward speed) is only implemented on the dense path.
_assembly_mesh(mesh::Mesh) = ChainRulesCore.ignore_derivatives() do
    _concrete_float_mesh(mesh) ? StaticArraysMesh(mesh) : mesh
end
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

# Solve single problem (one frequency and one radiating dof or wave direction).
# Problem structs are parameterized on the scalar type of ω/k (Float64, Dual, …).
function solve_problem(problem::LinearPotentialFlowProblem; direct::Bool=true, gf::String="Wu", greens_functions=nothing,
        static_SD=nothing, static_SD_normals=nothing, wave_gfs=nothing)
    omega, wavenumber = _problem_omega_k(problem)
    return _solve_problem(problem, omega, wavenumber; direct, gf, greens_functions, static_SD, static_SD_normals, wave_gfs)
end

@noinline function _solve_problem(problem, omega::T, wavenumber::T; direct, gf, greens_functions, static_SD, static_SD_normals, wave_gfs) where T
    # Apply boundary conditions
    bc_on_hull = vec(compute_bc(problem))
    gfs = isnothing(greens_functions) ? default_greens_functions(gf) : greens_functions
    mesh_including_lid, hull_mask = _mesh_including_lid(problem)
    bc = _pad_hull_bc(bc_on_hull, hull_mask)

    if isnothing(static_SD)
        S, D = assemble_matrices(gfs, _assembly_mesh(mesh_including_lid), wavenumber; direct=direct)
    else
        S, D = _add_wave_matrices(static_SD, wave_gfs, mesh_including_lid, wavenumber; direct)
    end

    # Solve linear system (include lid)
    potential, sources = solve(D, S, bc; direct=direct)

    # Compute pressure (include lid)
    pressure = 1im * SETTINGS.rho * omega * potential # uses encountered_omega

    # Slice pressure (excludes lid)
    pressure_on_hull = pressure[hull_mask]

    if problem.forward_speed!=0
        # change normals to all be unit vector in x direction
        if isnothing(static_SD_normals)
            S, D = assemble_matrices(gfs, mesh_including_lid, wavenumber; direct=direct, all_normals=[1,0,0])
        else
            S, D = _add_wave_matrices(static_SD_normals, wave_gfs, mesh_including_lid, wavenumber; direct, all_normals=[1,0,0])
        end
        if direct
            error("MarineHydro.jl has yet to be developed for nonzero forward speeds 
            with the direct method. Try changing direct to false.")
            # partial_phi_partial_x = S \ (D * potential)
            # sources = S \ potential
            # partial_phi_partial_x = K * sources
        else
            K = D
            partial_phi_partial_x = K * sources
        end
        additional_pressure = SETTINGS.rho * problem.forward_speed * partial_phi_partial_x
        pressure_on_hull = pressure_on_hull + additional_pressure[hull_mask]
    end

    forces = integrate_pressure(problem.floatingbody, problem.influenced_dofs, pressure_on_hull) # NamedTuple of complex forces, where each element corresponds to an influenced dof

    result = make_result(problem, forces)
    return result
end

function _same_geometry(problems)
    fb0 = problems[1].floatingbody
    return all(p -> p.floatingbody.mesh === fb0.mesh && p.floatingbody.lid_mesh === fb0.lid_mesh, problems)
end

# Solve multiple problems (multiple frequencies, radiating dofs, and/or wave directions)
# Equivalent to Capytaine's solve_all() function. Rankine is assembled once per mesh
# and reused across wavenumbers; the wave Green function is rebuilt each time.
# The cached Rankine matrices are read-only, so a later threaded map is safe.
function solve_all_problems(problems::AbstractVector{<:LinearPotentialFlowProblem}; direct::Bool=true, gf::String="Wu", greens_functions=nothing)
    isempty(problems) && return LinearPotentialFlowResult[]
    gfs = isnothing(greens_functions) ? default_greens_functions(gf) : greens_functions
    static_gfs, wave_gfs = partition_greens_functions(gfs)

    if isempty(static_gfs) || !_same_geometry(problems)
        return [solve_problem(problem; direct, gf, greens_functions=gfs) for problem in problems]
    end

    mesh, _ = _mesh_including_lid(problems[1])
    static_SD = assemble_matrices(static_gfs, _assembly_mesh(mesh), 1.0; direct)
    static_SD_normals = any(p -> p.forward_speed != 0, problems) ?
        assemble_matrices(static_gfs, mesh, 1.0; direct, all_normals=[1,0,0]) : nothing

    return [solve_problem(problem; direct, gf, greens_functions=gfs,
        static_SD, static_SD_normals, wave_gfs) for problem in problems]
end