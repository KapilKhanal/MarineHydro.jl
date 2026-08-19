# Old functions

function radiation_bc(mesh::Mesh, dof, omega)
    """
        radiation_bc(mesh::Mesh, dof, omega)

    Calculates the radiation boundary conditions for floating bodies at each panel.

    # Arguments
    - `mesh::Mesh`: The mesh of the floating body.
    - `dof`: The degrees of freedom (assumed same for each panel).
    - `omega`: The frequency of the incident ocean wave ~~~.

    # Returns
    - The (Neumann) radiation boundary condition values for each panel.
"""
    nfaces = mesh.nfaces
    d1 = dof[1]
    d2 = dof[2]
    d3 = dof[3]
    scale = -1im * omega
    return @inbounds [scale * (mesh.normals[i, 1] * d1 + mesh.normals[i, 2] * d2 + mesh.normals[i, 3] * d3) for i in 1:nfaces]
end

function radiation_bc(mesh::StaticArraysMesh, dof, omega)
    d1, d2, d3 = dof[1], dof[2], dof[3]
    scale = -1im * omega
    return [scale * (n[1] * d1 + n[2] * d2 + n[3] * d3) for n in mesh.normals]
end

function integrate_pressure(mesh::Mesh, pressure, dof)
    normal_dof_amp = -sum(transpose(dof) .* mesh.normals, dims=2)
    forces = sum(pressure .* normal_dof_amp .* mesh.areas)
    return forces
end

function integrate_pressure(mesh::StaticArraysMesh, pressure, dof)
    d1, d2, d3 = dof[1], dof[2], dof[3]
    s = zero(promote_type(eltype(pressure), eltype(mesh.areas)))
    @inbounds for i in 1:mesh.nfaces
        n = mesh.normals[i]
        s += pressure[i] * -(d1 * n[1] + d2 * n[2] + d3 * n[3]) * mesh.areas[i]
    end
    return s
end

function calculate_radiation_forces(mesh, dof, omega)
    k = omega^2 / SETTINGS.g
    S, D = assemble_matrix_wu(mesh, k)
    bc = radiation_bc(mesh, dof, omega)
    potential, sources = solve(D, S, bc)
    pressure = 1im * SETTINGS.rho * omega * potential
    forces = integrate_pressure(mesh, pressure, dof)
    return [real(forces)/omega^2, imag(forces)/omega]
end


function DiffractionForce(mesh::Mesh,ω,dof)
    green_functions = (
        Rankine(),
        RankineReflected(),
        GFWu(),
    )
    k = ω^2 / SETTINGS.g
    S, D = assemble_matrices(green_functions, mesh, k)
    bc = AiryBC(mesh, ω)
    potential, sources = solve(D, S, bc)
    forces = diffraction_force(potential,mesh, ω,dof)
    return forces
end

function airy_waves_potential(points, omega, beta=0)
    wavenumber = compute_wavenumber(omega)
    x, y, z = points[:, 1], points[:, 2], points[:, 3]
    wbar = x .* cos(beta) .+ y .* sin(beta)
    cih = exp.(wavenumber .* z)
    phi = -1im*SETTINGS.g/omega .* cih .* exp.(1im * wavenumber * wbar)
    return phi
end

function airy_waves_velocity(points, omega, beta=0, water_depth = Inf)
    """Compute the fluid velocity for Airy waves at a given point (or array of points)."""
    k = compute_wavenumber(omega)

    x, y, z = points[:, 1], points[:, 2], points[:, 3]

    wbar = x .* cos(beta) .+ y .* sin(beta)
    cih = exp.(k .* z)
    sih = exp.(k .* z)

    v = SETTINGS.g * k / omega .* exp.(1im * k .* wbar) .* 
        hcat(cos(beta) .* cih, sin(beta) .* cih, -1im .* sih)
    return v
end

function AiryBC(mesh,omega,beta=0)
    """Boundary condition for diffraction problem : the velocity on the floating body is the velocity of Airy wave field."""
    velocity = airy_waves_velocity(mesh.centers, omega, beta)
    nfaces = Int(mesh.nfaces)
    return @inbounds [-(velocity[i, 1] * mesh.normals[i, 1] + velocity[i, 2] * mesh.normals[i, 2] + velocity[i, 3] * mesh.normals[i, 3]) for i in 1:nfaces]
end

function airy_waves_pressure(points, omega, beta=0)
    """Compute the pressure for Airy waves."""  

    return 1im .* omega .* SETTINGS.rho .* airy_waves_potential(points, omega, beta)
end


function diffraction_force(potential,mesh, omega,dof)
    pressure = 1im*SETTINGS.rho* potential * omega 
    forces = integrate_pressure(mesh,pressure,dof) 
    return forces  
end

function FroudeKrylovForce(mesh::Mesh, ω,dof)
    """Compute the Froude-Krylov force."""
    pressure =  airy_waves_pressure(mesh.centers,  ω)
    return  integrate_pressure(mesh::Mesh, pressure, dof) 
end



# New functions

function compute_wavenumber(omega)
    return omega^2 / SETTINGS.g
end

function compute_encountered_values(omega, beta, forward_speed)
    k = compute_wavenumber(omega)    
    if forward_speed==0
        return omega, k, beta
    else
        doppler_omega = omega - k * forward_speed * cos(beta)
        encountered_omega = abs(doppler_omega)
        encountered_wavenumber = compute_wavenumber(encountered_omega)
        if doppler_omega >= 0
            encountered_wave_direction = beta
        else
            encountered_wave_direction = beta + pi
        end
        return encountered_omega, encountered_wavenumber, encountered_wave_direction
    end
end

function integrate_pressure(floatingbody::FloatingBody, influenced_dofs::Vector{Symbol}, pressure)
    mesh = floatingbody.mesh

    force_values = [
        _integrate_dof_pressure(mesh, floatingbody.dofs[dof_symbol], pressure)
        for dof_symbol in influenced_dofs
    ]
    forces = NamedTuple{Tuple(influenced_dofs)}(Tuple(force_values))
    return forces
end

function _integrate_dof_pressure(mesh::Mesh, dof_mat, pressure)
    normal_dof_amp_on_face = -sum(dof_mat .* mesh.normals, dims=2)
    return sum(pressure .* normal_dof_amp_on_face .* mesh.areas)
end

function _integrate_dof_pressure(mesh::StaticArraysMesh, dof_mat, pressure)
    s = zero(promote_type(eltype(pressure), eltype(mesh.areas), eltype(dof_mat)))
    @inbounds for i in 1:mesh.nfaces
        n = mesh.normals[i]
        disp = dof_mat[i, 1] * n[1] + dof_mat[i, 2] * n[2] + dof_mat[i, 3] * n[3]
        s += pressure[i] * -disp * mesh.areas[i]
    end
    return s
end

################################ Radiation methods #########################################


function radiation_bc(problem::RadiationProblem)
    """
    Calculates the radiation boundary conditions for floating bodies at each panel.

    # Returns
    - The (Neumann) radiation boundary condition values for each panel.
"""
    mesh = problem.floatingbody.mesh
    dof_mat = problem.floatingbody.dofs[problem.radiating_dof]
    displacement_on_face = _normal_displacement(mesh, dof_mat)
    bc =  -1im .* problem.encountered_omega .* displacement_on_face
    if problem.forward_speed!=0
        ddofdx = evaluate_gradient_of_motion(mesh, string(problem.radiating_dof))
        bc = bc .+ problem.forward_speed .* _normal_displacement(mesh, ddofdx)
    end
    return bc
end

_normal_displacement(mesh::Mesh, dof_mat) = sum(mesh.normals .* dof_mat, dims=2)

function _normal_displacement(mesh::StaticArraysMesh, dof_mat)
    T = promote_type(eltype(mesh.areas), eltype(dof_mat))
    out = Vector{T}(undef, mesh.nfaces)
    @inbounds for i in 1:mesh.nfaces
        n = mesh.normals[i]
        out[i] = n[1] * dof_mat[i, 1] + n[2] * dof_mat[i, 2] + n[3] * dof_mat[i, 3]
    end
    return out
end

function compute_bc(problem::RadiationProblem)
    return radiation_bc(problem)
end

################################ Diffraction and Excitation methods #########################################


function airy_waves_potential(points, problem::DiffractionProblem)
    x, y, z = points[:, 1], points[:, 2], points[:, 3]

    wavenumber = problem.wavenumber
    enc_beta = problem.encountered_beta
    wbar = x .* cos(enc_beta) .+ y .* sin(enc_beta)    
    
    cih = exp.(wavenumber .* z)
    phi = -1im*SETTINGS.g/problem.omega .* cih .* exp.(1im * wavenumber * wbar)
    return phi
end

function airy_waves_velocity(points, problem::DiffractionProblem)
    """Compute the fluid velocity for Airy waves at a given point (or array of points)."""

    x, y, z = points[:, 1], points[:, 2], points[:, 3]
    k = problem.wavenumber
    enc_beta = problem.encountered_beta

    wbar = x .* cos(enc_beta) .+ y .* sin(enc_beta)

    cih = exp.(k .* z)
    sih = exp.(k .* z)

    v = SETTINGS.g * k / problem.omega .* 
        exp.(1im * k .* wbar) .* 
        hcat(cos(problem.beta) .* cih, sin(problem.beta) .* cih, -1im .* sih)
    return v
end


#boundary conditions from airy wave for solving diffraction problem
function AiryBC(problem::DiffractionProblem)
    """Boundary condition for diffraction problem : the velocity on the floating body is the velocity of Airy wave field."""
    bcs = -sum(airy_waves_velocity(problem.floatingbody.mesh.centers, problem) .* problem.floatingbody.mesh.normals, dims = 2)
    return bcs
end

function compute_bc(problem::DiffractionProblem)
    return AiryBC(problem)
end


function airy_waves_pressure(points, problem::DiffractionProblem)
    """Compute the pressure for Airy waves."""  

    return 1im .* problem.omega .* SETTINGS.rho .* airy_waves_potential(points, problem)
end

function FroudeKrylovForce(problem::DiffractionProblem, influenced_dofs::Vector{Symbol})
    """Compute the Froude-Krylov force."""

    mesh = problem.floatingbody.mesh
    pressure =  airy_waves_pressure(mesh.centers,  problem)
    forces = integrate_pressure(problem.floatingbody, influenced_dofs, pressure) 
    return forces 
end
