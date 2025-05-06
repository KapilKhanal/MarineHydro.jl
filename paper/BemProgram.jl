

function added_mass_program(radius,omega,dof)  
    mesh = differentiableMesh(radius) #fd
    A = calculate_radiation_forces(mesh,dof,omega)[1]
    return A
end

function damping_program(radius,omega,dof)  
    mesh = differentiableMesh(radius) #fd
    B = calculate_radiation_forces(mesh,dof,omega)[2]
    return B
end

function diffraction_program(ω,radius,dof)  
    mesh = differentiableMesh(radius) #fd
    force = DiffractionForce(mesh,ω,dof)
    return force
end

function omega_added_mass_bem_program(omega,radius=1,dof = [0,0,1])  
    mesh = differentiableMesh(radius) #fd
    A = calculate_radiation_forces(mesh,dof,omega)[1]
    return A
end

function omega_damping_bem_program(omega,radius=1,dof = [0,0,1])  
    mesh = differentiableMesh(radius) #fd
    B = calculate_radiation_forces(mesh,dof,omega)[2]
    return B
end


function rankine_program(radius,dof = [1,0,0])  
    mesh = differentiableMesh(radius) #fd
    S,D = assemble_matrices([Rankine(), RankineReflected()], mesh, 1.0)
    BC = radiation_bc(mesh, dof, 1.0)
    ϕ = implicit_linear(D,S*BC)
    pressure = 1im * 1023 * omega * ϕ
    forces = integrate_pressure(mesh, pressure, dof)
    return real(forces)
end
