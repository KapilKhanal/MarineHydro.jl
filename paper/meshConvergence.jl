using PyCall, MarineHydro, Zygote, Plots, ColorTypes

orange = RGB(230/255,159/255,0/255)  
vermillion = RGB(213/255, 94/255, 0/255) 
bluishgreen = RGB(0/255, 158/255, 115/255) 

include("BemProgram.jl")
include("MeshGradients_singlebody.jl")


function radiation_forces(mesh::Mesh, dof, omega; direct = true)
    k = omega^2 / 9.8
    S, D = assemble_matrix_wu(mesh, k;direct)
    bc = radiation_bc(mesh, dof, omega)
    potential = solve(D, S, bc)
    pressure = 1im * 1023 * omega * potential
    forces = integrate_pressure(mesh, pressure, dof)
    return [real(forces)/omega^2, imag(forces)/omega]
end


using Plots

resolutions = [(6,6), (8,8), (10,10), (12,12),(14,14), (16,16), (18,18), (20,20)]
radius = 1.0
omega = 1.03
dof = [0.0, 0.0, 1.0]

panel_numbers = []
added_mass_values_direct = []
added_mass_values_indirect = []
for res in resolutions
    cptmesh = cpt.mesh_sphere(name="sphere", radius=radius, center=(0, 0, 0), resolution=res)
    cptmesh.keep_immersed_part(inplace=true)
    mesh = Mesh(cptmesh)
    push!(panel_numbers, mesh.nfaces)
    print(res)
    print(mesh.nfaces)
    direct = true
    push!(added_mass_values_direct, radiation_forces(mesh, dof,omega; direct)[1])
    direct  = false
    push!(added_mass_values_indirect, radiation_forces(mesh, dof,omega; direct)[1])
end

plot(panel_numbers, added_mass_values_direct, marker=:o, xlabel="Number of Panels", ylabel="Added Mass", linecolor = bluishgreen, label ="Direct BEM",markercolor = bluishgreen)
plot!(panel_numbers, added_mass_values_indirect, marker=:o, xlabel="Number of Panels", ylabel="Added Mass",linecolor = vermillion,label ="Indirect BEM",markercolor = vermillion)
savefig("MarineHydro.jl/paper/Plots/mesh_convergence_study.pdf")