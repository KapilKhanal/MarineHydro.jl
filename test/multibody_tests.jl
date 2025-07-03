using Test
using Revise
using MarineHydro
using PyCall

# Create two spheres as in paper/Multibody.jl
cpt = pyimport("capytaine")
cptmesh = cpt.mesh_sphere(radius=1.0, center=(0, 0, 0), resolution=(10, 10)).immersed_part()
mesh1 = Mesh(cptmesh)
mesh2 = Mesh(cptmesh.translated_x(5.0))


# Body 1: Heave
dofs1 = [DOF(:heave, [0.0, 0.0 ,1.0])]    
body1 = FloatingBody(mesh1, "sphere1", dofs1)

# Body 2: Heave only
dofs2 = [DOF(:heave, [0.0, 0.0, 1.0])]  # Heave only
body2 = FloatingBody(mesh2, "sphere2", dofs2)

mb = MultiBody([body1, body2])

omega = 2.03

A, B = compute_multibody_radiation(mb, omega)

println("Added mass matrix shape: ", size(A))
println("Damping matrix shape: ", size(B))

# Print some key elements
println("A ", A)
println("B ", B)

# Symmetry checks for identical spheres (heave-heave coupling)
@test isapprox(A[1,1,1,1], A[2,1,2,1]; atol=1e-6)  # Self-heave should be same for both bodies
@test isapprox(B[1,1,1,1], B[2,1,2,1]; atol=1e-6)
@test isapprox(A[1,1,2,1], A[2,1,1,1]; atol=1e-6)  # Cross-heave coupling should be symmetric
@test isapprox(B[1,1,2,1], B[2,1,1,1]; atol=1e-6)

println("Added mass matrix: ", A)
println("Damping matrix: ", B) 