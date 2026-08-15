using MarineHydro
using PyCall
using Test

mesh = MarineHydro.Mesh(MarineHydro.example_mesh_from_capytaine())
smesh = MarineHydro.StaticArraysMesh(MarineHydro.example_mesh_from_capytaine())
greens_functions = (Rankine(), RankineReflected(), GFWu())

@testset "Matrix shape" begin
    wavenumber = 1.0
    ω = √(wavenumber*9.8)
    dof = [0,0,1]
    green_functions = (
        Rankine(),
        RankineReflected(),
        #=GFWu()=#
    )
    S, D = assemble_matrices(green_functions, mesh, 1.0);
    S_, K = assemble_matrices(green_functions, mesh, 1.0; direct=false);
    @test size(S) == size(S_) == size(D) == size(K)
    @test S ≈ S_
    @test !(D ≈ K)
end

S, D = MarineHydro.assemble_matrices_comprehension(greens_functions, mesh, 1.0)
S_, D_ = MarineHydro.assemble_matrices_broadcasting(greens_functions, mesh, 1.0)
S__, D__ = MarineHydro.assemble_matrices_broadcasting(greens_functions, smesh, 1.0)
elements = [element(smesh, i) for i in 1:smesh.nfaces]
S_el, D_el = MarineHydro.assemble_matrices_broadcasting(greens_functions, elements, 1.0)
S_sa, D_sa = assemble_matrices(greens_functions, smesh, 1.0)
@test S ≈ S_ ≈ S__ ≈ S_el ≈ S_sa
@test D ≈ D_ ≈ D__ ≈ D_el ≈ D_sa


