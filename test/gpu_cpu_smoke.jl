# CPU-only smoke of the GPU-generic broadcasting / linsolve path (no CUDA, no Capytaine).
using Test
using LinearAlgebra
using StaticArrays
using ForwardDiff
using Zygote
using Enzyme
using DifferentiationInterface
using MarineHydro

const MH = MarineHydro

vertices = [
    SVector(-1.0, -1.0, -1.0),
    SVector( 0.0, -1.0, -1.0),
    SVector( 0.0,  0.0, -1.0),
    SVector(-1.0,  0.0, -1.0),
    SVector( 1.0, -1.0, -1.0),
    SVector( 1.0,  0.0, -1.0),
]
faces = [SVector(1, 2, 3, 4), SVector(2, 5, 6, 3)]
centers = [SVector(-0.5, -0.5, -1.0), SVector(0.5, -0.5, -1.0)]
normals = [SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 1.0)]
smesh = MH.StaticArraysMesh(vertices, faces, centers, normals, [1.0, 1.0], [sqrt(0.5), sqrt(0.5)], 6, 2)

gfs = (Rankine(), RankineReflected(), GFWu())
gfs_v = (VRankine(), VRankineReflected(), GFWu())
k = 1.0
ω = 1.0
bc = [-1im * ω * n[3] for n in smesh.normals]

@testset "CPU broadcasting autodispatch" begin
    S, D = MH.assemble_matrices_broadcasting(gfs, smesh, k)
    elements = [element(smesh, i) for i in 1:smesh.nfaces]
    S2, D2 = MH.assemble_matrices_broadcasting(gfs, elements, k)
    S3, D3 = assemble_matrices(gfs, smesh, k)
    @test S ≈ S2 ≈ S3
    @test D ≈ D2 ≈ D3

    Sv, Dv = assemble_matrices(gfs_v, smesh, k)
    @test size(Sv) == (2, 2)
    @test Sv ≈ S rtol=1e-2 atol=1e-3

    ϕ, _ = solve(D, S, bc; direct=true)
    @test length(ϕ) == 2
    @test all(isfinite, ϕ)

    ϕi, src = solve(D, S, bc; direct=false)
    @test length(ϕi) == 2
    @test src !== nothing

    fd = ForwardDiff.derivative(κ -> real(sum(assemble_matrices(gfs, smesh, κ)[1])), k)
    zy = Zygote.gradient(κ -> real(sum(assemble_matrices(gfs, smesh, κ)[1])), k)[1]
    @test fd ≈ zy rtol=1e-5 atol=1e-6

    enzyme_sum(κ) = real(sum(assemble_matrices(gfs, smesh, κ)[1]))
    ez = DifferentiationInterface.derivative(
        enzyme_sum, AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Reverse)), k)
    @test ez ≈ fd rtol=1e-5 atol=1e-6

    S_ka, D_ka = MH.assemble_matrices_ka(gfs, smesh, k)
    @test Array(S_ka) ≈ S rtol=1e-8 atol=1e-10
    @test Array(D_ka) ≈ D rtol=1e-8 atol=1e-10
end
