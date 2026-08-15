using Test
using MarineHydro
using PyCall
using ForwardDiff
using Zygote
using FiniteDifferences

const MH = MarineHydro

cpt = pyimport("capytaine")
cpt_mesh = cpt.mesh_sphere(name="sphere", radius=1.0, center=(0, 0, 0), resolution=(4, 4)).keep_immersed_part()
hull = Mesh(cpt_mesh)
# Slightly submerged lid so generate_lid stays finite (z=0 + lowest_lid_position can NaN).
cpt_lid = cpt_mesh.generate_lid(z=-0.05)
lid = Mesh(cpt_lid)

body_nolid = FloatingBody(hull, [:Heave], "nolid")
body_lid = FloatingBody(hull, lid, [:Heave], "withlid")
ω = 1.2

function heave_force(body, w)
    k = compute_wavenumber(w)
    prob = RadiationProblem(body, w, nothing, k, 0.0, :Heave, [:Heave])
    return real(solve_problem(prob; direct=true, gf="Wu").forces[:Heave])
end

@testset "solve_problem lid plumbing" begin
    combined = hull + lid
    @test combined.nfaces == hull.nfaces + lid.nfaces
    mask = eachrow(combined.faces) .∈ Ref(eachrow(hull.faces))
    @test count(mask) == hull.nfaces

    bc_hull = fill(1.0 + 2.0im, hull.nfaces)
    padded = MH._pad_hull_bc(bc_hull, mask)
    @test length(padded) == combined.nfaces
    @test padded[mask] == bc_hull
    @test all(iszero, padded[.!mask])
    @test MH._pad_hull_bc(bc_hull, trues(hull.nfaces)) == bc_hull
end

@testset "solve_problem primal with and without lid" begin
    F0 = heave_force(body_nolid, ω)
    F1 = heave_force(body_lid, ω)
    @test isfinite(F0)
    @test isfinite(F1)
end

@testset "solve_problem dF/dω no lid" begin
    fd = ForwardDiff.derivative(w -> heave_force(body_nolid, w), ω)
    zy = Zygote.gradient(w -> heave_force(body_nolid, w), ω)[1]
    cfd = FiniteDifferences.central_fdm(5, 1)(w -> heave_force(body_nolid, w), ω)
    @test zy ≈ fd rtol=1e-5 atol=1e-6
    @test zy ≈ cfd rtol=1e-4 atol=1e-4
end

@testset "solve_problem dF/dω with lid" begin
    fd = ForwardDiff.derivative(w -> heave_force(body_lid, w), ω)
    zy = Zygote.gradient(w -> heave_force(body_lid, w), ω)[1]
    cfd = FiniteDifferences.central_fdm(5, 1)(w -> heave_force(body_lid, w), ω)
    @test zy ≈ fd rtol=1e-5 atol=1e-6
    @test zy ≈ cfd rtol=1e-4 atol=1e-4
end

@testset "solve_all_problems Rankine reuse is differentiable in ω" begin
    function A_two_freq(w)
        parameters = (wave_frequencies=[w, w + 0.25], radiating_dofs=[:Heave], influenced_dofs=[:Heave])
        data = compute_hydrodynamic_coefficients(parameters, body_nolid)
        return sum(data.added_mass)
    end
    fd = ForwardDiff.derivative(A_two_freq, ω)
    zy = Zygote.gradient(A_two_freq, ω)[1]
    cfd = FiniteDifferences.central_fdm(5, 1)(A_two_freq, ω)
    @test zy ≈ fd rtol=1e-5 atol=1e-6
    @test zy ≈ cfd rtol=1e-4 atol=1e-4
end
