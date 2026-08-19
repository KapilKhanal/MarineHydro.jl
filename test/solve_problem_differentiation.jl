using Test
using MarineHydro
using PyCall
using ForwardDiff
using FiniteDifferences
using StaticArrays

const MH = MarineHydro

const HAS_ENZYME = try
    using Enzyme
    true
catch
    false
end

cpt = pyimport("capytaine")
cpt_mesh = cpt.mesh_sphere(name="sphere", radius=1.0, center=(0, 0, 0), resolution=(4, 4)).keep_immersed_part()
hull = Mesh(cpt_mesh)
# Slightly submerged lid so generate_lid stays finite (z=0 + lowest_lid_position can NaN).
cpt_lid = cpt_mesh.generate_lid(z=-0.05)
lid = Mesh(cpt_lid)

body_nolid = FloatingBody(hull, [:Heave], "nolid")
body_lid = FloatingBody(hull, lid, [:Heave], "withlid")
smesh = StaticArraysMesh(cpt_mesh)
body_smesh = FloatingBody(smesh, [:Heave], "smesh")
ω = 1.2

function heave_force(body, w)
    return real(solve_problem(RadiationProblem(body, w); direct=true, gf="Wu").forces[:Heave])
end

# Uniform scale of an existing smesh (not through Capytaine). Dual `r` promotes
# vertices/areas; Enzyme traces the same reconstruction.
function smesh_from_radius(r, smesh0)
    s = r
    z = zero(s)
    return StaticArraysMesh(
        [s * v for v in smesh0.vertices],
        smesh0.faces,
        [s * c for c in smesh0.centers],
        [n .+ z for n in smesh0.normals],
        [a * (s * s) for a in smesh0.areas],
        [ρ * s for ρ in smesh0.radii],
        smesh0.nvertices,
        smesh0.nfaces,
    )
end

function heave_force_from_radius(r, smesh0, w)
    body = FloatingBody(smesh_from_radius(r, smesh0), [:Heave], "scaled")
    return heave_force(body, w)
end

# Enzyme: active argument first; do not capture the body in a closure.
enzyme_heave_of_omega(w, body) = heave_force(body, w)

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
    res = solve_problem(RadiationProblem(body_smesh, ω))
    coef = radiation_coefficients(res)
    @test coef.added_mass.Heave ≈ real(res.forces[:Heave]) / ω^2
    @test coef.radiation_damping.Heave ≈ imag(res.forces[:Heave]) / ω
end

@testset "solve_problem dF/dω no lid" begin
    fd = ForwardDiff.derivative(w -> heave_force(body_nolid, w), ω)
    cfd = FiniteDifferences.central_fdm(5, 1)(w -> heave_force(body_nolid, w), ω)
    @test fd ≈ cfd rtol=1e-4 atol=1e-4
end

@testset "solve_problem dF/dω with lid" begin
    fd = ForwardDiff.derivative(w -> heave_force(body_lid, w), ω)
    cfd = FiniteDifferences.central_fdm(5, 1)(w -> heave_force(body_lid, w), ω)
    @test fd ≈ cfd rtol=1e-4 atol=1e-4
end

@testset "solve_all_problems Rankine reuse is differentiable in ω" begin
    function A_two_freq(w)
        parameters = (wave_frequencies=[w, w + 0.25], radiating_dofs=[:Heave], influenced_dofs=[:Heave])
        data = compute_hydrodynamic_coefficients(parameters, body_nolid)
        return sum(data.added_mass)
    end
    fd = ForwardDiff.derivative(A_two_freq, ω)
    cfd = FiniteDifferences.central_fdm(5, 1)(A_two_freq, ω)
    @test fd ≈ cfd rtol=1e-4 atol=1e-4
end

@testset "smesh FloatingBody through solve_problem" begin
    @test body_smesh isa FloatingBody{<:StaticArraysMesh}
    @test body_nolid isa FloatingBody{<:Mesh}
    F_dense = heave_force(body_nolid, ω)
    F_smesh = heave_force(body_smesh, ω)
    @test F_smesh ≈ F_dense rtol=1e-8 atol=1e-10

    fd_ω = ForwardDiff.derivative(w -> heave_force(body_smesh, w), ω)
    cfd_ω = FiniteDifferences.central_fdm(5, 1)(w -> heave_force(body_smesh, w), ω)
    @test fd_ω ≈ cfd_ω rtol=1e-4 atol=1e-4

    fd_r = ForwardDiff.derivative(r -> heave_force_from_radius(r, smesh, ω), 1.0)
    @test isfinite(fd_r)
    @test abs(fd_r) > 0

    function A_two_freq_smesh(w)
        parameters = (wave_frequencies=[w, w + 0.25], radiating_dofs=[:Heave], influenced_dofs=[:Heave])
        data = compute_hydrodynamic_coefficients(parameters, body_smesh)
        return sum(data.added_mass)
    end
    fd_all = ForwardDiff.derivative(A_two_freq_smesh, ω)
    cfd_all = FiniteDifferences.central_fdm(5, 1)(A_two_freq_smesh, ω)
    @test fd_all ≈ cfd_all rtol=1e-4 atol=1e-4

    if HAS_ENZYME
        mode = Enzyme.set_runtime_activity(Enzyme.Reverse)
        ez_ω = first(Enzyme.gradient(mode, enzyme_heave_of_omega, ω, Enzyme.Const(body_smesh)))
        @test ez_ω ≈ fd_ω rtol=1e-5 atol=1e-6

        ez_r = first(Enzyme.gradient(mode, heave_force_from_radius, 1.0,
            Enzyme.Const(smesh), Enzyme.Const(ω)))
        @test ez_r ≈ fd_r rtol=1e-5 atol=1e-6
    else
        @test_skip "Enzyme not installed"
    end
end
