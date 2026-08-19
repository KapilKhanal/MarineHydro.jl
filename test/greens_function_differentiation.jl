
using Test
using Zygote
using ForwardDiff
using MarineHydro
using PyCall
using LinearAlgebra
using StaticArrays

const HAS_ENZYME = try
    using Enzyme
    true
catch
    false
end

@testset "Greens Function Differentiability Tests" begin
    # Define elements
    e1 = (center=[0.0, 0.0, -1.0],)
    e2 = (
        center=[1.0, 1.0, -2.0],
        vertices=[
            -0.5 -0.5 0.0;
             0.5 -0.5 0.0;
             0.5  0.5 0.0;
            -0.5  0.5 0.0
        ] .+ [1.0, 1.0, -2.0]',
        normal=[0.0, 0.0, 1.0],
        radius=sqrt(2)/2,
        area=1.0,
    )

    @testset "Rankine Differentiability" begin
        # Define functions to test with only 1 parameter
        function greens_center1(center)
            e1_new = (center=center,)
            return greens(Rankine(), e1_new, e2)
        end

        function greens_center2(center)
            e2_new = (
                center=center,
                vertices=[
                    -0.5 -0.5 0.0;
                     0.5 -0.5 0.0;
                     0.5  0.5 0.0;
                    -0.5  0.5 0.0
                ] .+ center',
                normal=[0.0, 0.0, 1.0],
                radius=sqrt(2)/2,
                area=1.0,
            )
            return greens(Rankine(), e1, e2_new)
        end

        # Test differentiability of `greens` with respect to `e1.center`
        gradient1 = Zygote.gradient(greens_center1, e1.center)
        @test typeof(gradient1) == Tuple{Vector{Float64}}  # Ensure valid gradient type
        @test !any(isnan, gradient1[1])                  # Ensure gradient is valid (not NaN)

        # Test differentiability of `greens` with respect to `e2.center`
        gradient2 = Zygote.gradient(greens_center2, e2.center)
        @test typeof(gradient2) == Tuple{Vector{Float64}}
        @test !any(isnan, gradient2[1])

        # Define functions for `gradient_greens`
        function gradient_greens_center1(center)
            e1_new = (center=center,)
            return gradient_greens(Rankine(), e1_new, e2, with_respect_to_first_variable=true)
        end

        function gradient_greens_center2(center)
            e2_new = (
                center=center,
                vertices=[
                    -0.5 -0.5 0.0;
                     0.5 -0.5 0.0;
                     0.5  0.5 0.0;
                    -0.5  0.5 0.0
                ] .+ center',
                normal=[0.0, 0.0, 1.0],
                radius=sqrt(2)/2,
                area=1.0,
            )
            return gradient_greens(Rankine(), e1, e2_new, with_respect_to_first_variable=false)
        end

        # Test differentiability of `gradient_greens` with respect to `e1.center`
        grad_grad1 = Zygote.jacobian(gradient_greens_center1, e1.center)
        @test typeof(grad_grad1) == Tuple{Matrix{Float64}}  # Ensure the Jacobian is a matrix
        @test !any(isnan, grad_grad1[1])                   # Ensure no NaN values in the Jacobian


        # Test differentiability of `gradient_greens` with respect to `e2.center`
        grad_grad2 = Zygote.jacobian(gradient_greens_center2, e2.center)
        @test typeof(grad_grad1) == Tuple{Matrix{Float64}}
        @test !any(isnan, grad_grad2[1])
    end

    @testset "Rankine Integral Differentiability" begin
        function panel_at(center)
            return (
                center=center,
                vertices=[
                    -0.5 -0.5 0.0;
                     0.5 -0.5 0.0;
                     0.5  0.5 0.0;
                    -0.5  0.5 0.0
                ] .+ center',
                normal=[0.0, 0.0, 1.0],
                radius=sqrt(2)/2,
                area=1.0,
            )
        end

        @testset "primal agreement $name" for (name, gf_old, gf_new) in (
            ("Rankine", DelhommeauRankine(), Rankine()),
            ("RankineReflected", DelhommeauRankineReflected(), RankineReflected()),
        )
            @test integral(gf_new, e1, e2) ≈ integral(gf_old, e1, e2) rtol=1e-4 atol=1e-6
            @test collect(integral_gradient(gf_new, e1, e2; with_respect_to_first_variable=true)) ≈
                  collect(integral_gradient(gf_old, e1, e2; with_respect_to_first_variable=true)) rtol=1e-4 atol=1e-6
            @test collect(integral_gradient(gf_new, e1, e2; with_respect_to_first_variable=false)) ≈
                  collect(integral_gradient(gf_old, e1, e2; with_respect_to_first_variable=false)) rtol=1e-4 atol=1e-6
        end

        @testset "$name first-order AD" for (name, gf) in (
            ("DelhommeauRankine", DelhommeauRankine()),
            ("Rankine", Rankine()),
            ("DelhommeauRankineReflected", DelhommeauRankineReflected()),
            ("RankineReflected", RankineReflected()),
        )
            integral_center1(center) = integral(gf, (center=center,), e2)
            integral_center2(center) = integral(gf, e1, panel_at(center))

            g_zygote1 = Zygote.gradient(integral_center1, e1.center)[1]
            g_analytic1 = collect(integral_gradient(gf, e1, e2; with_respect_to_first_variable=true))
            g_fd1 = ForwardDiff.gradient(integral_center1, e1.center)
            @test g_zygote1 ≈ g_analytic1 rtol=1e-6 atol=1e-8
            @test g_zygote1 ≈ g_fd1 rtol=1e-5 atol=1e-7
            @test !any(isnan, g_zygote1)

            g_zygote2 = Zygote.gradient(integral_center2, e2.center)[1]
            g_analytic2 = collect(integral_gradient(gf, e1, e2; with_respect_to_first_variable=false))
            g_fd2 = ForwardDiff.gradient(integral_center2, e2.center)
            @test g_zygote2 ≈ g_analytic2 rtol=1e-5 atol=1e-6
            @test g_zygote2 ≈ g_fd2 rtol=1e-5 atol=1e-6
            @test !any(isnan, g_zygote2)
        end

        @testset "DelhommeauRankine vs Rankine Zygote agreement" begin
            for (gf_old, gf_new) in ((DelhommeauRankine(), Rankine()), (DelhommeauRankineReflected(), RankineReflected()))
                g_old1 = Zygote.gradient(c -> integral(gf_old, (center=c,), e2), e1.center)[1]
                g_new1 = Zygote.gradient(c -> integral(gf_new, (center=c,), e2), e1.center)[1]
                @test g_old1 ≈ g_new1 rtol=1e-4 atol=1e-6

                g_old2 = Zygote.gradient(c -> integral(gf_old, e1, panel_at(c)), e2.center)[1]
                g_new2 = Zygote.gradient(c -> integral(gf_new, e1, panel_at(c)), e2.center)[1]
                @test g_old2 ≈ g_new2 rtol=1e-4 atol=1e-6
            end
        end

        @testset "$name second-order AD" for (name, gf) in (
            ("DelhommeauRankine", DelhommeauRankine()),
            ("Rankine", Rankine()),
        )
            grad_center1(center) = collect(integral_gradient(gf, (center=center,), e2; with_respect_to_first_variable=true))
            J_zygote = Zygote.jacobian(grad_center1, e1.center)[1]
            J_fd = ForwardDiff.jacobian(grad_center1, e1.center)
            H_fd = ForwardDiff.hessian(c -> integral(gf, (center=c,), e2), e1.center)
            @test J_zygote ≈ J_fd rtol=1e-5 atol=1e-7
            @test J_zygote ≈ H_fd rtol=1e-5 atol=1e-7
            @test !any(isnan, J_zygote)
        end

        @testset "DelhommeauRankine vs Rankine Hessian agreement" begin
            J_old = Zygote.jacobian(c -> collect(integral_gradient(DelhommeauRankine(), (center=c,), e2; with_respect_to_first_variable=true)), e1.center)[1]
            J_new = Zygote.jacobian(c -> collect(integral_gradient(Rankine(), (center=c,), e2; with_respect_to_first_variable=true)), e1.center)[1]
            @test J_old ≈ J_new rtol=1e-4 atol=1e-6
        end

        @testset "Rankine ChainRules use analytic Birk derivatives" begin
            Tmat, qgc, local_corners, _ = MarineHydro.birk_panel_geometry(e2.vertices)
            p = Tmat' * (SVector{3}(e1.center) - SVector{3}(qgc))
            x, y, z = p[1], p[2], p[3]
            v = MarineHydro.velocity_derivatives(x, y, z, local_corners)
            H = MarineHydro.velocity_hessian(x, y, z, local_corners)

            gφ = Zygote.gradient((a, b, c) -> MarineHydro.velocity_potential(a, b, c, local_corners), x, y, z)
            @test [gφ[1], gφ[2], gφ[3]] ≈ collect(v) atol=1e-12

            Jv = Zygote.jacobian(q -> MarineHydro.velocity_derivatives(q[1], q[2], q[3], local_corners), [x, y, z])[1]
            @test Jv ≈ collect(H) atol=1e-12
            @test Jv ≈ ForwardDiff.jacobian(q -> MarineHydro.velocity_derivatives(q[1], q[2], q[3], local_corners), [x, y, z]) rtol=1e-8
        end

        if HAS_ENZYME
            @testset "Enzyme matches Zygote on Rankine geometry" begin
                mode = Enzyme.set_runtime_activity(Enzyme.Reverse)
                function enzyme_integral_c1(center, e2c, gf)
                    return integral(gf, (center=center,), e2c)
                end
                function enzyme_integral_verts(verts, c1, c2, n, r, a, gf)
                    e2n = (center=c2, vertices=verts, normal=n, radius=r, area=a)
                    return integral(gf, (center=c1,), e2n)
                end
                for gf in (Rankine(), RankineReflected())
                    g_zy = Zygote.gradient(c -> integral(gf, (center=c,), e2), e1.center)[1]
                    g_ez = first(Enzyme.gradient(mode, enzyme_integral_c1, copy(e1.center),
                        Enzyme.Const(e2), Enzyme.Const(gf)))
                    @test all(isfinite, g_ez)
                    @test g_ez ≈ g_zy rtol=1e-8 atol=1e-10

                    g_zy_v = Zygote.gradient(v -> integral(gf, e1, (center=e2.center, vertices=v,
                        normal=e2.normal, radius=e2.radius, area=e2.area)), e2.vertices)[1]
                    g_ez_v = first(Enzyme.gradient(mode, enzyme_integral_verts, copy(e2.vertices),
                        Enzyme.Const(e1.center), Enzyme.Const(e2.center), Enzyme.Const(e2.normal),
                        Enzyme.Const(e2.radius), Enzyme.Const(e2.area), Enzyme.Const(gf)))
                    @test all(isfinite, g_ez_v)
                    @test g_ez_v ≈ g_zy_v rtol=1e-8 atol=1e-10
                end

                Tmat, qgc, local_corners, _ = MarineHydro.birk_panel_geometry(e2.vertices)
                p = Tmat' * (SVector{3}(e1.center) - SVector{3}(qgc))
                x, y, z = p[1], p[2], p[3]
                v = MarineHydro.velocity_derivatives(x, y, z, local_corners)
                function enzyme_φ(a, b, c, lc)
                    return MarineHydro.velocity_potential(a, b, c, lc)
                end
                gφ = Enzyme.gradient(mode, enzyme_φ, x, y, z, Enzyme.Const(local_corners))
                @test [gφ[1], gφ[2], gφ[3]] ≈ collect(v) atol=1e-12
            end
        end
    end

    @testset "GFWu Differentiability" for k in [0.1, 1.0, 10.0]
        # Define functions to test real and imaginary parts separately
        function real_greens_gfwu_center1(center)
            e1_new = (center=center,)
            return real(greens(GFWu(), e1_new, e2, k))
        end

        function imag_greens_gfwu_center1(center)
            e1_new = (center=center,)
            return imag(greens(GFWu(), e1_new, e2, k))
        end

        function real_gradient_greens_gfwu_center1(center)
            e1_new = (center=center,)
            return real(gradient_greens(GFWu(), e1_new, e2, k))
        end

        function imag_gradient_greens_gfwu_center1(center)
            e1_new = (center=center,)
            return imag(gradient_greens(GFWu(), e1_new, e2, k))
        end

        # Test differentiability of `greens` with respect to `e1.center`
        real_grad1 = Zygote.gradient(real_greens_gfwu_center1, e1.center)
        imag_grad1 = Zygote.gradient(imag_greens_gfwu_center1, e1.center)
        @test typeof(real_grad1) == Tuple{Vector{Float64}}
        @test typeof(imag_grad1) == Tuple{Vector{Float64}}
        @test !any(isnan, real_grad1[1])
        @test !any(isnan, imag_grad1[1])

        # Test differentiability of `gradient_greens` with respect to `e1.center`
        real_grad_grad1 = Zygote.jacobian(real_gradient_greens_gfwu_center1, e1.center)
        imag_grad_grad1 = Zygote.jacobian(imag_gradient_greens_gfwu_center1, e1.center)
        @test typeof(real_grad_grad1) == Tuple{Matrix{Float64}}
        @test typeof(imag_grad_grad1) == Tuple{Matrix{Float64}}
        @test !any(isnan, real_grad_grad1[1])
        @test !any(isnan, imag_grad_grad1[1])
    end

    @testset "ExactGuevelDelhommeau Differentiability" for k in [0.1, 1.0, 10.0]
        function real_greens_gfwu_center1(center)
            e1_new = (center=center,)
            return real(greens(ExactGuevelDelhommeau(), e1_new, e2, k))
        end

        function imag_greens_gfwu_center1(center)
            e1_new = (center=center,)
            return imag(greens(ExactGuevelDelhommeau(), e1_new, e2, k))
        end

        function real_gradient_greens_gfwu_center1(center)
            e1_new = (center=center,)
            return real(gradient_greens(ExactGuevelDelhommeau(), e1_new, e2, k))
        end

        function imag_gradient_greens_gfwu_center1(center)
            e1_new = (center=center,)
            return imag(gradient_greens(ExactGuevelDelhommeau(), e1_new, e2, k))
        end

        # Test differentiability of `greens` with respect to `e1.center`
        real_grad1 = Zygote.gradient(real_greens_gfwu_center1, e1.center)
        imag_grad1 = Zygote.gradient(imag_greens_gfwu_center1, e1.center)
        @test typeof(real_grad1) == Tuple{Vector{Float64}}
        @test typeof(imag_grad1) == Tuple{Vector{Float64}}
        @test !any(isnan, real_grad1[1])
        @test !any(isnan, imag_grad1[1])

        # Test differentiability of `gradient_greens` with respect to `e1.center`
        real_grad_grad1 = Zygote.jacobian(real_gradient_greens_gfwu_center1, e1.center)
        imag_grad_grad1 = Zygote.jacobian(imag_gradient_greens_gfwu_center1, e1.center)
        @test typeof(real_grad_grad1) == Tuple{Matrix{Float64}}
        @test typeof(imag_grad_grad1) == Tuple{Matrix{Float64}}
        @test !any(isnan, real_grad_grad1[1])
        @test !any(isnan, imag_grad_grad1[1])
    end
end
