using Test

# Skip GPU files only when they fail to load (missing/optional packages).
# Real assertion failures must still fail the suite.
function include_optional_tests(path)
    try
        include(path)
    catch e
        root = e
        while root isa LoadError
            root = root.error
        end
        if root isa ArgumentError || root isa InitError
            @testset "$(basename(path))" begin
                @test_skip "Could not load $(basename(path)): $root"
            end
        else
            rethrow()
        end
    end
end

@testset "MarineHydro.jl test suite" begin
    include("./consistency_with_Capytaine.jl")
    include("./greens_function.jl")
    include("./rankine_vectorized.jl")
    include("./greens_function_differentiation.jl")
    include("./matrix_assembly.jl")
    include("./matrix_assembly_differentiation.jl")
    include("./ad_rules_geometry.jl")
    include("./consistency_with_analytical_solutions.jl")
    include("./outputs_differentiation.jl")
    include("./solve_problem_differentiation.jl")
    include("./forward_speed_tests.jl")
    include_optional_tests("./gpu_cpu_smoke.jl")
    include_optional_tests("./gpu.jl")
end
