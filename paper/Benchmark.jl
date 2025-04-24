# The goal of this script is to evaluate the time to compute a single term of
# the interaction matrices, for each term of the Green function and for each
# kind of integral.

using BenchmarkTools
using MarineHydro
using StaticArrays

wavenumber = 1.0

green_functions = ["Rankine" => Rankine(), "Wu" => GFWu(), "ExactDelhommeau" => ExactGuevelDelhommeau()]
integrals = ["S" => MarineHydro.integral, "D" => MarineHydro.integral_gradient, "both" => MarineHydro.both_integral_and_integral_gradient]


element_1 = (center=[0.0, 0.0, -1.0],)
element_2 = (
    center=[1.0, 1.0, -2.0],
    vertices= [-0.5 -0.5 0.0; 0.5 -0.5 0.0; 0.5 0.5 0.0; -0.5 0.5 0.0] .+ [1.0, 1.0, -2.0]',
    normal=[0.0, 0.0, 1.0],
    radius=sqrt(2)/2,
    area=1.0,
)

suite = BenchmarkGroup()
for (gf_name, gf) in green_functions
    for (term_name, term) in integrals
      #  suite["StaticElement"][gf_name][term_name] = @benchmarkable ($term)($(gf), $static_element_1, $static_element_2, $wavenumber)
        suite["NamedTuple"][gf_name][term_name] = @benchmarkable ($term)($(gf), $element_1, $element_2, $wavenumber)
    end
end

tune!(suite)
results = run(suite)
BenchmarkTools.save("latest_results.json", results)


for (gf_name, gf_benchmarks) in results["NamedTuple"]
    for (term_name, term_benchmark) in gf_benchmarks
        median_time = median(term_benchmark)
        println("Green Function: $gf_name, Integral: $term_name, Median Time: $median_time ns")
    end
end

using DataFrames


benchmark_df = DataFrame(
    GreenFunction = String[],
    IntegralType = String[],
    MedianTime_ns = Float64[]
)


for (gf_name, gf_benchmarks) in results["NamedTuple"]
    for (term_name, term_benchmark) in gf_benchmarks
        median_time_ns = median(term_benchmark).time  
        push!(benchmark_df, (
            GreenFunction = gf_name,
            IntegralType = term_name,
            MedianTime_ns = median_time_ns
        ))
    end
end


println(benchmark_df)


using DataFrames

gradient_df = DataFrame(
    GreenFunction = String[],
    GradientType = String[],
    MedianGradientTime_ns = Float64[]
)

# Define the gradient functions


