# Method choices live on the algorithm, not on kwargs / strings.
# `solve(prob, DirectBEM())` is the SciML `solve(prob, Tsit5())` pattern.

abstract type BEMAlgorithm end

struct DirectBEM{G} <: BEMAlgorithm
    greens::G
end
DirectBEM() = DirectBEM((Rankine(), RankineReflected(), GFWu()))
DirectBEM(gf::AbstractString) = DirectBEM(default_greens_functions(String(gf)))

struct IndirectBEM{G} <: BEMAlgorithm
    greens::G
end
IndirectBEM() = IndirectBEM((Rankine(), RankineReflected(), GFWu()))
IndirectBEM(gf::AbstractString) = IndirectBEM(default_greens_functions(String(gf)))

is_direct(::DirectBEM) = true
is_direct(::IndirectBEM) = false

function default_greens_functions(gf::String)
    if gf == "Wu"
        wave_gf = GFWu()
    elseif gf == "ExactGuevelDelhommeau"
        wave_gf = ExactGuevelDelhommeau()
    else
        error("Unknown Green's function \"$gf\". Use \"Wu\" or \"ExactGuevelDelhommeau\".")
    end
    return (Rankine(), RankineReflected(), wave_gf)
end

function bem_algorithm(; alg=nothing, direct::Bool=true, gf::String="Wu", greens_functions=nothing)
    alg !== nothing && return alg
    gfs = greens_functions === nothing ? default_greens_functions(gf) : greens_functions
    return direct ? DirectBEM(gfs) : IndirectBEM(gfs)
end
