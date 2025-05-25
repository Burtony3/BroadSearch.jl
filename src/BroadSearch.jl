module BroadSearch

# Package Imports
using LinearAlgebra, StaticArrays, Roots, Rotations
using ForwardDiff
using DataFrames
using Dates, Printf, ProgressBars, Requires

# Orbit mechanics scripts
include("orbits/kepler-prop.jl")
include("orbits/ephemeris.jl")
include("orbits/lambert.jl")
include("orbits/default-data.jl")

# Search setup scripts
include("search/abstracts.jl")
include("search/tree.jl")
include("search/costs.jl")
include("search/goals.jl")
include("search/policy.jl")
include("search/solutions.jl")

# Search algorithms scripts
include("search/algorithms/depth-first.jl")

# Needed for Requires.jl to work
function __init__()
    @require SPICE="5bab7191-041a-5c2e-a744-024b9c3a5062" begin
        @info "Loading SPICE.jl related extensions"
        include("orbits/spice.jl")
    end
end

end
