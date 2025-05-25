using Revise
# ENV["REQUIRES_DEBUG"] = "true";
# import Logging; Logging.global_logger(Logging.ConsoleLogger());
using BroadSearch
using SPICE, Downloads
using StaticArrays, LinearAlgebra
using BenchmarkTools
using Dates

# Retrieving and loading NAIF Kernels
kernels_files = [
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc",
]
kernels = [Downloads.download(f) for f in kernels_files]
for f in kernels; furnsh(f); end

# Testing
et = utc2et("2018-02-06T20:45:00")
spkpos("Earth", et, "ECLIPJ2000", "none", "Solar System Barycenter")

# Creating Ephemeris
ssb      = SpiceEphemeris(0, nothing; bodname="Sun")
venus    = SpiceEphemeris("Venus Barycenter", ssb; bodname="Venus")
earth    = SpiceEphemeris("Earth Barycenter", ssb; bodname="Earth")
mars     = SpiceEphemeris("Mars Barycenter", ssb; bodname="Mars")
jupiter  = SpiceEphemeris("Jupiter Barycenter", ssb; bodname="Jupiter")


# =====================================================================
# === Testing middleware

# TODO:
#   - something wrong with `outside` flag... maybe???

using BenchmarkTools
using Profile

# Getting Epochs
dt = [
    DateTime(1989, 10, 18, 00, 00, 00),
    DateTime(1990, 02, 10, 00, 00, 00),
    DateTime(1990, 12, 08, 00, 00, 00),
    DateTime(1992, 12, 30, 06, 00, 00),
    DateTime(1995, 12, 07, 00, 00, 00),
]

# Performing Lamberts (original method)
# r1,  v1,  r2, v2 = lambert(earth, venus, dt[1], dt[2])
# lam1 = BasicEphemeris(r1, v1, dt[1]; parent=ssb)
# r2′, v2′, r3, v3 = lambert(venus, earth, dt[2], dt[3])
# lam2 = BasicEphemeris(r2′, v2′, dt[2]; parent=ssb)
# r3′, v3′, r4, v4 = lambert(earth, earth, dt[3], dt[4], revs=1)
# lam3 = BasicEphemeris(r3′, v3′, dt[3]; parent=ssb)

# Performing Lamberts (new method)
r1,  v1,  r2, v2 = lambert(earth, venus, dt[1], dt[2])
lam1 = BasicEphemeris(r1, v1, dt[1]; parent=ssb)
r2′, v2′, r3, v3 = lambert(venus, earth, dt[2], dt[3])
lam2 = BasicEphemeris(r2′, v2′, dt[2]; parent=ssb)
r3′, v3′, r4, v4 = lambert(earth, earth, dt[3], dt[4], revs=1)
lam3 = BasicEphemeris(r3′, v3′, dt[3]; parent=ssb)
r4′, v4′, r5, v5 = lambert(earth, jupiter, dt[4], dt[5])
lam4 = BasicEphemeris(r4′, v4′, dt[4]; parent=ssb)

# Benchmarking difference
Profile.clear_malloc_data()
# @profview begin
    # for _ in 1:500
    #     lambert(earth, venus, dt[1], dt[2])
    # end
# end
@profview_allocs lambert(earth, venus, dt[1], dt[2]) sample_rate=0.1
# @profview_allocs lambert(earth, venus, dt[1], dt[2], izzo2=true) sample_rate=0.1
@benchmark r1,  v1,  r2, v2 = lambert(earth, venus, dt[1], dt[2])
# @benchmark r1,  v1,  r2, v2 = lambert(earth, venus, dt[1], dt[2], izzo2=true)
# @benchmark r3′, v3′, r4, v4 = lambert(earth, earth, dt[3], dt[4], revs=1)
# @benchmark r3′, v3′, r4, v4 = lambert(earth, earth, dt[3], dt[4], izzo2=true, revs=1)

# Finding vinf diff
rV, vV = state(venus, dt[2])
vinfin  = v2 - vV
vinfout = v2′ - vV
norm(vinfin) - norm(vinfout)

# Plotting
using CairoMakie
begin
    fig = Figure()
    ax  = Axis(fig[1, 1], aspect=DataAspect())

    lines!(ax, venus, color=:black, linestyle=:dash)
    lines!(ax, earth, color=:black, linestyle=:dash)
    lines!(ax, jupiter, color=:black, linestyle=:dash)

    lines!(ax, lam1; tspan=dt[1:2])
    lines!(ax, lam2; tspan=dt[2:3])
    lines!(ax, lam3; tspan=dt[3:4])
    lines!(ax, lam4; tspan=dt[4:5])

    fig
end


# =====================================================================
# === Tree Search

# Building Tree Search
N = 70
goal = ConstrainedArrival(jupiter, 7.0, 8*365*86400.0)
cost = BasicFlyby(3.0, 49.0)
policy = BasicPolicy([venus, earth, mars, jupiter], [LambertNode], N)
tspan = (DateTime(1989, 01, 01, 00, 00, 00), DateTime(1991, 01, 01, 00, 00))
tree = initialize_tree(goal, cost, policy, earth, tspan, N)

# Executing search
bf = DepthFirst(5)
sols = search!(bf, tree; multithread=false);


# =====================================================================
# === Plotting

begin
    # Setting up figure
    fig = Figure()
    ax  = Axis(fig[1, 1], aspect=DataAspect())

    # Plotting considered bodies
    for body in policy.bodies
        lines!(ax, body, color=:black, linestyle=:dash)
    end

    # Sorting solutions
    # metric = s -> norm(s.v∞)
    # metric = s -> norm(s.C₃)
    # metric = s -> (s._epochs[end] - s._epochs[1])
    # metric = s -> s._epochs[1]
    metric = s -> s.cost

    
    # Plotting solutions
    # sequences = nothing
    sequences = ["EVVEJ"]
    # lines!(ax, sols) # Plot everything
    lines!(ax, sols, only=sequences) # Plot specific sequences
    # lines!(ax, sort(sols, by = metric)[1:50]) # Plot top x
    Legend(fig[1, 2], ax)
    fig
end