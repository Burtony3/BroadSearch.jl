using Revise
using BroadSearch
using StaticArrays, LinearAlgebra
using BenchmarkTools
using Dates

# Setting up ephems
t0 = Dates.julian2datetime(2460806.5)
ssb = BasicEphemeris(1.3211e11; name="Solar System Barycenter")
earth = BasicEphemeris(
    SA[-9.739065086101781E+07, -1.168835899939733E+08, 3.215712293333560E+04], 
    SA[ 2.241415755885093E+01, -1.916764497189185E+01, 1.647478203533836E-03], 
    t0; 
    μ=398_600.435436, 
    parent=ssb, 
    name="Earth"
)
mars = BasicEphemeris(
    SA[-2.427376587813935E+08,  5.648226215753347E+07,  7.159669810561325E+06], 
    SA[-4.661523228103734E+00, -2.151121188496233E+01, -3.363192259658989E-01], 
    t0;
    μ=42_828.375214, 
    parent=ssb, 
    name="Mars"
)

# Propagaintg
r, v = state(earth, Dates.julian2datetime(2460807.5))

# Testing lambert
dt = -( Month(6) + Day(20) )
tof = Month(11)
r1, v1 = lambert(earth, mars, t0+dt, t0+tof+dt)
lam = BasicEphemeris(r1, v1, t0+dt;
    μ=0., 
    parent=ssb, 
    name="Transfer"
)
vinf = v1 - state(earth, t0+dt)[2]

# Plotting orbits
using CairoMakie
begin
    fig = Figure()
    ax  = Axis(fig[1, 1], aspect=DataAspect())

    # Using extension
    lines!(ax, earth, linestyle=:dash)
    lines!(ax, mars, linestyle=:dash)
    lines!(ax, lam, tspan=(t0+dt, t0+tof+dt), color=:black)
    fig
end