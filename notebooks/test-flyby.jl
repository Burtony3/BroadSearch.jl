using Revise
using BroadSearch
using StaticArrays, LinearAlgebra
using BenchmarkTools
using Dates

# Setup
dt = +( Month(2) )*0
tLaunch = DateTime(1999, 08, 18, 0, 0, 0) + dt
tFlyby  = DateTime(2000, 12, 30, 0, 0, 0) + dt
tArrive = DateTime(2004, 07, 01, 0, 0, 0) + dt
seq     = [3, 5, 6]

# Getting lamberts
Ephems = DefaultEphemeris.basic;
r1, v1, r2, v2 = lambert(Ephems.earth, Ephems.jupiter, tLaunch, tFlyby)
_, v2′, r3, v3 = lambert(Ephems.jupiter, Ephems.saturn, tFlyby, tArrive)

# Calculating flyby values
function flyby(ephem, t, v, v′; tol = 5.0)

    vP       = state(ephem, t)[2]
    v∞₊, v∞₋ = v - vP, v′ - vP
    Δv∞      = abs( norm(v∞₊) - norm(v∞₋) )
    v̅∞       = 0.5*(norm(v∞₊) + norm(v∞₋))

    δ    = acosd( ( v∞₊⋅v∞₋ ) / ( v̅∞^2 ) )
    δmax = 2*asind( ephem.μ / (ephem.μ + ephem.R*v̅∞^2) )

    valid = Δv∞ < tol && δ < δmax

    return Δv∞, valid, δ
end

Δv∞, valid, δ = flyby(Ephems.jupiter, tFlyby, v2, v2′)

# Plotting
using CairoMakie
begin
    fig = Figure()
    ax  = Axis(fig[1, 1])

    lines!(ax, Ephems.earth, linestyle=:dash)
    lines!(ax, Ephems.jupiter, tspan=(tLaunch, tArrive), linestyle=:dash)
    lines!(ax, Ephems.saturn, tspan=(tLaunch, tArrive), linestyle=:dash)

    lams = [
        BasicEphemeris(r1, v1, tLaunch;
            parent=Ephems.ssb, 
            name="Transfer"
        ),
        BasicEphemeris(r2, v2′, tFlyby;
            parent=Ephems.ssb, 
            name="Transfer"
        )
    ];
    lines!(ax, lams[1], tspan=(tLaunch, tFlyby), color=:black)
    lines!(ax, lams[2], tspan=(tFlyby, tArrive), color=:black)

    fig
end