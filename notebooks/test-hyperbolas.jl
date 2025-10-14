using Revise
using BroadSearch
using LambertsProblem
using StaticArrays, LinearAlgebra

# =====================================================================
# === Delete later

using ForwardDiff

func(r::SVector{3, <:Real}, μ::Real) = -μ/norm(r)^3 * r

μ = 1.0
r = SA[randn(3)...]
func(r, μ)

ForwardDiff.jacobian(x -> func(x, μ), r)

∇func(r::SVector{3, <:Real}, μ::Real) = begin
    r̂ = r/norm(r)

    2μ/norm(r)^4*r*r̂' - μ/norm(r)*(Matrix{Float64}(1.0I, 3, 3) - r̂*r̂')
end

∇func(r, μ)
ForwardDiff.jacobian(x -> func(x, μ), r)

# =====================================================================
# === Setting up function

function hyperbola(r⃗::SVector{3, T}, v⃗∞::SVector{3, T}, μ::T; incoming=false, retrograde=false) where T<:Real

    # Processing inputs
    r, v∞ = norm(r⃗), norm(v⃗∞) 
    io = incoming ? -1.0 : 1.0                 # Sign flip for incoming solutions
    h⃗ᵣ = SA[0.0, 0.0, retrograde ? -1.0 : 1.0] # Required ang momentum vec

    # Finding initial values
    h⃗ = r⃗×v⃗∞                            # Psuedo angular momentum
    q = io*sign(h⃗ᵣ⋅h⃗)*acos( r⃗⋅v⃗∞/r/v∞ ) # Angle between radius and vinf
    h⃗ = sign(h⃗ᵣ⋅h⃗)*h⃗                    # Final angular momentum direction

    # Orbit parameters
    a = -μ/v∞^2               # Semi-major axis
    u = -( r*sin(q) + sqrt(r^2*sin(q)^2 - 4a*(r - r*cos(q))) )/2/a
    e = sqrt(u^2 + 1)         # Eccentricity
    ν = io*( acos(-1/e) - q ) # True anomaly
    h = sqrt( μ*a*(1 - e^2) ) # Actual angular momentum scalar
    h⃗ = h⃗/norm(h⃗)*h           # Actual angular momentum vector

    # Returning velocity and FPA
    r̂ = r⃗/r
    θ̂ = (h⃗/h)×r̂
    v⃗ = μ/h*e*sin(ν)*r̂ + h/r*θ̂
    ϕ = acosd(h/r/norm(v⃗))

    # Getting minimum radius
    # If initial position is on opposite side of asymptote
    # from vinf, then periapsis, else target position radius
    rmin = io*sin(ν) < 0.0 ? a*(1-e) : r

    @info "Final values" a ϕ rmin

    return v⃗
end


# =====================================================================
# === Solving

# Target Vinf
v∞ = 4.0
ra, dec = 215.0, 0.0
v⃗∞ = v∞*SA[ cosd(ra)*cosd(dec), sind(ra)*cosd(dec), sind(dec) ]

# Creating ephemerides
center = BasicEphemeris(398_600.0; name="Earth Center")
R   = 35_786 + 6378.137
orb = BasicEphemeris(
    SA[ R, 0.0, 0.0 ],
    SA[ 0.0, sqrt(center.μ/R), 0.0 ],
    0; 
    parent=center, 
    name="Body",
)

# Solving for hyperbola
t = 2*3600.0 + 0*86400.0
v⃗ = hyperbola(state(orb, t)[1], v⃗∞, center.μ, incoming=false, retrograde=false)

# Creating ephemeris for new solution
escape = BasicEphemeris(
    state(orb, t)[1],
    v⃗,
    Int(t); 
    parent=center, 
    name="Escape",
)


# =====================================================================
# === Plotting solution

using CairoMakie

begin
    fig = Figure()
    ax  = Axis(fig[1, 1], aspect=DataAspect())

    # Plotting starting orbit
    lines!(ax, orb, color=:black, linestyle=:dash)

    # Plotting escape hyperbola
    tmax = t + 12*3600.0
    Rmax = norm( state(escape, tmax)[1] )
    lines!(ax, [0.0, Rmax]*cosd(ra)*cosd(dec), [0.0, Rmax]*sind(ra)*cosd(dec), color=:black)

    # Plotting escape trajectory
    lines!(ax, escape, tspan=(t, tmax))  # Outgoing
    # lines!(ax, escape, tspan=(-tmax, t)) # Incoming

    fig
end