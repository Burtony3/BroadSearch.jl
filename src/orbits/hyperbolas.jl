
function _hyperbolic_entry_escape(r⃗::SVector{3, <:Real}, v⃗∞::SVector{3, <:Real}, μ::Real; incoming::Bool=false, retrograde::Bool=false)

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

    # Getting minimum radius
    # If initial position is on opposite side of asymptote
    # from vinf, then periapsis, else target position radius
    rmin = io*sin(ν) < 0.0 ? a*(1-e) : r

    return v⃗, rmin
end