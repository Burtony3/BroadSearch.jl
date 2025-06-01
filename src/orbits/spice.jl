# Not a module, uses Requires instead...

# =====================================================================
# === Creating base structure

export SpiceEphemeris

"""
    SpiceEphemeris(name::String, parent::Union{AbstractEphemeris,Nothing, String}; frame::String="ECLIPJ2000", bodname::Union{String,Nothing}=nothing)
    SpiceEphemeris(id::Int, parent::Union{AbstractEphemeris,Nothing, String}; frame::String="ECLIPJ2000", bodname::Union{String,Nothing}=nothing)

Create a SPICE‐based ephemeris by specifying either the body name or NAIF ID.  
Automatically retrieves the body’s radius and gravitational parameter from SPICE kernels, and links to an optional parent ephemeris.
"""
struct SpiceEphemeris{T} <: AbstractEphemeris
    μ::T
    R::T
    parent::Union{AbstractEphemeris, Nothing}
    id::Int
    name::String
    frame::String

    # Primary Constructor (string name)
    function SpiceEphemeris(name::String, parent::Union{AbstractEphemeris, Nothing, String}; frame::String="ECLIPJ2000", bodname::Union{String, Nothing}=nothing)

        R, μ = 0.0, 0.0
        try
            R  = SPICE.bodvrd(isnothing(bodname) ? name : bodname, "RADII") 
            if length(R) > 1; R = sum(R)/length(R); end
        catch
            @warn "Could not find radius for $name, setting to 0"
            R  = 0.0
        end
        try
            μ  = SPICE.bodvrd(isnothing(bodname) ? name : bodname, "GM")[1]
        catch
            @warn "Could not find GM for $name, setting to 0"
            μ  = 0.0
        end
        id = SPICE.bodn2c(name) 
    
        return new(μ, R, parent, id, name, frame)
    end

    # Secondary constructor (id number)
    function SpiceEphemeris(id::Int, parent::Union{AbstractEphemeris, Nothing, String}; frame::String="ECLIPJ2000", bodname::Union{String, Nothing}=nothing)

        name = SPICE.bodc2n(id) 
        R, μ = 0.0, 0.0
        try
            R  = SPICE.bodvrd(isnothing(bodname) ? name : bodname, "RADII") 
            if length(R) > 1; R = sum(R)/length(R); end
        catch
            @warn "Could not find radius for $name, setting to 0"
            R  = 0.0
        end
        try
            μ  = SPICE.bodvrd(isnothing(bodname) ? name : bodname, "GM")[1]
        catch
            @warn "Could not find GM for $name, setting to 0"
            μ  = 0.0
        end
    
        return new(μ, R, parent, id, name, frame)
    end
end


# Methods
"""
    state(ephem::SpiceEphemeris, t::Real; abcorr::String="None")

Propagate the SPICE‐based ephemeris to time `t` (seconds since J2000 or ephemeris epoch) under SPICE spacecraft‐body dynamics.  
Returns the position and velocity vectors in Cartesian coordinates.
"""
function state(ephem::SpiceEphemeris, t::Real; abcorr::String="None")

    # utc = Dates.format( J2000 + Second(Int(round(t))), dateformat"yyyy-mm-ddTHH:MM:SS")
    # utc = "$(year)-$(month)-$(day)T$(hour):$(minute):$(second)"
    # et = SPICE.utc2et(utc)
    et = t

    x = SPICE.spkezr(ephem.name, et, ephem.frame, abcorr, ephem.parent.name)[1]

    return SA[x[1:3]...], SA[x[4:6]...]
end


"""
    period(ephem::SpiceEphemeris) -> Float64

Compute the orbital period of a SPICE‐based ephemeris assuming Keplerian motion about its parent.  
Uses a sample state at a fixed time to derive semimajor axis and calculate period.
"""
function period(ephem::SpiceEphemeris)
    rvec, vvec = state(ephem, 500_000_000)
    r, v = norm(rvec), norm(vvec)
    ϵ    = 0.5*v^2 - ephem.parent.μ / r
    a    = -ephem.parent.μ/2/ϵ
    return 2π*sqrt( a^3 / ephem.parent.μ )
end

"""
    eccentricity(ephem::SpiceEphemeris) -> Float64

Calculate the orbital eccentricity of a SPICE‐based ephemeris from a sample propagated state.  
Computes specific energy and angular momentum to return a scalar eccentricity.
"""
function eccentricity(ephem::SpiceEphemeris)
    rvec, vvec = state(ephem, 500_000_000)
    r, v = norm(rvec), norm(vvec)
    ϵ    = 0.5*v^2 - ephem.parent.μ / r

    # Angular momentum
    h = ephem.r⃗₀ × ephem.v⃗₀ |> norm

    return sqrt( 1 + 2*ϵ*h^2 / ephem.parent.μ^2 )
end

# Extending base
import Base: ==, isequal
"""
    ==(e1::SpiceEphemeris, e2::SpiceEphemeris) -> Bool

Check equality of two SPICE ephemerides by comparing μ, radius, and NAIF ID.
"""
==(e1::SpiceEphemeris, e2::SpiceEphemeris) = e1.μ == e2.μ && e1.R == e2.R && e1.id == e2.id
"""
    isequal(e1::SpiceEphemeris, e2::SpiceEphemeris) -> Bool

Alias for `==` to support hashing of SPICE ephemerides.
"""
isequal(e1::SpiceEphemeris, e2::SpiceEphemeris) = e1 == e2



# =====================================================================
# === Creating meshed structure

# Other Options
# RapidSpiceEphemeris
# InterpolatedSpiceEphemeris

# export MeshedSpiceEphemeris

"""
    MeshedSpiceEphemeris(μ::T, R::T, parent::Union{AbstractEphemeris,Nothing}, id::Int, name::String, frame::String, _coefficients::SMatrix{6,4,T}, tspan::Tuple{T,T}) where T<:Real

Represent an interpolated SPICE ephemeris using polynomial coefficients over a specified time span.  
Provides fast Keplerian state evaluation via mesh‐based interpolation of classical elements.
"""
struct MeshedSpiceEphemeris{T} <: AbstractEphemeris
    μ::T
    R::T
    parent::Union{AbstractEphemeris, Nothing}
    id::Int
    name::String
    frame::String
    _coefficients::SMatrix{6, 4, T}
    tspan::Tuple{T, T}
end


# Methods
"""
    state(ephem::MeshedSpiceEphemeris, t::Real) -> Tuple{SVector{3,T},SVector{3,T}}

Interpolate Keplerian elements for the given time `t` using stored polynomial coefficients, then convert to Cartesian position and velocity.
"""
function state(ephem::MeshedSpiceEphemeris, t::Real)

    # TODO: Need to double check this is correct still

    # Normalizing t
    t̅ = (t - ephem.tspan[1]) / diff(ephem.tspan)
    t⃗ = SA[0., t̅, t̅^2, t̅^3]

    # Calculating interpolated keplerian elements
    a, e, i, Ω, ω, M = ephem._coefficients*t⃗
    E = _mean2eccAnomaly(M, e)
    ν = 2*atan( sqrt(1 + e)*sin(E/2), sqrt(1 - e)*sin(E/2) )

    # Finding intermediate quantities
    r = a*(1 - e*cos(E))
    v = sqrt(ephem.μ*a)/r

    # Finding perifocal vectors
    r⃗ₚ = r * SA[cos(ν), sin(ν), 0.0]
    v⃗ₚ = v * SA[-sin(E), sqrt(1 - e^2)*cos(E), 0.0]

    # Creating rotation from perifocal to base frame
    # R = RotZXZ(-Ω, -i, -ω)
    R = RotZXZ(Ω, i, ω)

    return R*r⃗ₚ, R*v⃗ₚ
end

# Helper functions
"""
    _mean2eccAnomaly(M::Real, e::Real; tol::Float64=1e-10, maxiters::Int=20) -> Real

Solve Kepler’s equation for eccentric anomaly `E` given mean anomaly `M` and eccentricity `e`, using Newton’s method.
"""
function _mean2eccAnomaly(M::Real, e::Real; tol::Float64=1e-10, maxiters::Int=20)::Real
    # Creating initial guess for E
    E = M < π ? M + e/2 : M - e/2

    # Iterating
    for k in 1:maxiters
        ΔE = -(E - e*sin(E) - M)/(1 - e*cos(E))
        E += ΔE
        if abs(ΔE) < tol
            break
        end
    end

    return E
end