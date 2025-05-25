# Not a module, uses Requires instead...

# =====================================================================
# === Creating base structure

export SpiceEphemeris

struct SpiceEphemeris{T} <: AbstractEphemeris
    μ::T
    R::T
    parent::Union{AbstractEphemeris, Nothing}
    id::Int
    name::String
    frame::String
end

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

    return SpiceEphemeris(μ, R, parent, id, name, frame)
end

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

    return SpiceEphemeris(μ, R, parent, id, name, frame)
end


# Methods
function state(ephem::SpiceEphemeris, t::Real; abcorr::String="None")

    # utc = Dates.format( J2000 + Second(Int(round(t))), dateformat"yyyy-mm-ddTHH:MM:SS")
    # utc = "$(year)-$(month)-$(day)T$(hour):$(minute):$(second)"
    # et = SPICE.utc2et(utc)
    et = t

    x = SPICE.spkezr(ephem.name, et, ephem.frame, abcorr, ephem.parent.name)[1]

    return SA[x[1:3]...], SA[x[4:6]...]
end

# state(ephem::AbstractEphemeris, t::DateTime) = state(ephem, secondsPast(t))


function period(ephem::SpiceEphemeris)
    rvec, vvec = state(ephem, 500_000_000)
    r, v = norm(rvec), norm(vvec)
    ϵ    = 0.5*v^2 - ephem.parent.μ / r
    a    = -ephem.parent.μ/2/ϵ
    return 2π*sqrt( a^3 / ephem.parent.μ )
end


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
==(e1::SpiceEphemeris, e2::SpiceEphemeris) = e1.μ == e2.μ && e1.R == e2.R && e1.id == e2.id
isequal(e1::SpiceEphemeris, e2::SpiceEphemeris) = e1 == e2



# =====================================================================
# === Creating meshed structure

# Other Options
# RapidSpiceEphemeris
# InterpolatedSpiceEphemeris

# export MeshedSpiceEphemeris

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