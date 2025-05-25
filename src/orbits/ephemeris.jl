export AbstractEphemeris,
       BasicEphemeris

export state, period

const J2000 = DateTime( 2000, 01, 01, 12, 0, 0 )

secondsPast(t::DateTime, ref::DateTime=J2000) = Second( t - ref ).value

# =====================================================================
# === Type Declarations

abstract type AbstractEphemeris end

# === Basic Ephemeris
# Stucture type
struct BasicEphemeris{T} <: AbstractEphemeris
    r⃗₀::SVector{3, T}
    v⃗₀::SVector{3, T}
    t₀::Int
    μ::T
    R::T
    parent::Union{AbstractEphemeris, Nothing}
    id::Int
    name::String
end

# Constructors
function BasicEphemeris(r::SVector{3, T}, v::SVector{3, T}, t::Int; μ::T=0.0, R::T=0.0, parent=nothing, name="", id=0) where T<:Real
    return BasicEphemeris(r, v, t, μ, R, parent, id, name)
end

function BasicEphemeris(r::SVector{3, T}, v::SVector{3, T}, t::DateTime; kwargs...) where T<:Real

    return BasicEphemeris(r, v, secondsPast(t); kwargs...)
end

BasicEphemeris(μ::Real; name="") = BasicEphemeris(SA[0., 0., 0.], SA[0., 0., 0.], 0, μ, 0.0, nothing, 0, name)

# Methods
function state(ephem::BasicEphemeris, t::Real)

    # Assuming elliptical for now
    Δt = t - ephem.t₀
    if eccentricity(ephem) < 1.0
        return _propKepTE(ephem.r⃗₀, ephem.v⃗₀, Δt, ephem.parent.μ)
    else
        return _propKepTH(ephem.r⃗₀, ephem.v⃗₀, Δt, ephem.parent.μ)
    end
end

state(ephem::AbstractEphemeris, t::DateTime) = state(ephem, secondsPast(t))

function period(ephem::BasicEphemeris)
    r, v = norm(ephem.r⃗₀), norm(ephem.v⃗₀)
    ϵ    = 0.5*v^2 - ephem.parent.μ / r
    a    = -ephem.parent.μ/2/ϵ
    return 2π*sqrt( a^3 / ephem.parent.μ )
end

function eccentricity(ephem::BasicEphemeris)
    # Total energy
    r, v = norm(ephem.r⃗₀), norm(ephem.v⃗₀)
    ϵ    = 0.5*v^2 - ephem.parent.μ / r

    # Angular momentum
    h = ephem.r⃗₀ × ephem.v⃗₀ |> norm

    return sqrt( 1 + 2*ϵ*h^2 / ephem.parent.μ^2 )
end

# Extending base
import Base: ==, isequal
==(e1::BasicEphemeris, e2::BasicEphemeris) = e1.t₀ == e2.t₀ && e1.μ == e2.μ && e1.R == e2.R && e1.id == e2.id
isequal(e1::BasicEphemeris, e2::BasicEphemeris) = e1 == e2