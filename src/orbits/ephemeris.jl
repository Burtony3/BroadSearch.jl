export AbstractEphemeris,
       BasicEphemeris

export state, period

const J2000 = DateTime( 2000, 01, 01, 12, 0, 0 )

"""
    secondsPast(t::DateTime, ref::DateTime = J2000)

Compute the number of seconds elapsed from the reference epoch `ref` to time `t`.

Returns the elapsed time in seconds.
"""
secondsPast(t::DateTime, ref::DateTime=J2000) = Second( t - ref ).value

# =====================================================================
# === Type Declarations

"""
    AbstractEphemeris

Abstract supertype for ephemeris data providers.
"""
abstract type AbstractEphemeris end

"""
    BasicEphemeris(r⃗, v⃗, t₀::Int; μ=0.0, R=0.0, parent=nothing, name="", id=0)

    BasicEphemeris(r⃗, v⃗, t::DateTime; μ=0.0, R=0.0, parent=nothing, name="", id=0)

    BasicEphemeris(μ; name="")

Construct a Keplerian ephemeris with initial state `(r⃗, v⃗)` at epoch `t₀` (or `DateTime t`), or create a zero‐state ephemeris with just gravitational parameter `μ`. Additional metadata (`μ`, `R`, `parent`, `name`, `id`) may be supplied via keywords.
"""
struct BasicEphemeris{T} <: AbstractEphemeris
    r⃗₀::SVector{3, T}
    v⃗₀::SVector{3, T}
    t₀::Int
    μ::T
    R::T
    parent::Union{AbstractEphemeris, Nothing}
    id::Int
    name::String

    # Primary Constructor
    function BasicEphemeris(r::SVector{3, T}, v::SVector{3, T}, t::Int; μ::T=0.0, R::T=0.0, parent=nothing, name="", id=0) where T<:Real
        return new{T}(r, v, t, μ, R, parent, id, name)
    end

    # Datetime Constructor
    function BasicEphemeris(r::SVector{3, T}, v::SVector{3, T}, t::DateTime; kwargs...) where T<:Real
        return new{T}(r, v, secondsPast(t); kwargs...)
    end

    # Center Body Constructor
    function BasicEphemeris(μ::Real; name="")
        return new{T}(SA[0., 0., 0.], SA[0., 0., 0.], 0, μ, 0.0, nothing, 0, name)
    end
end

# Methods
"""
    state(ephem::BasicEphemeris, t::Real)

    state(ephem::AbstractEphemeris, t::DateTime)

Propagate the ephemeris `ephem` to time `t` (seconds since its epoch) under Keplerian motion.

Returns a tuple `(r⃗, v⃗)` for position and velocity at time `t`.
"""
function state(ephem::BasicEphemeris, t::Real)

    # Assuming elliptical for now
    Δt = t - ephem.t₀
    if eccentricity(ephem) < 1.0
        return LambertsProblem._propKepTE(ephem.r⃗₀, ephem.v⃗₀, Δt, ephem.parent.μ)
    else
        return LambertsProblem._propKepTH(ephem.r⃗₀, ephem.v⃗₀, Δt, ephem.parent.μ)
    end
end

state(ephem::AbstractEphemeris, t::DateTime) = state(ephem, secondsPast(t))

"""
    period(ephem::BasicEphemeris)

Compute the orbital period of `ephem` assuming an elliptical trajectory.

Returns the period in the same time unit as the ephemeris epoch (seconds).
"""
function period(ephem::BasicEphemeris)
    r, v = norm(ephem.r⃗₀), norm(ephem.v⃗₀)
    ϵ    = 0.5*v^2 - ephem.parent.μ / r
    a    = -ephem.parent.μ/2/ϵ
    return 2π*sqrt( a^3 / ephem.parent.μ )
end

"""
    eccentricity(ephem::BasicEphemeris)

Calculate the orbital eccentricity for `ephem` based on its specific energy and angular momentum.

Returns the scalar eccentricity value.
"""
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

"""
    ==(e1::BasicEphemeris, e2::BasicEphemeris) -> Bool

Compare two `BasicEphemeris` objects for equality.
"""
==(e1::BasicEphemeris, e2::BasicEphemeris) = e1.t₀ == e2.t₀ && e1.μ == e2.μ && e1.R == e2.R && e1.id == e2.id

"""
    isequal(e1::BasicEphemeris, e2::BasicEphemeris) -> Bool

Alias for `==` to satisfy hashing protocols.
"""
isequal(e1::BasicEphemeris, e2::BasicEphemeris) = e1 == e2