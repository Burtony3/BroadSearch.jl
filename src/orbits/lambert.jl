# =====================================================================
# === Extending LambertsProblem Interface

import LambertsProblem

"""
    LambertsProblem.BallisticLambertsProblem(e1::AbstractEphemeris, e2::AbstractEphemeris, t1::T, t2::T; kwargs...) where T<:Real

Construct a `BallisticLambertsProblem` from two ephemerides at real‐valued times, using their Keplerian states to compute the time of flight and boundary position vectors.
"""
function LambertsProblem.BallisticLambertsProblem(
    e1::AbstractEphemeris, e2::AbstractEphemeris, 
    t1::T, t2::T; kwargs...) where T<:Real

    # Pulling states
    tof = t2-t1 |> Float64
    r1, v1 = state(e1, t1)
    r2, v2 = state(e2, t2)

    # Returning Problem
    return BallisticLambertsProblem(r1, r2, tof, e1.parent.μ; kwargs...)
end

"""
    LambertsProblem.BallisticLambertsProblem(e1::AbstractEphemeris, e2::AbstractEphemeris, t1::T, t2::T; kwargs...) where T<:DateTime

Construct a `BallisticLambertsProblem` from two ephemerides at `DateTime` times, converting to seconds past the epoch to compute time of flight and boundary position vectors.
"""
function LambertsProblem.BallisticLambertsProblem(
    e1::AbstractEphemeris, e2::AbstractEphemeris, 
    t1::T, t2::T; kwargs...) where T<:DateTime

    # Pulling states
    tof = secondsPast(t2, t1) |> Float64
    r1, v1 = state(e1, t1)
    r2, v2 = state(e2, t2)

    # Returning Problem
    return BallisticLambertsProblem(r1, r2, tof, e1.parent.μ; kwargs...)
end


# === Old interface
#=
function lambert(
    e1::AbstractEphemeris, e2::AbstractEphemeris, 
    t1::T, t2::T; kwargs...
    ) where T<:Real

    tof = t2-t1
    r1, v1 = state(e1, t1)
    r2, v2 = state(e2, t2)

    _, v1′, _, v2′ = _lam_izzo2(r1, r2, tof, e1.parent.μ; kwargs...)
    # _, v1′, _, v2′ = _lam_izzo(r1, r2, tof, e1.parent.μ; kwargs...)

    return r1, v1′, r2, v2′
end

function lambert(
    e1::AbstractEphemeris, e2::AbstractEphemeris, 
    t1::T, t2::T; kwargs...
    ) where T<:DateTime

    tof = secondsPast(t2, t1)
    r1, v1 = state(e1, t1)
    r2, v2 = state(e2, t2)

    _, v1′, _, v2′ = _lam_izzo2(r1, r2, tof, e1.parent.μ; kwargs...)
    # _, v1′, _, v2′ = _lam_izzo(r1, r2, tof, e1.parent.μ; kwargs...)

    return r1, v1′, r2, v2′

end
=#