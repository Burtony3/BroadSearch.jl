# =====================================================================
# === Exports

export BasicFlyby

export evaluate!


# =====================================================================
# === Types

"""
    BasicFlyby(Δvmax::Float64, C₃max::Float64)

Define a flyby cost evaluator using maximum allowable Δv∞ and C₃ limits.
"""
struct BasicFlyby <: AbstractSearchCost
    Δvmax::Float64
    C₃max::Float64
end


# =====================================================================
# === Cost FUnctions

"""
    evaluate!(cost::AbstractSearchCost, node::LambertNode, parent::LaunchNode)

Evaluate C₃‐based launch cost to a `LambertNode`, updating its cost, incoming v∞, and status.
"""
function evaluate!(cost::AbstractSearchCost, node::LambertNode, parent::LaunchNode)::Nothing

    # Creating problem interface
    prob = BallisticLambertsProblem(parent.body, node.body, parent.epoch, node.epoch)
    sol  = LambertsProblem.solve(prob, Izzo())

    # Finding cost
    v1 = state(parent.body, parent.epoch)[2]
    C₃ = norm(sol.v⃗₁ - v1)^2

    # Updating node
    node.cost = max(C₃ - cost.C₃max, 0.0) |> sqrt
    v2 = state(node.body, node.epoch)[2]
    node.v∞in = sol.v⃗₂ - v2
    node.status = node.cost < cost.Δvmax ? :valid : :invalid

    nothing
end

"""
    evaluate!(cost::BasicFlyby, node::LambertNode, parent::LambertNode)

Evaluate flyby maneuver cost to a `LambertNode`, updating its cost, incoming v∞, and status.
"""
function evaluate!(cost::BasicFlyby, node::LambertNode, parent::LambertNode)::Nothing

    # Creating problem interface
    prob = BallisticLambertsProblem(parent.body, node.body, parent.epoch, node.epoch)
    sol  = LambertsProblem.solve(prob, Izzo())

    # Calculating flyby
    v1 = state(parent.body, parent.epoch)[2]
    v∞₊, v∞₋ = parent.v∞in, sol.v⃗₁ - v1
    Δv∞      = abs( norm(v∞₊) - norm(v∞₋) )
    v̅∞       = 0.5*(norm(v∞₊) + norm(v∞₋))
    dotprod  =  v∞₊⋅v∞₋ / ( v̅∞^2 )
    if dotprod >= -1.0 && dotprod <= 1.0 # Ensuring domain is valid
        δ    = acosd( dotprod )
    else                                 # Catching numerical issues
        δ    = dotprod < 0.0 ? 180.0 : 0.0
    end
    δmax = 2*asind( parent.body.μ / (parent.body.μ + parent.body.R*v̅∞^2) )

    # Updating node
    node.cost = Δv∞
    v2 = state(node.body, node.epoch)[2]
    node.v∞in = sol.v⃗₂ - v2
    node.status = node.cost < cost.Δvmax && δ < δmax ? :valid : :invalid

    nothing
end


# =====================================================================
# === Helper FUnctions

"""
    flyby(ephem, t, v, v′)

Compute the Δv∞ and validity flag for a planetary flyby at time `t` given ephemeris `ephem` and velocities `v` (incoming) and `v′` (outgoing).
"""
function flyby(ephem, t, v, v′)
    vP       = state(ephem, t)[2]
    v∞₊, v∞₋ = v - vP, v′ - vP
    Δv∞      = abs( norm(v∞₊) - norm(v∞₋) )
    v̅∞       = 0.5*(norm(v∞₊) + norm(v∞₋))

    δ    = acosd( ( v∞₊⋅v∞₋ ) / ( v̅∞^2 ) )
    δmax = 2*asind( ephem.μ / (ephem.μ + ephem.R*v̅∞^2) )

    valid = Δv∞ < tol && δ < δmax

    return Δv∞, valid
end