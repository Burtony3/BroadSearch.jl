# =====================================================================
# === Exports

export BasicFlyby

export evaluate!


# =====================================================================
# === Types

struct BasicFlyby <: AbstractSearchCost
    Δvmax::Float64
    C₃max::Float64
end


# =====================================================================
# === Cost FUnctions

function evaluate!(cost::AbstractSearchCost, node::LambertNode, parent::LaunchNode)::Nothing

    # TODO: Check if node is new?

    # Calculating lambert transfer
    v1 = state(parent.body, parent.epoch)[2]
    v2 = state(node.body, node.epoch)[2]
    _, v1′, _, v2′ = lambert(parent.body, node.body, parent.epoch, node.epoch)

    # Finding cost
    C₃ = norm(v1′ - v1)^2

    # Updating node
    node.cost = max(C₃ - cost.C₃max, 0.0) |> sqrt
    node.v∞in = v2′ - v2
    node.status = node.cost < cost.Δvmax ? :valid : :invalid

    nothing
end


function evaluate!(cost::BasicFlyby, node::LambertNode, parent::LambertNode)::Nothing

    # TODO: Check if node is new?

    # Calculating lambert transfer
    v1 = state(parent.body, parent.epoch)[2]
    v2 = state(node.body, node.epoch)[2]
    _, v1′, _, v2′ = lambert(parent.body, node.body, parent.epoch, node.epoch)

    # Calculating flyby
    v∞₊, v∞₋ = parent.v∞in, v1′ - v1
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
    node.v∞in = v2′ - v2
    node.status = node.cost < cost.Δvmax && δ < δmax ? :valid : :invalid

    nothing
end


# =====================================================================
# === Helper FUnctions

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