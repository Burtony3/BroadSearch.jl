# =====================================================================
# === Exports

export FreeArrival,
       ConstrainedArrival


# =====================================================================
# === Free Arrival

"""
    FreeArrival(body::AbstractEphemeris)

Goal for free arrival at a target body without constraints.
"""
struct FreeArrival <: AbstractSearchGoal
    body::AbstractEphemeris
    # TODO: t̅::Float64
end

"""
    atGoal(goal::FreeArrival, node::AbstractTreeNode) -> Bool

Check if `node` has reached the target body specified in `goal`.
"""
atGoal(goal::FreeArrival, node::AbstractTreeNode) = goal.body == node.body


# =====================================================================
# === Constrained Vinf Arrival

"""
    ConstrainedArrival(body::AbstractEphemeris, v̅∞::Float64, t̅::Float64)

Goal for arrival at a target body with constraints on excess velocity (`v̅∞`) and time of flight (`t̅`).
"""
struct ConstrainedArrival <: AbstractSearchGoal
    body::AbstractEphemeris
    v̅∞::Float64
    t̅::Float64
end

"""
    atGoal(goal::ConstrainedArrival, node::AbstractTreeNode) -> Bool

Check if `node` has reached the target body with excess velocity below `v̅∞` and flight time, in seconds, less than `t̅`.
"""
function atGoal(goal::ConstrainedArrival, node::AbstractTreeNode)

    if goal.body == node.body
        # Originally used deepcopy(), that ground speed to a halt
        parent = copy(node)
        # parent = n.parent
        while !any([isa(parent, LaunchNode), isa(parent, HyperbolaNode)])
            parent = parent.parent
        end

        tof = node.epoch - parent.epoch

        return norm(node.v∞in) < goal.v̅∞ && tof < goal.t̅

    else
        return false
    end
end


# =====================================================================
# === Vinf Exit Target

"""
    VinfEscape(v∞::SVector{3, Float64}, epoch::Int, Δv∞::Float64, Δepoch::Float64)

Goal for departing with a specified hyperbolic excess velocity at a given epoch within tolerances.
"""
struct VinfEscape <: AbstractSearchGoal
    v∞::SVector{3, Float64}
    epoch::Int
    Δv∞::Float64  # Allowable vector norm difference
    Δepoch::Float64 # Allowable time different
end

# Is there a valid flyby at a current node that will allow this escape vector?