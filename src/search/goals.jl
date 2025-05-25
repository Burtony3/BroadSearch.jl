# =====================================================================
# === Exports

export FreeArrival,
       ConstrainedArrival


# =====================================================================
# === Free Arrival

struct FreeArrival <: AbstractSearchGoal
    body::AbstractEphemeris
    # TODO: t̅::Float64
end

atGoal(goal::FreeArrival, node::AbstractTreeNode) = goal.body == node.body


# =====================================================================
# === Constrained Vinf Arrival

struct ConstrainedArrival <: AbstractSearchGoal
    body::AbstractEphemeris
    v̅∞::Float64
    t̅::Float64
end

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

struct VinfEscape <: AbstractSearchGoal
    v∞::SVector{3, Float64}
    epoch::Int
    Δv∞::Float64  # Allowable vector norm difference
    Δepoch::Float64 # Allowable time different
end

# Is there a valid flyby at a current node that will allow this escape vector?