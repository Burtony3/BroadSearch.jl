export BasicPolicy

export expand!,
       prune!

"""
    BasicPolicy(bodies::Vector{<:AbstractEphemeris}, nodes::Vector{DataType}, N::Int)

Define a basic expansion/pruning policy for search tree nodes over given `bodies`, using `nodes` types and sampling `N` time steps.
"""
struct BasicPolicy <: AbstractSearchPolicy
    bodies::Vector{<:AbstractEphemeris}
    nodes::Vector{DataType}
    N::Int
end


# =====================================================================
# === Constructing child array based on policy

"""
    expand!(parent::AbstractTreeNode, policy::BasicPolicy)

Generate children of `parent` according to `policy.bodies` and `policy.nodes`, sampling transfer times over `policy.N` steps.
"""
function expand!(parent::AbstractTreeNode, policy::BasicPolicy)::Nothing

    # TODO: Check if node has children already?
    
    children = AbstractTreeNode[]
    for body in policy.bodies
        # Building Tofs
        T̲, T̅ = tof_range(parent.body, body)
        tofs = [Int(round(tof)) for tof in LinRange(T̲, T̅, policy.N)]

        # Looping over values
        for tof in tofs
            # Looping over node types
            for type in policy.nodes
                # Appending to array
                push!(children, new_node(type, parent, body, tof))
            end
        end
    end

    # Setting children
    parent.children = children

    nothing
end


# =====================================================================
# === Building new children based on policy

# TODO: Need to handle hyperbola case later (no body field)

"""
    new_node(::Type{LambertNode}, parent::AbstractTreeNode, body::AbstractEphemeris, tof::Int) -> LambertNode

Create a new `LambertNode` child for `parent`, targeting `body` with time‐of‐flight `tof`.
"""
function new_node(::Type{LambertNode}, parent::AbstractTreeNode, body::AbstractEphemeris, tof::Int)

    return LambertNode(body, parent.epoch+tof, parent, nothing, SA[-1.0, -1, -1], :new, -1.0, 0)
end


# =====================================================================
# === Building new children based on policy

"""
    prune!(parent::AbstractTreeNode, policy::BasicPolicy; recursive::Bool=false)

Remove any invalid children from `parent.children`. If `recursive=true`, apply pruning to all descendants.
"""
function prune!(parent::AbstractTreeNode, policy::BasicPolicy; recursive::Bool=false)::Nothing
    filter!(node -> node.status != :invalid, parent.children)

    nothing
end


# =====================================================================
# === Helper function

"""
    tof_range(body1::AbstractEphemeris, body2::AbstractEphemeris; mults::Vector{Float64}=[0.1, 2.0]) -> Tuple{Float64,Float64}

Compute lower and upper bounds of time‐of‐flight between `body1` and `body2` using their orbital periods and optional multipliers.
"""
function tof_range(body1::AbstractEphemeris, body2::AbstractEphemeris; mults::Vector{Float64}=[0.1, 2.0])

    τ₁, τ₂ = period(body1), period(body2)
    if τ₁ > 1e8 && τ₂ > 1e8 # If both bodies are far out, reduce the lower end of the ranges
        mults .*= [0.5, 1.0]
    end
    T̲, T̅   = mults*(τ₁ + τ₂)

    return T̲, T̅
end