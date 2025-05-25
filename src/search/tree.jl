# =====================================================================
# === Exports

export Root,
       LaunchNode,
       HyperbolaNode,
       LambertNode

export initialize_tree

# =====================================================================
# === Setting up root of tree

# Root
mutable struct Root
    goal::AbstractSearchGoal
    cost::AbstractSearchCost
    policy::AbstractSearchPolicy
    children::Union{Vector{<:AbstractTreeNode}, Nothing}
    branches::Int
    status::Symbol
end

function Base.show(io::IO, ::MIME"text/plain", node::Root)
    println(io, "Search tree node root:")
    println(io, "   Children:    $(length(node.children))")
    println(io, "   Branches:    $(node.branches)")
    println(io, "   Goal Type:   $(string(typeof(node.goal)))")
    println(io, "   Cost Type:   $(string(typeof(node.cost)))")
    print(io,   "   Policy Type: $(string(typeof(node.policy)))")
    # println(io,)
end
Base.show(io::IO, node::Root) = show(io, MIME"text/plain"(), node)

function Base.show(io::IO, ::MIME"text/plain", node::AbstractTreeNode)
    if get(io, :compact, false)
        print(io, "$(string(typeof(node)))($(node.body.name), $(J2000 + Second(node.epoch)))")
    else
        isleaf = isnothing(node.children)
        println(io, "Search tree node $(string(typeof(node))):")
        println(io, "   Depth:    $(depth(node))")
        println(io, "   Epoch:    $(J2000 + Second(node.epoch))")
        println(io, "   Body:     $(node.body.name)")
        println(io, "   Cost:     $(node.cost)")
        println(io, "   Status:   $(node.status)")
        if !isleaf
            print(io, "   Children: $(length(node.children))")
        else
            print(io, "   Children: None (leaf)")
        end
    end
end
Base.show(io::IO, node::AbstractTreeNode) = show(io, MIME"text/plain"(), node)

# Valid Statuses
#   :valid    = Still explorable
#   :invalid  = Completely explored (nodes that exceed costs are removed?)
#   :complete = Goal achived

function root(node::AbstractTreeNode)::Root
    root = node.parent
    while !isa(root, Root)
        root = root.parent
    end
    return root
end

function depth(node::AbstractTreeNode)::Int
    k = 1
    root = node.parent
    while !isa(root, Root)
        root = root.parent
        k += 1
    end
    return k
end

# =====================================================================
# === Trajectory start nodes

mutable struct LaunchNode <: AbstractTreeNode
    body::AbstractEphemeris
    epoch::Int

    parent::Union{Root, Nothing}
    children::Union{Vector{<:AbstractTreeNode}, Nothing}

    status::Symbol
    cost::Float64
    visits::Int
end

mutable struct HyperbolaNode <: AbstractTreeNode
    v⃗∞::SVector{3, Float64}
    epoch::Int

    parent::Union{Root, Nothing}
    children::Union{Vector{<:AbstractTreeNode}, Nothing}

    status::Symbol
    cost::Float64
    visits::Int
end

# =====================================================================
# === Midpoint Nodes

mutable struct LambertNode <: AbstractTreeNode
    body::AbstractEphemeris
    epoch::Int
    
    parent::AbstractTreeNode
    children::Union{Vector{<:AbstractTreeNode}, Nothing}
    v∞in::SVector{3, Float64}

    status::Symbol
    cost::Float64
    visits::Int
end

# PlanetaryNode(body, epoch)

Base.copy(n::LambertNode) = LambertNode(n.body, n.epoch, n.parent, n.children, n.v∞in, n.status, n.cost, n.visits)

# =====================================================================
# === Initializing

function initialize_tree(
    goal::AbstractSearchGoal, cost::AbstractSearchCost, policy::AbstractSearchPolicy,
    ephem::AbstractEphemeris, tspan::Tuple{T, T}, 
    N::Int
    )::Root where T<:Union{DateTime, Int}

    # Creating children
    if isa(tspan[1], DateTime)
        tspan = (secondsPast(t) for t in tspan)
    end
    epochs = [Int(round(e)) for e in LinRange(tspan..., N)]
    children = [
        LaunchNode(ephem, epoch, nothing, nothing, :valid, 0.0, 0)
        for epoch in epochs
    ]

    # Initializing Root
    root = Root(goal, cost, policy, children, N, :valid)
    for child in root.children
        child.parent = root
    end

    return root
end

# function expand!(node::AbstractTreeNode)
#     root = root(node)
#     policy = root.policy
# end 