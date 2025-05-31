export BroadSearchSolution

export find_solutions,
       sequence

struct BroadSearchSolution
    sequence::String
    C₃::Float64
    v∞::SVector{3, Float64}
    epochs::Vector{DateTime}
    _epochs::Vector{Int}
    nodes::Vector{<:AbstractTreeNode}
    cost::Float64
    
    # TODO:
    # C3
    # Final vinf
    # TOF
end

function BroadSearchSolution(node::AbstractTreeNode)::BroadSearchSolution
    # Collecting all nodes
    nodes = AbstractTreeNode[]
    cost  = 0.0
    while ~isa(node, Root)
        # Adding to array and updating cost
        push!(nodes, node)
        cost += node.cost

        # Updating to parent
        node = node.parent
    end

    # Reversing list
    reverse!(nodes)

    # Collecting info from node list
    seq = join([node.body.name[1] for node in nodes])
    _epochs = [node.epoch for node in nodes]
    epochs = [J2000+Second(epoch) for epoch in _epochs]
    vinf = nodes[end].v∞in

    # C3
    # v1 = state(nodes[1].body, nodes[1].epoch)[2]
    # v2 = state(nodes[2].body, nodes[2].epoch)[2]
    # _, v1′, _, _ = lambert(nodes[1].body, nodes[2].body, nodes[1].epoch, nodes[2].epoch)
    prob = BallisticLambertsProblem(nodes[1].body, nodes[2].body, nodes[1].epoch, nodes[2].epoch)
    sol  = LambertsProblem.solve(prob, Izzo())

    # Finding cost
    v1 = state(nodes[1].body, nodes[1].epoch)[2]
    C₃ = norm(sol.v⃗₁ - v1)^2

    return BroadSearchSolution(seq, C₃, vinf, epochs, _epochs, nodes, cost)
end

function Base.show(io::IO, ::MIME"text/plain", sol::BroadSearchSolution)


    println(io, "Tree Search Solution:")
    println(io, "   Sequence:    $(sol.sequence)")
    println(io, "   Launch Date: $(sol.epochs[1])")
    println(io, "   Flight Time: $(round( (sol._epochs[end] - sol._epochs[1])/86400.0, digits=2)) days")
    println(io, "   Launch C₃:   $(round( sol.C₃, digits=2))")
    # println(io, "   Depth:    $(length(sol.nodes))")
    println(io, "   Cost:        $(round( sol.cost, digits=5))")
    # println(io,)
end
Base.show(io::IO, sol::BroadSearchSolution) = show(io, MIME"text/plain"(), sol)

# sequence(sol::BroadSearchSolution) = join([node.body.name[1] for node in sol.nodes])

# =====================================================================
# === Finding solutions in tree

function find_solutions(tree::Root)::Vector{BroadSearchSolution}
    solutions = BroadSearchSolution[]

    for launch in tree.children
        _find_solutions!(launch, tree, solutions)
    end

    return solutions
end

function _find_solutions!(node::AbstractTreeNode, tree::Root, solutions::Vector{BroadSearchSolution})::Nothing

    if atGoal(tree.goal, node)
        push!(solutions, BroadSearchSolution(node))
    elseif !isnothing(node.children) && length(node.children) > 0
        for child in node.children
            _find_solutions!(child, tree, solutions)
        end
    end

    nothing
end