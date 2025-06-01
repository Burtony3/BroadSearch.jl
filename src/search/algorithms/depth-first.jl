export DepthFirst

export search!

"""
    DepthFirst(maxdepth::Int)

Define a depth‐first search strategy with a maximum search depth of `maxdepth`.
"""
struct DepthFirst <: AbstractTreeSearch
    maxdepth::Int
end

# =====================================================================
# === Creating search

function search!(bf::DepthFirst, tree::Root, ::Val{true})::Vector{BroadSearchSolution}

    # Ensuring solution is possible
    if isdefined(Main, :SpiceEphemeris)
        @assert ~any(isa(x, SpiceEphemeris) for x in tree.policy.bodies) "SPICE.jl based ephemerides do not work in multithreaded environments, use "
    end

    solutions = BroadSearchSolution[]

    # Initial expansion of launch nodes
    pbar = ProgressBar(tree.children)
    Threads.@threads for launch in pbar

        # Expanding launch date
        _expand_and_prune!(launch, tree.policy, tree.cost)

        # Searching sub-nodes
        for node in launch.children
            _search!(node, tree, bf.maxdepth, solutions)
        end

        # Updating multi-line
        if length(solutions) > 0
            best_sol = argmin(s->s.cost, solutions)
            set_multiline_postfix(
                pbar, 
                "Num Solutions: $(length(solutions))\nBest Solution:\n   Sequence: $(best_sol.sequence)\n   Cost:     $(best_sol.cost)"
            )
        else
            set_multiline_postfix(
                pbar, 
                "Num Solutions: $(length(solutions))"
            )
        end
    end

    return solutions
end

function search!(bf::DepthFirst, tree::Root, ::Val{false})::Vector{BroadSearchSolution}

    # TODO: Ensure the numbers are possible, ie goal body is in policy

    solutions = BroadSearchSolution[]

    # Initial expansion of launch nodes
    pbar = ProgressBar(tree.children)
    for launch in pbar

        # Expanding launch date
        _expand_and_prune!(launch, tree.policy, tree.cost)

        # Searching sub-nodes
        for node in launch.children
            _search!(node, tree, bf.maxdepth, solutions)
        end

        # Updating multi-line
        if length(solutions) > 0
            best_sol = argmin(s->s.cost, solutions)
            set_multiline_postfix(
                pbar, 
                "Num Solutions: $(length(solutions))\nBest Solution:\n   Sequence: $(best_sol.sequence)\n   Cost:     $(best_sol.cost)"
            )
        else
            set_multiline_postfix(
                pbar, 
                "Num Solutions: $(length(solutions))"
            )
        end
    end

    return solutions
end

"""
    search!(bf::DepthFirst, tree::Root; multithread::Bool=false) -> Vector{BroadSearchSolution}

Perform a depth‐first traversal of the search tree rooted at `tree`, returning all found `BroadSearchSolution`s.  
If `multithread=true`, launch nodes are expanded in parallel using threads.
"""
search!(bf::DepthFirst, tree::Root; multithread::Bool=false) = search!(bf, tree, Val(multithread))

# =====================================================================
# === Sub-methods

"""
    _search!(node::AbstractTreeNode, tree::Root, maxdepth::Int, solutions::Vector{BroadSearchSolution}) -> Nothing

Recursively explore `node` and its descendants up to `maxdepth`.  
If a node satisfies the goal, mark it done and append its solution to `solutions`; otherwise expand or prune as appropriate.
"""
function _search!(node::AbstractTreeNode, tree::Root, maxdepth::Int, solutions::Vector{BroadSearchSolution})::Nothing

    # Do not expand if at goal
    if atGoal(tree.goal, node)
        node.status = :done
        push!( solutions, BroadSearchSolution(node) )
        # println(sol)

    # If below max depth, go one layer deeper
    elseif depth(node) < maxdepth
        # Expanding node
        _expand_and_prune!(node, tree.policy, tree.cost)

        # Checking if any subnodes remain
        if length(node.children) ≠ 0
            # Searching sub-nodes
            for child in node.children
                _search!(child, tree, maxdepth, solutions)
            end

        else # If empty, set invalid
            node.status = :invalid
            prune!(node.parent, tree.policy)

        end

    # If at max depth, set as invalid if necessary
    else
        node.status = :invalid

    end

    nothing
end

"""
    _expand_and_prune!(parent::AbstractTreeNode, policy::AbstractSearchPolicy, cost::AbstractSearchCost) -> Nothing

Expand `parent` according to `policy` and evaluate each child with `cost`, then prune invalid children.
"""
function _expand_and_prune!(parent, policy, cost)
    # Expanding 
    expand!(parent, policy)

    # Evaluating
    for node in parent.children
        evaluate!(cost, node, parent)
    end

    # Pruning
    # nstart = length(launch.children)
    prune!(parent, policy)
    # println("Pruned $(nstart - length(launch.children)) nodes")
end