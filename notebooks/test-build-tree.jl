using Revise
using BroadSearch
using StaticArrays, LinearAlgebra
using BenchmarkTools
using Dates
using ProgressBars

# Initializing Basic Data
defaultdata = DefaultEphemeris.basic;
N = 70

# Galileo Traj Search
goal = ConstrainedArrival(defaultdata.jupiter, 7.0, 8*365*86400.0)
cost = BasicFlyby(3.0, 49.0)
policy = BasicPolicy([defaultdata.venus, defaultdata.earth, defaultdata.mars, defaultdata.jupiter], [LambertNode], N)
tspan = (DateTime(1989, 01, 01, 00, 00, 00), DateTime(1991, 01, 01, 00, 00))

# Cassini Traj Search
# goal = ConstrainedArrival(defaultdata.saturn, 50.0, 8*365*86400.0)
# cost = BasicFlyby(3.0, 49.0)
# policy = BasicPolicy([defaultdata.venus, defaultdata.earth, defaultdata.mars, defaultdata.jupiter, defaultdata.saturn], [LambertNode], N)
# tspan = (DateTime(1997, 01, 01, 00, 00, 00), DateTime(1999, 01, 01, 00, 00))

# Trident Traj Search
# goal = ConstrainedArrival(defaultdata.neptune, 50.0, 15*365*86400.0)
# cost = BasicFlyby(3.0, 49.0)
# policy = BasicPolicy([defaultdata.venus, defaultdata.earth, defaultdata.mars, defaultdata.jupiter, defaultdata.saturn, defaultdata.neptune], [LambertNode], N)
# tspan = (DateTime(2026, 01, 01, 00, 00, 00), DateTime(2028, 01, 01, 00, 00))

# Building Tree
tree = initialize_tree(goal, cost, policy, defaultdata.earth, tspan, N)

# Using depth-first search
bf = DepthFirst(5)
sols = search!(bf, tree, multithread=true);

# Getting solutions
# sols = find_solutions(tree)
# sort!(sols, by=s->s.cost)


# =====================================================================
# === Plotting

using CairoMakie

begin
    # Setting up figure
    fig = Figure()
    ax  = Axis(fig[1, 1], aspect=DataAspect())

    # Plotting considered bodies
    for body in policy.bodies
        lines!(ax, body, color=:black, linestyle=:dash)
    end

    # Sorting solutions
    # metric = s -> norm(s.v∞)
    # metric = s -> norm(s.C₃)
    # metric = s -> (s._epochs[end] - s._epochs[1])
    # metric = s -> s._epochs[1]
    # metric = s -> s.cost

    # Filtering metrics
    # metric = s -> norm(s.v∞) < 5.75
    metric = s -> norm(s.C₃) < 15

    
    # Plotting solutions
    # sequences = nothing
    # sequences = ["EVVEJ"]
    # lines!(ax, sols) # Plot everything
    # lines!(ax, sols, only=sequences) # Plot specific sequences
    # lines!(ax, sort(sols, by = metric)[1:100]) # Plot top x
    lines!(ax, filter(metric, sols)) # All solutions under x
    Legend(fig[1, 2], ax)
    fig
end

begin
    # Setting up figure
    fig = Figure()
    ax  = Axis(fig[1, 1])

    # Scatterplot
    # sc = scatter!(ax, sols)
    # scatter!(ax, sols, ymetric=s->(s.epochs[end] - s.epochs[1]), colorby=x->norm(x.C₃), cmap=:viridis)
    scatter!(ax, sols, xmetric=s->s.epochs[1], ymetric=s->(s.epochs[end] - s.epochs[1]))
    # df = scatter!(ax, sols, xmetric=x->norm(x.C₃), ymetric=s -> (s._epochs[end] - s._epochs[1]))

    fig
end

# =====================================================================
# === Manual way

# Creating function
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

# Looping over all launch nodes
for launch in tree.children
    _expand_and_prune!(launch, policy, cost)
end

# Expanding best node
node = tree.children[1].children[1]
for launch in tree.children
    node′ = argmin(node->node.cost, [node for node in launch.children if node.body==defaultdata.venus])
    if node′.cost < node.cost
        node = node′
    end
end
_expand_and_prune!(node, policy, cost)