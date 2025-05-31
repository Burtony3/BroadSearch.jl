module CairoMakieExt

using BroadSearch, CairoMakie, LinearAlgebra, Dates, DataFrames, LambertsProblem

# =====================================================================
# === Plotting based on ephemerides

function CairoMakie.lines!(ax, ephem::AbstractEphemeris; tspan=nothing, kwargs...)

    # Getting states
    if isnothing(tspan) # Plot one period
        # Getting orbit period of planet
        T = period(ephem)
        t₀ = hasfield(typeof(ephem), :t₀) ? ephem.t₀ : 0
        epochs = LinRange(t₀, t₀ + T, 100)
    else # plot specified time span
        if typeof(tspan[1]) <: DateTime
            tspan = (
                BroadSearch.secondsPast(tspan[1]),
                BroadSearch.secondsPast(tspan[2])
            )
        end
        epochs = LinRange(tspan[1], tspan[2], 100)
    end
    R      = [ state(ephem, t)[1] for t in epochs ]
    x      = [r[1] for r in R]
    y      = [r[2] for r in R]

    CairoMakie.lines!(ax, x, y; kwargs...)

end


# =====================================================================
# === Plotting single broad search solutions

function CairoMakie.lines!(ax, sol::BroadSearchSolution; label=nothing, kwargs...)

    # Looping over pairs of nodes
    for (k, (n1, n2)) in enumerate(zip( sol.nodes, sol.nodes[2:end] ))
        # Running lambert
        # r1, v1, _, _ = lambert(n1.body, n2.body, n1.epoch, n2.epoch)
        prob = BallisticLambertsProblem(n1.body, n2.body, n1.epoch, n2.epoch)
        sol  = LambertsProblem.solve(prob, Izzo())
        lam  = BasicEphemeris(sol.r⃗₁, sol.v⃗₁, n1.epoch;
            μ=0., 
            parent=n1.body.parent, 
            name=n1.body.name*"-"*n2.body.name
        )

        # Plotting
        label = !isnothing(label) && k ≠ 1 ? nothing : label
        CairoMakie.lines!(ax, lam; tspan=(n1.epoch, n2.epoch), label=label, kwargs...)
    end

end


# =====================================================================
# === Plotting set of broad search solutions

function CairoMakie.lines!(ax, sols::Vector{BroadSearchSolution}; only=nothing, best=nothing, cmap=:seaborn_bright, kwargs...)

    # Setup
    colors = isa(cmap, Symbol) ? cgrad(cmap).colors : cmap
    k = 0
    d = Dict{String, typeof(colors[1])}()

    # Filtering results
    sols′ = BroadSearchSolution[]
    if ~isnothing(best) && ~isnothing(only)
        @assert isa(best, Function)

        # Getting best solution of a sequence by a metrics
        # lines!(ax, sort(sols, by = metric)[1:100]) # Plot top x
        # lines!(ax, filter(metric, sols)) # All solutions under x
        for sequence in only
            push!(sols′, sort(filter(s->s.sequence == sequence, sols), by=best)[1])
        end

        # Updating solutions
        sols = sols′
        only = nothing
    end

    # Plotting each solution
    for sol in sols
        # Getting color
        if !isnothing(only) && sol.sequence ∉ only; continue; end
        label = nothing
        if sol.sequence ∉ keys(d)
            k += 1
            if k > length(colors); k = 1; end
            d[sol.sequence] = colors[k]
            label = sol.sequence
        end

        # Plotting
        CairoMakie.lines!(ax, sol; color=d[sol.sequence], label=label, kwargs...)

    end
end


function CairoMakie.scatter!(ax, sols::Vector{BroadSearchSolution}; 
    xmetric::Function=s->s.epochs[1],
    ymetric::Function=s->s.cost,
    colorby::Function=s->s.sequence,
    cmap=:seaborn_bright
    )

    # Gathering data
    df = DataFrame(
        x = [xmetric(s) for s in sols],
        y = [ymetric(s) for s in sols],
        c = [colorby(s) for s in sols]
    )

    # Plotting
    if eltype(df.c) <: String
        colors = cgrad(cmap).colors
        color_idx = [mod(val, length(colors))+1 for val in groupindices(groupby(df, :c))]
        # println(length(colors))
        # println(maximum(color_idx))
        CairoMakie.scatter!(ax, df.x, df.y, color=colors[color_idx])
        cmap = cgrad(cmap, length(colors), categorical=true)
        ticks = unique(df.c)
        crange = nothing 
    else
        crange = (minimum(df.c), maximum(df.c))
        CairoMakie.scatter!(ax, df.x, df.y, color=df.c, colormap=cmap)
    end

    CairoMakie.Colorbar(ax.parent[1, 2], colormap=cmap, colorrange=crange)

    # return df
end

end