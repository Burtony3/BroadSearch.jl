

function _propKepTE(r⃗::T, v⃗::T, Δt::Real, μ::Float64)::Tuple{T, T} where T<:SVector

    # Defining Constants
    # @info "testing" x⃗₀ x⃗₀[1:3]
    # r⃗      = x⃗₀[1:3] |> SVector{3, eltype(x⃗₀)}; 
    r = norm(r⃗) 
    # v⃗      = x⃗₀[4:6] |> SVector{3, eltype(x⃗₀)};
    v = norm(v⃗)
    ε      = (0.5*v^2 - μ/r)
    a      = -0.5*μ/ε 
    ra, rμ = sqrt(a), sqrt(μ)
    σ      = dot(r⃗,v⃗)/rμ # Frequent Constant

    # Root Finding for Change in Eccentric Anomaly
    ΔEᵢ  = Δt*sqrt(μ/a^3)
    func = ΔE -> -ΔEᵢ + ΔE + (σ/ra)*(1 - cos(ΔE)) - (1 - r/a)*sin(ΔE)
    ΔE   = find_zero(func, ΔEᵢ)
    # ΔE = 0.3148

    # Creating Another Useful Constant
    cE, sE = cos(ΔE), sin(ΔE)
    ρ = a + (r - a)*cE + σ*(ra)*sE

    # Creating f, g, ḟ, ġ
    f = 1 - (a / r) * (1 - cE)
    g = a * σ / rμ * (1 - cE) + r * ra/rμ * sE
    ḟ = -rμ*ra / (ρ * r) * sE
    ġ = 1 - a / ρ * (1 - cE)

    # Converting to Cartesian States
    # r⃗ₖ = f*r⃗ + g*v⃗ # Most allocations
    # v⃗ₖ = ḟ*r⃗ + ġ*v⃗
    r⃗ₖ = @. f*r⃗ + g*v⃗ # Better
    v⃗ₖ = @. ḟ*r⃗ + ġ*v⃗
    
    return r⃗ₖ, v⃗ₖ
end

function _propKepTH(r⃗::T, v⃗::T, Δt::Real, μ::Float64)::Tuple{T, T} where T<:SVector

    # Defining Constants
    # r⃗      = x⃗₀[1:3] |> SVector{3, eltype(x⃗₀)}; 
    r = norm(r⃗) 
    # v⃗      = x⃗₀[4:6] |> SVector{3, eltype(x⃗₀)}; 
    v = norm(v⃗)
    ε      = (0.5*v^2 - μ/r)
    a      = -0.5*μ/ε 
    ra, rμ = sqrt(-a), sqrt(μ)
    σ      = dot(r⃗,v⃗)/rμ # Frequent Constant

    # Root Finding for Change in Hyperbolic Anomaly
    ΔHᵢ  = sign(Δt)
    func = ΔH -> -Δt*rμ/ra^3 - ΔH + (σ/ra)*(cosh(ΔH) - 1) + (1 - r/a)*sinh(ΔH)
    # ΔH   = find_zero(func, ΔHᵢ, Roots.Halley)
    ∇func = ΔH -> ForwardDiff.derivative(func, ΔH)
    δH = Inf
    ΔH = copy(ΔHᵢ)
    while abs(δH) > 1e-10
        δH = -func(ΔH)/∇func(ΔH)
        ΔH += δH
    end


    # Creating Another Useful Constant
    cH, sH = cosh(ΔH), sinh(ΔH)
    ρ = a + (r - a)*cH + σ*ra*sH

    # Creating f, g, ḟ, ġ
    f = 1 - (a / r) * (1 - cH)
    g = a * σ / rμ * (1 - cH) + r * ra/rμ * sH
    ḟ = -ra*rμ / (ρ * r) * sH
    ġ = 1 - a / ρ * (1 - cH)

    # Converting to Cartesian States
    r⃗ₖ = @. f*r⃗ + g*v⃗
    v⃗ₖ = @. ḟ*r⃗ + ġ*v⃗
    
    return r⃗ₖ, v⃗ₖ
end