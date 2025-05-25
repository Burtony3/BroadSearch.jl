"""
A  > 83
A- > 80
B+ > 77
B  > 75
B- > 73
C+ > 67
C  > 63
C- > 60
D+ > 55
D  > 50
"""


# =====================================================================
# === Testing new stuff

export lambert

function lambert(
    e1::AbstractEphemeris, e2::AbstractEphemeris, 
    t1::T, t2::T; kwargs...
    ) where T<:Real

    tof = t2-t1
    r1, v1 = state(e1, t1)
    r2, v2 = state(e2, t2)

    _, v1′, _, v2′ = _lam_izzo2(r1, r2, tof, e1.parent.μ; kwargs...)
    # _, v1′, _, v2′ = _lam_izzo(r1, r2, tof, e1.parent.μ; kwargs...)

    return r1, v1′, r2, v2′
end

function lambert(
    e1::AbstractEphemeris, e2::AbstractEphemeris, 
    t1::T, t2::T; kwargs...
    ) where T<:DateTime

    tof = secondsPast(t2, t1)
    r1, v1 = state(e1, t1)
    r2, v2 = state(e2, t2)

    _, v1′, _, v2′ = _lam_izzo2(r1, r2, tof, e1.parent.μ; kwargs...)
    # _, v1′, _, v2′ = _lam_izzo(r1, r2, tof, e1.parent.μ; kwargs...)

    return r1, v1′, r2, v2′

end

# =====================================================================
# === Exports

#=
export Lambert, 
       Izzo,
       lambert

abstract type Lambert end
struct Izzo <: Lambert end

Converts Vector -> SVector# Izzo
function lambert(
    lam::Lambert,
    x⃗1::Vector{<:Real}, 
    x⃗2::Vector{<:Real},
    tof::Real, 
    μ::Real;
    kwargs...
    )

    return lambert(lam, x⃗1|>SVector{6, eltype(x⃗1)}, x⃗2|>SVector{6, eltype(x⃗2)},
        tof, μ; kwargs...)
end
=#

# Handling State -> State
# # function lambert(
# #     lam::Lambert,
# #     sᵢ::State{D}, 
# #     sⱼ::State{D},
# #     μ::Real;
# #     kwargs...
# #     ) where D<:Epoch
# 
# #     # TODO: Compare Frames
# #     # assert sᵢ.frame ≠ sⱼ.frame "Frames are incompatible"
# 
# #     # Calculating Time of Flight
# #     tof = sⱼ.epoch - sᵢ.epoch
# #     if D <: TimePeriod
# #         tof = canonicalize(tof) |> compound2days # Converting to days
# #         tof *= 86400.                            # Converting to seconds
# #     end
# 
# #     # Calculating Lambert
# #     xᵢ, xⱼ = lambert(lam, sᵢ.state, sⱼ.state, tof, μ; kwargs...)
# 
# #     # Constructing Trajectory
# #     return DiscreteTrajectory([xᵢ, xⱼ], [sᵢ.epoch, sⱼ.epoch], sᵢ.frame)
# # end

# =====================================================================
# === Revamping

#=
export LambertProblem
export solve

struct LambertProblem
    xᵢ::SVector{6, <:Real}
    xⱼ::SVector{6, <:Real}
    tof::Real
    μ::Real
    revs::Int
    retrograde::Bool
end

@inline function LambertProblem(
    xᵢ::SVector{6, <:Real}, xⱼ::SVector{6, <:Real}, 
    tof::Real, μ::Real; 
    revs::Int=0, retrograde::Bool=false
    )

    return LambertProblem(xᵢ, xⱼ, tof, μ, revs, retrograde)

end

@inline function LambertProblem(
    xᵢ::Vector{T1}, xⱼ::Vector{T2}, 
    tof::Real, μ::Real; 
    revs::Int=0, retrograde::Bool=false
    ) where {T1<:Real, T2<:Real}

    return LambertProblem(SVector{6, T1}(xᵢ), SVector{6, T2}(xⱼ), tof, μ, revs, retrograde)

end

solve(prob::LambertProblem, ::Izzo; tol::Float64 = 1e-12, maxiters::Int=100) = _lam_izzo(
    prob.xᵢ, prob.xⱼ, prob.tof, prob.μ;
    retrograde=prob.retrograde, outside=false, revs=prob.revs, 
    iterMax=maxiters, tol=tol
)
=#

# =====================================================================
# === Specific Algorithms

#=
# Izzo
function lambert(::Izzo,
    x⃗1::SVector{6, <:Real}, 
    x⃗2::SVector{6, <:Real},
    tof::Real, 
    μ::Real;
    kwargs...
    )

    return _lam_izzo( x⃗1, x⃗2, tof, μ; kwargs...)
end

# No Lambert() specified
function lambert(
    x⃗1::SVector{6, <:Real}, 
    x⃗2::SVector{6, <:Real},
    tof::Real, 
    μ::Real;
    kwargs...
    )

    return lambert( Izzo(), x⃗1, x⃗2, tof, μ; kwargs... )
end
=#

# =====================================================================
# === Non-Exports

function _lam_izzo(
    r⃗1::SVector{3, <:Real}, 
    # v⃗1::SVector{3, <:Real}, 
    r⃗2::SVector{3, <:Real}, 
    # v⃗2::SVector{3, <:Real}, 
    tof::Real, 
    μ::Real;
    retrograde::Bool = false,
    outside::Bool = false,
    revs::Int = 0,
    iterMax::Int = 100,
    tol::Real = 1e-12)
    

    # Unwrapping Inputs
    # r⃗1 = x⃗1[1:3]
    # r⃗2 = x⃗2[1:3] 
    
    # Nondimensionalizing States
    LU = norm(r⃗1);   r⃗1  = r⃗1/LU
    VU = sqrt(μ/LU); r⃗2  = r⃗2/LU
    TU = LU/VU;      tof = tof/TU; 
    
    # Nondimensional Geometric Parameters
    r2 = norm(r⃗2)
    D, C⃗ = r⃗1⋅r⃗2, r⃗1×r⃗2
    C  = norm(C⃗) 
    δν = atan(C, D)
    if ~retrograde
        δν = C⃗[3] < 0 ? 2π - δν : δν
    else
        δν = C⃗[3] > 0 ? 2π - δν : δν
    end
    
    # Handling Branching
    longway = δν > π ? -1 : 1
    revs = abs(revs);   tof = abs(tof)
    code = 0
    
    # Handy Parameters
    c  = sqrt(1 + r2^2 - 2*r2*cos(δν)); # non-dimensional chord
    s  = (1 + r2 + c)/2;                # non-dimensional semi-perimeter
    aₘ = s/2;                           # minimum energy ellipse semi major axis()
    Λ  = sqrt(r2)*cos(δν/2)/s;          # Battin Λ parameter
    Ĉ  = C⃗/C;                  
    
    # Initializing root-finding values
    logt = log(tof);              # Commonly used value
    if revs == 0
        Xᵢ = [-0.5233, 0.5233]
        X  = log.(1 .+ Xᵢ)
    else
        # Xᵢ = outside ? [-0.5234, -0.2234] : [0.7234, 0.5234]
        Xᵢ = outside ? [-0.2234, -0.5234] : [0.7234, 0.5234]
        X  = tan.(Xᵢ*π/2)
    end
    A = @. aₘ/(1 - Xᵢ^2)
    β⃗ = @. longway * 2*asin(sqrt((s-c)/2/A))
    α⃗ = @. 2*acos(  max(-1.0, min(1.0, Xᵢ)) )
    Y = @. A*sqrt(A)*((α⃗ - sin(α⃗)) - (β⃗-sin(β⃗)) + 2π*revs)
    
    # Initial estimate for residual
    if revs == 0
        # Ensuring domain
        Y[Y .< 0.0] .= 0.0

        # Calculating
        Y[:] = @. log(Y) - logt
    else
        Y[:] = Y .- tof
    end
    
    # Iterating until converged
    err  = Inf  
    iter = 0
    tof′ = 0.0
    x′ = 0
    while err > tol
        iter += 1

        # Updating iteration variable
        x′ = (X[1]*Y[2] - Y[1]*X[2]) / (Y[2]-Y[1])

        # Handling multi-rev conditions
        x = revs == 0 ? exp(x′) - 1 : atan(x′)*2/π
        a = aₘ/(1 - x^2)

        # Handling high energy orbit conditions
        if x < 1 # ellipse
            β = longway * 2*asin(sqrt((s-c)/2/a))
            α = 2*acos( max(-1, min(1, x)) )
        else # hyperbola
            α = 2*acosh(x)
            β = longway * 2*asinh(sqrt((s-c)/(-2*a)))
        end

        # Evaluating for updated time of flight
        if a > 0
            tof′ =  a*sqrt(a)*((α - sin(α)) - (β-sin(β)) + 2π*revs)
        else
            tof′ = -a*sqrt(-a)*((sinh(α) - α) - (sinh(β) - β))
        end

        # Finding residual
        ynew = revs == 0 ? log(tof′) - logt : tof′ - tof

        # Saving new iterations
        X[:] = [X[2], x′]
        Y[:] = [Y[2], ynew]
        err = abs(X[2] - X[1])

        # Break condition
        if iter > iterMax
            code = 1
            break 
        end
    end
    
    # If the Newton-Raphson scheme failed; try to solve the problem.
    if code ≠ 0
        # NOTE: use the original; UN-normalized quantities
        # [V1, V2, extremal_distances, exitflag] = 
        #     lambert_LancasterBlanchard[r⃗1*LU, r⃗2*LU, longway*tof*TU, leftbranch*revs, μ]
        return NaN, NaN, NaN, NaN
    end
    
    # Finding transfer orbit SMA
    x = revs == 0 ? exp(x′) - 1 : atan(x′)*2/π
    a = aₘ/(1-x^2)
    
    # Calculating Psi (ψ)
    if x < 1 # Ellipse Case
        β = longway * 2*asin(sqrt((s-c)/2/a))
        α = 2*acos( max(-1, min(1, x)) )
        ψ  = (α-β)/2
        ηᵢ = 2*a*sin(ψ)^2/s
        η  = sqrt(ηᵢ)
    else     # Hyperbolic Case
        β = longway * 2*asinh(sqrt((c-s)/2/a))
        α = 2*acosh(x)
        ψ  = (α-β)/2
        ηᵢ = -2*a*sinh(ψ)^2/s
        η  = sqrt(ηᵢ)
    end
    
    # Useful unit vectors
    ĥ = longway * Ĉ
    r̂2 = r⃗2/r2

    # Radial and tangential velocities
    Vₜ = [sqrt(r2/aₘ/ηᵢ * sin(δν/2)^2)];      push!(Vₜ, Vₜ[1]/r2)
    Vᵣ = [1/η/sqrt(aₘ) * (2*Λ*aₘ - Λ - x*η)]; push!(Vᵣ, (Vₜ[1] - Vₜ[2])/tan(δν/2) - Vᵣ[1])
    
    # End point velocities
    V1 = SA[(Vᵣ[1]*r⃗1 + Vₜ[1]*(ĥ×r⃗1))*VU...]
    V2 = SA[(Vᵣ[2]*r̂2 + Vₜ[2]*(ĥ×r̂2))*VU...]
    
    # OUTPUTTING
    return r⃗1, V1, r⃗2, V2
end

# file:///Users/biii/Downloads/MScThesis_BrunoCorreia_69807.pdf
function _lam_izzo2(
    r⃗1::SVector{3, <:Real}, 
    r⃗2::SVector{3, <:Real}, 
    tof::Real, 
    μ::Real;
    retrograde::Bool=false,
    longway::Bool=false,
    # inside::Bool=false,
    maxiters::Int=35, tol::Real=1e-8,
    revs::Int=0)

    # Non-dimensionalizing
    # LU = norm(r⃗1);   r⃗1  = r⃗1/LU;   r1 = 1.0 # By definition
    # VU = sqrt(μ/LU); r⃗2  = r⃗2/LU;   r2 = norm(r⃗2)
    # TU = LU/VU;      tof = tof/TU
    r1 = norm(r⃗1)
    r2 = norm(r⃗2)
    
    # Finding transfer angle
    h⃗  = r⃗1×r⃗2; ĥ = h⃗/norm(h⃗)   # Orbit plane normal vector
    # δν = atan(norm(h⃗), r⃗1⋅r⃗2) # Interior transfer angle

    # # Handling angle definition
    # cond = retrograde ? h⃗[3]>0 : h⃗[3]<0 # Updating condition based on motion
    # δν   = cond ? 2π - δν : δν          # Converting to outside angle on condition

    # Finding useful constants
    # c  = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(δν)); # non-dimensional chord (law of cosines)
    c  = norm(r⃗2 - r⃗1);                       # non-dimensional chord (law of cosines)
    s  = (r1 + r2 + c)/2;                     # non-dimensional semi-perimeter
    # aₘ = s/2;                                 # minimum energy ellipse semi major axis()      
    λ  = sqrt(1 - min(1.0, c/s))     
    λ  = h⃗[3]<0 ? -λ : λ

    # Non-dimensional time of flight
    T  = sqrt( 2*μ/s^3 )*tof

    # Doing some checks
    # https://github.com/poliastro/poliastro/blob/21fd7719e89a7d22b4eac63141a60a7f1b01768c/src/poliastro/core/iod.py#L284
    @assert revs <= floor(T/π) "Input revolutions is more than the max possible: $(floor(T/π))"

    # Getting initial condition and root solving
    x₀ = _izzo_ic(T, λ, revs, longway)
    δx = Inf
    x  = copy(x₀)
    for k in 1:maxiters
        δx = _izzo_update(x, λ, T, revs)
        x += δx
        if abs(δx) < tol
            # @info "Successfully converged in $k iterations" δx
            break
        elseif k == maxiters
            # @warn "Could not reach tolerance in $maxiters iterations" δx
            if abs(δx) > tol*10
                out = SA[NaN, NaN, NaN]
                return out, out, out, out
            end
        end
    end
    # @info "This took $k iterations"
    y  = _izzo_y(x, λ)

    # Parameters for extracting solution
    γ = sqrt(μ*s/2)
    ρ = (r1 - r2)/c
    σ = sqrt( 1 - ρ^2 )

    # Reconstructing velocities norms
    v₁ᵣ =  γ*( (λ*y - x) - ρ*(λ*y + x) )/r1
    v₂ᵣ = -γ*( (λ*y - x) + ρ*(λ*y + x) )/r2
    v₁ₜ =  γ*σ*( y + λ*x )/r1
    v₂ₜ =  γ*σ*( y + λ*x )/r2

    # Creating velocity vectors
    sign   = retrograde ? -1 : 1
    r̂₁ = r⃗1/r1 
    t̂₁ = sign*( h⃗[3]<0 ? r̂₁×ĥ : ĥ×r̂₁ )
    r̂₂ = r⃗2/r2
    t̂₂ = sign*( h⃗[3]<0 ? r̂₂×ĥ : ĥ×r̂₂ )

    return r⃗1, v₁ᵣ*r̂₁ + v₁ₜ*t̂₁, r⃗2, v₂ᵣ*r̂₂ + v₂ₜ*t̂₂
end


function _izzo_update(x::Real, λ::Real, T::Real, revs::Int)

    # Getting constants w.r.t. x
    y = _izzo_y(x, λ)
    ψ = _izzo_ψ(x, y, λ)

    # Getting function values and derivatives
    f  = _izzo_root(x, λ, T, revs, y, ψ)
    f′ = _izzo_root′(x, λ, T, revs, y, ψ, f)
    f″ = _izzo_root″(x, λ, T, revs, y, ψ, f, f′)
    f‴ = _izzo_root‴(x, λ, T, revs, y, ψ, f, f′, f″)
    
    # Householder Third-Order update
    # https://en.wikipedia.org/wiki/Householder%27s_method
    return -f * ( f′^2 - 0.5*f*f″ )/( f′^3 - f*f′*f″ + f‴*f^2 / 6 )

end


function _izzo_root(
    x::Real, λ::Real, T::Real, 
    revs::Int,
    y::Real, ψ::Real)::Real
    
    return ( 1/(1-x^2) )*( (ψ + revs*π)/sqrt(abs(1 - x^2)) - x + λ*y ) - T
end


function _izzo_root′(
    x::Real, λ::Real, T::Real, 
    revs::Int,
    y::Real, ψ::Real, f::Real)::Real
    
    return ( 1/(1-x^2) )*( 3*x*(f + T) - 2 - 2*λ^3*x/y )
end


function _izzo_root″(
    x::Real, λ::Real, T::Real, 
    revs::Int,
    y::Real, ψ::Real, f::Real, f′::Real)::Real
    
    return ( 1/(1-x^2) )*( 3*(f + T) + 5*x*f′ + 2*λ^3*(1 - λ^2)/y^3 )
end


function _izzo_root‴(
    x::Real, λ::Real, T::Real, 
    revs::Int,
    y::Real, ψ::Real, f::Real, f′::Real, f″::Real)::Real
    
    return ( 1/(1-x^2) )*( 7*x*f″ + 8*f′ - 6*λ^5*(1 - λ^2)*x/y^5 )
end


function _izzo_ψ(x::Real, y::Real, λ::Real)::Real

    # Elliptic case
    if -1 ≤ x < 1
        ψ = acos( max(-1.0, min(1.0, x*y + λ*(1-x^2))) )

    # Hyperbolic case
    elseif x > 1
        ψ = asinh( (y-x*λ)*sqrt(x^2 - 1) )

    # Parabolic case
    else
        ψ = 0.0

    end

    return ψ
end


function _izzo_y(x::Real, λ::Real)::Real
    return sqrt( 1 - λ^2*(1-x^2))
end


function _izzo_ic(T::Real, λ::Real, revs::Int, longway::Bool)::Real

    # Single rev case
    if revs < 1
        T₀ = acos(λ) + λ*sqrt( 1 - λ^2 )
        T₁ = 2 * (1 - λ^3) / 3
        if T ≥ T₀
            x₀ = (T₀/T)^(2/3) - 1
        elseif T < T₁
            x₀ = (5/2)*(T₁/T)*( (T₁ - T)/(1 - λ^5) ) + 1
        else # T₁ < T < T₀
            x₀ = (T₀/T)^log2(T₁/T) - 1
        end

    # Multi rev case
    else
        # Finding long and short way boundaries
        x̲₀ = ( (( revs*π + π )/( 8*T ))^(2/3) - 1 ) / ( (( revs*π + π )/( 8*T ))^(2/3) + 1 )
        x̅₀ = ( (( 8*T )/( revs*π ))^(2/3) - 1 )     / ( (( 8*T )/( revs*π ))^(2/3) + 1 )

        x₀ = !longway ? max(x̲₀, x̅₀) : min(x̲₀, x̅₀)
    end

    return x₀
end