using Revise
using BroadSearch
using StaticArrays, LinearAlgebra
using BenchmarkTools
using Dates

# Setting up ephems
ssb = BasicEphemeris(1.3211e11; name="Solar System Barycenter")
earth = BasicEphemeris(
    SA[-9.739065086101781E+07, -1.168835899939733E+08, 3.215712293333560E+04], 
    SA[ 2.241415755885093E+01, -1.916764497189185E+01, 1.647478203533836E-03], 
    Dates.julian2datetime(2460806.5); 
    μ=398_600.435436, 
    parent=ssb, 
    name="Earth"
)
mars = BasicEphemeris(
    SA[-2.427376587813935E+08,  5.648226215753347E+07,  7.159669810561325E+06], 
    SA[-4.661523228103734E+00, -2.151121188496233E+01, -3.363192259658989E-01], 
    Dates.julian2datetime(2460806.5);
    μ=42_828.375214, 
    parent=ssb, 
    name="Mars"
)
jupiter = BasicEphemeris(
    SA[ 9.808456526822060E+06, 7.661066057682021E+08, -3.396683274403751E+06], 
    SA[-1.321463753771015E+01, 7.881160546847092E-01,  2.924111662178345E-01], 
    Dates.julian2datetime(2460806.5);
    μ=126_686_531.900, 
    parent=ssb, 
    name="Jupiter"
)

# Creating boundaries
N = 200
τ₁, τ₂ = period(earth), period(mars)
T̲, T̅   = 0.1*(τ₁ + τ₂), 2.0*(τ₁ + τ₂)
tof = LinRange(T̲, T̅, N)
t0  = DateTime(2024, 11, 02, 0, 0, 0)

# Running lambert and plotting boundary
using CairoMakie
begin
    fig = Figure()
    ax  = Axis(fig[1, 1], aspect=DataAspect())

    # Plotting planet arcs
    vinfs = Float64[]
    crange = (5.0, 15.0)
    val = Month(0)
    lines!(ax, earth, linestyle=:dash)
    lines!(ax, jupiter, linestyle=:dash)

    for dt in tof
        dt = Int(round(dt))
        # Calculating lambert and making trajectory
        r1, v1 = lambert(earth, jupiter, t0+val, t0+val+Second(dt))
        vinf = v1 - state(earth, t0+val)[2] |> norm
        lam  = BasicEphemeris(r1, v1, t0+val;
            μ=0., 
            parent=ssb, 
            name="Transfer"
        )
        # append!(lams, lam)
        append!(vinfs, vinf)

        # Calculating vinf

        # Plotting
        # if vinf < 5.0
        lines!(ax, lam, tspan=(t0+val, t0+val+Second(dt)), color=vinf*ones(100), colorrange=crange)
        # end
    end

    Colorbar(fig[1, 2], colorrange=crange, label="Launch vinf")

    fig
end