## --- Load packages
using Plots, Random, Distributions, StatGeochem

## --- Regular process: Brownian motion as a model for sediment accumulation/erosion

# Erosion/deposition once a year, bidirectional
function browniandistance(rng, N, T)
    @assert eachindex(T)==1:length(T)
    distance = zeros(N, length(T))
    rate = sqrt(pi)/sqrt(2)

    @inbounds for i in eachindex(T)
        for n in 1:N
            h = 0.
            for t in 1:T[i]
                r = rate*randn(rng)
                h += r
            end
            distance[n, i] = h
        end
    end
    return distance
end

# Erosion/deposition once a year, unidirectional
function unidirectionaldistance(rng, N, T)
    @assert eachindex(T)==1:length(T)
    distance = zeros(N, length(T))
    rate = sqrt(pi)/sqrt(2)

    @inbounds for i in eachindex(T)
        for n in 1:N
            h = 0.
            for t in 1:T[i]
                r = rate*abs(randn(rng))
                h += r
            end
            distance[n, i] = h
        end
    end
    return distance
end

## -- Run and plot erosion/deposition as a regular event

rng = MersenneTwister(0x12345)
N = 10_000
Nplot = 500
T = 10.0.^(0:0.5:6)

hr = plot(framestyle=:box,
    fontfamily=:Helvetica,
    xscale=:log10, 
    yscale=:log10,
    xticks=10.0.^(0:6),
    xlims=10.0.^(0,6),
    yticks=10.0.^(-4:1),
    ylims=10.0.^(-3.5,0.5),
    fg_color_legend=:white,
    ylabel="Rate [mm/yr]",
    xlabel="Time [yr]",
    size=(600,375),
    legend=:none,
)

@time bdistance = browniandistance(rng, N, T)
@time unidistance = unidirectionaldistance(rng, N, T)

rates = bdistance[1:Nplot, :]'./T
map!(x->x<=0 ? NaN : x, rates, rates)
plot!(T, rates, 
    label="",
    color=lines[4],
    seriestype=:scatter,
    markersize=1,
    mswidth=0,
    alpha=0.1,
)

rates = unidistance[1:Nplot, :]'./T
plot!(T, rates, 
    label="",
    color=lines[8],
    seriestype=:scatter,
    markersize=1,
    mswidth=0,
    alpha=0.1,
)

plot!(T, nanmean(abs.(bdistance), dim=1)./T,
    # label="Bidirectional (Brownian) process", 
    color=lines[4],
    seriestype=:scatter,
    markersize=5,
    mswidth=0,
)

plot!(T, nanmean(unidistance, dim=1)./T,
    # label="Unidirectional process", 
    color=lines[8],
    seriestype=:scatter,
    markersize=5,
    mswidth=0,
)

# Plot and annotate ideal constant rate line for unidirectional process
x = range(log10.(xlims())..., step=0.5)
plot!(10.0.^x, ones(size(x)), color=lines[8],  label="")
annotate!(10.0^3, 1.5*1, ("Unidirectional process: slope 0", 9, 0.0, :center, :bottom, lines[8]), fontfamily=:Helvetica)
# Plot and annotate ideal sqrt(t)/t rate line for bidirectional process
plot!(10.0.^x, 10.0.^(-0.5x), color=lines[4], label="", linestyle=:dash)
annotate!(10.0^3, 1.5*10.0.^(0 - 0.5*3), ("Bidirectional (Brownian) process: slope -0.5", 9, -45.0/2, :center, :bottom, lines[4]), fontfamily=:Helvetica)

savefig(hr, "regular_process.pdf")
display(hr)


## --- Consider erosion/deposition as a rare event

function distmean(x::Truncated{<:Pareto})
    H = if isnothing(x.upper)
        Inf
    else
        x.upper
    end
    L = if isnothing(x.lower)
        -Inf
    else
        x.lower
    end
    L = max(L, x.untruncated.θ)
    α = x.untruncated.α

    L^α/(1-(L/H)^α) * α/(α-1) * (1/L^(α-1) - 1/H^(α-1))
end

# Erosion/deposition at irregular intervals, bidirectional
function rarebrowniandistance(rng, dist, N, Nemax, T)
    distance = zeros(N, length(T))
    rateconst = sqrt(pi)/sqrt(2)*distmean(dist)

    hdist = zeros(Nemax+1)
    tdist = zeros(Nemax+1)
    @inbounds for n in 1:N
        h = t = 0.
        for j in 1:Nemax    
            h += rateconst*randn(rng) # Random height increment
            t += rand(rng, dist)  # Random time increment (hiatus duration)
            hdist[j+1] = h
            tdist[j+1] = t
        end
        for i in eachindex(T)
            t₀ = rand()*(last(tdist)-T[i])
            jₗ = jᵤ = 0
            for j in 1:Nemax
                if tdist[j+1] > t₀
                    jₗ = j
                    break
                end
            end
            for j in jₗ:Nemax
                if tdist[j+1] > t₀+T[i]
                    jᵤ = j
                    break
                end
            end
            
            distance[n, i] = hdist[jᵤ]-hdist[jₗ]
        end
    end

    return distance
end

# Erosion/deposition at irregular intervals, unidirectional
function rareunidirectionaldistance(rng, dist, N, Nemax, T)
    distance = zeros(N, length(T))
    rateconst = sqrt(pi)/sqrt(2)*distmean(dist)

    hdist = zeros(Nemax+1)
    tdist = zeros(Nemax+1)
    @inbounds for n in 1:N
        h = t = 0.
        for j in 1:Nemax    
            h += abs(rateconst*randn(rng)) # Random height increment
            t += rand(rng, dist)  # Random time increment (hiatus duration)
            # t += rand_pareto(rng, 0.5, 200_000)
            hdist[j+1] = h
            tdist[j+1] = t
        end
        for i in eachindex(T)
            t₀ = rand()*(last(tdist)-T[i])
            jₗ = jᵤ = 0
            for j in 1:Nemax
                if tdist[j+1] > t₀
                    jₗ = j
                    break
                end
            end
            for j in jₗ:Nemax
                if tdist[j+1] > t₀+T[i]
                    jᵤ = j
                    break
                end
            end
            
            distance[n, i] = hdist[jᵤ]-hdist[jₗ]
        end
    end

    return distance
end

## --- Run and plot erosion/deposition as a rare event (Pareto-distributed), bidirectional

alpha = 0.5 # Pareto shape parameter [unitless]
upper = 200_000 # Upper bound [years]
dist = truncated(Pareto(0.5); upper)

rng = MersenneTwister(0x12345)
N = 10_000
Nplot = 500
Nemax = 2^18
T = (10).^(0:0.5:8)

hp = plot(framestyle=:box,
    fontfamily=:Helvetica,
    xscale=:log10,
    yscale=:log10,
    xlabel="Time [yr]",
    ylabel="Rate [mm/yr]",
    xticks=10.0.^(0:8),
    xlims=10.0.^(0,8),
    yticks=10.0.^(-3:3),
    ylims=10.0.^(-3,3),
    size=(600, 420),
    legend=:bottomleft,
    fg_color_legend=:white,
)

@time distance = rarebrowniandistance(rng, dist, N, Nemax, T)
lightpurple = 0.4*lines[4]+0.6*parse(Color, "#ffffff")

# Excluding zero rates
map!(x->x==0 ? NaN : x, distance, distance)
plot!(hp, T, abs.(distance[1:Nplot,:]')./T, 
    label="", 
    color=lines[4],
    seriestype=:scatter,
    markersize=1,
    mswidth=0,
    alpha=0.1,
)
μₕ = nanmean(abs.(distance), dim=1)./T
plot!(hp, T, μₕ, 
    label="Excluding hiata", 
    color=lightpurple,
    seriestype=:scatter,
    markersize=5,
    mswidth=0,
)

# Fit annotated line to bidirectional process excluding hiata at short timescales 
# t = 0 .< T .<= 10^3
# x, y = log10.(T[t]), log10.(μₕ[t])
# a, b = linreg(x, y)
# plot!(10.0.^x, 10.0.^(a .+ b*x), color=lightpurple, label="", linestyle=:dash)
# annotate!(10.0^x[2], 0.8*10.0.^(a .+ b*x[2]), ("Slope: $(round(b, digits=3))", 9, atan(b)*180/pi, :left, :top, lightpurple))

# All rates
map!(x->isnan(x) ? 0 : x, distance, distance)
μₐ = nanmean(abs.(distance), dim=1)./T
plot!(hp, T, μₐ, 
    label="Including hiata", 
    color=lines[4],
    seriestype=:scatter,
    markersize=5,
    mswidth=0,
)

# Fit annotated line to bidirectional process at long timescales
t = upper .< T
x, y = log10.(T[t]), log10.(μₐ[t])
a, b = linreg(x, y)
x = range(log10(upper), log10(maximum(xlims(hp))), step=0.1)
plot!(hp, 10.0.^x, 10.0.^(a .+ b*x), color=lines[4], label="", linestyle=:dash)
annotate!(hp, 10.0^6.75, 2*10.0.^(a .+ b*6.75), ("Bidirectional: slope $(round(b, digits=2))", 9, -22.5, :center, :bottom, lines[4]), fontfamily=:Helvetica)

display(hp)


# --- Run and plot erosion/deposition as a rare event (Pareto-distributed), unidirectional

@time distance = rareunidirectionaldistance(rng, dist, N, Nemax, T)
lightblue = 0.4*lines[8]+0.6*parse(Color, "#ffffff")

# Excluding zero rates
map!(x->x==0 ? NaN : x, distance, distance)
rates = distance[1:Nplot,:]'./T
plot!(hp, T, rates, 
    label="", 
    color=lines[8],
    seriestype=:scatter,
    markersize=1,
    mswidth=0,
    alpha=0.1,
)
μₕ = nanmean(distance, dim=1)./T
plot!(hp, T, μₕ, 
    label="Excluding hiata", 
    color=lightblue,
    seriestype=:scatter,
    markersize=5,
    mswidth=0,
)

# Fit annotated line to unidirectional process excluding hiata at short timescales 
# t = 10 .< T .< upper/10
# x, y = log10.(T[t]), log10.(μₕ[t])
# a, b = linreg(x, y)
# plot!(10.0.^x, 10.0.^(a .+ b*x), color=lines[1], label="", linestyle=:dash)
# annotate!(10.0^x[2], 1.25*10.0.^(a .+ b*x[2]), ("Slope: $(round(b, digits=2))", 9, atan(b)*180/pi, :left, :bottom, lines[1]))

# All rates
map!(x->isnan(x) ? 0 : x, distance, distance)
plot!(hp, T, nanmean(distance, dim=1)./T, 
    label="Including hiata", 
    color=lines[8],
    seriestype=:scatter,
    markersize=5,
    mswidth=0,
)

# Plot and annotate ideal constant rate for unidirectional process
plot!(hp, T, ones(size(T)), color=lines[8], label="",)
annotate!(hp, 10.0^6.75, 2*1, ("Unidirectional: slope 0", 9, 0.0, :center, :bottom, lines[8]), fontfamily=:Helvetica)

# Illustrate the range of timescales larger than the maximum hiatus duration of the Pareto distribution 
xl = xlims(hp)
plot!(hp, [upper, maximum(xl)], [10^3, 10^3], fillto=[10^-3, 10^-3], alpha=0.1, color=lines[1], label="")
vline!(hp, [upper], color=lines[1], linestyle=:dot, label="")
xlims!(hp, xl)

savefig(hp, "pareto_process.pdf")
display(hp)


## ---

## ---

