## --- Load packages
using Plots, Random, Distributions, StatGeochem
using ProgressMeter: @showprogress

# Some colors we'll want later
lightpurple = 0.4*lines[4]+0.6*parse(Color, "#ffffff")
lightblue = 0.4*lines[8]+0.6*parse(Color, "#ffffff")

# Number of example simulations to plot as small dots
Nplot = 500

## --- Regular process: Brownian motion as a model for sediment accumulation/erosion

# Erosion/deposition once a year, bidirectional
function browniandistance(rng, N, T)
    @assert eachindex(T)==1:length(T)
    distance = zeros(N, length(T))
    rate = sqrt(pi)/sqrt(2)

    @inbounds @showprogress for n in 1:N
        for i in eachindex(T)
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

    @inbounds @showprogress for n in 1:N
        for i in eachindex(T)
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

# Number of simulations to run
N = 10_000
# Initialize PRNG
rng = MersenneTwister(0x1234)
# Set of time steps to run
T = 10.0.^(0:0.5:7.5)

hr = plot(framestyle=:box,
    fontfamily=:Helvetica,
    xscale=:log10, 
    yscale=:log10,
    xticks=10.0.^(0:8),
    xlims=10.0.^(0,8.01),
    yticks=10.0.^(-5:1),
    ylims=10.0.^(-4.5,0.5),
    fg_color_legend=:white,
    ylabel="Rate [mm/yr]",
    xlabel="Averaging timescale [yr]",
    size=(600,375),
    legend=:none,
    title="Regular events",
)

@time bdistance = browniandistance(rng, N, T)
@time unidistance = unidirectionaldistance(rng, N, T)

rates = bdistance[1:Nplot, :]'./T
map!(x->x<=0 ? NaN : x, rates, rates)
plot!(T, rates, 
    label="",
    color=lightpurple,
    seriestype=:scatter,
    markersize=1,
    mswidth=0,
)

rates = unidistance[1:Nplot, :]'./T
plot!(T, rates, 
    label="",
    color=lightblue,
    seriestype=:scatter,
    markersize=1,
    mswidth=0,
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
annotate!(10.0^4, 1.5*1, ("Unidirectional process: slope 0", 9, 0.0, :center, :bottom, lines[8]), fontfamily=:Helvetica)
# Plot and annotate ideal sqrt(t)/t rate line for bidirectional process
plot!(10.0.^x, 10.0.^(-0.5x), color=lines[4], label="", linestyle=:dash)
annotate!(10.0^4, 1.5*10.0.^(0 - 0.5*4), ("Bidirectional (Brownian) process: slope -0.5", 9, -45.0/2, :center, :bottom, lines[4]), fontfamily=:Helvetica)

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
function rarebrowniandistance(rng, dist, N, T)
    distance = zeros(N, length(T))
    rateconst = sqrt(pi)/sqrt(2)*distmean(dist)
    Nₜmax = round(Int, 2*maximum(T)/distmean(dist))
    bigenough = round(Int, 1.25*maximum(T))

    hdist = zeros(Nₜmax+1)
    tdist = zeros(Nₜmax+1)
    @inbounds @showprogress for n in 1:N
        Nₜlast = 1
        h = t = 0.
        for j in 1:Nₜmax    
            h += rateconst*randn(rng) # Random height increment
            t += rand(rng, dist)  # Random time increment (hiatus duration)
            Nₜlast = j+1
            hdist[Nₜlast] = h
            tdist[Nₜlast] = t
            if t > bigenough
                break
            end
        end

        for i in eachindex(T)
            t₀ = rand(rng)*(tdist[Nₜlast]-T[i])
            jₗ = jᵤ = 1
            for j in 1:Nₜmax
                if tdist[j+1] > t₀
                    jₗ = j
                    break
                end
            end
            for j in jₗ:Nₜmax
                if tdist[j+1] > (t₀+T[i])
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
function rareunidirectionaldistance(rng, dist, N, T)
    distance = zeros(N, length(T))
    rateconst = sqrt(pi)/sqrt(2)*distmean(dist)
    Nₜmax = round(Int, 2*maximum(T)/distmean(dist))
    bigenough = round(Int, 1.25*maximum(T))

    hdist = zeros(Nₜmax+1)
    tdist = zeros(Nₜmax+1)
    @inbounds @showprogress for n in 1:N
        Nₜlast = 1
        h = t = 0.
        for j in 1:Nₜmax    
            h += abs(rateconst*randn(rng)) # Random height increment
            t += rand(rng, dist)  # Random time increment (hiatus duration)
            Nₜlast = j+1
            hdist[Nₜlast] = h
            tdist[Nₜlast] = t
            if t > bigenough
                break
            end
        end

        for i in eachindex(T)
            t₀ = rand(rng)*(tdist[Nₜlast]-T[i])
            jₗ = jᵤ = 1
            for j in 1:Nₜmax
                if tdist[j+1] > t₀
                    jₗ = j
                    break
                end
            end
            for j in jₗ:Nₜmax
                if tdist[j+1] > (t₀+T[i])
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

# Number of simulations to run
N = 100_000

# Rare event distribution to draw times from
alpha = 0.5 # Pareto shape parameter [unitless]
upper = 200_000 # Upper bound [years]
dist = truncated(Pareto(0.5); upper)

# Initialize PRNG
rng = MersenneTwister(0x1234)
# Set of time steps to run
T = 10.0.^(0:0.5:8)


hp = plot(framestyle=:box,
    fontfamily=:Helvetica,
    xscale=:log10,
    yscale=:log10,
    xlabel="Averaging timescale [yr]",
    ylabel="Rate [mm/yr]",
    xticks=10.0.^(0:8),
    xlims=10.0.^(0,8.01),
    yticks=10.0.^(-3:3),
    ylims=10.0.^(-3,3),
    size=(600, 420),
    legend=:bottomleft,
    fg_color_legend=:white,
    title="Rare events",
)

@time distance = rarebrowniandistance(rng, dist, N, T)

# Excluding zero rates
map!(x->x==0 ? NaN : x, distance, distance)
plot!(hp, T, abs.(distance[1:Nplot,:]')./T, 
    label="", 
    color=lightpurple,
    seriestype=:scatter,
    markersize=1,
    mswidth=0,
)
μₕ = nanmean(abs.(distance), dim=1)./T
plot!(hp, T, μₕ, 
    label="Eventful intervals", 
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
    label="All intervals", 
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

@time distance = rareunidirectionaldistance(rng, dist, N, T)

# Excluding zero rates
map!(x->x==0 ? NaN : x, distance, distance)
rates = distance[1:Nplot,:]'./T
plot!(hp, T, rates, 
    label="", 
    color=lightblue,
    seriestype=:scatter,
    markersize=1,
    mswidth=0,
)
μₕ = nanmean(distance, dim=1)./T
plot!(hp, T, μₕ, 
    label="Eventful intervals", 
    color=lightblue,
    seriestype=:scatter,
    markersize=5,
    mswidth=0,
)

# # Fit annotated line to unidirectional process excluding hiata at short timescales 
# t = 10 .< T .< 100_000
# x, y = log10.(T[t]), log10.(μₕ[t])
# a, b = linreg(x, y)
# plot!(10.0.^x, 10.0.^(a .+ b*x), color=lightblue, label="", linestyle=:dot)
# annotate!(10.0^mean(x), 1.75*10.0.^(a .+ b*mean(x)), ("slope: $(round(b, digits=2))", 9, atan(b)*180/pi, :center, :bottom, lightblue))

# All rates
map!(x->isnan(x) ? 0 : x, distance, distance)
plot!(hp, T, nanmean(distance, dim=1)./T, 
    label="All intervals", 
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


## --- Plot fraction of intervals with an event

hc = plot(
    framestyle=:box,
    xlabel="Averaging timescale [yr]",
    ylabel="Percent of intervals",
    xticks=0:8,
    xlims=(0,8.01),
    size=(600,200),
    legend=:topleft,
    fg_color_legend=:white,
)

nodata = vec(sum(x->x==0, distance, dims=1))
hasdata = N .- nodata

plot!(hc, log10.(T), fill(100., length(T)),
    seriestype=:bar,
    color=:black,
    bar_width=0.25,
    linewidth=0,
    label="No event in interval",
)

plot!(hc, log10.(T), hasdata/N*100,
    seriestype=:bar,
    color=0.5*lines[8]+0.5*lines[1],
    bar_width=0.25,
    linewidth=0,
    label="Event(s) in interval",
)

# Illustrate the range of timescales larger than the maximum hiatus duration of the Pareto distribution 
xl = xlims(hc)
plot!(hc, [log10(upper), maximum(xl)], [100, 100], fillto=[0.,0.], alpha=0.15, color=lines[1], label="")
vline!(hc, [log10(upper)], color=lines[1], linestyle=:dot, label="")
xlims!(hc, xl)

savefig(hc, "fraction_event.pdf")
display(hc)


## --- Plot all together

h = plot(hr, hp, hc,
    layout=grid(3,1, heights=[0.4, 0.45, 0.15,]),
    size=(600,1000)
)
savefig(h, "numerical_model.pdf")
display(h)


## --- Plot example timeseries for regular event process

N = 100
hdist = randn(N)*sqrt(pi/2)
tdist = 1:N

h = plot(vec([tdist tdist tdist]'), vec([zeros(size(hdist)) abs.(hdist) zeros(size(hdist))]'), 
    color=lines[8], 
    lw=2,
    label="",
    # xticks=:none,
    # yticks=:none,
    # xlabel="Time [yr]",
    # ylabel="Erosion [mm]",
    # xlims=extrema(tdist),
    framestyle=:none,
    size=(600,100),
)

savefig(h, "regular_process_example.pdf")
display(h)

## --- Plot example timeseries for rare event process

N = 100
hdist = randn(N)*sqrt(pi/2)
tdist = cumsum(rand(truncated(Pareto(0.5), upper=200_000), N))

h = plot(vec([tdist tdist tdist]'), vec([zeros(size(hdist)) abs.(hdist) zeros(size(hdist))]'), 
    color=lines[8], 
    lw=2,
    label="",
    # xticks=:none,
    # yticks=:none,
    # xlabel="Time [yr]",
    # ylabel="Erosion [mm]",
    # xlims=extrema(tdist),
    framestyle=:none,
    size=(600,100),
)

savefig(h, "pareto_process_example.pdf")
display(h)

## --- End of File
