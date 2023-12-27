## --- Regular process: Brownian motion as a model for sediment accumulation/erosion

using Plots, Random, Distributions, StatGeochem
rng = Xoshiro()
N = 10_00
Nplot = 500
T = (10).^(0:0.5:6)

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

hr = plot(framestyle=:box,
    xscale=:log10, 
    yscale=:log10,
    legend=:bottomleft,
    xticks=10.0.^(0:1:6),
    yticks=10.0.^(-4:1),
    ylims=10.0.^(-3.5,0.5),
    fg_color_legend=:white,
    ylabel="Rate [mm/yr]",
    xlabel="Time [yr]",
    size=(600,350),
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
    label="Bidirectional (Brownian) process", 
    color=lines[4],
    seriestype=:scatter,
    markersize=5,
    mswidth=0,
)

plot!(T, nanmean(unidistance, dim=1)./T,
    label="Unidirectional process", 
    color=lines[8],
    seriestype=:scatter,
    markersize=5,
    mswidth=0,
)

x = 0:0.1:6
plot!(10.0.^x, ones(size(x)), color=:black,  label="", )
plot!(10.0.^x, 10.0.^(-0.5x), color=lines[4], label="", )
annotate!(10.0^2.5, 1.25*10.0.^(0 - 0.5*2.5), ("Slope: -0.5", 9, -45.0/2, :left, :bottom, lines[4]))

savefig(hr, "regular_process.pdf")
display(hr)


## --- Pareto process: Consider erosion/deposition as a rare event, bidirectional

alpha = 0.5 # Pareto shape parameter [unitless]
upper = 200_000 # Upper bound [years]
dist = truncated(Pareto(0.5); upper)

using Plots, Random
rng = Xoshiro()
N = 2^15
Nplot = 500
Nemax = 2^18
T = (10).^(0:0.5:8)

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

hp = plot(framestyle=:box,
    xscale=:log10,
    yscale=:log10,
    xlabel="Time [yr]",
    ylabel="Rate [mm/yr]",
    xticks=10.0.^(0:8),
    yticks=10.0.^(-5:3),
    ylims=10.0.^(-3,3),
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
    label="Bidirectional, excluding hiata", 
    color=lightpurple,
    seriestype=:scatter,
    markersize=5,
    mswidth=0,
)
t = 0 .< T .<= 10^3
x, y = log10.(T[t]), log10.(μₕ[t])
a, b = linreg(x, y)
plot!(10.0.^x, 10.0.^(a .+ b*x), color=lightpurple, label="", linestyle=:dash)
annotate!(10.0^x[2], 0.8*10.0.^(a .+ b*x[2]), ("Slope: $(round(b, digits=3))", 9, atan(b)*180/pi, :left, :top, lightpurple))

# All rates
map!(x->isnan(x) ? 0 : x, distance, distance)
μₐ = nanmean(abs.(distance), dim=1)./T
plot!(hp, T, μᵣ, 
    label="Bidirectional", 
    color=lines[4],
    seriestype=:scatter,
    markersize=5,
    mswidth=0,
)

t = upper .< T
x, y = log10.(T[t]), log10.(μₐ[t])
a, b = linreg(x, y)
plot!(hp, 10.0.^x, 10.0.^(a .+ b*x), color=lines[4], label="", linestyle=:dash)
annotate!(hp, 10.0^x[2], 1.25*10.0.^(a .+ b*x[2]), ("Slope: $(round(b, digits=2))", 9, -22.5, :left, :bottom, lines[4]))

display(hp)


## --- Pareto process: Consider erosion/deposition as a rare event, unidirectional

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

@time distance = rareunidirectionaldistance(rng, dist, N, Nemax, T)

# h = plot(framestyle=:box,
#     xscale=:log10,
#     yscale=:log10,
#     xticks=(10).^(0:8),
# )

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
    label="Unidirectional, excluding hiata", 
    color=lines[1],
    seriestype=:scatter,
    markersize=5,
    mswidth=0,
)

# Annotate unidirectional slope with hiata
# t = 10 .< T .< upper/10
# x, y = log10.(T[t]), log10.(μₕ[t])
# a, b = linreg(x, y)
# plot!(10.0.^x, 10.0.^(a .+ b*x), color=lines[1], label="", linestyle=:dash)
# annotate!(10.0^x[2], 1.25*10.0.^(a .+ b*x[2]), ("Slope: $(round(b, digits=2))", 9, atan(b)*180/pi, :left, :bottom, lines[1]))

# All rates
map!(x->isnan(x) ? 0 : x, distance, distance)
plot!(hp, T, nanmean(distance, dim=1)./T, 
    label="Unidirectional", 
    color=lines[8],
    seriestype=:scatter,
    markersize=5,
    mswidth=0,
)

plot!(hp, T, ones(size(T)), color=:black, label="",)

savefig(hp, "pareto_process.pdf")
display(hp)


## ---

## ---

# x = 0:0.01:8
# t = 10.0.^x
# y = cdf.(dist, t)
# plot(t, y, xscale=:log10)
