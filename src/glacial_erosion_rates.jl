using StatGeochem, Plots

datadir = "$(@__DIR__)/../data"
earth = importdataset("$datadir/glacial_erosion_Earth.tsv", '\t', importas=:Tuple);
mars = importdataset("$datadir/glacial_erosion_Mars.tsv", '\t', importas=:Tuple);

ngearth = importdataset("$datadir/nonglacial_erosion_Earth.tsv", '\t', importas=:Tuple)
ngmars = importdataset("$datadir/nonglacial_erosion_Mars.tsv", '\t', importas=:Tuple)


## --- Utility functions for plotting

# Scatter plot, but plot any points outside of ylims as up/down-ward pointing triangles
function scatternicely!(h, x, y; offset=1.08, label, kwargs...)
    plot!(h, x, y;
        seriestype=:scatter, label,
        kwargs...,
    )
    ymax = maximum(ylims(h))
    tu = y .> ymax
    plot!(h, x[tu], fill(ymax/offset, count(tu));
        seriestype=:scatter, markershape=:utriangle, label="",
        kwargs...,
    )
    ymin = minimum(ylims(h))
    td = y .< ymin
    plot!(h, x[td], fill(ymin*offset, count(td));
        seriestype=:scatter, markershape=:dtriangle, label="",
        kwargs...,
    )
    return h
end

# For custom histogram plots
stepify(x::AbstractVector) = vec(vcat(x', x'))
stepifyedges(x::AbstractVector) = vec(vcat(x[1:end-1]', x[2:end]'))

# # Linear Regression (already defined in StatGeochem)
# linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y

# Calculate the R-squared value for a line through x-y data
rsquared(x,y,line_y) = 1 - nansum((y .- line_y).^2)  / nansum((y .- nanmean(y)).^2)


## -- Area vs erosion rate

ha = plot(framestyle=:box,
    xlabel="Area [km²]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    xlims = (10^-2.5, 10^7.5),
    xticks = 10.0.^(-2:7),
    ylims = (10^-4, 10^3),
    yticks = 10.0.^(-4:3),
    fg_color_legend = :white,
)

t = (ngearth.Area_km2 .> 1e-6) .& (ngearth.Erosion_rate_mm_yr .>0)
scatternicely!(ha, ngearth.Area_km2[t], ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

t = (earth.Area_km2 .> 1e-6) .& (earth.Erosion_rate_mm_yr .>0)
scatternicely!(ha, earth.Area_km2[t], earth.Erosion_rate_mm_yr[t],
    label="Glacial",
    color=mineralcolors["azurite"],
    alpha=0.85,
    mswidth=0.25,
)

savefig(ha, "area_vs_erosion_rate.pdf")
display(ha)


## --- Latitude vs erosion rate, colored by glacier type (Continental, Alpine, etc.),

hl = plot(framestyle=:box,
    xlabel="Latitude [degrees, absolute]",
    ylabel="Erosion rate [mm/yr]",
    yscale=:log10,
    fontfamily=:Helvetica,
    fg_color_legend=:white,
    legend=:bottomleft,
    xlims=(0,90),
    xticks=0:15:90,
    ylims = (10^-5, 10^3),
    yticks = 10.0.^(-5:3),

)

t = .!isnan.(ngearth.Latitude) .& (ngearth.Erosion_rate_mm_yr .>0)
scatternicely!(hl, abs.(ngearth.Latitude[t]), ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

tdv = (earth.Type .== "Dry Valleys") .& contains.(earth.Methodology, "Cosmogenic")
scatternicely!(hl, abs.(earth.Latitude[tdv]), earth.Erosion_rate_mm_yr[tdv],
    label="Dry Valleys Cosmogenic",
    color=mineralcolors["spessartine"],
    alpha=0.85,
    mswidth=0.25,
)
scatternicely!(hl, abs.(earth.Latitude[.!tdv]), earth.Erosion_rate_mm_yr[.!tdv],
    label="Glacial",
    color=mineralcolors["azurite"],
    alpha=0.85,
    mswidth=0.25,
)

savefig(hl, "latitude_vs_erosion_rate.pdf")
display(hl)


## --- Slope vs erosion rate, colored by glacier type (Continental, Alpine, etc.),

hs = plot(framestyle=:box,
    xlabel="Regional slope [m/m]",
    ylabel="Erosion rate [mm/yr]",
    yscale=:log10,
    fontfamily=:Helvetica,
    fg_color_legend=:white,
    legend=:bottomright,
    xlims=(0,0.6),
    xticks=0:0.1:0.6,
    ylims = (10^-5, 10^3),
    yticks = 10.0.^(-5:3),

)

t = .!isnan.(ngearth.Slope_m_km) .& (ngearth.Erosion_rate_mm_yr .> 0)
slope = ngearth.Slope_m_km[t]./1000
scatternicely!(hs, slope, ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

tdv = (earth.Type .== "Dry Valleys") .& contains.(earth.Methodology, "Cosmogenic")
scatternicely!(hs, earth.Slope_m_km[tdv]/1000, earth.Erosion_rate_mm_yr[tdv],
    label="Dry Valleys Cosmogenic",
    color=mineralcolors["spessartine"],
    alpha=0.85,
    mswidth=0.25,
)
scatternicely!(hs, earth.Slope_m_km[.!tdv]/1000, earth.Erosion_rate_mm_yr[.!tdv],
    label="Glacial",
    color=mineralcolors["azurite"],
    alpha=0.85,
    mswidth=0.25,
)

aₙ,bₙ = linreg(slope, log10.(ngearth.Erosion_rate_mm_yr[t]))
xₙ = collect(xlims(hs))
yₙ = @. 10^(aₙ+bₙ*xₙ)
plot!(hs, xₙ, yₙ, color=parse(Color, "#777777"), linestyle=:dot, label="$(round(10^aₙ,digits=3))*10^($(round(bₙ,digits=3)) s)")

t = .!isnan.(earth.Slope_m_km) .& (earth.Erosion_rate_mm_yr .> 0) .& .!tdv
slope = earth.Slope_m_km[t]./1000
a,b = linreg(slope, log10.(earth.Erosion_rate_mm_yr[t]))
x = collect(xlims(hs))
y = @. 10^(a+b*x)

plot!(hs, x, y, color=:darkblue, linestyle=:dot, label="$(round(10^a,digits=3))*10^($(round(b,digits=3)) s)")

savefig(hs, "slope_vs_erosion_rate.pdf")
display(hs)

## --- Precipitation

hp = plot(framestyle=:box,
    xlabel="Precipitation [mm/yr]",
    ylabel="Erosion rate [mm/yr]",
    yscale=:log10,
    fontfamily=:Helvetica,
    fg_color_legend=:white,
    legend=:bottomright,
    xlims = (0, 9000),
    xticks = 0:1500:9000,
    ylims = (10^-5, 10^3),
    yticks = 10.0.^(-5:3),
)

t = ngearth.Erosion_rate_mm_yr .> 0
scatternicely!(hp, ngearth.Precipitation_mm_yr[t], ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color = parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

tdv = (earth.Type .== "Dry Valleys") .& contains.(earth.Methodology, "Cosmogenic")
scatternicely!(hp, earth.Precipitation_mm_yr[tdv], earth.Erosion_rate_mm_yr[tdv],
    label="Dry Valleys Cosmogenic",
    color=mineralcolors["spessartine"],
    alpha=0.85,
    mswidth=0.25,
)
scatternicely!(hp, earth.Precipitation_mm_yr[.!tdv], earth.Erosion_rate_mm_yr[.!tdv],
    label="Glacial",
    color=mineralcolors["azurite"],
    alpha=0.85,
    mswidth=0.25,
)

savefig(hp, "precipitation_vs_erosion_rate.pdf")
display(hp)

## --- Plot combined comparison

h = plot(ha, hl, hp, hs,
    size = (1600, 800/3),
    layout = (1,4),
)
savefig(h, "erosion_rate_comparison.pdf")
display(h)


## --- Timescale vs erosion rate, colored by glacier type (Continental, Alpine, etc.),

h = plot(framestyle=:box,
    xlabel="Timescale [yr]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    fg_color_legend=:white,
    legend=:bottomleft,
    xlims = (10^-2, 10^9),
    xticks = 10.0.^(-2:9),
    ylims = (10^-5, 10^3),
    yticks = 10.0.^(-5:3),
)

t = (ngearth.Time_interval_yr .> 0) .& (ngearth.Erosion_rate_mm_yr .>0)
scatternicely!(h, ngearth.Time_interval_yr[t], ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

type = ( "Alpine", "Alpine tidewater", "High-latitude", "Continental", "Outlet ice stream", "Dry Valleys",)
tlabel = ( "Alpine", "Alpine tidewater", "Non-cont. high-lat.", "Continental", "Outlet ice stream", "Dry Valleys",)
colors = [mineralcolors[m] for m in ("malachite", "zircon", "quartz", "azurite", "fluid", "spessartine", )]
for i in eachindex(type)
    t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Type .== type[i])
    plot!(h, earth.Time_interval_yr[t], earth.Erosion_rate_mm_yr[t],
        seriestype=:scatter,
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = tlabel[i],
    )
end

savefig(h, "timescale_vs_erosion_rate-type.pdf")
display(h)


## --- Histogram of erosion rates!

# Pick binning scheme
t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0)
μ = nanmean(log10.(earth.Erosion_rate_mm_yr[t]))
binedges = (-4:0.2:4) .+ μ
bincenters = cntr(binedges)
distx = first(binedges):0.01:last(binedges)

hist = plot(framestyle=:box,
    xlabel="Erosion Rate [mm/yr]",
    ylabel="N (nonglacial)",
    xflip=true,
    xlims=(10^-5,10^3),
    xticks=10.0.^(-5:3),
    xscale=:log10,
    size=(600,300),
    rotation=90.,
    xguidefontrotation=180.,
    yguidefontcolor=parse(Color, "#888888"),
    ytickfontcolor=parse(Color, "#888888"),
    fontfamily=:Helvetica,
)

# Nonglacial erosion
tng = (ngearth.Time_interval_yr .> 0) .& (ngearth.Erosion_rate_mm_yr .> 0)
logerosion = log10.(ngearth.Erosion_rate_mm_yr[tng])
Ns = histcounts(logerosion, binedges, T=Float64)
plot!(hist, 10.0.^stepifyedges(binedges), stepify(Ns),
    fill=true,
    color=parse(Color, "#dddddd"),
    label=""
)

ngμ, ngσ = nanmean(logerosion), nanstd(logerosion)
plot!(hist, 10.0.^distx, normpdf(ngμ,ngσ,distx) .* (sum(Ns) * step(binedges)),
    linestyle=:solid,
    color=parse(Color, "#888888"),
    label="",
)

# Glacial erosion
Ns = histcounts(log10.(earth.Erosion_rate_mm_yr[t]), binedges)
hist2 = twinx()
plot!(hist2,
    framestyle=:box,
    ylabel="N (glacial)",
    xflip=true,
    xlims=(10^-5,10^3),
    xticks=10.0.^(-5:3),
    xscale=:log10,
    size=(400,200),
    rotation=90.,
)

plot!(hist2, 10.0.^stepifyedges(binedges), stepify(Ns), fill=true, color=mineralcolors["spessartine"], label="")

t .&= .!((earth.Type .== "Dry Valleys") .& contains.(earth.Methodology, "Cosmogenic"))
logerosion = log10.(earth.Erosion_rate_mm_yr[t])
Ns = histcounts(logerosion, binedges, T=Float64)
plot!(hist2, 10.0.^stepifyedges(binedges), stepify(Ns), fill=true, color=mineralcolors["fluid"], label="")

μ, σ = nanmean(logerosion), nanstd(logerosion)
plot!(hist2, 10.0.^distx, normpdf(μ,σ,distx) .* (sum(Ns)*step(binedges)),
    linestyle=:solid,
    color=:darkblue,
    label="",
)

# Clean up y axis limits
ylims!(hist, 0, 800)
ylims!(hist2, 0, 70)

# Dry valleys
annotate!(hist2, [10^-4.3], [maximum(ylims(hist2))], text("Dry Valleys Cosmogenic ", 7, :right, :bottom, rotation=90, color=mineralcolors["spessartine"]))

# Glacial Gaussian
annotate!(hist2, [10.0.^μ], [0], text(" $(round(10^μ, sigdigits=2)) mm/yr", 7, :left, :bottom, rotation=90, color=:darkblue))
vline!(hist2, [10.0.^μ], color=:darkblue, linestyle=:dot, label="")
annotate!(hist2, [10^μ], [maximum(ylims(hist2))], text("glacial ", 7, :right, :bottom, rotation=90, color=:darkblue))

# Nonglacial Gaussian
annotate!(hist2, [10.0.^ngμ], [0], text(" $(round(10^ngμ, sigdigits=2)) mm/yr", 7, :left, :bottom, rotation=90, color=parse(Color, "#888888")))
vline!(hist2, [10.0.^ngμ], color=parse(Color, "#888888"), linestyle=:dot, label="")
annotate!(hist2, [10.0.^ngμ], [maximum(ylims(hist2))], text("nonglacial ", 7, :right, :bottom, rotation=90, color=parse(Color, "#888888")))


savefig(hist, "erosion_rate_histogram.pdf")
display(hist)


plot!(hist,
    xlims=(10^-8,10^3),
    xticks=10.0.^(-8:3),
)
plot!(hist2,
    xlims=(10^-8,10^3),
    xticks=10.0.^(-8:3),
)
savefig(hist, "erosion_rate_histogram_taller.pdf")
# display(hist)


## --- Timescale vs erosion rate, binned by area, wide version (glacial)

hg = plot(layout = (1,4),
    framestyle=:box,
    xlabel="Timescale [yr]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    size=(1200,300),
    xlims = (10^-3, 10^8),
    xticks = 10.0.^(-3:1:8),
    ylims = (10^-5.5, 10^5.5),
    yticks = 10.0.^(-5:1:5),
    legend=:none,
)

areas = (0, 0.1, 10, 1000, 10^7)
titles = ("<0.1 km²", "0.1-10 km²", "10-1000 km²", ">1000 km²")

# Add cosmogenic line of constant 600mm thickness
x = collect(xlims(hg[1]))
plot!(hg[1], x, 600.0./x, color=:red, label="", linestyle=:dash)
xl = minimum(x)
annotate!(hg[1], xl*10, 600/(xl*10), text("600 mm", 10, :bottom, :left, color=:red, rotation=-45.0))

# Add log-mean line for all others
t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Type .!= "Dry Valleys")
logμ = 10^nanmean(log10.(earth.Erosion_rate_mm_yr[t]))
for j in 1:length(areas)-1
    hline!(hg[j], [logμ], color=:darkblue, label="", linestyle=:dot)
end
annotate!(hg[1], xl*3, logμ, text("$(round(logμ, digits=2)) mm/yr", 10, :top, :left, color=:darkblue,))


method = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric", "Relief", )
mlabel = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric", "Relief", )
colors = [mineralcolors[m] for m in ( "spessartine", "rhodochrosite", "corundum", "sodalite", "sapphirine", )]

for j in 1:length(areas)-1
    t = (areas[j] .< earth.Area_km2 .< areas[j+1]) .& (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0)
    j==1 && (t .|= .!(earth.Area_km2 .> 0.1) .& (earth.Methodology .!= "Volumetric"))
    plot!(hg[j], title=titles[j])
    j > 1 && plot!(hg[j], ylabel="")

    for i in eachindex(method)
        tt = t .& (earth.Methodology .== method[i])
        plot!(hg[j], earth.Time_interval_yr[tt], earth.Erosion_rate_mm_yr[tt],
            seriestype=:scatter,
            color = colors[i],
            alpha = 0.85,
            mswidth = 0.1,
            label = mlabel[i],
        )
    end
end

savefig(hg, "timescale_vs_erosion_rate-area-glacial.pdf")
display(hg)

## ---

hng = plot(layout = (1,4),
    framestyle=:box,
    xlabel="Timescale [yr]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    size=(1200,300),
    xlims = (10^-3, 10^8),
    xticks = 10.0.^(-3:1:8),
    ylims = (10^-5.5, 10^5.5),
    yticks = 10.0.^(-5:1:5),
    legend=:none,
)

areas = (0, 0.1, 10, 1000, 10^7)
titles = ("<0.1 km²", "0.1-10 km²", "10-1000 km²", ">1000 km²")

# Add cosmogenic line of constant 600mm thickness
x = collect(xlims(hng[1]))
for j in 1:length(areas)-1
    plot!(hng[j], x, 600.0./x, color=:red, label="", linestyle=:dash)
end
xl = minimum(x)
annotate!(hng[1], xl*10, 600/(xl*10), text("600 mm", 10, :bottom, :left, color=:red, rotation=-45.0))

# Add log-mean line for all others
t = (ngearth.Time_interval_yr .> 0) .& (ngearth.Erosion_rate_mm_yr .>0) .& (ngearth.Type .!= "Dry Valleys")
logμ = 10^nanmean(log10.(ngearth.Erosion_rate_mm_yr[t]))
for j in 1:length(areas)-1
    hline!(hng[j], [logμ], color=parse(Color, "#666666"), label="", linestyle=:dot)
end
annotate!(hng[1], xl*3, logμ, text("$(round(logμ, sigdigits=2)) mm/yr", 10, :top, :left, color=parse(Color, "#666666"),))


method = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric", "Relief", )
mlabel = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric", "Relief", )
colors = [mineralcolors[m] for m in ( "spessartine", "rhodochrosite", "corundum", "sodalite", "sapphirine", )]

for j in 1:length(areas)-1
    t = (areas[j] .< ngearth.Area_km2 .< areas[j+1]) .& (ngearth.Time_interval_yr .> 0) .& (ngearth.Erosion_rate_mm_yr .>0)
    # j==1 && (t .|= .!(ngearth.Area_km2 .> 0.1))
    plot!(hng[j], title=titles[j])
    j > 1 && plot!(hng[j], ylabel="")

    for i in eachindex(method)
        tt = t .& (ngearth.Methodology .== method[i])
        plot!(hng[j], ngearth.Time_interval_yr[tt], ngearth.Erosion_rate_mm_yr[tt],
            seriestype=:scatter,
            color = colors[i]*0.25 + parse(Color, "#ffffff")*0.75,
            # alpha = 0.85,
            mswidth = 0,
            label = mlabel[i],
        )
    end
end

# Plot error bars, but on top
for j in 1:length(areas)-1
    t = (areas[j] .< ngearth.Area_km2 .< areas[j+1]) .& (ngearth.Time_interval_yr .> 0) .& (ngearth.Erosion_rate_mm_yr .>0)
    plot!(hng[j], title=titles[j])

    for i in eachindex(method)
        tt = t .& (ngearth.Methodology .== method[i])
        if count(tt) .> 50
            x = nanmean(log10.(ngearth.Time_interval_yr[tt]))
            y = nanmean(log10.(ngearth.Erosion_rate_mm_yr[tt]))
            y_sigma = nanstd(log10.(ngearth.Erosion_rate_mm_yr[tt]))

            plot!(hng[j], [10^x], [10^y],
                yerror = ([10^y - 10^(y-y_sigma)], [10^(y+y_sigma) - 10^y]),
                seriestype = :scatter,
                color = colors[i]*0.85,
                msc = colors[i]*0.85,
                label = "",
            )
        end
    end
end

savefig(hng, "timescale_vs_erosion_rate-area-nonglacial.pdf")
display(hng)


## --- Combined timescale vs erosion rate area

h = plot(hg, hng,
    size = (1200,600),
    layout = (2,1),
)
savefig(h, "timescale_vs_erosion_rate-area.pdf")
display(h)


## --- Continental icesheet erosion

tc = (earth.Type .== "Continental") .| (earth.Type .== "High-latitude") .| ((earth.Type .== "Dry Valleys") .& .!contains.(earth.Methodology, "Cosmogenic"))
h = plot(framestyle=:box,
    xlabel="Timescale [yr]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    xlims=(10^-2, 10^9),
    xticks=10.0.^(-2:9),
    ylims=(10^-2, 10^0.5),
    yticks=10.0.^(-2:0.5:0.5),
    fg_color_legend=:white,
    legend=:none,
)

# plot!(h, earth.Time_interval_yr, earth.Erosion_rate_mm_yr, color=:black, alpha=0.3, seriestype=:scatter, label="")
regions = ("Greenland", "Fennoscandia", "Antarctica", )

colors = [mineralcolors[m] for m in ("fluid", "kyanite", "azurite", )]

for i in eachindex(regions)
    region = regions[i]
    t = (earth.Region .== region) .& tc

    scatternicely!(h, earth.Time_interval_yr[t], earth.Erosion_rate_mm_yr[t];
        color=colors[i],
        label=region,
        alpha=0.85,
        mswidth=0.25,
    )
    # Draw lines for unweighted means of each ice sheet
    mean_erosion = nanmean(earth.Erosion_rate_mm_yr[t])
    erosion_sigma = nanstd(earth.Erosion_rate_mm_yr[t])./sqrt(count(t))
    hline!(h, [mean_erosion],
        linestyle=:dot,
        color=colors[i],
        label="",
    )
    annotate!(0.9e9, mean_erosion*1.05, text("$region ", colors[i], :right, :bottom, 9, fontfamily=:Helvetica))
    # # (optional) error bars on the means
    # plot!(h, [5e8], [mean_erosion], yerr=[erosion_sigma], seriestype=:scatter, msc=colors[i], color=colors[i], label="")
end

# Optional: Plot 3-5 km cryogenian erosion, for comparison
cryocolor=mineralcolors["magnetite"]
plot!(h, [64e6], [4000*1000/64e6], yerr=[1000*1000/64e6], msc=cryocolor, color=cryocolor, seriestype=:scatter, label="")
annotate!(h, 2.5e8, 4000*1000/64e6, text("3-5 km\nCryogen.\nerosion", cryocolor, :center, :8))
savefig(h, "continental_erosion_rates.pdf")

# # Add Mars!
# marscolor = mineralcolors["almandine"] # parse(Colorant, "#ff0000")
# scatternicely!(h, mars.Time_interval_yr, mars.Erosion_rate_mm_yr;
#     color=marscolor,
#     label="Mars",
#     alpha=0.85,
#     mswidth=0.25,
# )
#
# # Draw line for unweighted mean
# mean_erosion = nanmean(mars.Erosion_rate_mm_yr)
# erosion_sigma = nanstd(earth.Erosion_rate_mm_yr)./sqrt(length(mars.Erosion_rate_mm_yr))
# hline!(h, [mean_erosion],
#     linestyle=:dot,
#     color=marscolor,
#     label="",
# )
# annotate!(0.9e9, mean_erosion, text("Mars ", marscolor, :right, :top, 9, fontfamily=:Helvetica))
# savefig(h, "continental_erosion_rates_withmars.pdf")

display(h)


## --- Histogram of continental erosion rates

hist = plot(framestyle=:box,
    xlabel="Erosion Rate [mm/yr]",
    ylabel="N",
    xflip=true,
    xlims=(10^-8,10^3),
    xticks=10.0.^(-8:3),
    xscale=:log10,
    size=(400,200),
    rotation=90,
    xguidefontrotation=180.,
    fontfamily=:Helvetica,
    fg_color_legend=:white,
)

tc = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0)
tc .&= (earth.Type .== "Continental") .| (earth.Type .== "High-latitude") .| ((earth.Type .== "Dry Valleys") .& .!contains.(earth.Methodology, "Cosmogenic"))
tc .&= (earth.Region .== "Greenland") .| (earth.Region .== "Fennoscandia") .| (earth.Region .== "Antarctica")

logerosion = log10.(earth.Erosion_rate_mm_yr[tc])
μ, σ = nanmean(logerosion), nanstd(logerosion)
binedges = (-4:0.2:4) .+ μ
bincenters = cntr(binedges)

# Ns = histcounts(log10.(earth.Erosion_rate_mm_yr[t]), binedges)
# plot!(hist, 10.0.^stepifyedges(binedges), stepify(Ns), fill=true, color=mineralcolors["spessartine"], label="")

regions = ("Greenland", "Fennoscandia", "Antarctica", )
colors = [mineralcolors[m] for m in ("fluid", "kyanite", "azurite", )]
t = copy(tc)
for i in eachindex(regions)
    logerosion = log10.(earth.Erosion_rate_mm_yr[t])
    Ns = histcounts(logerosion, binedges)
    plot!(hist, 10.0.^stepifyedges(binedges), stepify(Ns),
        fill=true,
        color=colors[i],
        label=regions[i],
    )

    t .&= earth.Region .!= regions[i]
end

x = first(binedges):0.01:last(binedges)
plot!(hist, 10.0.^x, normpdf(μ,σ,x).*(count(tc)*step(binedges)),
    linestyle=:solid,
    color=:black,
    label="",
)
ylims!(hist, 0, maximum(ylims(hist)))

vline!(hist, [10.0.^μ], color=:black, linestyle=:dot, label="")
annotate!(hist, [10.0.^μ], [0], text(" $(round(10^μ, sigdigits=2)) mm/yr ", 7, :left, :bottom, rotation=90))
annotate!(hist, [10.0.^μ], [maximum(ylims(hist))], text(" continental ", 7, :right, :bottom, rotation=90))

savefig(hist, "continental_erosion_rate_histogram.pdf")
display(hist)

## -- Area vs erosion rate, Mars

h = plot(framestyle=:box,
    xlabel="Area [km²]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    xlims = (10^-2.5, 10^7.5),
    xticks = 10.0.^(-2:7),
    fg_color_legend = :white,
)

t = (ngmars.Area_km2 .> 1e-6) .& (ngmars.Erosion_rate_mm_yr .>0)
scatternicely!(h, ngmars.Area_km2[t], ngmars.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=mineralcolors["muscovite"],
    alpha=1,
    mswidth=0,
)

t = (mars.Area_km2 .> 1e-6) .& (mars.Erosion_rate_mm_yr .>0)
scatternicely!(h, mars.Area_km2[t], mars.Erosion_rate_mm_yr[t],
    label="Glacial",
    color=mineralcolors["rhodochrosite"],
    alpha=0.85,
    mswidth=0.25,
)

# savefig(h, "Mars_area_vs_erosion_rate.pdf")
display(h)

## -- Timescale vs erosion rate, Mars

h = plot(framestyle=:box,
    xlabel="Timescale [yr]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    fg_color_legend=:white,
    xlims = (10^-0.5, 10^10.5),
    xticks = 10.0.^(-0:10),
    ylims = (10^-8, 10^3),
    yticks = 10.0.^(-8:3),
    size = (400,400),
    legend = :bottomleft,
)

t = (ngmars.Time_interval_yr .> 0) .& (ngmars.Erosion_rate_mm_yr .>0)
scatternicely!(h, ngmars.Time_interval_yr[t], ngmars.Erosion_rate_mm_yr[t],
    label="Nonglacial (Mars)",
    color=parse(Color, "#f9dfe9"),
    alpha=1,
    mswidth=0,
)

t = (mars.Time_interval_yr .> 0) .& (mars.Erosion_rate_mm_yr .>0)
scatternicely!(h, mars.Time_interval_yr[t], mars.Erosion_rate_mm_yr[t],
    label="Glacial (Mars)",
    color=mineralcolors["rhodochrosite"],
    alpha=0.85,
    mswidth=0.25,
)

x = collect(xlims(h))
plot!(h, x, (10^5)./x, color=:darkred, label="", linestyle=:dash)
xl = minimum(x)
annotate!(h, xl*1000, (10^5)/(xl*1000), text("100 m", 10, :bottom, :left, color=:darkred, rotation=-45.0))

savefig(h, "Mars_timescale_vs_erosion_rate.pdf")
display(h)


## --- Histogram of erosion rates... on Mars!

# Pick binning scheme
t = (ngmars.Time_interval_yr .> 0) .& (ngmars.Erosion_rate_mm_yr .>0)
μ = nanmean(log10.(ngmars.Erosion_rate_mm_yr[t]))
binedges = (-6:0.2:6) .+ μ
bincenters = cntr(binedges)
distx = first(binedges):0.01:last(binedges)

hist = plot(framestyle=:box,
    xlabel="Erosion Rate [mm/yr]",
    ylabel="N (nonglacial)",
    xflip=true,
    xlims=(10^-8,10^3),
    xticks=10.0.^(-8:3),
    xscale=:log10,
    size=(400,200),
    rotation=90.,
    xguidefontrotation=180.,
    yguidefontcolor=parse(Color, "#888888"),
    ytickfontcolor=parse(Color, "#888888"),
    fontfamily=:Helvetica,
)

# Nonglacial erosion
tng = (ngmars.Time_interval_yr .> 0) .& (ngmars.Erosion_rate_mm_yr .> 0)
logerosion = log10.(ngmars.Erosion_rate_mm_yr[tng])
Ns = histcounts(logerosion, binedges, T=Float64)
plot!(hist, 10.0.^stepifyedges(binedges), stepify(Ns),
    fill=true,
    color=parse(Color, "#d8cad4"),
    label=""
)

ngμ, ngσ = nanmean(logerosion), nanstd(logerosion)
plot!(hist, 10.0.^distx, normpdf(ngμ,ngσ,distx) .* (sum(Ns) * step(binedges)),
    linestyle=:solid,
    color=parse(Color, "#764a54"),
    label="",
)

# Glacial erosion
t = (mars.Time_interval_yr .> 0) .& (mars.Erosion_rate_mm_yr .>0)
logerosion = log10.(mars.Erosion_rate_mm_yr[t])
Ns = histcounts(logerosion, binedges)
hist2 = twinx()
plot!(hist2,
    framestyle=:box,
    ylabel="N (glacial)",
    xflip=true,
    xlims=(10^-8,10^3),
    xticks=10.0.^(-8:3),
    xscale=:log10,
    size=(400,200),
    rotation=90.,
)

plot!(hist2, 10.0.^stepifyedges(binedges), stepify(Ns), fill=true, color=mineralcolors["rhodochrosite"], label="")

μ, σ = nanmean(logerosion), nanstd(logerosion)
plot!(hist2, 10.0.^distx, normpdf(μ,σ,distx) .* (sum(Ns)*step(binedges)),
    linestyle=:solid,
    color=:darkred,
    label="",
)

# Clean up y axis limits
ylims!(hist, 0, 10)
ylims!(hist2, 0, 4)


# Glacial Gaussian
annotate!(hist2, [10.0.^μ], [0], text(" $(round(10^(μ+3), digits=2))e-3 mm/yr", 7, :left, :bottom, rotation=90, color=:darkred))
vline!(hist2, [10.0.^μ], color=:darkred, linestyle=:dot, label="")
annotate!(hist2, [10^μ], [maximum(ylims(hist2))], text("glacial ", 7, :right, :bottom, rotation=90, color=:darkred))

# Nonglacial Gaussian
annotate!(hist2, [10.0.^ngμ], [0], text(" $(round(10^(ngμ+3), digits=2))e-3 mm/yr", 7, :left, :bottom, rotation=90, color=parse(Color, "#764a54")))
vline!(hist2, [10.0.^ngμ], color=parse(Color, "#764a54"), linestyle=:dot, label="")
annotate!(hist2, [10.0.^ngμ], [maximum(ylims(hist2))], text("nonglacial ", 7, :right, :bottom, rotation=90, color=parse(Color, "#764a54")))


savefig(hist, "Mars_erosion_rate_histogram.pdf")
display(hist)


## --- End of File
