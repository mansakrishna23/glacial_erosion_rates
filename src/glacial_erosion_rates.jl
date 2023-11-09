using StatGeochem, Plots

datadir = "$(@__DIR__)/../data"
earth = importdataset("$datadir/glacial_erosion_Earth.tsv", '\t', importas=:Tuple);
mars = importdataset("$datadir/glacial_erosion_Mars.tsv", '\t', importas=:Tuple);

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

# h = plot(framestyle=:box,
#     xlabel="Area [km²]",
#     ylabel="Erosion rate [mm/yr]",
#     xscale=:log10,
#     yscale=:log10,
#     fontfamily=:Helvetica,
# )
#
# t = (earth.Area_km2 .> 0) .& (earth.Erosion_rate_mm_yr .>0)
#
# plot!(h, earth.Area_km2[t], earth.Erosion_rate_mm_yr[t],
#     seriestype=:scatter,
#     label="",
#     color=mineralcolors["corundum"],
#     alpha=0.85,
#     mswidth=0.25,
# )
#
# savefig(h, "erosion_rate_vs_area_all.pdf")

h = plot(framestyle=:box,
    xlabel="Area [km²]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    xlims = (10^-1, 10^7),
    xticks = 10.0.^(-1:7),
    ylims = (10^-3, 10^3),
    yticks = 10.0.^(-3:3),
)

t = (earth.Area_km2 .> 0) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Methodology .!== "Cosmogenic surface")

scatternicely!(h, earth.Area_km2[t], earth.Erosion_rate_mm_yr[t],
    label="",
    color=mineralcolors["corundum"],
    alpha=0.85,
    mswidth=0.25,
)

savefig(h, "area_vs_erosion_rate.pdf")
display(h)


## --- Latitude vs erosion rate, colored by glacier type (Continental, Alpine, etc.),

h = plot(framestyle=:box,
    xlabel="Latitude",
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

type = ( "Outlet ice stream", "Alpine tidewater", "Alpine", "High-latitude", "Continental", "Dry Valleys",)
tlabel = ( "Outlet ice stream", "Alpine tidewater", "Alpine", "High-latitude", "Continental", "Dry Valleys",)
colors = [mineralcolors[m] for m in ("muscovite", "quartz", "azurite", "apatite", "fluid", "rhodochrosite", )]
for i in eachindex(type)
    t = .!isnan.(earth.Latitude) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Type .== type[i])
    scatternicely!(h, abs.(earth.Latitude[t]), earth.Erosion_rate_mm_yr[t],
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = tlabel[i],
    )
end

savefig(h, "latitude_vs_erosion_rate-type.pdf")
display(h)

## --- Latitude vs erosion rate, colored by measurement method (Cosmogenic, Volumetric, etc.)

h = plot(framestyle=:box,
    xlabel="Latitude",
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

method = ("Cosmogenic", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
mlabel = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
colors = [mineralcolors[m] for m in ( "spessartine", "rhodochrosite", "corundum", "sodalite", "glaucophane")]
for i in eachindex(method)
    t = .!isnan.(earth.Latitude) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Methodology .== method[i])
    scatternicely!(h, abs.(earth.Latitude[t]), earth.Erosion_rate_mm_yr[t],
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = mlabel[i],
    )
end

savefig(h, "latitude_vs_erosion_rate-method.pdf")
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

type = ( "Outlet ice stream", "Alpine tidewater", "Alpine", "High-latitude", "Continental", "Dry Valleys",)
tlabel = ( "Outlet ice stream", "Alpine tidewater", "Alpine", "High-latitude", "Continental", "Dry Valleys",)
colors = [mineralcolors[m] for m in ("muscovite", "quartz", "azurite", "apatite", "fluid", "rhodochrosite", )]
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


## --- Timescale vs erosion rate, colored by measurement method (Cosmogenic, Volumetric, etc.)

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

method = ("Cosmogenic", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
mlabel = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
colors = [mineralcolors[m] for m in ( "spessartine", "rhodochrosite", "corundum", "sodalite", "glaucophane")]
for i in eachindex(method)
    t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Methodology .== method[i])
    plot!(h, earth.Time_interval_yr[t], earth.Erosion_rate_mm_yr[t],
        seriestype=:scatter,
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = mlabel[i],
    )
end

# t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0) .& (findmatches(earth.Methodology, method) .== 0)
# plot!(h, earth.Time_interval_yr[t], earth.Erosion_rate_mm_yr[t],
#     seriestype=:scatter,
#     color = mineralcolors["magnetite"],
#     alpha=0.85,
#     mswidth=0.25,
#     label = "All other methods",
# )

savefig(h, "timescale_vs_erosion_rate-method.pdf")
display(h)

## --- Histogram of erosion rates!

t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0)
μ = nanmean(log10.(earth.Erosion_rate_mm_yr[t]))
binedges = (-4:0.2:4) .+ μ
bincenters = cntr(binedges)
hist = plot(framestyle=:box,
    xlabel="Erosion Rate [mm/yr]",
    ylabel="N",
    xflip=true,
    fg_color_legend=:white,
    xlims=(10^-5,10^3),
    xticks=10.0.^(-5:3),
    xscale=:log10,
)
Ns = histcounts(log10.(earth.Erosion_rate_mm_yr[t]), binedges)
plot!(10.0.^stepifyedges(binedges), stepify(Ns), fill=true, color=mineralcolors["rhodochrosite"], label="")

t .&= earth.Type .!= "Dry Valleys"
logerosion = log10.(earth.Erosion_rate_mm_yr[t])
Ns = histcounts(logerosion, binedges)
plot!(10.0.^stepifyedges(binedges), stepify(Ns), fill=true, color=mineralcolors["fluid"], label="")

x = first(binedges):0.01:last(binedges)
μ, σ = nanmean(logerosion), nanstd(logerosion)
plot!(10.0.^x, normpdf(μ,σ,x).*(sum(Ns)*step(binedges)),
    linestyle=:dot,
    color=:black,
    label="",
    )

vline!([10.0.^μ], color=:black, label="")
annotate!([10.0.^μ], [0], text(" $(round(10^μ, digits=2)) mm/yr", 10, :left, :bottom, rotation=90))

annotate!([10^-4.9], [maximum(ylims())], text("Dry Valleys, East Antarctica ", 10, :right, :bottom, rotation=90, color=mineralcolors["rhodochrosite"]))
annotate!([10^2.9], [maximum(ylims())], text("All other rates ", 10, :right, :top, rotation=90, color=mineralcolors["fluid"]))

ylims!(0, maximum(ylims()))

## --- Timescale vs erosion rate, binned by area and colored by measurement method

h = plot(layout = (2,2),
    framestyle=:box,
    xlabel="Timescale [yr]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    size=(800,800),
    xlims = (10^-2, 10^8),
    xticks = 10.0.^(-2:1:8),
    ylims = (10^-5, 10^5),
    yticks = 10.0.^(-5:1:5),
)

areas = (0, 0.1, 10, 1000, 10^7)
titles = ("0-0.1 km²", "0.1-10 km²", "10-1000 km²", ">1000 km²")

# Add cosmogenic line of constant 600mm thickness
plot!(h[1], collect(xlims(h[1])), 600.0./collect(xlims(h[1])), color=:red, label="")
x = minimum(xlims(h[1]))*90
annotate!(h[1], x, 600/x, text("600 mm", 10, :bottom, color=:red, rotation=-45.0))

# Add log-mean line for all others
t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Type .!= "Dry Valleys")
logμ = 10^nanmean(log10.(earth.Erosion_rate_mm_yr[t]))
for j in 1:4
    hline!(h[j], [logμ], color=mineralcolors["fluid"], label="")
end
annotate!(h[1], x, logμ, text("$(round(logμ, digits=2)) mm/yr", 10, :top, color=mineralcolors["fluid"],))


method = ("Cosmogenic", "Cosmogenic surface", "Thermochronometric", "Volumetric",)
mlabel = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
colors = [mineralcolors[m] for m in ( "spessartine", "rhodochrosite", "corundum", "sodalite", "glaucophane")]

for j in 1:length(areas)-1
    t = (areas[j] .< earth.Area_km2 .< areas[j+1]) .& (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0)
    plot!(h[j], title=titles[j])

    for i in eachindex(method)
        tt = t .& (earth.Methodology .== method[i])
        plot!(h[j], earth.Time_interval_yr[tt], earth.Erosion_rate_mm_yr[tt],
            seriestype=:scatter,
            color = colors[i],
            alpha = 0.85,
            mswidth = 0.1,
            label = mlabel[i],
        )
    end
end

savefig(h, "timescale_vs_erosion_rate-area.pdf")
display(h)


## --- Continental icesheet erosion

tc = (earth.Type .== "Continental") .| (earth.Type .== "High-latitude") # .| (earth.Type .== "Dry Valleys")  .| (earth.Type .== "Outlet Ice Stream")
h = plot(framestyle=:box,
    xlabel="Timescale [yr]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    xlims=(10^-1, 10^9),
    xticks=10.0.^(-1:9),
    ylims=(10^-2.1, 10^0.6),
    yticks=10.0.^(-2:0.5:0.5),
    fg_color_legend=:white,
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
    annotate!(0.9e9, mean_erosion*1.3, text("$region,\nunweighted ave.", colors[i], :right, 8, fontfamily=:Helvetica))
    # # (optional) error bars on the means
    # plot!(h, [5e8], [mean_erosion], yerr=[erosion_sigma], seriestype=:scatter, msc=colors[i], color=colors[i], label="")
end

# Optional: Plot 3-5 km cryogenian erosion, for comparison
cryocolor=mineralcolors["magnetite"]
plot!(h, [64e6], [4000*1000/64e6], yerr=[1000*1000/64e6], msc=cryocolor, color=cryocolor, seriestype=:scatter, label="")
annotate!(h, 2.5e8, 4000*1000/64e6, text("3-5 km\nCryogen.\nerosion", cryocolor, :center, :8))

# savefig(h, "continental_erosion_rates.pdf")

# Add Mars!
marscolor = mineralcolors["almandine"] # parse(Colorant, "#ff0000")
scatternicely!(h, mars.Time_interval_yr, mars.Erosion_rate_mm_yr;
    color=marscolor,
    label="Mars",
    alpha=0.85,
    mswidth=0.25,
    legend=:bottomleft,
)

# Draw lines for unweighted means of each ice sheet
mean_erosion = nanmean(mars.Erosion_rate_mm_yr)
erosion_sigma = nanstd(earth.Erosion_rate_mm_yr)./sqrt(length(mars.Erosion_rate_mm_yr))
hline!(h, [mean_erosion],
    linestyle=:dot,
    color=marscolor,
    label="",
)
annotate!(0.9e9, mean_erosion/1.3, text("Mars,\nunweighted ave.", marscolor, :right, 8, fontfamily=:Helvetica))

savefig(h, "continental_erosion_rates_withmars.pdf")

display(h)

## ---
