using StatGeochem, Plots

datadir = "$(@__DIR__)/../data"
earth = importdataset("$datadir/glacial_erosion_Earth.tsv", '\t', importas=:Tuple);
mars = importdataset("$datadir/glacial_erosion_Mars.tsv", '\t', importas=:Tuple);

## --- Utility functions for plotting

# Scatter, but plot any points outside of ylims as up/down-ward pointing triangles
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


## --- All data by method, log erosion rate vs log time

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

method = ("Cosmogenic surface", "Cosmogenic", "Thermochronometric", "Volumetric",  "Relief")
mlabel = ("Cosmogenic surface", "Cosmogenic detrital", "Thermochronometric", "Volumetric",  "Relief")

colors = [mineralcolors[m] for m in ("fluid", "zircon", "kyanite", "azurite", "glaucophane")]
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

savefig(h, "erosion_rate_vs_timescale.pdf")
display(h)

## --- Everything but cosmogenic surface, log erosion rate vs log time
#
# h = plot(framestyle=:box,
#     xlabel="Timescale [yr]",
#     ylabel="Erosion rate [mm/yr]",
#     xscale=:log10,
#     yscale=:log10,
#     fontfamily=:Helvetica,
#     xlims = (10^-2, 10^9),
#     xticks = 10.0.^(-2:9),
#     ylims = (10^-3, 10^3),
#     yticks = 10.0.^(-3:3),
# )
#
# t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Methodology .!= "Cosmogenic surface")
# plot!(h, earth.Time_interval_yr[t], earth.Erosion_rate_mm_yr[t],
#     seriestype = :scatter,
#     color = mineralcolors["azurite"],
#     alpha = 0.85,
#     mswidth = 0.25,
#     label="",
# )
#
# savefig(h, "erosion_rate_vs_timescale_noncosmogenic.pdf")
# display(h)

## --- Everything but Dry Valleys, log erosion rate vs log time
#
# h = plot(framestyle=:box,
#     xlabel="Timescale [yr]",
#     ylabel="Erosion rate [mm/yr]",
#     xscale=:log10,
#     yscale=:log10,
#     fontfamily=:Helvetica,
#     xlims = (10^-2, 10^9),
#     xticks = 10.0.^(-2:9),
#     ylims = (10^-3, 10^3),
#     yticks = 10.0.^(-3:3),
# )
#
# t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0)
# t .&=  earth.Locality .!= "Dry Valleys, East Antarctica"
# plot!(h, earth.Time_interval_yr[t], earth.Erosion_rate_mm_yr[t],
#     seriestype = :scatter,
#     color = mineralcolors["azurite"],
#     alpha = 0.85,
#     mswidth = 0.25,
#     label="",
# )
#
# savefig(h, "erosion_rate_vs_timescale_nondryvalleys.pdf")
# display(h)

## Calculate some means

t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Area_km2 .>0)
t .&=  earth.Locality .!= "Dry Valleys, East Antarctica"

μ = nanmean(earth.Erosion_rate_mm_yr[t], earth.Area_km2[t])

logμ = nanmean(log10.(earth.Erosion_rate_mm_yr[t]))
logσ = nanstd(log10.(earth.Erosion_rate_mm_yr[t]))


## --- All data by glacier type (Continental, Alpine, etc.), log erosion rate vs log time

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

type = ("Alpine", "High-latitude", "Continental", "Outlet ice stream", "Dry Valleys",)

colors = [mineralcolors[m] for m in ("azurite", "zircon", "kyanite", "quartz", "fluid", )]
for i in eachindex(type)
    t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Type .== type[i])
    plot!(h, earth.Time_interval_yr[t], earth.Erosion_rate_mm_yr[t],
        seriestype=:scatter,
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = type[i],
    )
end

savefig(h, "glacier_type_erosion_rate_vs_timescale.pdf")
display(h)


## -- Area vs erosion rate


h = plot(framestyle=:box,
    xlabel="Area [km²]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
)

t = (earth.Area_km2 .> 0) .& (earth.Erosion_rate_mm_yr .>0)

plot!(h, earth.Area_km2[t], earth.Erosion_rate_mm_yr[t],
    seriestype=:scatter,
    label="",
    color=mineralcolors["corundum"],
    alpha=0.85,
    mswidth=0.25,
)

savefig(h, "erosion_rate_vs_area_all.pdf")


h = plot(framestyle=:box,
    xlabel="Area [km²]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    xlims = (10^-1, 10^7),
    xticks = 10.0.^(-1:7),
)

t = (earth.Area_km2 .> 0) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Methodology .!== "Cosmogenic surface")

plot!(h, earth.Area_km2[t], earth.Erosion_rate_mm_yr[t],
    seriestype=:scatter,
    label="",
    color=mineralcolors["corundum"],
    alpha=0.85,
    mswidth=0.25,
)

savefig(h, "erosion_rate_vs_area.pdf")
display(h)

## --- Erosion rate vs duration, binned by area

areas = (0, 0.1, 10, 1000, 10^7)
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

for i in 1:length(areas)-1
    t = areas[i] .< earth.Area_km2 .< areas[i+1]
    plot!(h[i], earth.Time_interval_yr[t], earth.Erosion_rate_mm_yr[t],
        seriestype=:scatter,
        title = "$(areas[i]) - $(areas[i+1]) km²",
        label="",
    )
end
savefig(h, "Sadlerianness_vs_area.pdf")
display(h)

## --- Sadler Slope vs area

# N.B. linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y
rsquared(x,y,line_y) = 1 - nansum((y .- line_y).^2)  / nansum((y .- nanmean(y)).^2) # Function for rsquared value

# areas = 10.0.^(-8:2:7)
areas = [0, 0.1, 10, 1000, 10^7]
labels = ["0-0.1", "0.1-10", "10-1000", ">1000"]

Ns = fill(0, length(areas)-1)
slopes = fill(NaN, length(areas)-1)
rsquareds = fill(NaN, length(areas)-1)


for i in eachindex(slopes)
    t = areas[i] .<= earth.Area_km2 .< areas[i+1]
    if any(t)
        Ns[i] = count(t)

        x = log10.(earth.Time_interval_yr[t])
        y = log10.(earth.Erosion_rate_mm_yr[t])

        a, b = linreg(x, y)
        slopes[i] = b
        ŷ = a .+ b.*x

        rsquareds[i] = rsquared(x,y,ŷ)
    end
end

h = plot(framestyle=:box,
    xlabel="Area [km²]",
    ylabel="Sadler Slope",
    xticks=(1:length(areas)-1,labels),
    ylims=(-1, 0.1)
)
plot!(h, 1:length(areas)-1, slopes,
    label="",
    markersize=5,
    seriestype=:scatter,
    color=:black,
)
plot!(h, 1:length(areas)-1, slopes,
    markersize=sqrt.(rsquareds.*Ns),
    label="",
    seriestype=:scatter,
    color=mineralcolors["crocoite"],
)
hline!(h, [0], color=:black, label="")

savefig(h, "Sadler_slope_vs_area.pdf")
display(h)


## --- Continental icesheet erosion

minarea = 0

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

colors = (mineralcolors["fluid"], mineralcolors["chrome diopside"], mineralcolors["azurite"])
for i in eachindex(regions)
    region = regions[i]
    t = (earth.Region .== region) .& tc
    if minarea>0
        t .&= earth.Area_km2 .> minarea
    end
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

# Optional: Add 3-5 km cryogenian erosion, for comparison
cryocolor=mineralcolors["magnetite"]
plot!(h, [64e6], [4000*1000/64e6], yerr=[1000*1000/64e6], msc=cryocolor, color=cryocolor, seriestype=:scatter, label="")
annotate!(h, 2.5e8, 4000*1000/64e6, text("3-5 km\nCryogen.\nerosion", cryocolor, :center, :8))

savefig(h, "Continental_glaciation_rates.pdf")

# Add Mars!
marscolor = mineralcolors["proustite"]
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

savefig(h, "Continental_glaciation_rates_withmars.pdf")

display(h)

## ---
