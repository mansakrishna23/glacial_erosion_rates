using StatGeochem, Plots
cd(@__DIR__)
earth = importdataset("glacial_erosion_Earth.tsv", '\t', importas=:Tuple);
mars = importdataset("glacial_erosion_Mars.tsv", '\t', importas=:Tuple);


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
plotmethod = ("Cosmogenic surface", "Cosmogenic detrital", "Thermochronometric", "Volumetric",  "Relief")

colors = [mineralcolors[m] for m in ("fluid", "zircon", "kyanite", "azurite", "glaucophane")]
for i in eachindex(method)
    t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Methodology .== method[i])
    plot!(h, earth.Time_interval_yr[t], earth.Erosion_rate_mm_yr[t],
        seriestype=:scatter,
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = plotmethod[i],
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

h = plot(framestyle=:box,
    xlabel="Timescale [yr]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    xlims = (10^-2, 10^9),
    xticks = 10.0.^(-2:9),
    ylims = (10^-3, 10^3),
    yticks = 10.0.^(-3:3),
)

t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Methodology .!= "Cosmogenic surface")
plot!(h, earth.Time_interval_yr[t], earth.Erosion_rate_mm_yr[t],
    seriestype = :scatter,
    color = mineralcolors["azurite"],
    alpha = 0.85,
    mswidth = 0.25,
    label="",
)

savefig(h, "erosion_rate_vs_timescale_noncosmogenic.pdf")
display(h)


## -- Area vs erosion rate

h = plot(framestyle=:box,
    xlabel="Area [km2]",
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

## --- Continental icesheet erosion

minarea = 0

tc = (earth.Type .== "Continental") .| (earth.Type .== "High-latitude")
h = plot(framestyle=:box,
    xlabel="Timescale [yr]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    xlims=(10^-1, 10^9),
    ylims=(10^-2, 10^1),
    xticks=10.0.^(-1:9),
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
    plot!(h, earth.Time_interval_yr[t], earth.Erosion_rate_mm_yr[t],
        seriestype=:scatter,
        # msc=colors[i],
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

# Add 3-5 km cryogenian erosion
cryocolor=mineralcolors["magnetite"]
plot!([64e6], [4000*1000/64e6], yerr=[1000*1000/64e6], msc=cryocolor, color=cryocolor, seriestype=:scatter, label="")
annotate!(2.5e8, 4000*1000/64e6, text("3-5 km\nCryogen.\nerosion", cryocolor, :center, :8))

savefig("Continental_glaciation_rates.pdf")

# Add Mars!
marscolor = mineralcolors["proustite"]
plot!(h, mars.Time_interval_yr, mars.Erosion_rate_mm_yr,
    seriestype=:scatter,
    color=marscolor,
    label="Mars",
    alpha=0.85,
    mswidth=0.25,
    legend=:bottomleft,
    xlims=(10^-1, 10^9),
    ylims=(10^-2.5, 10^0.5),
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


savefig("Continental_glaciation_rates_withmars.pdf")

ylims!(10^-4, 10^1)
savefig("Continental_glaciation_rates_withmars_all.pdf")

display(h)

## ---
