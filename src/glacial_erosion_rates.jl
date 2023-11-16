using StatGeochem, Plots

datadir = "$(@__DIR__)/../data"
earth = importdataset("$datadir/glacial_erosion_Earth.tsv", '\t', importas=:Tuple);
mars = importdataset("$datadir/glacial_erosion_Mars.tsv", '\t', importas=:Tuple);

nonglacial = importdataset("$datadir/nonglacial_erosion_Earth.tsv", '\t', importas=:Tuple)


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
    xlims = (10^-2.5, 10^7.5),
    xticks = 10.0.^(-2:7),
    ylims = (10^-4, 10^3),
    yticks = 10.0.^(-4:3),
    fg_color_legend = :white,
)

t = (nonglacial.Area_km2 .> 1e-6) .& (nonglacial.Erosion_rate_mm_yr .>0)
scatternicely!(h, nonglacial.Area_km2[t], nonglacial.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

t = (earth.Area_km2 .> 1e-6) .& (earth.Erosion_rate_mm_yr .>0)
scatternicely!(h, earth.Area_km2[t], earth.Erosion_rate_mm_yr[t],
    label="Glacial",
    color=mineralcolors["azurite"],
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

t = .!isnan.(nonglacial.Latitude) .& (nonglacial.Erosion_rate_mm_yr .>0)
scatternicely!(h, abs.(nonglacial.Latitude[t]), nonglacial.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

type = ( "Alpine", "Alpine tidewater", "High-latitude", "Continental", "Outlet ice stream", "Dry Valleys",)
tlabel = ( "Alpine", "Alpine tidewater", "Non-cont. high-lat.", "Continental", "Outlet ice stream", "Dry Valleys",)
colors = [mineralcolors[m] for m in ("malachite", "zircon", "quartz", "azurite", "fluid", "rhodochrosite", )]
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

t = .!isnan.(nonglacial.Latitude) .& (nonglacial.Erosion_rate_mm_yr .>0)
scatternicely!(h, abs.(nonglacial.Latitude[t]), nonglacial.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

method = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
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


## --- Slope vs erosion rate

h = plot(framestyle=:box,
    xlabel="Slope [m/m]",
    ylabel="Erosion rate [mm/yr]",
    yscale=:log10,
    fontfamily=:Helvetica,
    fg_color_legend=:white,
    legend=:bottomright,
    xlims=(0,0.6),
    # xticks=0:15:90,
    ylims = (10^-5, 10^3),
    yticks = 10.0.^(-5:3),

)

t = .!isnan.(nonglacial.Slope_m_km) .& (nonglacial.Erosion_rate_mm_yr .> 0)
slope = nonglacial.Slope_m_km[t]./1000
scatternicely!(h, slope, nonglacial.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)
aₙ,bₙ = linreg(slope, log10.(nonglacial.Erosion_rate_mm_yr[t]))
xₙ = collect(xlims(h))
yₙ = @. 10^(aₙ+bₙ*xₙ)
plot!(h, xₙ, yₙ, color=parse(Color, "#777777"), linestyle=:dash, label="$(round(10^aₙ,digits=3))*10^($(round(bₙ,digits=3)) s)")


type = ( "Alpine", "Alpine tidewater", "High-latitude", "Continental", "Outlet ice stream", "Dry Valleys",)
tlabel = ( "Alpine", "Alpine tidewater", "Non-cont. high-lat.", "Continental", "Outlet ice stream", "Dry Valleys",)
colors = [mineralcolors[m] for m in ("malachite", "zircon", "quartz", "azurite", "fluid", "rhodochrosite", )]
for i in eachindex(type)
    t = .!isnan.(earth.Slope_m_km) .& (earth.Erosion_rate_mm_yr .> 0) .& (earth.Type .== type[i])
    scatternicely!(h, earth.Slope_m_km[t]./1000, earth.Erosion_rate_mm_yr[t],
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = tlabel[i],
    )
end

t = .!isnan.(earth.Slope_m_km) .& (earth.Erosion_rate_mm_yr .> 0) .& (earth.Type .!= "Dry Valleys")
slope = earth.Slope_m_km[t]./1000
a,b = linreg(slope, log10.(earth.Erosion_rate_mm_yr[t]))
x = collect(xlims(h))
y = @. 10^(a+b*x)

plot!(h, x, y, color=:darkblue, linestyle=:dash, label="$(round(10^a,digits=3))*10^($(round(b,digits=3)) s)")

savefig(h, "slope_vs_erosion_rate-type.pdf")
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

t = (nonglacial.Time_interval_yr .> 0) .& (nonglacial.Erosion_rate_mm_yr .>0)
scatternicely!(h, nonglacial.Time_interval_yr[t], nonglacial.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

type = ( "Alpine", "Alpine tidewater", "High-latitude", "Continental", "Outlet ice stream", "Dry Valleys",)
tlabel = ( "Alpine", "Alpine tidewater", "Non-cont. high-lat.", "Continental", "Outlet ice stream", "Dry Valleys",)
colors = [mineralcolors[m] for m in ("malachite", "zircon", "quartz", "azurite", "fluid", "rhodochrosite", )]
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

t = (nonglacial.Time_interval_yr .> 0) .& (nonglacial.Erosion_rate_mm_yr .>0)
scatternicely!(h, nonglacial.Time_interval_yr[t], nonglacial.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

method = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
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
tng = (nonglacial.Time_interval_yr .> 0) .& (nonglacial.Erosion_rate_mm_yr .> 0)
logerosion = log10.(nonglacial.Erosion_rate_mm_yr[tng])
Ns = histcounts(logerosion, binedges, T=Float64)
plot!(hist, 10.0.^stepifyedges(binedges), stepify(Ns),
    fill=true,
    color=parse(Color, "#dddddd"),
    label=""
)

ngμ, ngσ = nanmean(logerosion), nanstd(logerosion)
plot!(hist, 10.0.^distx, normpdf(ngμ,ngσ,distx) .* (sum(Ns) * step(binedges)),
    linestyle=:dot,
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
    size=(600,300),
    rotation=90.,
)

plot!(hist2, 10.0.^stepifyedges(binedges), stepify(Ns), fill=true, color=mineralcolors["rhodochrosite"], label="")

t .&= earth.Type .!= "Dry Valleys"
logerosion = log10.(earth.Erosion_rate_mm_yr[t])
Ns = histcounts(logerosion, binedges, T=Float64)
plot!(hist2, 10.0.^stepifyedges(binedges), stepify(Ns), fill=true, color=mineralcolors["fluid"], label="")

μ, σ = nanmean(logerosion), nanstd(logerosion)
plot!(hist2, 10.0.^distx, normpdf(μ,σ,distx) .* (sum(Ns)*step(binedges)),
    linestyle=:dot,
    color=:darkblue,
    label="",
)

# Clean up y axis limits
ylims!(hist, 0, 800)
ylims!(hist2, 0, 70)

# Dry valleys
annotate!(hist2, [10^-4.3], [maximum(ylims(hist2))], text("Dry Valleys, East Antarctica ", 9, :right, :bottom, rotation=90, color=mineralcolors["rhodochrosite"]))

# Glacial gaussian
annotate!(hist2, [10.0.^μ], [0], text(" $(round(10^μ, digits=2)) mm/yr", 9, :left, :bottom, rotation=90, color=:darkblue))
vline!(hist2, [10.0.^μ], color=:darkblue, linestyle=:dash, label="")
annotate!(hist2, [10^μ], [maximum(ylims(hist2))], text("glacial mean ", 9, :right, :bottom, rotation=90, color=:darkblue))

# Nonglacial gaussian
annotate!(hist2, [10.0.^ngμ], [0], text(" $(round(10^ngμ, digits=2)) mm/yr", 9, :left, :bottom, rotation=90, color=parse(Color, "#888888")))
vline!(hist2, [10.0.^ngμ], color=parse(Color, "#888888"), linestyle=:dash, label="")
annotate!(hist2, [10.0.^ngμ], [maximum(ylims(hist2))], text("nonglacial mean ", 9, :right, :bottom, rotation=90, color=parse(Color, "#888888")))


savefig(hist, "erosion_rate_histogram.pdf")
display(hist)


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


method = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",)
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
    legend=:bottomleft,
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
# annotate!(0.9e9, mean_erosion/1.3, text("Mars,\nunweighted ave.", marscolor, :right, 8, fontfamily=:Helvetica))
# savefig(h, "continental_erosion_rates_withmars.pdf")

display(h)


## --- Histogram of continental erosion rates

hist = plot(framestyle=:box,
    xlabel="Erosion Rate [mm/yr]",
    ylabel="N",
    xflip=true,
    xlims=(10^-5,10^3),
    xticks=10.0.^(-5:3),
    xscale=:log10,
    size=(600,300),
    rotation=90,
    xguidefontrotation=180.,
    fontfamily=:Helvetica,
)

tc = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0)
tc .&= (earth.Type .== "Continental") .| (earth.Type .== "High-latitude")
tc .&= (earth.Region .== "Greenland") .| (earth.Region .== "Fennoscandia") .| (earth.Region .== "Antarctica")

logerosion = log10.(earth.Erosion_rate_mm_yr[tc])
μ, σ = nanmean(logerosion), nanstd(logerosion)
binedges = (-4:0.2:4) .+ μ
bincenters = cntr(binedges)

# Ns = histcounts(log10.(earth.Erosion_rate_mm_yr[t]), binedges)
# plot!(hist, 10.0.^stepifyedges(binedges), stepify(Ns), fill=true, color=mineralcolors["rhodochrosite"], label="")

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
    linestyle=:dot,
    color=:black,
    label="",
)
vline!(hist, [10.0.^μ], color=:black, linestyle=:dash, label="")
annotate!(hist, [10.0.^μ], [0], text(" $(round(10^μ, digits=2)) mm/yr", 10, :left, :bottom, rotation=90))

ylims!(hist, 0, 25)
savefig(hist, "continental_erosion_rate_histogram.pdf")
display(hist)

## ---
