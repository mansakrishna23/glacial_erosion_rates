using StatGeochem, Plots

datadir = "$(@__DIR__)/../data"
earth = importdataset("$datadir/glacial_erosion_Earth.tsv", '\t', importas=:Tuple)
mars = importdataset("$datadir/glacial_erosion_Mars.tsv", '\t', importas=:Tuple)

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

## --- Slope vs erosion rate, colored by glacier type (Continental, Alpine, etc.),

h = plot(framestyle=:box,
    xlabel="Regional slope [m/m]",
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

t = .!isnan.(ngearth.Slope_m_km) .& (ngearth.Erosion_rate_mm_yr .> 0)
slope = ngearth.Slope_m_km[t]./1000
scatternicely!(h, slope, ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)
aₙ,bₙ = linreg(slope, log10.(ngearth.Erosion_rate_mm_yr[t]))
xₙ = collect(xlims(h))
yₙ = @. 10^(aₙ+bₙ*xₙ)
plot!(h, xₙ, yₙ, color=parse(Color, "#777777"), linestyle=:dot, label="$(round(10^aₙ,digits=3))*10^($(round(bₙ,digits=3)) s)")

type = ( "Alpine", "Alpine tidewater", "High-latitude", "Continental", "Outlet ice stream", )
tlabel = ( "Alpine", "Alpine tidewater", "Non-cont. high-lat.", "Continental", "Outlet ice stream", )
colors = [mineralcolors[m] for m in ("malachite", "zircon", "azurite", "quartz", "fluid", )]
for i in eachindex(type)
    t = .!isnan.(earth.Slope_m_km) .& (earth.Erosion_rate_mm_yr .> 0) .& (earth.Type .== type[i])
    scatternicely!(h, earth.Slope_m_km[t]./1000, earth.Erosion_rate_mm_yr[t],
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = tlabel[i],
    )
end

t = .!isnan.(earth.Slope_m_km) .& (earth.Erosion_rate_mm_yr .> 0) 
slope = earth.Slope_m_km[t]./1000
a,b = linreg(slope, log10.(earth.Erosion_rate_mm_yr[t]))
x = collect(xlims(h))
y = @. 10^(a+b*x)

plot!(h, x, y, color=:darkblue, linestyle=:dot, label="$(round(10^a,digits=3))*10^($(round(b,digits=3)) s)")

savefig(h, "slope_vs_erosion_rate-type.pdf")
display(h)

## --- Latitude vs erosion rate, colored by glacier type (Continental, Alpine, etc.),

h = plot(framestyle=:box,
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
scatternicely!(h, abs.(ngearth.Latitude[t]), ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

type = ( "Alpine", "Alpine tidewater", "High-latitude", "Continental", "Outlet ice stream",)
tlabel = ( "Alpine", "Alpine tidewater", "Non-cont. high-lat.", "Continental", "Outlet ice stream", )
colors = [mineralcolors[m] for m in ("malachite", "zircon", "azurite", "quartz", "fluid", )]
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
scatternicely!(h, abs.(ngearth.Latitude[t]), ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

method = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
mlabel = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
colors = [mineralcolors[m] for m in ( "spessartine", "rhodochrosite", "corundum", "sodalite", "sapphirine")]
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

t = (ngearth.Time_interval_yr .> 0) .& (ngearth.Erosion_rate_mm_yr .>0)
scatternicely!(h, ngearth.Time_interval_yr[t], ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

method = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
mlabel = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
colors = [mineralcolors[m] for m in ( "spessartine", "rhodochrosite", "corundum", "sodalite", "sapphirine")]
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
t = (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0) 
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
minsamples = 50 # Minimum number of samples for a method to plot an error bar

# Add cosmogenic line of constant 600mm thickness
x = collect(xlims(hng[1]))
for j in 1:length(areas)-1
    plot!(hng[j], x, 600.0./x, color=:red, label="", linestyle=:dash)
end
xl = minimum(x)
annotate!(hng[1], xl*10, 600/(xl*10), text("600 mm", 10, :bottom, :left, color=:red, rotation=-45.0))

# Add log-mean line for all others
t = (ngearth.Time_interval_yr .> 0) .& (ngearth.Erosion_rate_mm_yr .>0) .& (ngearth.Type .== "Fluvial")
logμ = 10^nanmean(log10.(ngearth.Erosion_rate_mm_yr[t]))
for j in 1:length(areas)-1
    hline!(hng[j], [logμ], color=parse(Color, "#666666"), label="", linestyle=:dot)
end
annotate!(hng[1], xl*3, logμ, text("$(round(logμ, sigdigits=2)) mm/yr", 10, :top, :left, color=parse(Color, "#666666"),))


method = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric", "Relief", )
mlabel = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric", "Relief", )
colors = [mineralcolors[m] for m in ( "spessartine", "rhodochrosite", "corundum", "sodalite", "sapphirine", )]

for j in 1:length(areas)-1
    t = (areas[j] .< ngearth.Area_km2 .< areas[j+1]) .& (ngearth.Time_interval_yr .> 0) .& (ngearth.Erosion_rate_mm_yr .>0) .& (ngearth.Type .== "Fluvial")
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
    t = (areas[j] .< ngearth.Area_km2 .< areas[j+1]) .& (ngearth.Time_interval_yr .> 0) .& (ngearth.Erosion_rate_mm_yr .>0) .& (ngearth.Type .== "Fluvial")
    plot!(hng[j], title=titles[j])

    for i in eachindex(method)
        tt = t .& (ngearth.Methodology .== method[i])
        if count(tt) .> minsamples
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

savefig(hng, "timescale_vs_erosion_rate-area-fluvial.pdf")
display(hng)


## --- Combined timescale vs erosion rate area

h = plot(hg, hng,
    size = (1200,600),
    layout = (2,1),
)
savefig(h, "timescale_vs_erosion_rate-area.pdf")
display(h)


## --- Plot all data by method
method = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric", "Relief", )

hm = plot(layout = (1,length(method)),
    framestyle=:box,
    xlabel="Timescale [yr]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    size=(300*length(method),300),
    xlims = (10^-3, 10^8),
    xticks = 10.0.^(-3:1:8),
    ylims = (10^-5.5, 10^5.5),
    yticks = 10.0.^(-5:1:5),
    legend=:none,
)

# Add cosmogenic line of constant 600mm thickness
x = collect(xlims(hm[1]))
for j in 1:2
    plot!(hm[j], x, 600.0./x, color=:red, label="", linestyle=:dash)
    xl = minimum(x)
    annotate!(hm[j], xl*10, 600/(xl*10), text("600 mm", 10, :bottom, :left, color=:red, rotation=-45.0))
end

plot!(hm[3], x, 2e6./x, color=:blue, label="", linestyle=:dash)
plot!(hm[5], x, 2e6./x, color=:blue, label="", linestyle=:dash)


for j in eachindex(method)
    t = (ngearth.Methodology .== method[j])

    tt = t .& (ngearth.Type .== "Fluvial")
    plot!(hm[j], ngearth.Time_interval_yr[tt], ngearth.Erosion_rate_mm_yr[tt],
        seriestype=:scatter,
        color = parse(Color, "#dddddd"),
        mswidth = 0,
        label = "Fluvial",
    )

    tt = t .& (ngearth.Type .== "Subaerial")
    plot!(hm[j], ngearth.Time_interval_yr[tt], ngearth.Erosion_rate_mm_yr[tt],
        seriestype=:scatter,
        color = parse(Color, "#aaaaaa"),
        mswidth = 0,
        label = "Subaerial",
    )

    t = (earth.Methodology .== method[j])
    plot!(hm[j], earth.Time_interval_yr[t], earth.Erosion_rate_mm_yr[t],
        seriestype=:scatter,
        color = :darkblue,
        alpha = 0.85,
        mswidth = 0,
        label = "Glacial",
        title = method[j],
    )
end

display(hm)

## -- Area-weighted averages

t = earth.Area_km2 .> 0
area_weighted_glacial = nanmean(earth.Erosion_rate_mm_yr[t], earth.Area_km2[t])

## ---
t = ngearth.Area_km2 .> 0
area_weighted_nonglacial = nanmean(ngearth.Erosion_rate_mm_yr[t], ngearth.Area_km2[t])
