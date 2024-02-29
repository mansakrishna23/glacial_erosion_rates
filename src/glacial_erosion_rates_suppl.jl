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

hst = plot(framestyle=:box,
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
scatternicely!(hst, slope, ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)
aₙ,bₙ = linreg(slope, log10.(ngearth.Erosion_rate_mm_yr[t]))
xₙ = collect(xlims(hst))
yₙ = @. 10^(aₙ+bₙ*xₙ)
plot!(hst, xₙ, yₙ, color=parse(Color, "#777777"), linestyle=:dot, label="$(round(10^aₙ,digits=3))*10^($(round(bₙ,digits=3)) s)")

type = ( "Alpine", "Alpine tidewater", "High-latitude", "Continental", "Outlet ice stream", )
tlabel = ( "Alpine", "Alpine tidewater", "Non-cont. high-lat.", "Continental", "Outlet ice stream", )
colors = [mineralcolors[m] for m in ("zircon", "malachite", "azurite", "quartz", "fluid", )]
for i in eachindex(type)
    t = .!isnan.(earth.Slope_m_km) .& (earth.Erosion_rate_mm_yr .> 0) .& (earth.Type .== type[i])
    scatternicely!(hst, earth.Slope_m_km[t]./1000, earth.Erosion_rate_mm_yr[t],
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = tlabel[i],
    )
end

t = .!isnan.(earth.Slope_m_km) .& (earth.Erosion_rate_mm_yr .> 0) 
slope = earth.Slope_m_km[t]./1000
a,b = linreg(slope, log10.(earth.Erosion_rate_mm_yr[t]))
x = collect(xlims(hst))
y = @. 10^(a+b*x)

plot!(hst, x, y, color=:darkblue, linestyle=:dot, label="$(round(10^a,digits=3))*10^($(round(b,digits=3)) s)")

savefig(hst, "slope_vs_erosion_rate-type.pdf")
display(hst)

## --- Precipitation vs erosion rate, colored by glacier type (Continental, Alpine, etc.),

hpt = plot(framestyle=:box,
    xlabel="Latitude [degrees, absolute]",
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

t = .!isnan.(ngearth.Precipitation_mm_yr) .& (ngearth.Erosion_rate_mm_yr .>0)
scatternicely!(hpt, ngearth.Precipitation_mm_yr[t], ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

type = ( "Alpine", "Alpine tidewater", "High-latitude", "Continental", "Outlet ice stream",)
tlabel = ( "Alpine", "Alpine tidewater", "Non-cont. high-lat.", "Continental", "Outlet ice stream", )
colors = [mineralcolors[m] for m in ("zircon", "malachite", "azurite", "quartz", "fluid", )]
for i in eachindex(type)
    t = .!isnan.(earth.Precipitation_mm_yr) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Type .== type[i])
    scatternicely!(hpt, earth.Precipitation_mm_yr[t], earth.Erosion_rate_mm_yr[t],
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = tlabel[i],
    )
end

savefig(hpt, "precipitation_vs_erosion_rate-type.pdf")
display(hpt)


## --- Latitude vs erosion rate, colored by glacier type (Continental, Alpine, etc.),

hlt = plot(framestyle=:box,
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
scatternicely!(hlt, abs.(ngearth.Latitude[t]), ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

type = ( "Alpine", "Alpine tidewater", "High-latitude", "Continental", "Outlet ice stream",)
tlabel = ( "Alpine", "Alpine tidewater", "Non-cont. high-lat.", "Continental", "Outlet ice stream", )
colors = [mineralcolors[m] for m in ("zircon", "malachite", "azurite", "quartz", "fluid", )]
for i in eachindex(type)
    t = .!isnan.(earth.Latitude) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Type .== type[i])
    scatternicely!(hlt, abs.(earth.Latitude[t]), earth.Erosion_rate_mm_yr[t],
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = tlabel[i],
    )
end

savefig(hlt, "latitude_vs_erosion_rate-type.pdf")
display(hlt)


## --- Slope vs erosion rate, colored by measurement method (Cosmogenic, Volumetric, etc.)

hsm = plot(framestyle=:box,
    xlabel="Regional slope [m/m]",
    ylabel="Erosion rate [mm/yr]",
    yscale=:log10,
    fontfamily=:Helvetica,
    fg_color_legend=:white,
    legend=:bottomright,
    xlims=(0,0.6),
    ylims = (10^-5, 10^3),
    yticks = 10.0.^(-5:3),

)

t = .!isnan.(ngearth.Slope_m_km) .& (ngearth.Erosion_rate_mm_yr .> 0)
slope = ngearth.Slope_m_km[t]./1000
scatternicely!(hsm, slope, ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)
aₙ,bₙ = linreg(slope, log10.(ngearth.Erosion_rate_mm_yr[t]))
xₙ = collect(xlims(hsm))
yₙ = @. 10^(aₙ+bₙ*xₙ)
plot!(hsm, xₙ, yₙ, color=parse(Color, "#777777"), linestyle=:dot, label="$(round(10^aₙ,digits=3))*10^($(round(bₙ,digits=3)) s)")

method = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
mlabel = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
colors = [mineralcolors[m] for m in ( "spessartine", "rhodochrosite", "corundum", "sodalite", "sapphirine")]
for i in eachindex(type)
    t = .!isnan.(earth.Slope_m_km) .& (earth.Erosion_rate_mm_yr .> 0) .& (earth.Methodology .== method[i])
    scatternicely!(hsm, earth.Slope_m_km[t]./1000, earth.Erosion_rate_mm_yr[t],
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = mlabel[i],
    )
end

t = .!isnan.(earth.Slope_m_km) .& (earth.Erosion_rate_mm_yr .> 0) 
slope = earth.Slope_m_km[t]./1000
a,b = linreg(slope, log10.(earth.Erosion_rate_mm_yr[t]))
x = collect(xlims(hsm))
y = @. 10^(a+b*x)

plot!(hsm, x, y, color=:darkblue, linestyle=:dot, label="$(round(10^a,digits=3))*10^($(round(b,digits=3)) s)")

savefig(hsm, "slope_vs_erosion_rate-method.pdf")
display(hsm)

## --- Precipitation vs erosion rate, colored by measurement method (Cosmogenic, Volumetric, etc.)

hpm = plot(framestyle=:box,
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

t = .!isnan.(ngearth.Precipitation_mm_yr) .& (ngearth.Erosion_rate_mm_yr .>0)
scatternicely!(hpm, ngearth.Precipitation_mm_yr[t], ngearth.Erosion_rate_mm_yr[t],
    label="Nonglacial",
    color=parse(Color, "#dddddd"),
    alpha=1,
    mswidth=0,
)

method = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
mlabel = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
colors = [mineralcolors[m] for m in ( "spessartine", "rhodochrosite", "corundum", "sodalite", "sapphirine")]
for i in eachindex(method)
    t = .!isnan.(earth.Precipitation_mm_yr) .& (earth.Erosion_rate_mm_yr .>0) .& (earth.Methodology .== method[i])
    scatternicely!(hpm, earth.Precipitation_mm_yr[t], earth.Erosion_rate_mm_yr[t],
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = mlabel[i],
    )
end

savefig(hpm, "precipitation_vs_erosion_rate-method.pdf")
display(hpm)


## --- Latitude vs erosion rate, colored by measurement method (Cosmogenic, Volumetric, etc.)

hlm = plot(framestyle=:box,
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
scatternicely!(hlm, abs.(ngearth.Latitude[t]), ngearth.Erosion_rate_mm_yr[t],
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
    scatternicely!(hlm, abs.(earth.Latitude[t]), earth.Erosion_rate_mm_yr[t],
        color = colors[i],
        alpha = 0.85,
        mswidth = 0.25,
        label = mlabel[i],
    )
end

savefig(hlm, "latitude_vs_erosion_rate-method.pdf")
display(hlm)

## --- Combined comparisons

h = plot(hst, hsm, hpt, hpm, hlt, hlm, 
    layout=(3,2),
    size=(400*3, 600*2),
    mswidth=0,
)
savefig(h, "erosion_rate_comparison_method_type.pdf")
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

## --- Plot all data, glacial and nonglacial, by method (grouped)

method = ("Cosmogenic", "Thermochronometric", "Volumetric", "Relief", )

hm = plot(layout = (2,2),
    framestyle=:box,
    xlabel="Timescale [yr]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    size=(300*2,300*2),
    xlims = (10^-3, 10^8),
    xticks = 10.0.^(-3:1:8),
    ylims = (10^-5.5, 10^5.5),
    yticks = 10.0.^(-5:1:5),
    legend=:none,
)

for j in eachindex(method)
    t = contains.(ngearth.Methodology, method[j])

    tt = t .& (ngearth.Type .== "Fluvial")
    plot!(hm[j], ngearth.Time_interval_yr[tt], ngearth.Erosion_rate_mm_yr[tt],
        seriestype=:scatter,
        color = parse(Color, "#bbbbbb"),
        mswidth = 0,
        label = "Fluvial",
    )

    tt = t .& (ngearth.Type .== "Subaerial")
    plot!(hm[j], ngearth.Time_interval_yr[tt], ngearth.Erosion_rate_mm_yr[tt],
        seriestype=:scatter,
        color = parse(Color, "#888888"),
        mswidth = 0,
        label = "Subaerial",
    )

    t = contains.(earth.Methodology, method[j])
    plot!(hm[j], earth.Time_interval_yr[t], earth.Erosion_rate_mm_yr[t],
        seriestype=:scatter,
        color = :darkblue,
        alpha = 0.85,
        mswidth = 0,
        label = "Glacial",
        title = method[j],
    )
end

# Add cosmogenic line of constant 600mm thickness
x = collect(xlims(hm[1]))
for j in 1:1
    plot!(hm[j], x, 600.0./x, color=:red, label="", linestyle=:dash)
    xl = minimum(x)
    annotate!(hm[j], xl*10, 600/(xl*10), text("600 mm", 10, :bottom, :left, color=:red, rotation=-45.0))
end

plot!(hm[2], x, 10e6./x, color=:blue, label="", linestyle=:dash)
xl = 3e2
annotate!(hm[2], xl, 10e6/(xl), text("10 km", 10, :bottom, :left, color=:blue, rotation=-45.0))

# plot!(hm[4], x, 10e6./x, color=:blue, label="", linestyle=:dash)
# annotate!(hm[4], xl, 10e6/(xl), text("10 km", 10, :bottom, :left, color=:blue, rotation=-45.0))

plot!(hm[1], legend=:bottomleft)
savefig("sadlerianness_by_method.pdf")
display(hm)

## -- Area-weighted averages

t = earth.Area_km2 .> 0
area_weighted_glacial = nanmean(earth.Erosion_rate_mm_yr[t], earth.Area_km2[t])

## ---
t = ngearth.Area_km2 .> 0
area_weighted_nonglacial = nanmean(ngearth.Erosion_rate_mm_yr[t], ngearth.Area_km2[t])
