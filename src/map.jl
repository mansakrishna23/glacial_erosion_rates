
## --- Load packages and read data
using GeoMakie, CairoMakie, Colors
using StatGeochem

# ETOPO gobal elevation dataset
etopo = get_etopo()

# Country outlines, from GeoMakie examples
countries_file = download("https://datahub.io/core/geo-countries/r/countries.geojson")
countries = GeoMakie.GeoJSON.read(read(countries_file, String))

# Erosion rate data
datadir = "$(@__DIR__)/../data"
earth = importdataset("$datadir/glacial_erosion_Earth.tsv", '\t', importas=:Tuple);
ngearth = importdataset("$datadir/nonglacial_erosion_Earth.tsv", '\t', importas=:Tuple)


## --- Plot background

# Our Makie scene
# We'll use the Kavrayskiy VII projection
fig = Figure()
ax = GeoAxis(fig[1,1],
    dest="+proj=kav7",
)

scalefactor = 20
lon = downsample(etopo["x_lon_cntr"], scalefactor)
lat = downsample(etopo["y_lat_cntr"], scalefactor)
elev = downsample(etopo["elevation"]', scalefactor)

# Could NaN things out, but elevations here don't include ice, so Antarctica gets patchy
# elev[elev .< 0] .= NaN
# Instead, let's try making the colormap symmetrical about zero
# (so we're not wasting too much of the colormap on trenches)
emax = nanpctile(elev, 99.9)
elev[elev .> emax] .= emax
elev[elev .< -emax] .= -emax

sf = surface!(ax, lon, lat, elev;
    shading = false
)
cb1 = Colorbar(fig[2,1], sf;
    label = "Elevation [m]",
    vertical = false,
    width = Relative(0.7)
)

n = length(countries)
hm = poly!(ax, countries;
    color=(:black, 0.0),
    strokecolor = :black,
    strokewidth = 0.5,
)

fig

## --- Plot sample locations

t = ngearth.Erosion_rate_mm_yr .> 0
scatter!(ax, ngearth.Longitude[t], ngearth.Latitude[t];
    color = parse(Color, "#cccccc"),
    label = "Nonglacial ($(count(t)))",
    markersize = 10,
)
# t = earth.Erosion_rate_mm_yr .> 0
# scatter!(ax, earth.Longitude[t], earth.Latitude[t];
#     color = mineralcolors["fluid"],
#     label = "Glacial (N=$(count(t)))",
#     markersize = 15,
#     strokewidth = 0.5,
#     strokecolor = :white,
# )

type = ( "Outlet ice stream", "Continental", "High-latitude", "Dry Valleys", "Alpine", )
tlabel = ( "Ice stream", "Continental", "Non-cont. high-lat.", "Dry Valleys", "Alpine", )
colors = [mineralcolors[m] for m in ("fluid", "azurite", "quartz", "spessartine", "zircon", )]
for i in eachindex(type)
    t = contains.(earth.Type, type[i]) .& (earth.Erosion_rate_mm_yr .> 0)
    scatter!(ax, earth.Longitude[t], earth.Latitude[t];
        color = colors[i],
        label = "$(tlabel[i]) ($(count(t)))",
        markersize = 16,
        strokewidth = 0.5,
        strokecolor = :white,
    )
end

# # Legend
# fig[1, 1] = Legend(fig, ax, "Location",
#     framevisible = false,
# )

fig

## -- Save results

save("map.png", fig, px_per_unit=2)


## --- End of File
