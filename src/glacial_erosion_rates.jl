using StatGeochem, Plots
cd(@__DIR__)
ds = importdataset("glacial_erosion.csv", '\t', importas=:Tuple);

## --- Continental icesheet erosion

minarea = 0 # km^2

tc = (ds.Type .== "Continental") .| (ds.Type .== "High-latitude")
h = plot(framestyle=:box,
    xlabel="Timescale [yr]",
    ylabel="Erosion rate [mm/yr]",
    xscale=:log10,
    yscale=:log10,
    fontfamily=:Helvetica,
    xlims=(0^0, 10^9),
    ylims=(10^-2, 10^1)
)

# plot!(h, ds.Time_interval_yr, ds.Erosion_rate_mm_yr, color=:black, alpha=0.3, seriestype=:scatter, label="")
regions = ("Greenland", "Fennoscandia", "Antarctica", )

colors = (mineralcolors["fluid"], mineralcolors["chrome diopside"], mineralcolors["azurite"])
for i in eachindex(regions)
    region = regions[i]
    t = (ds.Region .== region) .& tc
    if minarea>0
        t .&= ds.Area_km2 .> minarea
    end
    plot!(h, ds.Time_interval_yr[t], ds.Erosion_rate_mm_yr[t], seriestype=:scatter, msc=colors[i], color=colors[i], label=region, alpha=0.85)
    mean_erosion = nanmean(ds.Erosion_rate_mm_yr[t])
    erosion_sigma = nanstd(ds.Erosion_rate_mm_yr[t])./sqrt(count(t))
    hline!(h, [mean_erosion],
        linestyle=:dot,
        color=colors[i],
        label="",
    )
    annotate!(0.9e9, mean_erosion*1.3, text("$region,\nunweighted ave.", colors[i], :right, 8, fontfamily=:Helvetica))
    # plot!(h, [5e8], [mean_erosion], yerr=[erosion_sigma], seriestype=:scatter, msc=colors[i], color=colors[i], label="")
end

# x = [64e6,64e6]
# y = [3000., 5000.].* 1000 ./ x
# plot!(x, y, color=mineralcolors["rhodochrosite"], ms=4, lw=4, marker=:circle, label="")#, label="3-5 km Cryogenian")

cryocolor=mineralcolors["proustite"]
plot!([64e6], [4000*1000/64e6], yerr=[1000*1000/64e6], msc=cryocolor, color=cryocolor, seriestype=:scatter, label="")
annotate!(2.5e8, 4000*1000/64e6, text("3-5 km\nCryogen.\nerosion", cryocolor, :center, :8))

savefig("Continental_glaciation_rates.pdf")
display(h)

## ---
