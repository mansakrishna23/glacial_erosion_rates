## -- Beryllium solubility

using StatGeochem, Plots

h = plot(framestyle=:box,
    xlabel="pH",
    ylabel="Be [mol/L]",
    yscale=:log10,
    ylims=10.0.^(-12, -2),
    xlims=(4,9),
)

Ksp_BeOH2 = 10^-21.16 # Lambert and Cleverly 1992, ≈ 6.92e−22

pH = 4:0.01:9
H = 10.0.^-pH
OH = 10.0.^-(14.0.-pH)

# Beₐ * OH^2 == Ksp_BeOH2
Beₐ = Ksp_BeOH2./OH.^2

# BeOH2 saturation
plot!(pH, Beₐ, color=:black, lw=2, label="")
plot!(pH, Beₐ, color=:black, lw=2, fillto=1.0, alpha=0.15, label="")
annotate!(pH[200:200],Beₐ[200:200], text("BeOH₂ saturation", 10, :left, :bottom, color=:black, rotation=-33.33333))

# Ocean Be
hline!([2.33E-02*1e-9,], ribbon=(10e-12, 40e-12), fillalpha=0.15, color=:blue, label="")
annotate!([4.1,],[2.33E-02*1e-9,],text("Seawater Be", 10, :left, :bottom, color=:darkblue, ))
# Ocean pH
yl = ylims()
vline!([8.2,], color=:blue, linestyle=:dot, label="")
plot!([7.8, 8.4], fill(max(yl...), 2), fillto=min(yl...), color=:blue, alpha=0.15, label="")
annotate!([8.2,],[sqrt(prod(yl)),],text("Seawater pH", 10, :center, :bottom, color=:darkblue, rotation=90))


## --- manually compiled Be

water = importdataset("data/beryllium.tsv", '\t', importas=:Tuple)
plot!(water.pH, water.Be_dissolved_nmol_kg.*1e-9, color=:darkblue, seriestype=:scatter, label="Observed surface water Beaq")
# plot!(water.pH, water.Be_nmol_kg.*1e-9, seriestype=:scatter, label="")
# plot!(water.pH, water.Be_colloidal_nmol_kg.*1e-9, seriestype=:scatter, label="")

## ---
savefig("beryllium.pdf")



## --- Net denudation / sediment production

sedflux = 19.1e15 # grams per year, Milliman and Farnsworth 2013
denudation = sedflux / 1000 / 2700 / (148940000*1000^2) * 1000 # mm/yr


## --- Beryllium 10/9 ratios

beprodflux = 3.9e-2 # atoms/cm^2/s, McHargue and Damon 1991
beprod_at_yr = beprodflux * 365.25*24*60*60 * 510072000 * 1000^2 * 100^2
beprod = beprod_at_yr/6.022e23

sedflux = 19.1e15 # grams per year
besed = sedflux * 2.1e-6 / 9.0122 # moles/yr

beprod/besed

## --- Fraction of sedimentary be that actually dissolves
besedreq = beprod/1e-7
besedreq/besed

