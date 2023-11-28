using StatGeochem, MAT

bulk = matread("/Users/cbkeller/Documents/Projects/EarthChem/Bulk/bulk.mat")["bulk"]
mcbulk = matread("/Users/cbkeller/Documents/Projects/EarthChem/Bulk/mcbulk.mat")["mcbulk"]

t = mcbulk["Elevation"] .> 0
be_ave = nanmean(mcbulk["Be"][t])

# Himalaya
t = (26.5 .< mcbulk["Latitude"] .< 38.5) .& (70 .< mcbulk["Longitude"] .< 102)
be_himalaya = nanmean(mcbulk["Be"][t])

# SE Asian Islands
t = (-11 .< mcbulk["Latitude"] .< 11) .& (94 .< mcbulk["Longitude"] .< 150)
be_seai = nanmean(mcbulk["Be"][t])


# Igneous
t = floor.(mcbulk["Type"]).==3
be_ign = nanmean(mcbulk["Be"][t])

# Metamorphic
t = floor.(mcbulk["Type"]).==2
be_met = nanmean(mcbulk["Be"][t])

# Sedimentary
t = floor.(mcbulk["Type"]).==1
be_sed = nanmean(mcbulk["Be"][t])

## --- Igneous differentiation

t = (floor.(mcbulk["Type"]) .== 3) .& (mcbulk["Be"] .> 0)
c,m,e = binmeans(mcbulk["SiO2"][t], mcbulk["Be"][t], 40, 80, 20)
Plots.plot(c, m, yerror=2e)
