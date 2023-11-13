# This code will get lithological information at each point from the Burwell geologic
# map database (https://macrostrat.org/burwell).

    # API Packages
    using HTTP 
    using JSON 

    # File / General Packages
    using DelimitedFiles
    using StatGeochem
    using ProgressMeter
    using Plots

    # Local utilites
    using Static
    include("lithology_utilities.jl")


## --- Load data
    earth = importdataset("data/glacial_erosion_Earth.tsv", '\t', importas=:Tuple);

    # Get locations, only query each unique location once
    loc = [(earth.Latitude[i], earth.Longitude[i]) for i in eachindex(earth.Latitude)]
    loc = unique(loc)

    npoints = length(loc)
    lats = [loc[i][1] for i=1:npoints]
    lons = [loc[i][2] for i=1:npoints]
    

## --- API request (commented out after we get the data)
    # Definitions / Preallocate 
    zoom = 1
    responses = Array{Any}(undef, npoints, 1)

    p = Progress(npoints, desc="Querying Macrostrat...")
    for i = 1:npoints
        try
            responses[i] = query_macrostrat(lats[i], lons[i], zoom)
        catch
            try
                # Wait and try again
                sleep(10)
                responses[i] = query_macrostrat(lats[i], lons[i], zoom)
            catch
                # If still nothing, move on
                responses[i] = "No response"
            end
        end
        sleep(0.05)
        next!(p)
    end

    # Parse output
    parsed = parse_burwell_responses(responses, npoints)

    # Save data and references
    writedlm("data/lithology_unparsed.tsv", unelementify(parsed))
    writedlm("data/lithology_reference.tsv", unique(parsed.refstrings))


## --- Parse responses into useable data
    # Load responses if you already have them!
    parsed = importdataset("data/lithology_unparsed.tsv", '\t', importas=:Tuple)

    # Match each to a rock name
    rock_cats = match_rocktype(parsed.rocktype, parsed.rockname, parsed.rockdescrip, 
        major=false)
    rock_cats.cover .= false

    # Remove diorite false positive from granodiorite
    rock_cats.diorite .&= .!rock_cats.granodiorite

    # Convert into a human-readable format
    rocktypes = Array{String}(undef, npoints)
    for i in eachindex(rock_cats.sed)
        rocks = get_type(rock_cats, i, all_keys=true)
        out = "" 

        if rocks !==nothing
            for s in rocks
                out *= (string(s) * ", ")
            end
        end
        rocktypes[i] = out
    end

    # Get additional data about the rocks
    unparsed = parsed.rocktype .* " | " .* parsed.rockname .* " | " .* parsed.rockdescrip

    age_uncert = @. (parsed.agemax - parsed.agemin) / 2
    t = .!isnan.(age_uncert)
    age_uncert[t] .= round.(Int, age_uncert[t])

    # Create a hash string for lookups.
    hash = Array{String}(undef, npoints)
    for i in eachindex(hash)
        latᵢ = lats[i]
        lonᵢ = lons[i]

        # Round whole numbers to integers, e.g. 20.0 becomes 20
        if latᵢ % 1 == 0
            latᵢ = Int(latᵢ)
        end
        if lonᵢ % 1 == 0
            lonᵢ = Int(lonᵢ)
        end

        # Save
        hash[i] = string(latᵢ) * ";" * string(lonᵢ)
    end
    
    
    string.(lats) .* ";" .* string.(lons)
    # hash = string.(trunc.(lats, digits=4)) .* ";" .* string.(trunc.(lons, digits=4))

    # Write to file
    header = ["Latitude" "Longitude" "Hash" "Lithology" "Age" "Age Uncert" "Unparsed Lithology"]
    writedlm("data/lithology_parsed.tsv", vcat(header, hcat(lats, lons, hash, rocktypes,
        parsed.age, age_uncert, unparsed))
    )
    

## --- Match lithologies to erosion rates based on location lookup
    iloc = [(earth.Latitude[i], earth.Longitude[i]) for i in eachindex(earth.Latitude)]
    point_lith = Array{String}(undef, length(iloc))

    for i in eachindex(point_lith)
        # Find the index of parsed / loc that corresponds to the sample
        k = 0
        for j in eachindex(loc)
            if iloc[i] == loc[j]
                k = j
                break
            end
        end
        
        if k == 0
            # If there aren't any matches (e.g. NaN lat / lon) then put in an empty string
            point_lith[i] = "" 
        else
            # Use that index to get the lithology
            point_lith[i] = rocktypes[k]
        end
    end

    # Translate that back into a BitVector. There has to be a better way to do this
    cats = match_rocktype(point_lith)
    # minorsed, minorign, = get_minor_types()
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    for type in minorsed
        cats.sed .|= cats[type]
    end
    for type in minorvolc
        cats.volc .|= cats[type]
    end
    for type in minorplut
        cats.plut .|= cats[type]
    end
    for type in minorign
        cats.ign .|= cats[type]
    end


## --- Plot erosion by major rock type (sed, ign, met)
    t = @. !isnan(earth.Time_interval_yr) & !isnan(earth.Erosion_rate_mm_yr);
    t .&= (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .> 0);
    h = plot(framestyle=:box, xlabel="Time Interval [yr]", ylabel="Erosion Rate [mm/yr]")
    
    plot!(h, earth.Time_interval_yr[cats.ign .& t], earth.Erosion_rate_mm_yr[cats.ign .& t], 
        seriestype=:scatter, label="Igneous",)
    plot!(h, earth.Time_interval_yr[cats.sed .& t], earth.Erosion_rate_mm_yr[cats.sed .& t], 
        seriestype=:scatter, label="Sedimentary",)
    plot!(h, earth.Time_interval_yr[cats.met .& t], earth.Erosion_rate_mm_yr[cats.met .& t], 
        seriestype=:scatter, label="Metamorphic",
        yaxis=:log10, xaxis=:log10, legend=:bottomleft
    )


## --- Analyze + Plot
    # Exclude areas without data and Antarctic Dry Valleys
    t = (earth.Area_km2 .> 0) .& (earth.Erosion_rate_mm_yr .>0);
    t .&=  earth.Locality .!= "Dry Valleys, East Antarctica";
    
    # Average erosion rate by rock type
    ignμ = nanmean(log10.(earth.Erosion_rate_mm_yr[cats.ign .& t]))
    volcμ = nanmean(log10.(earth.Erosion_rate_mm_yr[cats.volc .& t]))
    plutμ = nanmean(log10.(earth.Erosion_rate_mm_yr[cats.plut .& t]))
    
    sedμ = nanmean(log10.(earth.Erosion_rate_mm_yr[cats.sed .& t]))
    metμ = nanmean(log10.(earth.Erosion_rate_mm_yr[cats.met .& t]))

    # Weighted by area
    ignμ = nanmean(log10.(earth.Erosion_rate_mm_yr[cats.ign .& t]), earth.Area_km2[cats.ign .& t])
    volcμ = nanmean(log10.(earth.Erosion_rate_mm_yr[cats.volc .& t]), earth.Area_km2[cats.volc .& t])
    plutμ = nanmean(log10.(earth.Erosion_rate_mm_yr[cats.plut .& t]), earth.Area_km2[cats.plut .& t])
    
    sedμ = nanmean(log10.(earth.Erosion_rate_mm_yr[cats.sed .& t]), earth.Area_km2[cats.sed .& t])
    metμ = nanmean(log10.(earth.Erosion_rate_mm_yr[cats.met .& t]), earth.Area_km2[cats.met .& t])

## --- Functions
    stepify(x::AbstractVector) = vec(vcat(x', x'))
    stepifyedges(x::AbstractVector) = vec(vcat(x[1:end-1]', x[2:end]'))


## --- Histograms of erosion rate by lithology
    # Exclude areas without data and Antarctic Dry Valleys
    t = (earth.Area_km2 .> 0) .& (earth.Erosion_rate_mm_yr .>0);
    # dv = (earth.Locality .== "Dry Valleys, East Antarctica") .& t;
    t .&= (earth.Locality .!= "Dry Valleys, East Antarctica");

    # All erosion rates, including dry valleys
    μ = nanmean(log10.(earth.Erosion_rate_mm_yr[t]))
    binedges = (-4:0.2:4) .+ μ
    bincenters = cntr(binedges)
    stepc = step(binedges)

    # Calculate mean and standard deviation for each lithology
    x = first(binedges):0.01:last(binedges)
    logerosion = log10.(earth.Erosion_rate_mm_yr[cats.sed .& t])
    sedμ, sedσ = nanmean(logerosion), nanstd(logerosion)

    logerosion = log10.(earth.Erosion_rate_mm_yr[cats.ign .& t])
    ignμ, ignσ = nanmean(logerosion), nanstd(logerosion)

    logerosion = log10.(earth.Erosion_rate_mm_yr[cats.met .& t])
    metμ, metσ = nanmean(logerosion), nanstd(logerosion)

    # And for dry valleys
    # logerosion = log10.(earth.Erosion_rate_mm_yr[dv])
    # dvμ, dvσ = nanmean(logerosion), nanstd(logerosion)

    # Plot histograms
    hist = plot(framestyle=:box,
        xlabel="Erosion Rate [mm/yr]",
        ylabel="N",
        xflip=true,
        fg_color_legend=:white,
        xlims=(10^-5,10^3),
        xticks=10.0.^(-5:3),
        xscale=:log10,
        size=(600,300)
    )

    nsed = histcounts(log10.(earth.Erosion_rate_mm_yr[cats.sed .& t]), binedges);
    nsed = float(nsed) ./ nansum(float(nsed) .* stepc)
    plot!(hist, 10.0.^stepifyedges(binedges), stepify(nsed), fill=true, alpha=0.5,
        color=mineralcolors["quartz"], label=""
    )

    nign = histcounts(log10.(earth.Erosion_rate_mm_yr[cats.ign .& t]), binedges);
    nign = float(nign) ./ nansum(float(nign) .* stepc)
    plot!(hist, 10.0.^stepifyedges(binedges), stepify(nign), fill=true, alpha=0.5,
        color=mineralcolors["melt"], label=""
    )

    nmet = histcounts(log10.(earth.Erosion_rate_mm_yr[cats.met .& t]), binedges);
    nmet = float(nmet) ./ nansum(float(nmet) .* stepc)
    plot!(hist, 10.0.^stepifyedges(binedges), stepify(nmet), fill=true, alpha=0.5,
        color=mineralcolors["feldspar"], label=""
    )

    # Dry valleys?
    # ndv = histcounts(log10.(earth.Erosion_rate_mm_yr[dv]), binedges);
    # ndv = float(nign) ./ nansum(float(nign) .* stepc)
    # plot!(hist, 10.0.^stepifyedges(binedges), stepify(nmet), fill=true, alpha=0.5,
    #     color=mineralcolors["crocite"], label=""
    # )

    # Plot PDFs    
    plot!(hist, 10.0.^x, normpdf(sedμ,sedσ,x).*(sum(nsed)*step(binedges)),
        linestyle=:solid, linewidth=2,
        color=mineralcolors["quartz"],
        label="Sedimentary",
    )
    plot!(hist, 10.0.^x, normpdf(ignμ,ignσ,x).*(sum(nign)*step(binedges)),
        linestyle=:solid, linewidth=2,
        color=mineralcolors["melt"],
        label="Igneous",
    )
    plot!(hist, 10.0.^x, normpdf(metμ,metσ,x).*(sum(nmet)*step(binedges)),
        linestyle=:solid, linewidth=2,
        color=mineralcolors["feldspar"],
        label="Metamorphic",
    )
    display(hist)
    savefig(hist, "erosion_rate_vs_lithology.pdf")


## --- Repeat erosion rate vs. timescale for area bins by lithology rather than method
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


    lithology = (:sed, :ign, :met)
    lith_label = ("Sedimentary", "Igneous", "Metamorphic")

    # method = ("Cosmogenic", "Cosmogenic surface", "Thermochronometric", "Volumetric",)
    # mlabel = ("Cosmogenic detrital", "Cosmogenic surface", "Thermochronometric", "Volumetric",  "Relief")
    colors = [mineralcolors[m] for m in ( "quartz", "melt", "feldspar")]

    for j in 1:length(areas)-1
    t = (areas[j] .< earth.Area_km2 .< areas[j+1]) .& (earth.Time_interval_yr .> 0) .& (earth.Erosion_rate_mm_yr .>0)
    plot!(h[j], title=titles[j])

    for i in eachindex(lithology)
        tt = t .& cats[lithology[i]]
        plot!(h[j], earth.Time_interval_yr[tt], earth.Erosion_rate_mm_yr[tt],
            seriestype=:scatter,
            color = colors[i],
            alpha = 0.85,
            mswidth = 0.1,
            label = lith_label[i],
        )
    end
    end

    savefig(h, "timescale_vs_erosion_rate-area_lithology.pdf")
    display(h)


## --- End of file