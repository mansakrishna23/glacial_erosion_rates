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

    # Functions 
    stepify(x::AbstractVector) = vec(vcat(x', x'))
    stepifyedges(x::AbstractVector) = vec(vcat(x[1:end-1]', x[2:end]'))


## --- [GLACIAL] API request
    # Load data 
    glacial_data = importdataset("data/glacial_erosion_Earth.tsv", '\t', importas=:Tuple);
    g_location = [(glacial_data.Latitude[i], glacial_data.Longitude[i]) 
        for i in eachindex(glacial_data.Latitude)]
    unique_g_location = unique(g_location)
    g_npoints = length(unique_g_location)
    g_lats = [unique_g_location[i][1] for i=1:g_npoints]
    g_lons = [unique_g_location[i][2] for i=1:g_npoints]

    # Definitions / Preallocate 
    zoom = 1
    responses = Array{Any}(undef, g_npoints, 1)

    p = Progress(g_npoints, desc="Querying Macrostrat...")
    for i = 1:g_npoints
        try
            responses[i] = query_macrostrat(g_lats[i], g_lons[i], zoom)
        catch
            try
                # Wait and try again
                sleep(10)
                responses[i] = query_macrostrat(g_lats[i], g_lons[i], zoom)
            catch
                # If still nothing, move on
                responses[i] = "No response"
            end
        end
        sleep(0.05)
        next!(p)
    end

    # Parse output
    parsed = parse_burwell_responses(responses, g_npoints)

    # Save data and references
    writedlm("data/lithology_unparsed.tsv", hcat(
        unelementify(parsed), vcat("Latitude", g_lats), vcat("Longitude", g_lons)
    ))
    writedlm("data/lithology_reference.tsv", unique(parsed.refstrings))


## --- [NONGLACIAL] API request
    # Load data 
    nonglacial_data = importdataset("data/nonglacial_erosion_Earth.tsv", '\t', importas=:Tuple);
    ng_location = [(nonglacial_data.Latitude[i], nonglacial_data.Longitude[i]) 
        for i in eachindex(nonglacial_data.Latitude)]
    unique_ng_location = unique(ng_location)
    ng_npoints = length(unique_ng_location)
    ng_lats = [unique_ng_location[i][1] for i=1:ng_npoints]
    ng_lons = [unique_ng_location[i][2] for i=1:ng_npoints]

    # Definitions / Preallocate 
    zoom = 1
    responses = Array{Any}(undef, ng_npoints, 1)

    p = Progress(ng_npoints, desc="Querying Macrostrat...")
    for i = 1:ng_npoints
        try
            responses[i] = query_macrostrat(ng_lats[i], ng_lons[i], zoom)
        catch
            try
                # Wait and try again
                sleep(10)
                responses[i] = query_macrostrat(ng_lats[i], ng_lons[i], zoom)
            catch
                # If still nothing, move on
                responses[i] = "No response"
            end
        end
        sleep(0.05)
        next!(p)
    end

    # Parse output
    parsed = parse_burwell_responses(responses, ng_npoints)

    # Save data and references
    writedlm("data/lithology_unparsed_nonglacial.tsv", hcat(
        unelementify(parsed), vcat("Latitude", ng_lats), vcat("Longitude", ng_lons)
    ))
    writedlm("data/lithology_reference_nonglacial.tsv", unique(parsed.refstrings))

## --- Some notes on methodology
    # We don't want to doule-count samples, but also some of the super broad locations
    # (i.e., country scale) may have more than one glacier.
    
    # Therefore, I want to match each erosion rate to each location. Then look at the
    # individual erosion rates for each glacier. Samples with the same location, lithology, 
    # and that are taken from the same glacier...

## --- Get lithology for each erosion rate measurement
    # Load all data again
    glacial_data = importdataset("data/glacial_erosion_Earth.tsv", '\t', importas=:Tuple);
    glacial_parsed = importdataset("data/lithology_unparsed.tsv", '\t', importas=:Tuple);

    nonglacial_data = importdataset("data/nonglacial_erosion_Earth.tsv", '\t', importas=:Tuple);
    nonglacial_parsed = importdataset("data/lithology_unparsed_nonglacial.tsv", '\t', importas=:Tuple);

    # Get lithology for each unique location
    glacial_cats = match_rocktype(glacial_parsed.rocktype, glacial_parsed.rockname, 
        glacial_parsed.rockdescrip, major=false)

    nonglacial_cats = match_rocktype(nonglacial_parsed.rocktype, nonglacial_parsed.rockname, 
        nonglacial_parsed.rockdescrip, major=false) 

    # Fix false positives and make sure major types are inclusive of minor types
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();

    # Cover can be sedimentary, why not
    # glacial_cats.sed .|= glacial_cats.cover
    # nonglacial_cats.sed .|= nonglacial_cats.cover
    
    glacial_cats.cover .= false
    nonglacial_cats.cover .= false

    # Diorite is not granodiorite, not that it really matters here
    glacial_cats.diorite .&= .!glacial_cats.granodiorite
    nonglacial_cats.diorite .&= .!nonglacial_cats.granodiorite
    
    for type in minorsed
        glacial_cats.sed .|= glacial_cats[type]
        nonglacial_cats.sed .|= nonglacial_cats[type]
    end
    for type in minorvolc
        glacial_cats.volc .|= glacial_cats[type]
        nonglacial_cats.volc .|= nonglacial_cats[type]
    end
    for type in minorplut
        glacial_cats.plut .|= glacial_cats[type]
        nonglacial_cats.plut .|= nonglacial_cats[type]
    end
    for type in minorign
        glacial_cats.ign .|= glacial_cats[type]
        nonglacial_cats.ign .|= nonglacial_cats[type]
    end

    # Match glacial lithology to individual points
    npoints = length(glacial_data.Latitude)
    gcats = NamedTuple{keys(glacial_cats)}([falses(npoints) for _ in keys(glacial_cats)])

    g_location = [(glacial_data.Latitude[i], glacial_data.Longitude[i]) 
        for i in eachindex(glacial_data.Latitude)]
    unique_g_location = [(glacial_parsed.Latitude[i], glacial_parsed.Longitude[i])
        for i in eachindex(glacial_parsed.Latitude)]
    
    for i = 1:npoints
        # Find the index of the parsed location corresponding to sample i 
        k = 0
        for j in eachindex(unique_g_location)
            if unique_g_location[j] == g_location[i]
                k = j
                break
            end
        end

        # If there is an index, get the lithology and put it in the new BitVector
        if k != 0
            rocks = get_type(glacial_cats, k, all_keys=true)
            if rocks !== nothing
                for j in rocks
                    gcats[j][i] = true
                end
            end
        end
    end

    # Match nonglacial lithology to individual points
    npoints = length(nonglacial_data.Latitude)
    ncats = NamedTuple{keys(nonglacial_cats)}([falses(npoints) for _ in keys(nonglacial_cats)])

    ng_location = [(nonglacial_data.Latitude[i], nonglacial_data.Longitude[i]) 
        for i in eachindex(nonglacial_data.Latitude)]
    unique_ng_location = [(nonglacial_parsed.Latitude[i], nonglacial_parsed.Longitude[i])
        for i in eachindex(nonglacial_parsed.Latitude)]
    
    for i = 1:npoints
        # Find the index of the parsed location corresponding to sample i 
        k = 0
        for j in eachindex(unique_ng_location)
            if unique_ng_location[j] == ng_location[i]
                k = j
                break
            end
        end

        # If there is an index, get the lithology and put it in the new BitVector
        if k != 0
            rocks = get_type(nonglacial_cats, k, all_keys=true)
            if rocks !== nothing
                for j in rocks
                    ncats[j][i] = true
                end
            end
        end
    end


## --- Plot erosion rate by lithology
    # Limit glacial data
    # g = (glacial_data.Area_km2 .> 0) .& (glacial_data.Erosion_rate_mm_yr .>0);
    g = (glacial_data.Erosion_rate_mm_yr .>0);
    g .&= (glacial_data.Locality .!= "Dry Valleys, East Antarctica");
    # g .&= (glacial_data.Methodology .!= "Cosmogenic surface")

    # Limit nonglacial data
    n = (nonglacial_data.Erosion_rate_mm_yr .> 0);
    # n .&= (nonglacial_data.Methodology .!= "Cosmogenic surface")

    # Set up for histograms in a loop
    lim = [g, n]
    data = [glacial_data, nonglacial_data]
    cats = [gcats, ncats]
    fig_names = ["Glacial", "Nonglacial"]

    lith = [:sed, :ign, :met]
    colors = [mineralcolors[m] for m in ("apatite", "melt", "quartz", )]
    labels = ["Sedimentary", "Igneous", "Metamorphic"]

    # Make histograms
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, 2)
    for f in eachindex(fig)
        t = lim[f]

        # Set up histogram
        μ = nanmean(log10.(data[f].Erosion_rate_mm_yr[t]))
        binedges = (-4:0.2:4) .+ μ
        bincenters = cntr(binedges)
        stepc = step(bincenters)

        hist = plot(framestyle=:box,
            # xlabel="Erosion Rate [mm/yr]",
            ylabel="Abundance",
            xflip=false,
            grid=false,
            fg_color_legend=:white,
            xlims=(10^-5,10^3),
            # xticks=10.0.^(-5:3),
            xticks = false,
            xscale=:log10,
            size=(600,300),
            bottom_margin=(15,:px), left_margin=(15,:px)
        )

        # Plot each lithology
        for i in eachindex(lith)
            n = histcounts(log10.(data[f].Erosion_rate_mm_yr[cats[f][lith[i]] .& t]), 
                binedges, T=Float64);
            n ./= nansum(n) * stepc
            plot!(hist, 10.0.^stepifyedges(binedges), stepify(n), fill=true, alpha=0.3,
                color=colors[i], label=""
            )
        end

        # Plot each PDF
        for i in eachindex(lith)
            logerosion = log10.(data[f].Erosion_rate_mm_yr[cats[f][lith[i]] .& t])
            μ, σ = nanmean(logerosion), nanstd(logerosion)

            plot!(hist, 10.0.^x, normpdf(μ,σ,x).*(sum(nsed)*step(binedges)),
            linestyle=:solid, linewidth=2,
            color=colors[i],
            label=labels[i],
        )
        end

        # Label
        annotate!(((0.05, 0.95), (fig_names[f] * "\nn = $(count(t))", 12, :left, :top)))
        # annotate!((fig_names[f], 14, :left, :top))

        # display(hist)
        fig[f] = hist
    end

    # Assemble subplots into one plot
    xticks!(fig[2], 10.0.^(-5:3))
    xlabel!(fig[2], "Erosion Rate [mm/yr]")
    fig[2] = plot(fig[2], legend=false)

    h = plot(fig..., layout=(2,1), size=(600,600))
    display(h)
    savefig(h, "erosion_rate_vs_lithology_all.pdf")


## --- Parse responses into useable data
    # Load responses if you already have them!
    parsed = importdataset("data/lithology_unparsed.tsv", '\t', importas=:Tuple)
    # parsed = importdataset("data/lithology_unparsed_nonglacial.tsv", '\t', importas=:Tuple)

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
    # writedlm("data/lithology_parsed.tsv", vcat(header, hcat(lats, lons, hash, rocktypes,
    #     parsed.age, age_uncert, unparsed))
    # )
    writedlm("data/lithology_parsed_nonglacial.tsv", vcat(header, hcat(lats, lons, hash, rocktypes,
        parsed.age, age_uncert, unparsed))
    )

## --- Match lithologies to erosion rates based on location lookup
    iloc = [(data.Latitude[i], data.Longitude[i]) for i in eachindex(data.Latitude)]
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
    t = (data.Time_interval_yr .> 0) .& (data.Erosion_rate_mm_yr .>0) .& (data.Type .!= "Dry Valleys")
    logμ = 10^nanmean(log10.(data.Erosion_rate_mm_yr[t]))
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
    t = (areas[j] .< data.Area_km2 .< areas[j+1]) .& (data.Time_interval_yr .> 0) .& (data.Erosion_rate_mm_yr .>0)
    plot!(h[j], title=titles[j])

    for i in eachindex(lithology)
        tt = t .& cats[lithology[i]]
        plot!(h[j], data.Time_interval_yr[tt], data.Erosion_rate_mm_yr[tt],
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