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
    minorsed, minorign, = get_minor_types()
    for type in minorsed
        cats.sed .|= cats[type]
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


## --- 
    # Means for each lithology, potentially weighted by area (second argument to nanmean)
## --- End of file