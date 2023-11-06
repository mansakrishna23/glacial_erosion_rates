# This code will get lithological information at each point from the Burwell geologic
# map database (https://macrostrat.org/burwell).

    # API Packages
    using HTTP 
    using JSON 

    # File / General Packages
    using DelimitedFiles
    using StatGeochem
    using ProgressMeter

    # Local utilites
    using Static
    include("lithology_utilities.jl")


## --- Load data
    data = importdataset("data/glacial_erosion.tsv", '\t', importas=:Tuple)

    # Get locations, only query each unique location once
    loc = [(data.Latitude[i], data.Longitude[i]) for i in eachindex(data.Latitude)]
    loc = unique(loc)

    npoints = length(loc)
    lats = [loc[i][1] for i=1:npoints]
    lons = [loc[i][2] for i=1:npoints]
    

## --- API request (commented out after we get the data)
    # # Definitions / Preallocate 
    # zoom = 1
    # responses = Array{Any}(undef, npoints, 1)

    # p = Progress(npoints, desc="Querying Macrostrat...")
    # for i = 1:npoints
    #     try
    #         responses[i] = query_macrostrat(lats[i], lons[i], zoom)
    #     catch
    #         try
    #             # Wait and try again
    #             sleep(10)
    #             responses[i] = query_macrostrat(lats[i], lons[i], zoom)
    #         catch
    #             # If still nothing, move on
    #             responses[i] = "No response"
    #         end
    #     end
    #     sleep(0.05)
    #     next!(p)
    # end

    # # Parse output
    # parsed = parse_burwell_responses(responses, npoints)

    # # Save data as a .csv
    # writedlm("data/lithology_unparsed.tsv", unelementify(parsed))


## --- Parse responses into useable data
    # Load responses if you already have them!
    parsed = importdataset("data/lithology_unparsed.tsv", '\t', importas=:Tuple)

    # Match each to a rock name
    rock_cats = match_rocktype(parsed.rocktype, parsed.rockname, parsed.rockdescrip, 
        major=false)

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

    writedlm("data/lithology_parsed.tsv", vcat(["lat" "lon" "type"], hcat(lats, lons, rocktypes)))


## --- End of file