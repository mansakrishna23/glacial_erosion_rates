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
    include("utilities_burwell.jl")


## --- Load data
    data = importdataset("data/glacial_erosion.tsv", '\t', importas=:Tuple)

    # Get locations, only query each unique location once
    loc = [(data.Latitude[i], data.Longitude[i]) for i in eachindex(data.Latitude)]
    loc = unique(loc)

    npoints = length(loc)
    lats = [loc[i][1] for i=1:npoints]
    lons = [loc[i][2] for i=1:npoints]
    

## --- API request 
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

    # Save data as a .csv
    writedlm("data/lithology_unparsed.tsv", unelementify(parsed))

## --- Parse responses into useable data
    

## --- End of file