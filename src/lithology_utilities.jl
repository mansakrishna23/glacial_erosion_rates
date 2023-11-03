# Functions for querying the Macrostrat / Burwell API for lithologic information

## --- Ping Macrostrat API
    """
    ```julia
    query_macrostrat(lat, lon, zoom::Number=11)
    ```
    Get lithological data for rocks at `lat`, `lon` coordinate from the Macrostrat API.

    Argument `zoom` controls precision; default is approximately 5km. Automatically retry with
    less precise window if initial query does not return data.
    """
    function query_macrostrat(lat, lon, zoom::Number=1)
        resp = HTTP.get("https://macrostrat.org/api/mobile/map_query?lat=$lat&lng=$lon&z=$zoom")
        str = String(resp.body)
        return JSON.Parser.parse(str)
    end


## --- Get data from the unparsed response dictionary
    function get_burwell_min_age(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["t_int_age"]::Number
        catch error
            return NaN
        end
    end

    function get_burwell_max_age(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["b_int_age"]::Number
        catch error
            return NaN
        end
    end

    function get_burwell_map_id(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["map_id"]::Number
        catch error
            return NaN
        end
    end

    function get_burwell_lith(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["lith"]
        catch error
            return "NA"
        end
    end

    function get_burwell_descrip(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["descrip"]
        catch error
            return "NA"
        end
    end

    function get_burwell_name(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["name"]
        catch error
            return "NA"
        end
    end

    function get_burwell_strat_name(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["strat_name"]
        catch error
            return "NA"
        end
    end

    function get_burwell_comments(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["comments"]
        catch error
            return "NA"
        end
    end

    function get_burwell_ref_title(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["ref"]["ref_title"]
        catch error
            return "NA"
        end
    end

    function get_burwell_ref_authors(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["ref"]["authors"]
        catch error
            return "NA"
        end
    end

    function get_burwell_ref_year(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["ref"]["ref_year"]
        catch error
            return "NA"
        end
    end

    function get_burwell_ref_doi(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["ref"]["isbn_doi"]
        catch error
            return "NA"
        end
    end


## --- Parse responses
    function parse_burwell_responses(responses, stop)
        # Preallocate
        rocktype = Array{String}(undef, stop, 1)
        rockdescrip = Array{String}(undef, stop, 1)
        rockname = Array{String}(undef, stop, 1)
        rockstratname = Array{String}(undef, stop, 1)
        rockcomments = Array{String}(undef, stop, 1)
        agemax = Array{Float64}(undef, stop, 1)
        agemin = Array{Float64}(undef, stop, 1)
        age = Array{Float64}(undef, stop, 1)
        authors = Array{String}(undef, stop, 1)
        years = Array{String}(undef, stop, 1)
        titles = Array{String}(undef, stop, 1)
        dois = Array{String}(undef, stop, 1)

        # Parse responses into preallocated arrays
        for i in 1:stop
            # Catch undefined references
            rocktype[i] = get_burwell_lith(responses[i])
            rockdescrip[i] = get_burwell_descrip(responses[i])
            rockname[i] = get_burwell_name(responses[i])
            rockstratname[i] = get_burwell_strat_name(responses[i])
            rockcomments[i] = get_burwell_comments(responses[i])
            agemax[i] = get_burwell_max_age(responses[i])
            agemin[i] = get_burwell_min_age(responses[i])

            authors[i] = get_burwell_ref_authors(responses[i])
            years[i] =get_burwell_ref_year(responses[i])
            titles[i] = get_burwell_ref_title(responses[i])
            dois[i] = get_burwell_ref_doi(responses[i])
        end

        # Filter ages younger than 0 and older than 4000 and compute age from min / max
        # Put age bounds in the correct order
        for i in 1:stop
            agemax[i] = ifelse(0<agemax[i]<4000, agemax[i], NaN)
            agemin[i] = ifelse(0<agemin[i]<4000, agemin[i], NaN)
            age[i] = nanmean([agemax[i], agemin[i]])

            if agemin[i] > agemax[i]
                agemin[i] = agemin[i] + agemax[i]
                agemax[i] = agemin[i] - agemax[i]
                agemin[i] = agemin[i] - agemax[i]
            end
        end

        # Convert strings to lowercase so they can be matched to known names of rock types
        rocktype = lowercase.(rocktype)
        rockdescrip = lowercase.(rockdescrip)
        rockname = lowercase.(rockname)
        rockstratname = lowercase.(rockstratname)
        rockcomments = lowercase.(rockcomments)

        # Replace tabs with spaces so they will not be confused with the \t delimitator
        rocktype = replace.(rocktype, "    " => " ")
        rockdescrip = replace.(rockdescrip, "    " => " ")
        rockname = replace.(rockname, "    " => " ")
        rockstratname = replace.(rockstratname, "    " => " ")
        rockcomments = replace.(rockcomments, "    " => " ")

        # Condense references
        refstrings = @. authors * " | " * years * " | " * titles * " | " * dois

        # Return as a tuple
        return (
            agemax = agemax,
            agemin = agemin,
            age = age,
            rocktype = rocktype,
            rockname = rockname,
            rockdescrip = rockdescrip,
            rockstratname = rockstratname,
            rockcomments = rockcomments,
            refstrings = refstrings
        )
    end

    
## --- Match Macrostrat responses to usable rock names
    """
    ```julia
    match_rocktype(rocktype, rockname, rockdescrip; 
        source::Symbol, 
        [major::Bool], 
        [unmultimatch::Bool])
    ```

    Classify rock samples as sedimentary, igneous, or metamorphic (and associated subtypes)
    based on `rocktype`, `rockname`, and `rockdescrip` sample metadata. Use `rocktype` 
    metadata first; if no matches are made, attempt to classify sample using `rockname` 
    metadata, etc. 

    ### Optional kwarg `major`
    `true` returns: `sed, ign, met`

    `false` returns: `siliciclast, shale, carb, chert, evaporite, coal, sed, volc, plut, 
    ign, metased, metaign, met, cover`

    Major rock types include subclasses; i.e. `ign` includes volcanic and plutonic samples.

    ### Optional kwarg `unmultimatch`
    Setting `unmultimatch=false` will not remove multiply-matched samples. Defaults to `true`.

    # Example
    ```julia
    cats = match_rocktype(rocktype, rockname, rockdescrip, source=:macrostrat, major=true)
    NamedTuple with 4 elements:
    sed    = BitVector(50000,)    [true ... true]
    ign    = BitVector(50000,)    [false ... false]
    met    = BitVector(50000,)    [false ... false]
    cover  = BitVector(50000,)    [false ... false]
    ```
    """
    function match_rocktype(rocktype, rockname, rockdescrip;
        major::Bool=false, 
        unmultimatch::Bool=true,
        inclusive::Bool=true,
        )

        # Get rock type classifications and initialized BitVector
        typelist, cats = get_cats(major, length(rocktype))
        p = Progress(length(typelist)*4, desc="Finding Macrostrat rock types...")

        # Check major lithology 
        for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    cats[j][k] |= match(r"major.*?{(.*?)}", rocktype[k]) |> x -> 
                    isa(x,RegexMatch) ? containsi(x[1], typelist[j][i]) : false
                end
            end
            next!(p)
        end

        # Check the rest of rocktype
        not_matched = find_unmatched(cats)
        @inbounds for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    not_matched[k] && (cats[j][k] |= containsi(rocktype[k], typelist[j][i]))
                end
            end
            next!(p)
        end

        # Then rockname
        not_matched = find_unmatched(cats)
        @inbounds for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    not_matched[k] && (cats[j][k] |= containsi(rockname[k], typelist[j][i]))
                end
            end
            next!(p)
        end

        # Then rockdescrip
        not_matched = find_unmatched(cats)
        @inbounds for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    not_matched[k] && (cats[j][k] |= containsi(rockdescrip[k], typelist[j][i]))
                end
            end
            next!(p)
        end

        # If subtypes are true, major types must also be true
        if !major && inclusive
            minorsed, minorign, minormet = get_minor_types()
            for type in minorsed
                cats.sed .|= cats[type]
            end
            for type in minorign
                cats.ign .|= cats[type]
            end
            for type in minormet
                cats.met .|= cats[type]
            end
        end

        unmultimatch && return un_multimatch!(cats, major)
        return cats
    end


## --- End of file