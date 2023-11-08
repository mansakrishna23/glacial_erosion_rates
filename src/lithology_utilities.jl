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

    
## --- Definitions that let the lithology matching functions work
    """
    ```julia
    get_rock_class([major::Bool], [inclusive::Bool])
    ```

    Define sedimentary, igneous, and metamorphic rock types and subtypes. Metamorphic rocks
    are grouped with their protoliths when possible.

    ### Possible subtypes
    *  Sedimentary: siliciclastic, shale, carbonate, evaporite, chert, phosphorite, coal,
        sedimentary (uncategorized)
    *  Igneous: 
        *  Volcanic: komatiite, basalt, andesite, dacite, rhyolite, alkaline, 
            volcaniclastic, volcanic (uncategorized)
        *  Plutonic: peridotite, pyroxenite, gabbro, diorite, trondhjemite, tonalite, 
            granodiorite, granite, alkaline, plutonic (uncategorized)
        *  Carbonatite
        *  Igneous (uncategorized)
    * Metamorphic: metamorphic (uncategorized)

    ### Optional kwargs
        
        major

    Return only Tuples for sedimentary, metamorphic, and igneous types. Boolean; defaults 
    to `false`

        inclusive

    Sedimentary, igneous, and metamorphic lists include terms for all included subtypes. 
    Boolean; defaults to `false` unless `major` is `true`.

    Note that volcanic and plutonic rocks will include all nested volcanic / plutonic 
    subtypes.
    """
    function get_rock_class(; major::Bool=false)
        # Sedimentary
        siliciclast = ("siliciclast", "conglo", "sand", "psamm", "arenit", "arkos", "silt",
            "breccia", "quartzite")
        shale = ("lutite", "mud", "clay", "shale", "wacke", "argillite", "argillaceous", 
            "flysch", "pelit", "turbidite", "tasmanite", "slate", "phyllite",)
        carb = ("carbonate", "limestone", "dolo", "marl", "chalk", "coquina", "biogenic", 
            "travertine", "tavertine", "tufa", "calcarenite", "teravertine", "marble", 
            "calc silicate", "calcsilicate", "skarn", )
        evap = ("evaporite", "anhydrite", "gypsum", "trona", "halite", "sylvite", 
            "salt flat", "caliche", "exhalite")
        chert = ("chert", "opal", "porcellanite", "diatomite", "novaculite", "iron", 
            "taconite", "banded iron")
        phosphorite = ("phosphorite", "phosphate")
        coal = ("coal", "anthracite", "lignite", "bitumen")
        sed = ("sediment", "clast", "diamict", "tillite", "stream deposits", 
            "beach deposits", "terrace",  "marine deposits",  "paleosol", "spiculite", 
            "glauconite", "meta-sed", "metased", "paragneiss", "para")

        # Volcanic
        komatiite = ("komatiite", "meimechite", "ultramafitite")
        basalt = ("basalt", "pillow", "scoria", "picrite", "anamesite", "hawaiite", 
            "mafite", "mugearite", "oceanite", "palagonite", "tachylyte", "tholeiite", 
            "mafic", "melaphyre", "greenstone", "spilite", "greenschist", "blueschist", 
            "basite", "metabasite",)
        andesite = ("andesit", "boninite", "icelandite", "marianite", "adakite", 
            "propylite",)
        dacite = ("dacit", "santorinite", "ignimbrite",)
        rhyolite = ("rhyolit", "felsite", "liparite", "felsic", "silicic", "pumice", 
            "obsidian", "dellenite", "rhyodacite", "ignimbrite", "lenticulite", 
            "halleflinta", "leptite",)
        alk_volc = ( "polzenite", "hauynite", "arsoite", "benmoreite", "camptonite", 
            "ciminite", "damkjernite", "damtjernite", "domite", "fortunite", "gauteite",
            "kenyte", "keratophyre", "kersantite", "kivite", "lampro", "madupite",
            "minette", "monchiquite", "mondhaldeite", "orendite", "phonolite", "sannaite", 
            "trachyte", "wyomingite", "fonolito", "tinguaite",  "ordanchite", "melilit", 
            "katungite", "vsbergite", "analcimite", "ankaratrite", "etindite", "foidite", 
            "grazinite", "hauynophyre", "kalsilit", "leucitite", "mafurite", "melafoidite",
            "nephelinite","ugandite", "ottajanite", "melnoite", "pantellerite", "comendite", 
            "latite", "tristanite", "augitite", "absarokite", "shoshonite", "linosaite",
            "bergalite", "alnoite", "kimberlite",  "orangeite", "diatreme", "pipe", )
        volcaniclast = ("tonstein", "peperite", "volcaniclastic", "lahar",)
        volc = ("volcanic", "extrusive", "lava", "eutaxite", "vitrophyre", "volcan", 
            "ash", "ashfall", "tuff",  "tephra", "cinder", "porphyrite", 
            "vulsinite", "glass", "vitrophere",)
            
        # Plutonic
        peridotite = ("olivinite", "dunit", "lherzolite", "peridot", "harzburg", 
            "wehrlite", "serpentin", "soapstone", "talc", "alkremite", )
        pyroxenite = ("bronzitite", "pyroxenite", "enstatitite", "websterite",
            "hornblendite", "cortlandite",)
        gabbro = ("gabbro", "mafraite", "allivalite", "anorthosite", "diabase", "dolerit", 
            "leucophyre", "glenmuirite", "jotunite", "labradorite", "luscladite", 
            "theralite", "norite", "troctolite", "sebastianite", "eclogite", "amphibolit", 
            "rodingite", "corganite", "corgaspinite",)
        diorite = ("diorit", "jotunite", "marscoite", "sanukite",)
        trondhjemite =  ("trondhjemite", "trond",)
        tonalite = ("tonalit", "adamellite", "enderbite", "enderbite",)
        granodiorite = ("granodiorite",)
        granite = ("granit", "microgranite", "adamellite", "aplite", "charnockite", 
            "granophyre", "rapakivi", "monzonit", "mangerite", "greisen", "pegmat",)
        alk_plut = ("syenit", "alaskite", "borolanite", "bostonite", "durbachite", 
            "foyaite", "jacupirangite", "juvite", "kentallenite", "larvikite", "lujavrite",
            "nordmarkite", "orthoclasite", "shonkinite", "sommaite", "kaersutitite",
            "lestiwarite", "puglianite", "vaugnerite","fergusite", "ijolite", 
            "melteigite", "missourite", "tannbuschite", "buchonite", "campanite", 
            "murambite", "tephrite", "tahitite", "vicoite", "urtite", "ankaramite", 
            "basanit", "limburgite", "biotitite", "riedenite", "glimmerite", "kamafugite",
            "turjaite", "essexite", "yamaskite", "teschenite", "crinanite", "vibetoite",  
            "uncompahgrite", "apatitite", "nelsonite", "phoscorite", "kullaite", )
        plut = ("plutonic", "pluton", "intrusive", "intrus", "sill", "dike", "stock", 
            "laccolith", "lopolith", "batholith", "porphyry", "megacryst",
            "hypabyssal", "chromitite", "topazite", )
            
        # Undefined igneous
        carbonatite = ("alvikite", "carbonatite", "beforsite", "rauhaugite", "sovite",
            "breunneritite",)
        ign = ("igneous", "metaign", "orthogneiss", "ortho", "meta-ign", "zeolite",)

        # Undefined metamorphic
        met = ("crystalline", "migma", "alter", "hydrothermal", "basement", 
            "high grade metamorphic", "meta", "granulit", "granofels", "gneiss", "schist", 
            "hornfels", "garnet", "buchite", "epidot", "fenite", "albitite", "chloritite", 
            "phlogopitite", "sericitite", "tactite", "tourmalinite", "unakite", 
            "vogesite", "gossan", "palagonite", "sanidinite",)

        # Cover
        cover = ("lluv", "fluv", "boulder", "gravel", "aleurite", "glaci", "till", "loess", 
            "regolith", "debris", "fill", "slide", "unconsolidated", "talus", "stream", 
            "beach", "terrace", "placer", "paleosol", "mass-wasting", "pebble", "cover", 
            "quaternary", "soil", "laterite", "surficial deposits", "scree", "peat", 
            "swamp", "marsh", "water", "ice")

        # Define major and minor rock types
        minorsed = (:siliciclast, :shale, :carb, :evap, :chert, :phosphorite, :coal,)
        minorvolc = (:komatiite, :basalt, :andesite, :dacite, :rhyolite, :alk_volc, 
            :volcaniclast,)
        minorplut = (:peridotite, :pyroxenite, :gabbro, :diorite, :trondhjemite, :tonalite, 
            :tonalite, :granodiorite, :granite, :alk_plut,)
        minorign = (:volc, :plut, :carbonatite)

        if major
            typelist = (sed=(minorsed, sed...,), ign=(minorign, ign...,), met=met, cover=cover)
            minors = ()
        else
            typelist = (
                # Sedimentary
                siliciclast=siliciclast, shale=shale, carb=carb, evap=evap, chert=chert, 
                phosphorite=phosphorite, coal=coal, sed=sed,

                # Volcanic
                komatiite=komatiite, basalt=basalt, andesite=andesite, dacite=dacite, 
                rhyolite=rhyolite, alk_volc=alk_volc, volcaniclast=volcaniclast, volc=volc, 
                
                # Plutonic
                peridotite=peridotite, pyroxenite=pyroxenite, gabbro=gabbro, diorite=diorite, 
                trondhjemite=trondhjemite, tonalite=tonalite, granodiorite=granodiorite, 
                granite=granite, alk_plut=alk_plut, plut=plut, 

                # Igneous
                carbonatite=carbonatite,
                ign=ign,

                # Metamorphic
                met=met,

                # Cover
                cover=cover,
            )
            minors = (minorsed, minorvolc, minorplut, minorign)
        end

        return typelist, minors...
    end

    """
    ```julia
    get_cats(major::Bool, npoints::Int64)
    ```

    Meow! 🐈 Initialize a NamedTuple of `npoints`-element BitVectors for each defined rock 
    type.

    If `major` is `true`, only major types are defined. If `major` is `false`, all subtypes
    are defined.

    See also: `get_minor_types` and `get_rock_class`.

    # Example
    ```julia
    typelist, cats = get_cats(true, 10)
    ```
    """
    function get_cats(major::Bool, npoints::Int64)
        typelist = get_rock_class(major=major)[1]

        return typelist, NamedTuple{keys(typelist)}([falses(npoints) for _ in 1:length(typelist)]) 
    end

    """
    ```julia
    get_minor_types()
    ```

    Return types nested under the sed, ign, and met "major" types.

    **Important: this function will break if major types are listed before minor types in
    the `get_cats` function, and if sedimentary types are not listed first!**

    ### Minor Types:
    * Sed: siliciclast, shale, carb, chert, evaporite, coal, phosphorite, volcaniclast
    * Ign: volc, plut
    * Met: metased, metaign

    # Example
    ```julia
    minorsed, minorign, minormet = get_minor_types()
    ```
    """
    function get_minor_types()
        types = get_cats(false, 1)[2]
        allkeys = collect(keys(types))

        sed = findfirst(==(:sed), allkeys)
        ign = findfirst(==(:ign), allkeys)
        met = findfirst(==(:met), allkeys)

        return Tuple(allkeys[1:sed-1]), Tuple(allkeys[sed+1:ign-1]), Tuple(allkeys[ign+1:met-1])
    end

    """
    ```julia
    un_multimatch!(cats, major::Bool)
    ```

    Exclude rock type matches from each other so each sample is only classified as one
    rock type.

    If `major` is `true`: 
    * Cover is excluded from all rock types.
    * Rocks classified as both metamorphic and sedimentary / igneous (i.e., metasedimentary
        and metaigneous rocks) are respectively re-classified as sedimentary and igneous rocks.
    * Rocks classified as sedimentary and igneous are re-classified as igneous rocks.

    If `major` is false:
    * Cover is excluded from all rock types.
    * Rocks classified as both metamorphic and sedimentary / igneous (i.e., metasedimentary
        and metaigneous rocks) are **excluded** from sedimentary and igneous rocks, and 
        re-classified as metasedimentary and metaigneous rocks.
    * Rocks classified as more than one subtype of sedimentary rocks are excluded from
        each other, in arbitrary order.
    * Rocks classified as both volcanic and plutonic, or both metasedimentary and 
        metaigneous, are respectively re-classified as undifferentiated igneous and metamorphic
        rocks.
    * Rocks classified as sedimentary and igneous are re-classified as igneous rocks.
    """
    un_multimatch!(cats, major::Bool) = _un_multimatch!(cats, static(major))
    function _un_multimatch!(cats, major::True)
        # Exclude cover
        cats.sed .&= .! cats.cover
        cats.ign .&= .! cats.cover
        cats.met .&= .! cats.cover

        # Classify metased as sed and metaign and ign
        cats.met .&= .! cats.sed
        cats.ign .&= .! cats.ign

        # Sed / ign rocks are classified as ign
        cats.sed .&= .! cats.ign

        return cats
    end

    function _un_multimatch!(cats, major::False)
        # Define types
        minorsed, minorign, minormet = get_minor_types()
        minortypes = (minorsed..., minorign..., minormet...)

        # Exclude cover from all major and minor rock types
        cats.sed .&= .! cats.cover
        cats.ign .&= .! cats.cover
        cats.met .&= .! cats.cover
        for i in minortypes
            cats[i] .&= .! cats.cover
        end

        # Exclude metamorphic rocks from sed and igns. Class as metased and metaign
        cats.metased .|= (cats.sed .& cats.met)
        cats.metaign .|= (cats.ign .& cats.met)

        cats.sed .&= .! cats.met
        for i in minorsed
            cats[i] .&= .! cats.met
        end

        cats.ign .&= .! cats.met
        for i in minorign
            cats[i] .&= .! cats.met
        end

        # Exclude sed subtypes from other sed subtypes
        for i in 1:(length(minorsed)-1)
            for j in minorsed[i+1:end]
                cats[minorsed[i]] .&= .! cats[j]
            end
        end

        # Both volcanic and plutonic is undifferentiated igneous
        cats.volc .&= .!(cats.volc .& cats.plut)
        cats.plut .&= .!(cats.volc .& cats.plut)

        # Both metaigneous and metasedimentary is undifferentiated metamorphic
        cats.metased .&= .!(cats.metased .& cats.metaign)
        cats.metaign .&= .!(cats.metased .& cats.metaign)

        # Sed / ign rocks are classified as ign
        cats.sed .&= .! cats.ign            
        for i in minorsed
            cats[i] .&= .! cats.ign 
        end

        return cats
    end


## --- Match Macrostrat responses to usable rock names
    """
    ```julia
    find_unmatched(cats)
    ```

    Given a `Tuple` of `BitVectors`, return a `BitVector` that is `true` at index `i` iff 
    all elements of the `Tuple` are `false` at index `i`.

    If `cats` is a `NamedTuple` of rock types defined by `get_cats`, specify `major` 
    as `true` or `false` to decrease runtime. `major` is `true` if `cats` contains only 
    `sed`, `ign`, `met`, and `cover`.

    # Example
    ```julia-repl
    julia> find_unmatched(cats)
    500-element BitVector:
    0
    1
    1
    1
    ⋮
    1
    1
    1
    1
    ```
    """
    function find_unmatched(cats)
        matched = falses(length(cats[1]))
        @inbounds for i in eachindex(cats)
            matched .|= cats[i]
        end
        return .!matched
    end

    """
    ```julia
    match_rocktype(rocktype, rockname, rockdescrip; 
        [major::Bool], 
        [unmultimatch::Bool]
        [inclusive::Bool=true])
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

    ### Optional kwarg `inclusive`
    Defines if major types should include their mapped minor subtypes; e.g. does `sedimentary`
    include `siliciclastic`? Defaults to `true`.

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
        if !major
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

        return cats
    end

    """
    ```julia
    match_rocktype(writtentype::AbstractArray{String})
    ```

    Return a `NamedTuple` of `BitVector`s catagorizing Macrostrat samples as sedimentary, 
    igneous, metamorphic, and associated subtypes, or cover from types stored as strings 
    in `writtentype`.

    Major types do not include minor subtypes.
    """
    function match_rocktype(writtentype::AbstractArray{String})
        cats = get_cats(false, length(writtentype))[2]

        # Parse all of the written types into cats
        for i in eachindex(writtentype)
            types = split(writtentype[i], ",")
            for t in eachindex(types)
                try
                    cats[Symbol(types[t])][i] = true
                catch
                    continue
                end
            end
        end

        return cats
    end


## --- Additional metadata for rock types
    """
    ```julia
    get_type(cats::NamedTuple, i; [all_keys=false])
    ```

    Return the first key of `cats` where `i` is `true`. Optionally specify `allkeys`=`true`
    to return _all_ keys where `i` is true. Assumes equal length elements of `cats`.

    # Examples
    ```julia-repl
    julia> get_type(macro_cats, 254)
    :sed

    julia> get_type(name_cats, 3, all_keys=true)
    (:gravel, :sand, :silt, :clay)
    ```
    """
    get_type(cats::NamedTuple, i::Int64; all_keys::Bool=false) = _get_type(cats, i, static(all_keys))

    function _get_type(cats, i, all_keys::False)
        @assert 0 < i <= length(cats[1]) "Index $i out of bounds."

        @inbounds for k in keys(cats)
            cats[k][i] && return Symbol(k)
        end

        return nothing
    end

    function _get_type(cats, i, all_keys::True)
        @assert 0 < i <= length(cats[1]) "Index $i out of bounds."

        catkeys = keys(cats)
        keymatches = falses(length(catkeys))

        @inbounds for k in eachindex(catkeys)
            cats[k][i] && (keymatches[k]=true)
        end

        count(keymatches)==0 && return nothing
        return catkeys[keymatches]
    end

    
## --- End of file