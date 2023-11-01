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
    """
    function get_rock_class()
        # Sedimentary
        siliciclast = ("siliciclast", "conglo", "sand", "psamm", "arenit", "arkos", "silt",
            "breccia", "quartzite", "leptite")
        shale = ("lutite", "mud", "clay", "shale", "wacke", "argillite", "argillaceous", 
            "flysch", "pelit", "turbidite", "tasmanite", "slate", "phyllite",)
        carb = ("carbonate", "limestone", "dolo", "marl", "chalk", "coquina", "biogenic", 
            "travertine", "tavertine", "tufa", "calcarenite", "teravertine", "marble", 
            "calc silicate", "calcsilicate", "skarn", )
        evap = ("evaporite", "anhydrite", "gypsum", "trona", "halite", "sylvite", 
            "salt flat", "caliche", "exhalite")
        chert = ("chert", "opal", "porcellanite", "diatomite", "novaculite", "iron", 
            "taconite", "banded iron")
        sed = ("sediment", "clast", "diamict", "tillite", "stream deposits", 
            "beach deposits", "terrace",  "marine deposits",  "paleosol", "spiculite", 
            "glauconite", "meta-sed", "metased", "paragneiss", "para")

        # Volcanic
        komatiite = ("komatiite", "bergalite", "kimberlite", "meimechite", "picrite", 
            "polzenite", "ultramafitite")
        basalt = ("basalt", "pillow", "scoria", "absarokite", "anamesite", "hawaiite", 
            "linosaite", "mafite", "mugearite", "oceanite", "picrite", "palagonite",
            "shoshonite", "tachylyte", "tholeiite", "mafic", "hauynite", "porphyrite",
            "ottajanite",)
        andesite = ("andesit", "boninite", "icelandite", "marianite", "tristanite", 
            "vulsinite",)
        dacite = ("dacit", "santorinite", )
        rhyolite = ("rhyolit", "dellenite", "rhyodacite", "felsite", "liparite", 
            "pantellerite", "felsic", "silicic", "comendite", "latite", "vulsinite",
            "pumice", "obsidian", "glass",)
        alk_volc = ("adakite", "alnoite", "arsoite", "augitite", "benmoreite", "camptonite", 
            "ciminite", "damkjernite", "damtjernite", "domite", "fortunite", "gauteite",
            "kenyte", "keratophyre", "kersantite", "kivite", "lampro", "madupite",
            "melnoite", "minette", "monchiquite", "mondhaldeite", "orangeite",
            "orendite", "phonolite", "sannaite", "trachyte", "wyomingite", "fonolito",
            "tinguaite", "vsbergite", "ordanchite", "melilit", "katungite", "analcimite", 
            "ankaratrite", "etindite", "foidite", "grazinite", "hauynophyre", "kalsilit", 
            "leucitite", "mafurite", "melafoidite", "nephelinite","ugandite", )
        volc = ("volcanic", "extrusive", "lava", "eutaxite", "vitrophyre", "volcan", 
            "diatreme", "pipe", "ash", "ashfall", "tuff",  "tephra", "cinder", )

        # Plutonic
        peridotite = ("olivinite", "dunit", "lherzolite", "peridot", "harzburg", 
            "bronzitite", "enstatitite", "hornblendite", "pyroxenite", "turjaite",
            "uncompahgrite", "websterite", "wehrlite", "cortlandite", "vibetoite",
            "yamaskite",)
        gabbro = ("gabbro", "mafraite", "allivalite", "anorthosite", "crinanite", 
           "teschenite", "diabase", "dolerit", "essexite", "glenmuirite", "jotunite",
           "labradorite", "luscladite", "theralite", "norite", "troctolite",
           "sebastianite", "amphibolit",)
        diorite = ("diorit", "jotunite", "marscoite", "sanukite", "trondhjemite", "trond",)
        tonalite = ("tonalit", "adamellite", "enderbite", "enderbite",)
        granodiorite = ("granodiorite")
        granite = ("granit", "microgranite", "adamellite", "aplite", "charnockite", 
            "granophyre", "rapakivi", "monzonit", "mangerite", "greisen",)
        alk_plut = ("syenit", "alaskite", "borolanite", "bostonite", "durbachite", 
            "foyaite", "jacupirangite", "juvite", "kentallenite", "larvikite", "lujavrite",
            "nordmarkite", "orthoclasite", "shonkinite", "sommaite", "kaersutitite",
            "lestiwarite", "puglianite", "vaugnerite","fergusite", "ijolite", 
            "melteigite", "missourite", "tannbuschite", "buchonite", "campanite", 
            "murambite", "tephrite", "tahitite","vicoite", "urtite", "ankaramite", 
            "basanit", "limburgite", "biotitite", "riedenite", "glimmerite", "kamafugite",)
        plut = ("plutonic", "pluton", "intrusive", "intrus", "sill", "dike", "stock", 
            "laccolith", "lopolith", "batholith", "pegmatite", "porphyry", "megacryst",
            "hypabyssal", "pegmat",)

        # Igneous
        carbonatite = ("alvikite", "carbonatite", "beforsite", "rauhaugite", "sovite",)
        ign = ("igneous", "metaign", "orthogneiss", "ortho", "meta-ign",)

        # Unknown metamorphic
        met = ("crystalline", "migma", "alter", "hydrothermal", "basement", 
            "high grade metamorphic", "meta",)

        # Cover
        cover = ("lluv", "fluv", "boulder", "gravel", "aleurite", "glaci", "till", "loess", 
            "regolith", "debris", "fill", "slide", "unconsolidated", "talus", "stream", 
            "beach", "terrace", "placer", "paleosol", "mass-wasting", "pebble", "cover", 
            "quaternary", "soil", "laterite", "surficial deposits", "scree", "peat", 
            "swamp", "marsh", "water", "ice")

        # Sedimentary
        phosphorite = ("phosphorite", "phosphate")
        coal = ("coal", "anthracite", "peat", "lignite", "bitumen")
        volcaniclast = ("tonstein", "peperite", "volcaniclastic")

        # Igneous
        volc = (
            "lahar", # Mudflow--associated with volcanism, but are sediments volcanic?
            "ignimbrite", "lenticulite", # Associated with ignimbrites?  
            "breunneritite", # Magnesite?
            "vitrophere",   # Flow breccia?  
            )
        plut = ("corganite", "corgaspinite", # Something to do with garnet paragenesis? Maybe should be moved to met?
            "kullaite", # Maybe hypabyssal? Not a lot of information 
            "chromitite", "apatitite",   # Monomineralic..? Not really any alkalis
            "topazite", # hypabyssal, just topaz and quartz
            )
        ign = ("basite", "metabasite", # idk
            "phoscorite",   # magnetite, apatite and olivine, usually associated with carbonatites
            )

        # Metamorphic
        metased = ("schist",  "hornfels",) # Always?
        metaign = ( 
            "soapstone", "talc", # Ultramafic OR silicious dolomite protolith??
            "serpentin",  # Mafic to ultramafic serpentinization
            "greenstone", # Granofels with green minerals 
            "eclogite",   # Gabbro? Mostly garnet and pyroxene...
            "halleflinta", "leucophyre", "melaphyre", "propylite", "spilite", 
            "alkremite",  
            "greenschist", "blueschist", "zeolite", "rodingite",
            )

        met = (( "garnet", "buchite", "epidot", "fenite", "albitite", "chloritite", 
            "phlogopitite",  "sericitite", "tactite",  "tourmalinite", "unakite", 
            "vogesite", "gossan",  "palagonite", 
            "gneiss", "granulit", "granofels", "sanidinite",   
            )...,
            cataclastic...,
        )
    end


## --- End of file