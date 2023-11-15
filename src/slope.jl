# Calculate the maximum slope at each point using the SRTM15+ dataset

## --- Set up
    using StatGeochem
    using HDF5
    using DelimitedFiles


## --- Find the average value of slope over an area
    """
    ```julia
    function movingwindow(data::AbstractArray, lat::AbstractArray, lon::AbstractArray, 
        sf::Number=240; 
        n::Number=5, maxpossible::Number=0xffff
    )
    ```

    Find the average value of geospatial `data` in a `n` * `n` km window centered at 
    `lat`ᵢ, `lon`ᵢ.

    ### Optional Kwargs:
    * `sf::Number=240`: Scale factor, or cells per degree, for the geospatial `data`. 
        For 15 arc-second resolution, the scale factor is 240, because 15 arc-seconds go 
        into 1 arc-degree 60 (arc minutes / degree) * 4 (cells / arc minute) = 240 times.

    * `n::Number=5`: Defines the size of the window, in kilometers.

    * `maxpossible::Number=0xffff`: The maximum possible value for the variable of interest; 
        variables with values greater than this are ignored.


    """
    function movingwindow(data::AbstractArray, lat::AbstractArray, lon::AbstractArray, 
            sf::Number=240; n::Number=5, maxpossible::Number=0xffff
        )
        # Interpret user input - make sure DEM is the correct size
        (nrows, ncols) = size(data)
        @assert nrows == 180 * sf + 1   # Latitude
        @assert ncols == 360 * sf + 1   # Longitude

        # Preallocate
        out = fill(NaN ± NaN, length(lat))
        # NShalfwidth = Array{Int64}(undef, nrows, 1)   # Constant (for now?)
        EWhalfwidth = Array{Int64}(undef, nrows, 1)   # Depends on latitude

        # Precalculate the number of grid cells per n×n window at each latitude
        km_per_cell = 1852 * 60 / sf / 1000    # Kilometers per cell (lon) at the equator
        target = n/2
        NShalfwidth = Int(round(target / km_per_cell))

        for i = 1:nrows
            latᵢ = 90 - (1/sf*(i-1))
            gridᵢ = cos(deg2rad(latᵢ)) * km_per_cell
            EWhalfwidth[i] = min(Int(round(target / gridᵢ)), ncols)
        end

        # Find result by indexing into the varname matrix
        for i in eachindex(lat, lon)
            if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
                # Result is NaN if either input is NaN or out of bounds
                continue

            else
                # Convert latitude and longitude into indicies 
                row = 1 + round(Int,(90 + lat[i]) * sf)
                col = 1 + round(Int,(180 + lon[i]) * sf)

                row = (row-EWhalfwidth[row]):(row+EWhalfwidth[row])
                col = (col-NShalfwidth):(col+NShalfwidth)

                # Preallocate
                s = fill(NaN, length(row)*length(col))  # Hold the data values we'll use
                k = 0                                   # Counter

                # Index into the array
                for r in row
                    for c in col
                        # Only do the computation if we are in bounds
                        if 1 <= r <= nrows
                            res = data[r, mod(c-1,ncols-1)+1]

                            # Ignore values that are larger than the max possible value
                            if res < maxpossible
                                k +=1
                                s[k] = res
                            end
                        end
                    end
                end

                # Save the average value
                out[i] = (nansum(s) / k) ± (nanstd(s))
            end
        end

        return out
    end


## --- Calculate the maximum slope at each point
    # # We prefer maximum slope to average slope, because average slope will be biased 
    # lower by "sidehilling" along a slope.

    # # # Download SRTM15+ if necessary, but note that this doesn't work on Windows OS :(
    # # try
    # #     srtm = h5read("data/srtm15plus.h5", "vars/")
    # # catch
    # #     # Note that this will download the SRTM15+ file to \resources, not \data
    # #     srtm = get_srtm15plus()
    # # end

    # srtm = h5read("data/srtm15plus.h5", "vars/")

    # # Calculate slope and save to a file
    # slope = maxslope(srtm["elevation"], srtm["x_lon_cntr"], srtm["y_lat_cntr"], 
    #     srtm["cellsize"], minmatval=-12000
    # )

    # fid = h5open("output/srtm15plus_maxslope.h5","w")
    # g = create_group(fid, "vars")

    # # Copy over SRTM15+ location data
    # g["y_lat_cntr"] = srtm["y_lat_cntr"]
    # g["x_lon_cntr"] = srtm["x_lon_cntr"]
    # g["cellsize"] = srtm["cellsize"]
    # g["scalefactor"] = srtm["scalefactor"]

    # # Add a data set for slope and compress data (Takes about 2.5 minutes)
    # g["slope", compress=3] = slope
    # close(fid)

    
## --- Calculate the average slope (from our maximum slope dataset) for each sample point
    # Load SRTM15+ dataset
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Load erosion data
    glacial = importdataset("data/glacial_erosion_Earth.tsv", '\t', importas=:Tuple)
    nonglacial = importdataset("data/nonglacial_erosion_Earth.tsv", '\t', importas=:Tuple)

    # Get slope at each coordinate point
    # We use a 18 x 18 km moving window, which is the median basin area in our dataset
    slope_glacial = movingwindow(srtm15_slope, glacial.Latitude, glacial.Longitude, 
        srtm15_sf, n=18
    )
    slope_nonglacial = movingwindow(srtm15_slope, nonglacial.Latitude, nonglacial.Longitude, 
        srtm15_sf, n=18
    )

    # Save to file
    writedlm("data/slope_glacial.tsv", slope_glacial, '\t')
    writedlm("data/slope_nonglacial.tsv", slope_nonglacial, '\t')


## --- End of file