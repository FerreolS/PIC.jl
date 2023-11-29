"""
    fit_wavelamps_specpos(arguments...; keywords...) -> Result

For each lenslet, fit the parameters of the model for the data of Wave,Lamp and Spec,Pos.

# Arguments
- `lenses_positions ::AbstractMatrix{Float64}`: absolute position of the center of each lens on
  the detector for both Wave,Lamp and Spec,Pos files. First axis contains the x and y position.
- `wavelamps_λlasers ::AbstractVector{Float64}`: wavelength for each laser of the Wave,Lamp file.
  In YJ mode 3 lasers are expected and in YJH 4 lasers are expected.
- `wavelamps_data ::AbstractMatrix{T}`: Matrix containing the Wave,Lamp data.
- `wavelamps_weights ::AbstractMatrix{T}`: Matrix containing the Wave,Lamp weights.
- `wavelamps_fwhm_init ::AbstractVector{Float64}`: initial fwhm value for each Wave,Lamp laser.
- `wavelamps_cx_init ::AbstractVector{Float64}`: initial coefficient for the polynome of x center
  positions of the lasers.
- `wavelamps_cy_init ::AbstractVector{Float64}`: same as `wavelamps_cx_init` for y position.
- `specpos_data ::AbstractMatrix{T}`: Matrix containing the Spec,Pos data.
- `specpos_weights ::AbstractMatrix{T}`: Matrix containing the Spec,Pos weights.

The WaveLamp order (for polynomes of x and y positions) is defined to be the number of
lasers minus one.

# Keywords Arguments
- `specpos_order ::Int = 2`: SpecPos order for the polynomes of x center position and fwhm.
- `λ0 ::Float64 = mean(wavelamps_λlasers)`: reference wavelength for both Wave,Lamp and Spec,Pos.
- `λrange ::AbstractRange{Float64} = IFS_λRANGE`: wavelength range used for Spec,Pos.
- `init_box ::BoundingBox{Int} = INIT_BOX`: box surrounding each lens center.
- `valid_lenses_map ::AbstractVector{Bool} = trues(size(lenses_positions,2)`:
  if `valid_lenses_map[i]` is false, the fit of the lenslet `i` will be skipped.
"""
function fit_wavelamps_specpos(
    lenses_positions ::AbstractMatrix{Float64},
    
    wavelamps_λlasers ::AbstractVector{Float64},
    wavelamps_data    ::AbstractMatrix{T},
    wavelamps_weights ::AbstractMatrix{T},
    wavelamps_fwhm_init ::AbstractVector{Float64},
    wavelamps_cx_init   ::AbstractVector{Float64},
    wavelamps_cy_init   ::AbstractVector{Float64},
    
    specpos_data      ::AbstractMatrix{T},
    specpos_weights   ::AbstractMatrix{T},
    
    ; specpos_order ::Int = 2,
      λ0        ::Float64 = mean(wavelamps_λlasers),
      λrange    ::AbstractRange{Float64} = IFS_λRANGE,
      init_box ::BoundingBox{Int} = INIT_BOX,
      valid_lenses_map ::AbstractVector{Bool} = trues(size(lenses_positions,2))
      
) ::FitResult where {T<:Real}

    nb_lenses = size(lenses_positions,2)
    nλ = length(wavelamps_λlasers)
    wavelamp_order = nλ - 1
    
    size(wavelamps_data) == (2048,2048) || throw(DimensionMismatch(
        "size of `wavelamps_data` must be (2048,2048)"))
    
    size(wavelamps_weights) == (2048,2048) || throw(DimensionMismatch(
        "size of `wavelamps_weights` must be (2048,2048)"))
    
    nλ == length(wavelamps_fwhm_init) || throw(DimensionMismatch(
        "length of `wavelamps_fwhm_init` must be equal to length of `wavelamps_λlasers`"))
    
    wavelamp_order == length(wavelamps_cx_init) || throw(DimensionMismatch(
        "length of `wavelamps_cx_init` must be equal to length of `wavelamps_λlasers - 1`"))
    
    wavelamp_order == length(wavelamps_cy_init) || throw(DimensionMismatch(
        "length of `wavelamps_cy_init` must be equal to length of `wavelamps_λlasers - 1`"))
        
    size(specpos_data) == (2048,2048) || throw(DimensionMismatch(
        "size of `specpos_data` must be (2048,2048)"))
    
    size(specpos_weights) == (2048,2048) || throw(DimensionMismatch(
        "size of `specpos_weights` must be (2048,2048)"))
        
    nb_lenses == length(valid_lenses_map) || throw(DimensionMismatch(
        "`size(lenses_positions,2)` must be equal to length of `valid_lenses_map`"))
    
    res = FitResult(nb_lenses, nλ, λ0, wavelamp_order, specpos_order, size(init_box,2))
    
    valid_boxes = compute_lenses_boxes!(res, init_box, lenses_positions)
    
    progress = Progress(nb_lenses; showspeed=true)
    Threads.@threads for i in findall(valid_lenses_map .& valid_boxes)[12000]

        box = res.lenses_boxes[i]
        
        # Fit WaveLamp
        
        lens_wavelamp_data    = view(wavelamps_data,    box)
        lens_wavelamp_weights = view(wavelamps_weights, box)
        
        # if more than 50% bad pixels, we skip
        count(==(0), lens_wavelamp_weights) / length(lens_wavelamp_weights) ≥ 0.5 && begin
            next!(progress)
            continue
        end
        
        wavelamp_lkl = WaveLampLikelihood(
            λ0, wavelamps_λlasers, box, lens_wavelamp_data, lens_wavelamp_weights)
        
        wavelamp_fit_params ::Matrix{Float64} =
            [ wavelamps_fwhm_init                      ;; # first column
              lenses_positions[1,i]; wavelamps_cx_init ;; # second column
              lenses_positions[2,i]; wavelamps_cy_init  ] # third column
        
        try vmlmb!(wavelamp_lkl, wavelamp_fit_params
                  ; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
            # call one last time, to ensure that `last_amps` is produced from the final params
            wavelamp_lkl(wavelamp_fit_params)
        catch e
            @debug "Error on lenslet $i" exception=(e,catch_backtrace())
            next!(progress)
            continue
        end
        
        res.wavelamps_fits_fwhm[:,i]     .= wavelamp_fit_params[:,1]
        res.wavelamps_fits_cx[:,i]       .= wavelamp_fit_params[:,2]
        res.wavelamps_fits_cy[:,i]       .= wavelamp_fit_params[:,3]
        res.wavelamps_fits_background[i]  = wavelamp_lkl.last_background.x
        res.wavelamps_fits_amps[:,i]     .= wavelamp_lkl.last_amps
        
        (lens_centers_dists, lens_λvals) = compute_distance_map(
            box, λrange, λ0, res.wavelamps_fits_cx[:,i], res.wavelamps_fits_cy[:,i])
            
        res.wavelamps_centers_dists[box] .= lens_centers_dists
        res.wavelamps_λvals[box]         .= lens_λvals

        # Fit SpecPos

        lens_specpos_data    = view(specpos_data,    box)
        lens_specpos_weights = view(specpos_weights, box)

        # if more than 50% bad pixels, we skip
        count(==(0), lens_specpos_weights) / length(lens_specpos_weights) ≥ 0.5 && begin
            next!(progress)
            continue
        end

        specpos_fit_params = zeros(Float64, specpos_order + 1, 2)
        specpos_fit_params[1,1] = res.wavelamps_fits_cx[1,i]
        specpos_fit_params[1,2] = mean(res.wavelamps_fits_fwhm[:,i])

        specpos_lkl = SpecPosLikelihood(
            λ0, specpos_order, box, lens_specpos_data, lens_specpos_weights, lens_λvals)
            
        try vmlmb!(specpos_lkl, specpos_fit_params
                  ; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
            # call one last time, to ensure that `last_amps` is updated from final params
            specpos_lkl(specpos_fit_params)
        catch e
            @debug "Error on lenslet $i" exception=(e,catch_backtrace())
            next!(progress)
            continue
        end
        
        res.specpos_fits_cx[:,i]       .= specpos_fit_params[:,1]
        res.specpos_fits_cλ[:,i]       .= specpos_fit_params[:,2]
        res.specpos_fits_background[i]  = specpos_lkl.last_background.x
        res.sepcpos_fits_amps[:,i]     .= specpos_lkl.last_amps
        
        # if we are here the computation was complete
        res.computed_lenses_map[i] = true
        
        next!(progress)
    end
    finish!(progress)
   
    res ::FitResult
end

"""
    compute_lenses_boxes!(::FitResult, init_box::BoundingBox{Int}, lenses_positions) -> BitVector

Compute the box for each lenslet, from the given centers and using the given box shape.

Stores them in the given `FitResult`.

Returns a `BitVector` with `false` when the box is out of range of the data
"""
function compute_lenses_boxes!(
    res::FitResult, init_box::BoundingBox{Int}, lenses_positions::AbstractMatrix{Float64}
) ::BitVector

    valid_boxes = trues(size(lenses_positions,2))

    for (i, (x,y)) in enumerate(eachcol(lenses_positions))

        # shifting the std box to the given lens center
        boxR ::BoundingBox{Float64} = init_box + Point(x, y)
        
        # RoundNearestTiesUp to enforce box size (special case of `.5` floats)
        boxI ::BoundingBox{Int}     = round(Int, boxR, RoundNearestTiesUp)
        
        # failure if size is different
        size(boxI) == size(init_box) || error("incorrect box size")
        
        # register if box is out of data range
        (1 ≤ boxI.xmin < boxI.xmax ≤ 2048) && (1 ≤ boxI.ymin < boxI.ymax ≤ 2048) || begin
            valid_boxes[i] = false
        end
        
        res.lenses_boxes[i] = boxI
    end
    
    valid_boxes
end

"""
    compute_distance_map(::BoundingBox{Int}, wavelengthrange::AbstractRange, λ0, \
        cx::Vector{Float64}, cy::Vector{Float64}) -> NTuple{2,Matrix{Float64}}

Compute the distance to the gaussian center and the wavelength for each pixel of the box.

The given `cx` and `cy` are the coefficients that give the gaussian centers for each wavelength
of the `wavelengthrange`.
"""
function compute_distance_map(
    box::BoundingBox{Int}, wavelengthrange::AbstractRange, λ0,
    cx::Vector{Float64}, cy::Vector{Float64}) ::NTuple{2,Matrix{Float64}}
    
    (ax, ay) = axes(box)
    polynome_cx = polynome_with_reference(λ0, cx)
    polynome_cy = polynome_with_reference(λ0, cy)
    
    dist = fill(1e3, size(box))
    pixλ = fill(1e0, size(box))
    
    previous_index = 0
    for I in CartesianIndices(dist)
        previous_index = max(1, previous_index - 5)
        for (index,λ) in enumerate(wavelengthrange[previous_index:end])
            mx, my = polynome_cx(λ), polynome_cy(λ)
            rx = ax[I[1]] - mx
            ry = ay[I[2]] - my
            r = sign(rx) * sqrt(rx^2 + ry^2)
            if abs(r) < abs(dist[I[1],I[2]])
                dist[I[1],I[2]] = r;
                pixλ[I[1],I[2]] = λ;
            else
                previous_index += index - 1
                break
            end
        end
    end
    dist, pixλ
end


