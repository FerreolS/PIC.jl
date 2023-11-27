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
      box_frame ::NTuple{4,Int} = BOX_FRAME,
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
    
    # map of the computed lenses (i.e the non-catch-error ones)
    computed_lenses_map ::Vector{Bool} = falses(nb_lenses)

    (dxmin, dxmax, dymin, dymax) = box_frame
    box_size = (dxmax + dxmin + 1, dymin + dymax + 1)
    box_size[1] ≥ 1 && box_size[2] ≥ 1 || throw(DimensionMismatch())
    
    lenses_boxes = Vector{BoundingBox{Int}}(undef, nb_lenses)
    for (i, (x,y)) in enumerate(eachcol(lenses_positions))
        box = BoundingBox(x - dxmin, x + dxmax, y - dymin, y + dymax)
        # RoundNearestTiesUp to enforce box size (special case of `.5` floats)
        lenses_boxes[i] = round(Int, box, RoundNearestTiesUp)
        size(lenses_boxes[i]) == box_size || error("incorrect box size")
    end
    
    wavelamps_fits_cx         = fill(NaN64, wavelamp_order + 1, nb_lenses)
    wavelamps_fits_cy         = fill(NaN64, wavelamp_order + 1, nb_lenses)
    wavelamps_fits_fwhm       = fill(NaN64, nλ, nb_lenses)
    wavelamps_fits_background = fill(NaN64, nb_lenses)
    wavelamps_fits_amp        = fill(NaN64, nλ, nb_lenses)
    
    wavelamps_centers_dists = fill(NaN64, 2048, 2048)
    wavelamps_λvals         = fill(NaN64, 2048, 2048)
    
    specpos_fits_cx         = fill(NaN64, specpos_order + 1, nb_lenses)
    specpos_fits_cλ         = fill(NaN64, specpos_order + 1, nb_lenses)
    specpos_fits_background = fill(NaN64, nb_lenses)
    sepcpos_fits_amps       = fill(NaN64, box_size[2], nb_lenses)
    
    progress = Progress(nb_lenses; showspeed=true)
    Threads.@threads for i in findall(valid_lenses_map)

        box = lenses_boxes[i]
        
        # skip box if out of data range
        (1 ≤ box.xmin < box.xmax ≤ 2048) && (1 ≤ box.ymin < box.ymax ≤ 2048) || begin
            next!(progress)
            continue
        end
        
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
            # call one last time, to ensure that `last_amp` is produced from the final params
            wavelamp_lkl(wavelamp_fit_params)
        catch e
            @debug "Error on lenslet $i" exception=(e,catch_backtrace())
            next!(progress)
            continue
        end
        
        wavelamps_fits_fwhm[:,i]     .= wavelamp_fit_params[:,1]
        wavelamps_fits_cx[:,i]       .= wavelamp_fit_params[:,2]
        wavelamps_fits_cy[:,i]       .= wavelamp_fit_params[:,3]
        wavelamps_fits_background[i]  = wavelamp_lkl.last_background.x
        wavelamps_fits_amp[:,i]      .= wavelamp_lkl.last_amp
        
        (lens_centers_dists, lens_λvals) = compute_distance_map(
            box, λrange, λ0, wavelamps_fits_cx[:,i], wavelamps_fits_cy[:,i])
            
        wavelamps_centers_dists[box] .= lens_centers_dists
        wavelamps_λvals[box]         .= lens_λvals

        # Fit SpecPos

        lens_specpos_data    = view(specpos_data,    box)
        lens_specpos_weights = view(specpos_weights, box)

        # if more than 50% bad pixels, we skip
        count(==(0), lens_specpos_weights) / length(lens_specpos_weights) ≥ 0.5 && begin
            next!(progress)
            continue
        end

        specpos_fit_params = zeros(Float64, specpos_order + 1, 2)
        specpos_fit_params[1,1] = wavelamps_fits_cx[1,i]
        specpos_fit_params[1,2] = mean(wavelamps_fits_fwhm[:,i])

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
        
        specpos_fits_cx[:,i]       .= specpos_fit_params[:,1]
        specpos_fits_cλ[:,i]       .= specpos_fit_params[:,2]
        specpos_fits_background[i]  = specpos_lkl.last_background.x
        sepcpos_fits_amps[:,i]     .= specpos_lkl.last_amps
        
        # if we are here the computation was complete
        computed_lenses_map[i] = true
        
        next!(progress)
    end
    finish!(progress)
   
    FitResult(
        lenses_boxes, λ0,
        wavelamp_order, wavelamps_fits_cx, wavelamps_fits_cy,
        wavelamps_fits_fwhm, wavelamps_fits_background, wavelamps_fits_amp,
        wavelamps_centers_dists, wavelamps_λvals,
        specpos_order, specpos_fits_cx, specpos_fits_cλ,
        specpos_fits_background, sepcpos_fits_amps,
        computed_lenses_map)
end


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

struct FitResult

    lenses_boxes ::Vector{BoundingBox{Int}}
    λ0       ::Float64

    wavelamp_order      ::Int
    wavelamps_fits_cx   ::Matrix{Float64}
    wavelamps_fits_cy   ::Matrix{Float64}
    wavelamps_fits_fwhm ::Matrix{Float64}
    wavelamps_fits_background ::Vector{Float64}
    wavelamps_fits_amps ::Matrix{Float64}

    wavelamps_centers_dists ::Matrix{Float64}
    wavelamps_λval          ::Matrix{Float64}

    specpos_order           ::Int
    specpos_fits_cx         ::Matrix{Float64}
    specpos_fits_cλ         ::Matrix{Float64}
    specpos_fits_background ::Vector{Float64}
    sepcpos_fits_amps       ::Matrix{Float64}

    computed_lenses_map ::Vector{Bool}
end
