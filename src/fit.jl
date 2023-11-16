function fit_wavelamps_specpos(
    lenses_positions ::Matrix{Float64},
    
    wavelamps_λlasers ::Vector{Float64},
    wavelamps_data    ::Matrix{T},
    wavelamps_weights ::Matrix{T},
    wavelamps_fwhm_init ::Vector{Float64},
    wavelamps_cx_init   ::Vector{Float64},
    wavelamps_cy_init   ::Vector{Float64},
    
    specpos_data      ::Matrix{T},
    specpos_weights   ::Matrix{T},
    
    ; specpos_order ::Int = 2,
      λ0 ::Float64 = mean(wavelamps_λlasers),
      λrange ::AbstractRange{Float64} = IFS_λRANGE,
      box_frame ::NTuple{4,Int} = BOX_FRAME,
      valid_lenses_map ::AbstractVector{Bool} = trues(size(lenses_positions,1))
) where {T<:Real}

    nb_lenses = size(lenses_positions,1)
    
    nb_lenses == length(valid_lenses_map) || throw(DimensionMismatch(
        "size of `position` incompatible with size of `valid_lenses_map`"))
    
    length(wavelamps_λlasers) == length(wavelamps_fwhm_init) || throw(DimensionMismatch(
        "size of `wavelamps_λlasers` incompatible with size of `fwhminit`"))
    
    # map of the computed lenses (i.e the non-catch-error ones)
    computedlensmap ::Vector{Bool} = falses(nb_lenses)
    
    nλ = length(wavelamps_λlasers)
    wavelamp_order = nλ - 1

    (dxmin, dxmax, dymin, dymax) = box_frame
    box_size = (dxmax + dxmin + 1, dymin + dymax + 1)
    
    lensboxs = Vector{BoundingBox{Int}}(undef, nb_lenses)
    for (i, (x,y)) in enumerate(eachrow(lenses_positions))
        box = BoundingBox(x - dxmin, x + dxmax, y - dymin, y + dymax)
        # RoundNearestTiesUp to enforce box size (special case of `.5` floats)
        lensboxs[i] = round(Int, box, RoundNearestTiesUp)
        size(lensboxs[i]) == box_size || error("incorrect box size")
    end
    
    wavelamps_fits_cx = fill(NaN64, wavelamp_order + 1, nb_lenses)
    wavelamps_fits_cy = fill(NaN64, wavelamp_order + 1, nb_lenses)
    
    specpos_fits_cx = fill(NaN64, specpos_order + 1, nb_lenses)
    specpos_fits_cλ = fill(NaN64, specpos_order + 1, nb_lenses)
    
    laserAmplitude = Matrix{Float64}(undef, nλ, nb_lenses)
    lampBackground = Vector{Float64}(undef, nb_lenses)
    lampAmplitude  = Matrix{Float64}(undef, box_size[2], nb_lenses)
    laserfwhm = Matrix{Float64}(undef, nλ, nb_lenses)
    laserdist = Matrix{Float64}(undef, 2048, 2048)
    λMap = Matrix{Float64}(undef, 2048, 2048)
    
    progress = Progress(nb_lenses; showspeed=true)
    Threads.@threads for i in findall(valid_lenses_map)#[1:100:end]

        box = lensboxs[i]

        # Fit spectral law
        
        lens_wavelamp_data    = view(wavelamps_data,    box)
        lens_wavelamp_weights = view(wavelamps_weights, box)
        
        wavelamp_lkl = WaveLampLikelihood(
            λ0, wavelamps_λlasers, box, lens_wavelamp_data, lens_wavelamp_weights)
        
        wavelamp_fit_params ::Matrix{Float64} =
            [ wavelamps_fwhm_init                      ;; # first column
              lenses_positions[i,1]; wavelamps_cx_init ;; # second column
              lenses_positions[i,2]; wavelamps_cy_init  ] # third column
        
        try vmlmb!(wavelamp_lkl, wavelamp_fit_params
                  ; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug "Error on lenslet $i" exception=(e,catch_backtrace())
            continue
        end
        # call one last time, to ensure that `last_amp` is produced from the final `fit_params`
        wavelamp_lkl(wavelamp_fit_params)
        
        laserfwhm[:,i]         .= wavelamp_fit_params[:,1]
        wavelamps_fits_cx[:,i] .= wavelamp_fit_params[:,2]
        wavelamps_fits_cy[:,i] .= wavelamp_fit_params[:,3]
        laserAmplitude[:,i]    .= wavelamp_lkl.last_amp
        
        (lensletdist, lensletpixλ) = compute_distance_map(
            box, λrange, λ0, wavelamps_fits_cx[:,i], wavelamps_fits_cy[:,i])
            
        laserdist[box] .= lensletdist
        λMap[box] .= lensletpixλ

        # Fit profile

        lens_specpos_data    = view(specpos_data,    box)
        lens_specpos_weights = view(specpos_weights, box)

        specpos_fit_params = zeros(Float64, specpos_order + 1, 2)
        specpos_fit_params[1,1] = wavelamps_fits_cx[1,i]
        specpos_fit_params[1,2] = mean(laserfwhm[:,i])

        specpos_lkl = SpecPosLikelihood(
            λ0, specpos_order, box, lens_specpos_data, lens_specpos_weights, lensletpixλ)
        try
            vmlmb!(specpos_lkl, specpos_fit_params
                  ; verb=false, ftol =(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug "Error on lenslet  $i" exception=(e,catch_backtrace())
            continue
        end
        # call one last time, to ensure that `last_amps` is updated from the final `profilecoefs`
        specpos_lkl(specpos_fit_params)
        
        specpos_fits_cx[:,i] .= specpos_fit_params[:,1]
        specpos_fits_cλ[:,i] .= specpos_fit_params[:,2]
        lampAmplitude[:,i]   .= specpos_lkl.last_amps
        lampBackground[i]     = specpos_lkl.last_background.x
        
        # if we are here the computation was complete
        computedlensmap[i] = true
        
        next!(progress)
    end
    finish!(progress)
   
    (lensboxs, λ0,
     wavelamp_order, wavelamps_fits_cx, wavelamps_fits_cy,
     specpos_order, specpos_fits_cx, specpos_fits_cλ,
     laserAmplitude, lampBackground, lampAmplitude, laserfwhm,laserdist, λMap, computedlensmap)
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

