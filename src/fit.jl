




function fitSpectralLawAndProfile(laserdata::Matrix{T},
    laserweights ::Matrix{T},
    lampdata    ::Matrix{T},
    lampweights ::Matrix{T},
    λlaser      ::Vector{Float64},
    lensletsize ::NTuple{4,Int},
    position ::Matrix{Float64},
    cxinit   ::Vector{Float64},
    cyinit   ::Vector{Float64},
    fwhminit ::Vector{Float64},
    wavelengthrange ::AbstractVector{Float64}
    ; validlensmap  ::AbstractVector{Bool} = trues(size(position, 1)),
      specpos_order  ::Int = 2
    ) where T<:Real

    numberoflenslet = size(position,1)
    
    numberoflenslet == length(validlensmap) || throw(DimensionMismatch(
        "size of `position` incompatible with size of `validlensmap`"))
    
    length(λlaser) == length(fwhminit) || throw(DimensionMismatch(
        "size of `λlaser` incompatible with size of `fwhminit`"))
    
    # map of the computed lenses (i.e the non-catch-error ones)
    computedlensmap ::Vector{Bool} = falses(numberoflenslet)
    
    nλ = length(λlaser)
    wavelamp_order = nλ - 1
    λ0 = mean(λlaser) # reference
    refed_λs ::Vector{Float64}  = map(Base.Fix1(apply_reference, λ0), λlaser)
    poweredλs ::Matrix{Float64} = [ λ^o for o in 1:wavelamp_order, λ in refed_λs ]

    (dxmin, dxmax,dymin,dymax) = lensletsize
    
    lensboxs = Vector{BoundingBox{Int}}(undef, numberoflenslet)
    for (i, (x,y)) in enumerate(eachrow(position))
        box = BoundingBox(x - dxmin, x + dxmax, y - dymin, y + dymax)
        # RoundNearestTiesUp to enforce box size (special case of `.5` floats)
        lensboxs[i] = round(Int, box, RoundNearestTiesUp)
    end
    
    wavelamplkltab = Vector{WaveLampLikelihood}(undef, numberoflenslet)
    specposlkltab  = Vector{SpecPosLikelihood}( undef, numberoflenslet)
    
    wavelamps_fits_cx = fill(NaN64, wavelamp_order + 1, numberoflenslet)
    wavelamps_fits_cy = fill(NaN64, wavelamp_order + 1, numberoflenslet)
    
    specpos_fits_cx = fill(NaN64, specpos_order + 1, numberoflenslet)
    specpos_fits_cλ = fill(NaN64, specpos_order + 1, numberoflenslet)
    
    laserAmplitude = Matrix{Float64}(undef, nλ, numberoflenslet)
    lampAmplitude  = Matrix{Float64}(undef, 41, numberoflenslet)
    laserfwhm = Matrix{Float64}(undef, nλ, numberoflenslet)
    laserdist = Matrix{Float64}(undef, 2048, 2048)
    λMap = Matrix{Float64}(undef, 2048, 2048)
    
    progress = Progress(numberoflenslet; showspeed=true)
    Threads.@threads for i in findall(validlensmap)[1:100:end]

        box = lensboxs[i]

        # Fit spectral law
        
        laserDataView   = view(laserdata,    box)
        laserWeightView = view(laserweights, box)
        wavelamplkltab[i] = WaveLampLikelihood(
            wavelamp_order, nλ, poweredλs, box, laserDataView, laserWeightView)
        
        fit_params = [ fwhminit..., position[i,1], cxinit..., position[i,2], cyinit... ]
        
        try vmlmb!(wavelamplkltab[i], fit_params
                  ; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug "Error on lenslet $i" exception=(e,catch_backtrace())
            continue
        end
        # call one last time, to ensure that `last_amp` is produced from the final `fit_params`
        wavelamplkltab[i](fit_params)
        fit_amp = wavelamplkltab[i].last_amp
        (fit_fwhm, fit_cx, fit_cy) = decode_WaveLampLikelihood_input(nλ,wavelamp_order, fit_params)
        
        wavelamps_fits_cx[:,i] .= fit_cx
        wavelamps_fits_cy[:,i] .= fit_cy
        
        laserAmplitude[:,i] .= fit_amp
        laserfwhm[:,i] .= fit_fwhm
        
        (lensletdist, lensletpixλ) = compute_distance_map(box, wavelengthrange, λ0, fit_cx, fit_cy)
            
        laserdist[box] .= lensletdist
        λMap[box] .= lensletpixλ

        # Fit profile

        lampDataView = view(lampdata, box);
        lampWeightView = view(lampweights, box);

        profilecoefs = zeros(Float64, specpos_order + 1, 2)
        profilecoefs[1,1] = fit_cx[1]
        profilecoefs[1,2] = mean(fit_fwhm)

        specposlkltab[i] = SpecPosLikelihood(λ0, specpos_order, box, lampDataView, lampWeightView, lensletpixλ)
        try
            vmlmb!(specposlkltab[i], profilecoefs
                  ; verb=false,ftol = (0.0,1e-8),maxeval=500,autodiff=true);
        catch e
            @debug "Error on lenslet  $i" exception=(e,catch_backtrace())
            continue
        end
        # call one last time, to ensure that `last_amps` is updated from the final `profilecoefs`
        specposlkltab[i](profilecoefs)
        (fit_cx, fit_cλ) = decode_SpecPosLikelihood_input(profilecoefs)
        fit_amps = specposlkltab[i].last_amps
        
        specpos_fits_cx[:,i] .= fit_cx
        specpos_fits_cλ[:,i] .= fit_cλ
        lampAmplitude[:,i]   .= fit_amps
        
        # if we are here the computation was complete
        computedlensmap[i] = true
        
        next!(progress)
    end
    finish!(progress)
   
    (lensboxs, λ0,
     wavelamp_order, wavelamps_fits_cx, wavelamps_fits_cy,
     specpos_order, specpos_fits_cx, specpos_fits_cλ,
     laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap, computedlensmap)
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

