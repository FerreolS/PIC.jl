"""
Compute the value at position sqrt(r) 1D centered Gaussian
* `x`:  squared sampled position
* `fwhm` : full-width at half maximum
"""
function GaussianModel2(x::T, fwhm::T) where {T<:Real}
    local fwhm2sigma =T(1 / (2 * sqrt(2 * log( 2))))
    return exp(-x / (2 * (fwhm * fwhm2sigma )^2));
end

GaussianModel2(tpl::Tuple{T,T}) where {T<:Real} = GaussianModel2(tpl[1], tpl[2])




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
      profileorder  ::Int = 2
    ) where T<:Real

    numberoflenslet = size(position,1)
    
    numberoflenslet == length(validlensmap) || throw(DimensionMismatch(
        "size of `position` incompatible with size of `validlensmap`"))
    
    length(λlaser) == length(fwhminit) || throw(DimensionMismatch(
        "size of `λlaser` incompatible with size of `fwhminit`"))
    
    # map of the computed lenses (i.e the non-catch-error ones)
    computedlensmap ::Vector{Bool} = falses(numberoflenslet)
    
    nλ = length(λlaser)
    order = nλ - 1
    λ0 = mean(λlaser) # reference
    refed_λs ::Vector{Float64}  = map(Base.Fix1(apply_reference, λ0), λlaser)
    poweredλs ::Matrix{Float64} = [ λ^o for o in 1:order, λ in refed_λs ]

    (dxmin, dxmax,dymin,dymax) = lensletsize
    
    wavelamplkltab = Vector{WaveLampLikelihood}(undef, numberoflenslet)
    specposlkltab  = Vector{SpecPosLikelihood}( undef, numberoflenslet)
    
    wavelamps_fits_cx = Matrix{Float64}(undef, order + 1, numberoflenslet)
    wavelamps_fits_cy = Matrix{Float64}(undef, order + 1, numberoflenslet)
    
    laserAmplitude = Matrix{Float64}(undef, nλ, numberoflenslet)
    lampAmplitude  = Matrix{Float64}(undef, 41, numberoflenslet)
    laserfwhm = Matrix{Float64}(undef, nλ, numberoflenslet)
    laserdist = Matrix{Float64}(undef, 2048, 2048)
    λMap = Matrix{Float64}(undef, 2048, 2048)
    
    progress = Progress(numberoflenslet; showspeed=true)
    Threads.@threads for i in findall(validlensmap)[1:100:end]

        box = BoundingBox(position[i,1]-dxmin, position[i,1]+dxmax,
                                 position[i,2]-dymin, position[i,2]+dymax)
        # we use RoundNearestTiesUp to enforce box size (special case of `.5` floats)
        box = round(Int, box, RoundNearestTiesUp);

        # Fit spectral law
        
        laserDataView   = view(laserdata,    box)
        laserWeightView = view(laserweights, box)
        wavelamplkltab[i] = WaveLampLikelihood(
            order, nλ, poweredλs, box, laserDataView, laserWeightView)
        
        fit_params = [ fwhminit..., position[i,1], cxinit..., position[i,2], cyinit... ]
        
        try vmlmb!(wavelamplkltab[i], fit_params
                  ; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug "Error on lenslet $i" exception=(e, catch_backtrace())
            continue
        end
        # call one last time, to ensure that `last_amp` is produced from the final `fit_params`
        wavelamplkltab[i](fit_params)
        fit_amp = wavelamplkltab[i].last_amp
        (fit_fwhm, fit_cx, fit_cy) = decode_WaveLampLikelihood_input(nλ, order, fit_params)
        
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

        profilecoefs = zeros(Float64, 2, profileorder + 1)
        profilecoefs[2,:] .= [2.3, 2.5, 2.9] # maximum(fwhm)
        profilecoefs[1,1] = fit_cx[1]

        pmodel = ProfileModel(λ0, profileorder, profilecoefs[1,:], profilecoefs[2,:])
        specposlkltab[i] = SpecPosLikelihood(pmodel,box, lampDataView, lampWeightView, lensletpixλ)
        try
            vmlmb(specposlkltab[i], profilecoefs
                  ; verb=false,ftol = (0.0,1e-8),maxeval=500,autodiff=true);
        catch e
            @debug "Error on lenslet  $i" exception=(e, catch_backtrace())
            continue
        end
        #FIXME: are we sure that the last run contains the best values ?

        profile = @. GaussianModel2(specposlkltab[i].pmodel(lensletpixλ,($(axes(box,1)))))
        profile = profile ./ sum(profile,dims=1)
        try
            lampAmplitude[:,i] .= updateAmplitudeAndBackground(profile,lampDataView,lampWeightView)
        catch e
            @debug "Error on lenslet  $i" exception=(e, catch_backtrace())
            continue
        end
        
        # if we are here the computation was complete
        computedlensmap[i] = true
        
        next!(progress)
    end
    finish!(progress)
   
    (λ0, wavelamplkltab, wavelamps_fits_cx, wavelamps_fits_cy, specposlkltab, laserAmplitude,
     lampAmplitude, laserfwhm,laserdist, λMap, computedlensmap)
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

#function updateAmplitude(profile,data::Matrix{T},weight::Matrix{T}) where T<:AbstractFloat
#    A = similar(data)
#    b = similar(data)
#
#    @. b = profile * data * weight
#    @. A = profile^2 * weight
#    A = sum(A,dims=1)
#    b = sum(b,dims=1)
#    zA = (A .== T(0)).||(b.<=T(0))
#    if any(zA)
#        A[zA] .=1
#        b[zA] .=0
#    end
#
#    return b ./ A
#end

function updateAmplitudeAndBackground(profile,data::MA,weight::MB) where {T<:AbstractFloat,MA<:AbstractMatrix{T},MB<:AbstractMatrix{T}}

    c = @. profile *  weight
    b = @. profile * data * weight
    a = @. profile^2 * weight
    a = sum(a,dims=1)[:]
    b = sum(b,dims=1)[:]
    c = sum(c,dims=1)[:]
    za = (a .== T(0)).||(b.<=T(0))
    if any(za)
        a[za] .=T(1)
        b[za] .=T(0)
        c[za] .=T(0)
    end


    N = length(a)
    A = Matrix{T}(undef,N+1,N+1)
    A[1,1] = sum(weight)
    A[1,2:end] .= A[2:end,1] .= c[:]
    A[2:end,2:end] .= diagm(a)

    b =  vcat(sum(data .* weight),b[:])

    return  inv(A)*b
end

