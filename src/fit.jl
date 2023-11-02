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


 """
        updateAmplitude(nλ,m,d,W)

    return the `nλ` amplitudes `a` according the the model `m`, the data and the precision `W`
    such that
    `a = argmin_a || a*m - D||^2_W`
    where
    * `nλ` : is the number of spots in the model
    * `m`:  is the model composed of `nλ` images of spots
    * `d`:  is the data
    * `W`: is the precision (inverse variance) of the data
"""
function updateAmplitude(N::Int,spots::AbstractArray{T},data::AbstractArray{T},weight::AbstractArray{T}) where T<:Real
    A = @MMatrix zeros(Float64,N,N)
    b = @MVector zeros(Float64,N)
    mw = similar(spots);
    @inbounds for index=1:N
        mw[:,:,index] .=  spots[:,:,index].* weight ;
        b[index] = sum(mw[:,:,index].* data );
        A[index,index] = sum(mw[:,:,index].* spots[:,:,index]);
        for i=1:index-1
            A[i,index] = A[index,i] = sum(mw[:,:,index].* spots[:,:,i])
        end
    end
    return inv(A)*b
end

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
    λ0 = mean(λlaser) # reference

    (dxmin, dxmax,dymin,dymax) = lensletsize
    lenslettab = Vector{LensletModel}(undef,numberoflenslet)
    laserAmplitude = Matrix{Float64}(undef, nλ, numberoflenslet)
    lampAmplitude  = Matrix{Float64}(undef, 41, numberoflenslet)
    laserfwhm = Matrix{Float64}(undef, nλ, numberoflenslet)
    laserdist = Matrix{Float64}(undef, 2048, 2048)
    λMap = Matrix{Float64}(undef, 2048, 2048)
    
    progress = Progress(numberoflenslet; showspeed=true)
    Threads.@threads for i in findall(validlensmap)[1:100:end]

        lensletbox = BoundingBox(position[i,1]-dxmin, position[i,1]+dxmax,
                                 position[i,2]-dymin, position[i,2]+dymax)
        # we use RoundNearestTiesUp to enforce box size (special case of `.5` floats)
        lensletbox = round(Int, lensletbox, RoundNearestTiesUp);

        lenslettab[i] = LensletModel(λ0,nλ-1, profileorder,lensletbox);

        # Fit spectral law
        xinit = [ fwhminit..., position[i,1], cxinit..., position[i,2], cyinit... ]
        laserDataView   = view(laserdata,    lensletbox)
        laserWeightView = view(laserweights, lensletbox)
        spectrallkl = LikelihoodDisp(lenslettab[i], λlaser, laserDataView, laserWeightView)
        xopt =
            try vmlmb(spectrallkl, xinit; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
            catch e
                @debug "Error on lenslet  $i" exception=(e, catch_backtrace())
                continue
            end
        fwhm = xopt[1:nλ]
        laserAmplitude[:,i] = spectrallkl.amplitude;
        laserfwhm[:,i] = fwhm
        (lensletdist, lensletpixλ) = distanceMap(wavelengthrange, lenslettab[i])
        laserdist[lensletbox] .= lensletdist
        λMap[lensletbox] .= lensletpixλ

        # Fit profile

        lampDataView = view(lampdata, lensletbox);
        lampWeightView = view(lampweights,lensletbox);

        profilecoefs = zeros(Float64, 2, profileorder+1)
        profilecoefs[2,:] .= [2.3, 2.5, 2.9] # maximum(fwhm)
        profilecoefs[1,1] = lenslettab[i].dmodel.cx[1]

        pmodel = ProfileModel(λ0, profileorder, profilecoefs[1,:], profilecoefs[2,:])
        profilelkl = LikelihoodProfile(pmodel, lampDataView, lampWeightView, lensletpixλ, lensletbox)
        costpr(x::Matrix{Float64}) = profilelkl(x);
        try
            vmlmb!(costpr, profilecoefs; verb=false,ftol = (0.0,1e-8),maxeval=500,autodiff=true);
        catch e
            @debug "Error on lenslet  $i" exception=(e, catch_backtrace())
            continue
        end
        pmodel = ProfileModel(λ0, profileorder, profilecoefs[1,:], profilecoefs[2,:])
        lenslettab[i] = LensletModel(lensletbox, lenslettab[i].dmodel, pmodel)

        profile = @. GaussianModel2(pmodel(lensletpixλ,($(axes(lensletbox,1)))))
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
   
    (lenslettab, laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap, computedlensmap)
end


function distanceMap(wavelengthrange::AbstractArray{Float64,1},
                    lenslet::LensletModel
                    )
    bbox = lenslet.bbox;
    dist = ones(Float64,size(round(bbox))).*1000;
    pixλ = ones(Float64,size(round(bbox)));
    (ax,ay) = axes(bbox)
    previous_index = 0;
    for I in CartesianIndices(dist)
        previous_index = max(1,previous_index-5);
        for  (index,λ) in enumerate(wavelengthrange[previous_index:end])
            (mx, my)  = lenslet.dmodel(λ)
            rx = ax[I[1]]-mx;
            ry = ay[I[2]]-my;
            r = sign(rx) * sqrt(rx^2 + ry^2);
            if abs(r) < abs(dist[I[1],I[2]])
                dist[I[1],I[2]] = r;
                pixλ[I[1],I[2]] = λ;
            else
                previous_index = previous_index + index-1;
                break
            end
        end
    end
    return (dist,pixλ)
end

function updateAmplitude(profile,data::Matrix{T},weight::Matrix{T}) where T<:AbstractFloat
    A = similar(data)
    b = similar(data)

    @. b = profile * data * weight
    @. A = profile^2 * weight
    A = sum(A,dims=1)
    b = sum(b,dims=1)
    zA = (A .== T(0)).||(b.<=T(0))
    if any(zA)
        A[zA] .=1
        b[zA] .=0
    end

    return b ./ A
end

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

