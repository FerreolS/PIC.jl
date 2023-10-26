"""
    GaussianModel(A,fwhm,x,y)

Compute the value  at position (x,y) of a 2D centered Gaussian
* `fwhm` : full-width at half maximum
* `A` : amplitude at (x,y) = (0,0)
* `x`, `y`: sampled postions
"""
function GaussianModel(A::Float64, fwhm::Float64, x::Float64, y::Float64)
    local fwhm2sigma = 1 / (2 * sqrt(2 * log(2.)))::Float64
    return A * exp(-(x^2 + y^2) / (2 * (fwhm * fwhm2sigma )^2));
end

"""
    GaussianModel(A::Float64, fwhm::Float64, x::AbstractArray)

Compute the value at position x 1D centered Gaussian
* `A` : amplitude at x = 0
* `fwhm` : full-width at half maximum
* `x`: array of the sampled position
"""
function GaussianModel(A::Float64, fwhm::Float64, x::AbstractArray)
    local fwhm2sigma = 1 / (2 * sqrt(2 * log(2.)))::Float64
    return A .* exp.(.-(x.^2) ./ (2 * (fwhm * fwhm2sigma )^2));
end

"""
    GaussianModel2(A::Float64, fwhm::Float64, x::AbstractArray)

Compute the value at position sqrt(r) 1D centered Gaussian
* `A` : amplitude at x = 0
* `fwhm` : full-width at half maximum
* `x`: array of the squared sampled position

Equivalent to `GaussianModel(A, fwhm, sqrt.(x))`
"""
function GaussianModel2(A::Float64, fwhm::Float64, x::AbstractArray)
    local fwhm2sigma = Float64(1) / (2 * sqrt(2 * log(2.)))
    return A .* exp.(-x ./ (2 * (fwhm * fwhm2sigma )^2));
end

"""
    GaussianModel2(A::Float64, fwhm::Float64, x::Real)

Compute the value at position sqrt(x) 1D centered Gaussian
* `A` : amplitude at x = 0
* `fwhm` : full-width at half maximum
* `x`: sampled position

Equivalent to `GaussianModel(A, fwhm, sqrt(x))`
"""
function GaussianModel2(A::T, fwhm::T, x::T) where (T<:Real)
    local fwhm2sigma = T(1 / (2 * sqrt(2 * log(2.))))
    return A * exp(-x / (2 * (fwhm * fwhm2sigma )^2));
end

"""
    GaussianModel2(fwhm::Float64, x::AbstractArray)

Compute the value at position sqrt(r) 1D centered Gaussian
* `fwhm` : full-width at half maximum
* `x`: array of the squared sampled position

Equivalent to `GaussianModel(1.,fwhm, sqrt.(x))`
"""
function GaussianModel2(fwhm::T, x::AbstractArray{T})where (T<:Real)
    local fwhm2sigma = T(1 / (2 * sqrt(2 * log(2.))))
    return exp.(-x ./ (2 * (fwhm * fwhm2sigma )^2));
end

"""
    GaussianModel2(fwhm::Float64, x::AbstractArray)

Compute the value at position sqrt(r) 1D centered Gaussian
* `fwhm` : full-width at half maximum
* `x`:  squared sampled position

Equivalent to `GaussianModel(1.,fwhm, sqrt(x))`
"""
function GaussianModel2(fwhm::T, x::T) where (T<:Real)
    local fwhm2sigma =T(1 / (2 * sqrt(2 * log( 2))))
    return exp(-x / (2 * (fwhm * fwhm2sigma )^2));
end

GaussianModel2(tpl::Tuple{T, T}) where (T<:Real) = GaussianModel2(tpl...)

"""
    GaussianModel2!(ret::AbstractArray{T},fwhm::Float64, x::AbstractArray)

Compute inplace the value at position sqrt(r) 1D centered Gaussian
* `ret` : output array
* `fwhm` : full-width at half maximum
* `x`:  squared sampled position

Equivalent to `GaussianModel(1.,fwhm, sqrt(x))`
"""
function GaussianModel2!(ret::AbstractArray{T},fwhm, x::AbstractArray{T}) where (T<:Real)
        @. ret = exp(-x / T(2 * (fwhm * 1 / (2 * sqrt(2 * log(2.))) )^2));
        nothing
end

"""
    GaussianSpotsCost(data::Array{Float64,2}, weight::Array{Float64,2}, lmodel::LensletModel,  laser::LaserModel,A::Array{Float64,1}, fwhm::Array{Float64,1}, C::Array{Float64,2})

Compute a weighted quadratic cost of a lenslet model :
cost = weight .* || data - model ||^2
* `lmodel`:  model of the lenslet
* `A` : 1D array containing the amplitude  of all Gaussian spots
* `fwhm` : 1D array containing the fwhm  of all Gaussian spots
* `C` : 2D array containing the chromatic law coefficients.
"""
function GaussianSpotsCost(data::Array{Float64,2}, weight::Array{Float64,2}, lmodel::LensletModel,  laser::LaserModel,A::Array{Float64,1}, fwhm::Array{Float64,1}, C::Array{Float64,2})
    UpdateDispModel(lmodel.dmodel, C);
    UpdateLaserModel(laser,A,fwhm);
    s = 0.;
    for I in CartesianIndices(lmodel.bbox)
        spotsmodel = 0;
        for (index, λ) in enumerate(laser.λlaser)
            (mx, my)  = lmodel.dmodel(λ);
            spotsmodel += GaussianModel(laser.amplitude[index], laser.fwhm[index], I[1] - mx, I[2] - my)
        end
        s += weight[I] * ( data[I] - spotsmodel)^2;
    end
    return s;
end


"""
    GaussianSpotsModel(lmodel::LensletModel,laser::LaserModel, A::Array{Float64,1}, fwhm::Array{Float64,1}, C::Array{Float64,2})

Build the model of a lenslet
* `lmodel`:  model of the lenslet
* `A` : 1D array containing the amplitude  of all Gaussian spots
* `fwhm` : 1D array containing the fwhm  of all Gaussian spots
* `C` : 2D array containing the chromatic law coefficients.
"""
function GaussianSpotsModel(lmodel::LensletModel,laser::LaserModel, A::Array{Float64,1}, fwhm::Array{Float64,1}, C::Array{Float64,2})
    UpdateDispModel(lmodel.dmodel, C);
    UpdateLaserModel(laser,A,fwhm);
    bbox = lmodel.dmodel;
    model = zeros(Float64,round(bbox).xmax-round(bbox).xmin+1,round(bbox).ymax-round(bbox).ymin+1);
    t = Zygote.Buffer(model);
    t[:] = model[:];
    for I in CartesianIndices(bbox)
        for (index, λ) in enumerate(laser.λlaser)
            (mx, my)  = lmodel.dmodel(λ);
            t[I[1],I[2]] += GaussianModel(laser.amplitude[index], laser.fwhm[index], I[1] - mx, I[2] - my)
        end
    end
    model =  Zygote.copy(t);
end




 SimpleGauss(x,center::Float64,fwhm::Float64) = exp(-(x-center)^2 / (2 * (fwhm * Float64(1) / (2 * sqrt(2 * log(2.))) )^2));


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

"""
        (lenslettab, atab, fwhmtab,ctab) = fitSpectralLaw(laserData,weights,λlaser,lensletsize,cx0,cy0,cinit,fwhminit;validlenslets=true);

    fits the spectral of all lenslet identified as valid in the `valid` vector.
    * `laserData` : is the laser calibration data
    * `weight`:  is the precision (inverse variance) of the data
    * `valid`:  a boolean vector indicating valid lenslet
    * `λlaser`: is the vector of the wavelength of the lasers
    * `lensletsize` :  is a 4 Int tuple giving the size of lenslet on the detector
    * `position`: is the a priori position of the lenslets
    * `cxinit`: is the initial vector of the spectral law coefficients along x
    * `cxinit`: is the initial vector of the spectral law coefficients along y
    * `fwhminit`: is the initial vector of full width at half maximum of the spots
    * `validlenslets`: is an optionnal vector indicating the already known invalid lenslets
"""
function fitSpectralLaw(laserdata::Matrix{T},
                        weights::Matrix{T},
                        λlaser::Array{Float64,1},
                        lensletsize::NTuple{4, Int},
                        position::Matrix{Float64},
                        cxinit::Vector{Float64},
                        cyinit::Vector{Float64},
                        fwhminit::Array{Float64,1};
                        validlenslets::AbstractArray{Bool,1}=[true]
                        ) where T<:Real

    numberoflenslet = size(position)[1]
    if length(validlenslets)==1
        validlenslets = true(numberoflenslet)
    else
        numberoflenslet = min(length(validlenslets) ,numberoflenslet)
    end
    profileorder= 2
    nλ = length(λlaser)
    λ0 = mean(λlaser)# reference
    @assert length(fwhminit) == nλ

    (dxmin, dxmax,dymin,dymax) = lensletsize
    lenslettab = Array{Union{LensletModel,Missing}}(missing,numberoflenslet);
    atab = Array{Union{Float64,Missing}}(missing,nλ,numberoflenslet);
    fwhmtab = Array{Union{Float64,Missing}}(missing,nλ,numberoflenslet);
    ctab = Array{Union{Float64,Missing}}(missing,2,nλ,numberoflenslet);
    p = Progress(numberoflenslet; showspeed=true)
    Threads.@threads for i in findall(validlenslets)
        lensletbox = round(Int, BoundingBox(position[i,1]-dxmin, position[i,1]+dxmax, position[i,2]-dymin, position[i,2]+dymax));

        lenslettab[i] = LensletModel(λ0,nλ-1,profileorder, lensletbox);
        Cinit= [ [position[i,1] cxinit...]; [position[i,2] cyinit...] ];
        xinit = vcat([fwhminit[:],Cinit[:]]...);
        laserDataView = view(laserdata, lensletbox);
        weightView = view(weights,lensletbox);
        lkl = LikelihoodDisp(lenslettab[i],λlaser, laserDataView,weightView);
        cost(x::Vector{Float64}) = lkl(x);
        local xopt
        try
            xopt = vmlmb(cost, xinit; verb=false,ftol = (0.0,1e-8),maxeval=500,autodiff=true);
        catch e
            @debug showerror(stdout, e)
            @debug "Error on lenslet  $i"
            continue
        end
        (fwhmopt,copt) = (xopt[1:(nλ)],reshape(xopt[(nλ+1):(3*nλ)],2,:));
        atab[:,i] = lkl.amplitude;
        fwhmtab[:,i] = fwhmopt;
        ctab[:,:,i] = copt;
        next!(p);
    end
    ProgressMeter.finish!(p);
    return (lenslettab, atab, fwhmtab,ctab);
end

function fitSpectralLawAndProfile(laserdata::Matrix{T},
    laserweights::Matrix{T},
    lampdata::Matrix{T},
    lampweights::Matrix{T},
    λlaser::Array{Float64,1},
    lensletsize::NTuple{4, Int},
    position::Matrix{Float64},
    cxinit::Vector{Float64},
    cyinit::Vector{Float64},
    fwhminit::Array{Float64,1},
    wavelengthrange::AbstractArray{Float64,1}
    ; validlensmap ::AbstractVector{Bool} = trues(size(position, 1)),
      profileorder::Int = 2
    ) where T<:Real

    numberoflenslet = size(position,1)
    
    numberoflenslet == length(validlensmap) || throw(DimensionMismatch(
        "size of `position` incompatible with size of `validlensmap`"))
    
    nλ = length(λlaser)
    λ0 = mean(λlaser)# reference
    @assert length(fwhminit) == nλ

    (dxmin, dxmax,dymin,dymax) = lensletsize
    lenslettab = Array{LensletModel}(undef,numberoflenslet);
    laserAmplitude = Array{Float64,2}(undef,nλ,numberoflenslet);
    lampAmplitude = Array{Float64,2}(undef,41,numberoflenslet);
    laserfwhm = Array{Float64,2}(undef,nλ,numberoflenslet);
    laserdist = Array{Float64,2}(undef,2048,2048);
    λMap =  Array{Float64,2}(undef,2048,2048);
    p = Progress(numberoflenslet; showspeed=true)
    Threads.@threads for i in findall(validlensmap)[1:100:end]

        lensletbox = BoundingBox(position[i,1]-dxmin, position[i,1]+dxmax,
                                 position[i,2]-dymin, position[i,2]+dymax)
        # we use RoundNearestTiesUp to enforce box size (special case of `.5` floats)
        lensletbox = round(Int, lensletbox, RoundNearestTiesUp);

        lenslettab[i] = LensletModel(λ0,nλ-1, profileorder,lensletbox);

        # Fit spectral law

        Cinit= [ [position[i,1] cxinit...]; [position[i,2] cyinit...] ];
        xinit = vcat([fwhminit[:],Cinit[:]]...);
        laserDataView = view(laserdata, lensletbox);
        laserWeightView = view(laserweights,lensletbox);
        spectrallkl = LikelihoodDisp(lenslettab[i],λlaser, laserDataView,laserWeightView);
        cost(x::Vector{Float64}) = spectrallkl(x);
        local xopt
        try
            xopt = vmlmb(cost, xinit; verb=false,ftol = (0.0,1e-8),maxeval=500,autodiff=true);
        catch e
            @debug "Error on lenslet  $i" exception=(e, catch_backtrace())
            continue
        end
        fwhm = xopt[1:nλ]
        laserAmplitude[:,i] = spectrallkl.amplitude;
        laserfwhm[:,i] = fwhm
        (dist, pixλ) = distanceMap(wavelengthrange,lenslettab[i]);
        view(laserdist,lensletbox) .= dist;
        view(λMap,lensletbox) .= pixλ;

        # Fit profile

        lampDataView = view(lampdata, lensletbox);
        lampWeightView = view(lampweights,lensletbox);
        profilecoefs = zeros(Float64,2,profileorder+1)
        profilecoefs[1,1:3] .= [2.3, 2.5, 2.9] # maximum(fwhm)
        #profilecoefs[1,1] = 3
        profilecoefs[2,1] = lenslettab[i].dmodel.cx[1]
        pmodel  =  ProfileModel(λ0,profilecoefs)
        profilelkl = LikelihoodProfile(pmodel,lampDataView,lampWeightView,pixλ,lensletbox)
        costpr(x::Matrix{Float64}) = profilelkl(x);
        try
            vmlmb!(costpr, profilecoefs; verb=false,ftol = (0.0,1e-8),maxeval=500,autodiff=true);
        catch e
            @debug showerror(stdout, e)
            @debug "Error on lenslet  $i"
            continue
        end
        pmodel = ProfileModel(λ0, profilecoefs)
        lenslettab[i] = LensletModel(lensletbox,lenslettab[i].dmodel, pmodel)

        profile = @. GaussianModel2(pmodel(pixλ,($(axes(lensletbox,1)))))
        profile = profile ./ sum(profile,dims=1)
        try
            lampAmplitude[:,i] .= updateAmplitudeAndBackground(profile,lampDataView,lampWeightView)
        catch e
            @debug "Error on lenslet  $i" exception=(e, catch_backtrace())
            continue
        end
        next!(p);
    end
    ProgressMeter.finish!(p);
    return (lenslettab, laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap);
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

