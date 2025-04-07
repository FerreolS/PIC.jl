"""
    GaussianModel2(fwhm::Float64, x::AbstractArray)

Compute the value at position sqrt(r) 1D centered Gaussian
* `fwhm` : full-width at half maximum
* `x`:  squared sampled position

Equivalent to `GaussianModel(1.,fwhm, sqrt(x))`
"""
function GaussianModel2(fwhm::T, x::T) ::T where {T<:Real}
    fwhm2sigma = 1 / (2 * sqrt(2 * log(2)))
    exp(-x / (2 * (fwhm * fwhm2sigma)^2))
end

function GaussianModel2(t::NTuple{2,T}) ::T where {T<:Real} return GaussianModel2(t[1], t[2]) end

"""
    LikelihoodIFS(model::LensletModel,wavelengths::AbstractArray{<:Real,1},data::AbstractArray,weight::AbstractArray)

Build the likelihood function for a given lenslet
* `lmodel`: model of the lenslet
* `laser`: wavelengths of the illumination lasers
* `data` : data
* `weight`: precision (ie inverse variance) of the data
"""
struct LikelihoodIFS{T<:Real}
    nλ::Int
    lens::LensletModel
    wavelengths::Vector{T}
    data::Matrix{T}
    weight::Matrix{T}
    spots::Array{T,3}
    amplitude::Vector{T}
    function LikelihoodIFS{T}(nλ, lens, wavelengths, data, weight, spots, amplitude) where {T}
        nλ > lens.dmodel.order || error(" order of the law must be less than number of laser")
        size(data) == size(weight)          || throw(ArgumentError)
        size(spots)[1:2] == size(lens.bbox) || throw(ArgumentError)
        size(spots,3) == nλ                 || throw(ArgumentError)
        length(amplitude) == nλ             || throw(ArgumentError)
        new{T}(nλ, lens, wavelengths, data, weight, spots, amplitude)
    end
end

function LikelihoodIFS(
    lens::LensletModel, wavelengths::AbstractVector,
    data::AbstractMatrix{T}, weight::AbstractMatrix{T}
) where {T<:Real}
    nλ = length(wavelengths)
    spots = zeros(T, size(lens.bbox)..., nλ)
    amplitude = zeros(T, nλ)
    LikelihoodIFS{T}(nλ, lens, wavelengths, data, weight, spots, amplitude)
end

"""
    (self::LikelihoodIFS)(x::Vector{Float64})
compute the likelihood for a given lenslet for the parameters `x`
"""
function  (self::LikelihoodIFS)(x::Vector{T}) ::T where {T<:Real}
    (fwhm::Vector{T},c::Matrix{T}) = (x[1:(self.nλ)],reshape(x[(self.nλ+1):(3*self.nλ)],2,:));
    self(fwhm,c)
end

function  (self::LikelihoodIFS)(fwhm::Vector{T},C::Matrix{T}) ::T where {T<:Real}
    UpdateDispModel(self.lens.dmodel, C);
    bbox = self.lens.bbox;
    (rx,ry) = axes(bbox) # extracting bounding box range
    m = Zygote.Buffer(self.spots);
    @inbounds for (index, λ) in enumerate(self.wavelengths)  # For all laser
        (mx, my)  = self.lens.dmodel(λ);  # center of the index-th Gaussian spot
        r = ((rx.-mx).^2) .+ ((ry.-my).^2)';
        m[:,:,index] = GaussianModel2.(fwhm[index], r);
    end
    spots = copy(m)
    Zygote.@ignore self.amplitude .= updateAmplitude(self.nλ,spots,self.data,self.weight)
    sumspot = zeros(T, size(bbox))
    @inbounds for i =1:self.nλ
        sumspot += self.amplitude[i] *spots[:,:,i]
    end
    return sum(self.weight .* (self.data .-sumspot).^2)
end

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
function updateAmplitude(
    N::Int, spots::AbstractArray{T,3},
    data::AbstractMatrix{T}, weight::AbstractMatrix{T}
) ::Vector{T} where {T<:Real}
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

function fitSpectralLawAndProfile(
    laserdata::Matrix{T},
    laserweights::Matrix{T},
    lampdata::Matrix{T},
    lampweights::Matrix{T},
    λlaser::Vector{Float64},
    λ0::Float64,
    lensletsize::NTuple{4,Int},
    position::Matrix{Float64},
    cxinit::Vector{Float64},
    cyinit::Vector{Float64},
    fwhminit::Vector{Float64},
    wavelengthrange::AbstractVector{Float64};
    validlenslets::AbstractVector{Bool}=trues(size(position,1)),
    profileorder::Int = 2,
    smalltest ::Bool = false
) where {T<:Real}

    numberoflenslet = size(position,1)
    nλ = length(λlaser)

    length(fwhminit) == nλ || throw(ArgumentError)

    (dxmin, dxmax,dymin,dymax) = lensletsize
    nrows_lampAmplitude = (1 + dymin + dymax) + 1 # nrows(lenslet box) + 1 additional cell
    
    lenslettab = Vector{LensletModel}(undef,numberoflenslet);
    laserAmplitude = Matrix{Float64}(undef,nλ,numberoflenslet);
    lampAmplitude = Matrix{Float64}(undef,nrows_lampAmplitude,numberoflenslet);
    laserfwhm = Matrix{Float64}(undef,nλ,numberoflenslet);
    laserdist = Matrix{Float64}(undef,2048,2048);
    λMap =  Matrix{Float64}(undef,2048,2048);
    p = Progress(numberoflenslet; showspeed=true)
    indices = findall(validlenslets)
    indices = smalltest ? rand(MersenneTwister(1234), indices, 300) : indices
    Threads.@threads for i in indices
        lensletbox = round(Int, BoundingBox(
            position[i,1]-dxmin, position[i,1]+dxmax,
            position[i,2]-dymin, position[i,2]+dymax))

        lenslettab[i] = LensletModel(lensletbox, λ0, nλ-1, profileorder);

        # Fit spectral law

        Cinit= [ [position[i,1] cxinit...]; [position[i,2] cyinit...] ];
        xinit = vcat([fwhminit[:],Cinit[:]]...);
        laserDataView = view(laserdata, lensletbox);
        laserWeightView = view(laserweights,lensletbox);
        spectrallkl = LikelihoodIFS(lenslettab[i],λlaser, laserDataView,laserWeightView);
        cost(x::Vector{Float64}) = spectrallkl(x);
        local xopt
        try
            xopt = vmlmb(cost, xinit; verb=false,ftol = (0.0,1e-8),maxeval=500,autodiff=true);
        catch e
            @debug showerror(stdout, e)
            @debug "Error on lenslet  $i"
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
        lampAmplitude[:,i] .= updateAmplitudeAndBackground(profile,lampDataView,lampWeightView)
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
