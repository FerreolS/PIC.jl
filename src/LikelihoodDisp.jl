"""
    LikelihoodDisp(model::LensletModel,wavelengths::AbstractArray{<:Real,1},data::AbstractArray,weight::AbstractArray)

Build the likelihood function for a given lenslet
* `lmodel`: model of the lenslet
* `laser`: wavelengths of the illumination lasers
* `data` : data
* `weight`: precision (ie inverse variance) of the data
"""
struct LikelihoodDisp{T<:Real}
    nλ          ::Int
    lmodel      ::LensletModel
    wavelengths ::Vector{T}
    data        ::Matrix{T}
    weight      ::Matrix{T}
    spots       ::Array{T,3}
    amplitude   ::Vector{T}
    
    function LikelihoodDisp{T}(
        nλ, lmodel, wavelengths, data, weight, spots, amplitude
    ) where {T<:Real}
    
        nλ > lmodel.dmodel.order || throw(ArgumentError(
            "the order of the law must be less than the number of laser"))
        nλ == length(wavelengths) || throw(DimensionMismatch(
            "`nλ` is incompatible with size of `wavelengths`"))
        size(data) == size(weight) || throw(DimensionMismatch(
            "size of `data` must match size of `weight`"))
        size(spots) == (size(lmodel.bbox)..., nλ) || throw(DimensionMismatch(
            "incorrect size of `spots`"))
        size(amplitude) == (nλ,) || throw(DimensionMismatch(
            "incorrect size of `amplitude`"))
            
        new{T}(nλ, lmodel, wavelengths, data, weight, spots, amplitude)
    end
end

"""
    LikelihoodDisp(nλ, model, wavelengths, data, weight, spots, amplitude)

Version with automatic setting of the type parameter.
"""
function LikelihoodDisp(nλ, lmodel, wavelengths, data, weight, spots, amplitude)
    T = float(promote_type(
        eltype(wavelengths), eltype(data), eltype(weight), eltype(spots), eltype(amplitude) ))
    LikelihoodDisp{T}(nλ, lmodel, wavelengths, data, weight, spots, amplitude)
end

"""
    LikelihoodDisp(model, wavelengths, data, weight)

Version with automatic setting of `nλ`, `spots`, `amplitude`, and of the type parameter.
"""
function LikelihoodDisp(lmodel, wavelengths, data, weight)
    nλ = length(wavelengths);
    spots = zeros(Float64, size(lmodel.bbox)..., nλ)
    amplitude = zeros(Float64, nλ)
    LikelihoodDisp(nλ, lmodel, wavelengths, data, weight, spots, amplitude)
end


"""
    (::LikelihoodDisp)(x::Vector{Float64}) -> Float64
    compute the likelihood for a given lenslet for the parameters `x`
"""
function (self::LikelihoodDisp)(V::Vector{T}) ::Float64 where {T<:Real}
    
    startfwhm = 1
    endfwhm   = self.nλ
    startcx   = endfwhm + 1
    endcx     = startcx + self.lmodel.dmodel.order
    startcy   = endcx   + 1
    endcy     = startcy + self.lmodel.dmodel.order
    
    fwhm = V[startfwhm : endfwhm]
    cx   = V[startcx   : endcx]
    cy   = V[startcy   : endcy]
    
    self(fwhm, cx, cy)
end

function (self::LikelihoodDisp)(
    fwhm::Vector{T}, cx::Vector{T}, cy::Vector{T}
) ::Float64 where {T<:Real}

    updateDispModel(self.lmodel.dmodel, cx, cy)
    
    (rx, ry) = axes(self.lmodel.bbox) # extracting bounding box range
    m = Zygote.Buffer(self.spots);
    @inbounds for (i,λ) in enumerate(self.wavelengths) # for each laser
        (mx, my)  = self.lmodel.dmodel(λ)  # center of the index-th Gaussian spot
        r = ((rx.-mx).^2) .+ ((ry.-my).^2)'
        m[:,:,i] = GaussianModel2.(r, fwhm[i])
    end
    spots = copy(m)
    
    Zygote.@ignore self.amplitude .= updateAmplitude(self.nλ, spots, self.data, self.weight)
    any(isnan, self.amplitude) && error("NaN amplitude: \"$(self.amplitude)\"")
    
    sumspot = zeros(Float64, size(self.lmodel.bbox))
    @inbounds for i in 1:self.nλ
        sumspot += self.amplitude[i] * spots[:,:,i]
    end
    
    return Float64.(sum(self.weight .* (self.data .- sumspot).^2))
end
 
struct WaveLampLikelihood
    order ::Int
    nλ ::Int
    λs ::Vector{Float64}
    λ0 ::Float64
    λpowers ::Matrix{Float64}
    box ::BoundingBox{Int}
    data ::Matrix{Float64}

    last_fhwm ::Vector{Float64}
    last_cx   ::Vector{Float64}
    last_cy   ::Vector{Float64}
    last_amp  ::Vector{Float64}
end


λpowers = ( (λ - self.λ0) / self.λ0 ).^(1:self.order)


function (self:WaveLampLikelihood)(
    fwhm::Vector{Float64}, cx::Vector{Float64}, cy::Vector{Float64}
) ::Matrix{Float64}
    
    # register new inputs for convenience
    Zygote.@ignore begin
        self.last_fwhm .= fwhm
        self.last_cx   .= cx
        self.last_cy   .= cy
        self.last_amp  .= NaN # to be clean, in case of failure in the next computations
    end
    
    eachspot  = Zygote.Buffer(Array{Float64}(undef, size(box)))
    
    sumspot = Zygote.Buffer(fill(0e0, size(box)))
    
    spot_amps = Zygote.Buffer(Vector{Float64}(undef, 3))
    
    (xs, ys) = axes(self.box)
    
    for i in 1:self.nλ
    
        # center of the laser spot
        spot_ctr_x = cx[1] + sum(cx[2:end] .* self.λpowers[:,i])
        spot_ctr_y = cy[1] + sum(cy[2:end] .* self.λpowers[:,i])
        
        # gaussian spot
        matrixradius = (xs .- spot_ctr_x).^2 .+ ((ys .- spot_ctr_y).^2)'
        eachspot .= GaussianModel2.(matrixradius, fwhm[i])
    end
    
    # amplitude of gaussian spot
    Zygote.@ignore spot_amps .= updateAmplitude(self.nλ, spots, self.data, self.weight)
    # register for convenience
    Zygote.@ignore self.last_amp .= spot_amps
    onespot .*= spot_amp
    
    for i in 1:self.nλ
        allspots .+= onespot[:,:,i] .* spots_amps[i]
    end
    
    copy(allspots) # extract Array from Zygote Buffer
end

