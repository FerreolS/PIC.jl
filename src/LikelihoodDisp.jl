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
function (self::LikelihoodDisp)(x::Vector{T}) ::Float64 where {T<:Real}
    fwhm ::Vector{T} = x[1:self.nλ]
    vC   ::Vector{T} = x[self.nλ+1:end]
    C    ::Matrix{T} = reshape(vC, (2,:))
    self(fwhm, C)
end

function (self::LikelihoodDisp)(fwhm::Vector{T}, C::Matrix{T}) ::Float64 where {T<:Real}
    UpdateDispModel(self.lmodel.dmodel, C);
    bbox = self.lmodel.bbox;
    (rx,ry) = axes(bbox) # extracting bounding box range
    m = Zygote.Buffer(self.spots);
    @inbounds for (index, λ) in enumerate(self.wavelengths)  # For all laser
        (mx, my)  = self.lmodel.dmodel(λ);  # center of the index-th Gaussian spot
        r = ((rx.-mx).^2) .+ ((ry.-my).^2)';
        m[:,:,index] = GaussianModel2.( fwhm[index], r);
        #m[:,:,index] = SimpleGauss.(rx, mx, fwhm[index]) .* SimpleGauss.(ry, my, fwhm[index])';
    end
    spots = copy(m)
    Zygote.@ignore  self.amplitude .= updateAmplitude(self.nλ,spots,self.data,self.weight)
    sumspot =   zeros(Float64,size(round(bbox)));
    @inbounds for i =1:self.nλ
        sumspot += self.amplitude[i] *spots[:,:,i]
    end
    return Float64.(sum(self.weight .* (self.data .-sumspot).^2))
end
 
 