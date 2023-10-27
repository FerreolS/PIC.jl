struct LikelihoodProfile{T<:AbstractFloat, A<:AbstractMatrix{T}, B<:AbstractMatrix{T}}
    pmodel    ::ProfileModel
    data      ::A
    weight    ::B
    λMap      ::Matrix{T}
    bbox      ::BoundingBox{Int}
    amplitude ::Vector{T}
    # Inner constructor provided to force using outer constructors.
    function LikelihoodProfile{T,A,B}(pmodel, data, weight, λMap, bbox, amplitude) where {T,A,B}
        size(data) == size(weight) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `weight`"))
        size(data) == size(λMap) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `λMap`"))
        size(data) == size(bbox) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `bbox`"))
        size(data,2) + 1 == length(amplitude) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `amplitude`"))
        new{T,A,B}(pmodel, data, weight, λMap, bbox, amplitude)
    end
end

function LikelihoodProfile(
    pmodel::ProfileModel, data::A, weight::B,
    λMap::AbstractMatrix{T}, bbox::BoundingBox{Int}
) where {T<:AbstractFloat, A<:AbstractMatrix{T}, B<:AbstractMatrix{T}}

    amplitude = zeros(T, size(data,2) + 1)
    LikelihoodProfile{T,A,B}(pmodel, data, weight, λMap, bbox, amplitude)
end

 function  (self::LikelihoodProfile)(C::Matrix{T})::T where (T<:AbstractFloat)
    cx = C[1,:]
    cλ = C[2,:]
    updateProfileModel(self.pmodel, cx, cλ)
    p = @. GaussianModel2(self.pmodel(self.λMap,($(axes(self.bbox,1)))))
    profile = p ./ sum(p,dims=1)
    amp = Zygote.@ignore  updateAmplitudeAndBackground(profile,self.data,self.weight)
    Zygote.@ignore self.amplitude .= amp[:]
    return (sum(abs2,@. self.weight * (self.data - amp[1] - $(reshape(amp[2:end],1,:)) * profile)))
 end
 