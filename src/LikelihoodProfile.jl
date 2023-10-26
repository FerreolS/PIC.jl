struct LikelihoodProfile{T<:AbstractFloat,A<:AbstractMatrix{T},B<:AbstractMatrix{T}}
    model::ProfileModel
    data::A
    weight::B
    λMap::Matrix{T}
    bbox::BoundingBox{Int64}
    amplitude::Vector{T}#MVector{N, Float64}
    # Inner constructor provided to force using outer constructors.
    function LikelihoodProfile(model::ProfileModel,
                                    data::A,
                                    weight::B,
                                    λMap::Matrix{T},
                                    bbox::BoundingBox{Int64}) where {T<:AbstractFloat,A<:AbstractMatrix{T},B<:AbstractMatrix{T}}
        @assert size(data) == size(weight)
        @assert size(data) == size(λMap)
        @assert size(data) == size(bbox)
        amplitude =  zeros(T,size(data,2)+1)
        return new{T,A,B}(model,data, weight,λMap, bbox,amplitude)
    end
end


 function  (self::LikelihoodProfile)(coefs::Matrix{T})::T where (T<:AbstractFloat)
    UpdateProfileModel(self.model,coefs)
    p = @. GaussianModel2(self.model(self.λMap,($(axes(self.bbox,1)))))
    profile = p ./ sum(p,dims=1)
    amp = Zygote.@ignore  updateAmplitudeAndBackground(profile,self.data,self.weight)
    Zygote.@ignore self.amplitude .= amp[:]
    return (sum(abs2,@. self.weight * (self.data - amp[1] - $(reshape(amp[2:end],1,:)) * profile)))
 end
 