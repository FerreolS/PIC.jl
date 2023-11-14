struct SpecPosLikelihood{
        Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64}, Mm<:AbstractMatrix{Float64} }
    pmodel  ::ProfileModel
    box     ::BoundingBox{Int}
    data    ::Md
    weights ::Mw
    λMap    ::Mm
    
    function SpecPosLikelihood{Md,Mw,Mm}(pmodel, box, data, weights, λMap) where {Md,Mw,Mm}
        size(data) == size(box) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `box`"))
        size(data) == size(weights) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `weights`"))
        size(data) == size(λMap) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `λMap`"))
        new{Md,Mw,Mm}(pmodel, box, data, weights, λMap)
    end
end

function SpecPosLikelihood(
    pmodel::ProfileModel, box::BoundingBox{Int}, data::Md, weights::Mw, λMap::Mm
) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64}, Mm<:AbstractMatrix{Float64} }

    SpecPosLikelihood{Md,Mw,Mm}(pmodel, box, data, weights, λMap)
end

function (self::SpecPosLikelihood)(C::Matrix{Float64})::Float64
    cx = C[1,:]
    cλ = C[2,:]
    updateProfileModel(self.pmodel, cx, cλ)
    p = @. GaussianModel2(self.pmodel(self.λMap,($(axes(self.box,1)))))
    profile = p ./ sum(p,dims=1)
    amp = Zygote.@ignore  updateAmplitudeAndBackground(profile,self.data,self.weights)
    return (sum(abs2,@. self.weights * (self.data - amp[1] - $(reshape(amp[2:end],1,:)) * profile)))
end

