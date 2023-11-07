struct SpecPosLikelihood{
        Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64}, Mm<:AbstractMatrix{Float64} }
    pmodel  ::ProfileModel
    box     ::BoundingBox{Int}
    data    ::Md
    weights ::Mw
    λMap    ::Mm
    
    last_amp ::Vector{Float64}
    
    function SpecPosLikelihood{Md,Mw,Mm}(pmodel, box, data, weights, λMap, last_amp) where {Md,Mw,Mm}
        size(data) == size(box) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `box`"))
        size(data) == size(weights) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `weights`"))
        size(data) == size(λMap) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `λMap`"))
        size(data,2) + 1 == length(last_amp) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `last_amp`"))
        new{Md,Mw,Mm}(pmodel, box, data, weights, λMap, last_amp)
    end
end

function SpecPosLikelihood(
    pmodel::ProfileModel, box::BoundingBox{Int}, data::Md, weights::Mw, λMap::Mm
) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64}, Mm<:AbstractMatrix{Float64} }

    last_amp = Vector{Float64}(undef, size(data,2) + 1)
    SpecPosLikelihood{Md,Mw,Mm}(pmodel, box, data, weights, λMap, last_amp)
end

function (self::SpecPosLikelihood)(C::Matrix{Float64})::Float64
    cx = C[1,:]
    cλ = C[2,:]
    updateProfileModel(self.pmodel, cx, cλ)
    p = GaussianModel2.(cp(self).(self.λMap, axes(self.box,1)))
    profile = p ./ sum(p,dims=1)
    amp = Zygote.@ignore  updateAmplitudeAndBackground(profile,self.data,self.weights)
    Zygote.@ignore self.last_amp .= amp[:]
    return (sum(abs2,@. self.weights * (self.data - amp[1] - $(reshape(amp[2:end],1,:)) * profile)))
end

function cp(self)
    function (λ, a)
        λpowers = ( (λ-self.pmodel.λ0) / self.pmodel.λ0 ).^(1:self.pmodel.order)
        x = self.pmodel.cx[1] + sum(self.pmodel.cx[2:end] .* λpowers)
        w = self.pmodel.cλ[1] + sum(self.pmodel.cλ[2:end] .* λpowers)
        
        return ((x - a)^2, w)
    end
end
