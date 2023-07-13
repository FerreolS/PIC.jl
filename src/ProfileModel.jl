

mutable struct ProfileModel    
    λ0::Float64   # reference wavelength
    order::Int64  # order of the polynomial
    cλ::Vector{Float64} # coefficients of the polynomial along the wavelength axis
    cy::Vector{Float64} # coefficients of the polynomial along the y axis
    
    function ProfileModel(λ0,order::Int64,cλ::Vector{Float64},cy::Vector{Float64})
        @assert length(cλ)==(order+1) "coefficients size does not match the order"
        @assert length(cy)==(order+1) "coefficients size does not match the order"
        new(λ0,order,cλ,cy)
    end
end

function ProfileModel(λ0,order::Int64)
    cλ = zeros(order+1)
    cλ[1]=1
    cy = zeros(order+1)
    cy[1]=1
    ProfileModel(λ0,order,cλ,cy)
end

function ProfileModel(λ0,C::Matrix{Float64})
	order = size(C,2)-1
    cλ = C[1,:]
    cy = C[2,:]
    ProfileModel(λ0,order,cλ,cy)
end


function (self::ProfileModel)(λ::Float64)
    w = self.cλ[1];
    y = self.cy[1];
    @inbounds for o in 1:self.order
        λpo = (( λ - self.λ0)/self.λ0 )^(o)
        w += self.cλ[o + 1]  * λpo;
        y += self.cy[o + 1]  * λpo;
    end
    return (w,y)
end

function (self::ProfileModel)(λ::Float64,x)
    w = self.cλ[1];
    y = self.cy[1];
    @inbounds for o in 1:self.order
        λpo = (( λ - self.λ0)/self.λ0 )^(o)
        w += self.cλ[o + 1]  * λpo;
        y += self.cy[o + 1]  * λpo;
    end
    return (w,(y-x)^2)
end

function UpdateProfileModel(self::ProfileModel, C::Matrix{Float64})
    @assert size(C)==(2,self.order+1)  "coefficients size does not match the order"
    self.cλ = C[1,:];
    self.cy = C[2,:];
    return self
end


# # GaussianModel2!(ret::AbstractArray{T},fwhm, x::AbstractArray{T}) where (T<:Real)
# function getProfile(pmodel::ProfileModel,pixλ)  
#     p = Zygote.Buffer(pixλ);
#     #p = similar(dist)
#     for (index, λ) in enumerate(pixλ) 
#         p[index] = GaussianModel2(pmodel(λ)...)
#     end
#     cp = copy(p)
#     return cp./sum(cp,dims=1)
# end

# function getProfile!(p,pmodel::ProfileModel,pixλ)  
#     for (index, λ) in enumerate(pixλ)
#         p[index] = GaussianModel2(pmodel(λ)...)
#     end
#     p .= p./sum(p,dims=1)
#     #return p
# end



struct LikelihoodProfile{T<:AbstractFloat}
    model::ProfileModel
    data::Matrix{T}
    weight::Matrix{T}
    λMap::Matrix{T}
    bbox::BoundingBox{Int64}
    amplitude::Vector{T}#MVector{N, Float64}
    # Inner constructor provided to force using outer constructors.
    function LikelihoodProfile{T}(model::ProfileModel,
                                    data::Matrix{T},
                                    weight::Matrix{T},
                                    λMap::Matrix{T},
                                    bbox::BoundingBox{Int64}) where {T<:AbstractFloat}
        @assert size(data) == size(weight)
        @assert size(data) == size(λMap)
        @assert size(data) == size(bbox)
        amplitude =  zeros(T,size(data,2)+1)
        return new{T}(model,data, weight,λMap, bbox,amplitude)
    end
end


function  (self::LikelihoodProfile)(coefs::Vector{T})::T where (T<:AbstractFloat)
    UpdateProfileModel(self.model,coefs)
    #profile =Zygote.Buffer(self.distMap)
    #getProfile!(profile,self.model,self.λMap, self.distMap)
    #profile = copy(profile)
    #@show profile = getProfile(self.model,self.λMap, self.distMap)
    p = @. GaussianModel2(self.model.(self.λMap)...)
    profile = p ./ sum(p,dims=1)
    amp = Zygote.@ignore  updateAmplitude(profile,self.data,self.weight)
    Zygote.@ignore self.amplitude .= amp[:]
    return (sum(abs2,@. self.weight * (self.data - amp .* profile)))
 end

 function  (self::LikelihoodProfile)(coefs::Matrix{T})::T where (T<:AbstractFloat)
    UpdateProfileModel(self.model,coefs)
    #profile =Zygote.Buffer(self.distMap)
    #getProfile!(profile,self.model,self.λMap, self.distMap)
    #profile = copy(profile)
    #@show profile = getProfile(self.model,self.λMap, self.distMap)
    profile = @. GaussianModel2(self.model(self.λMap,($(axes(self.bbox,1)))))
    #profile = p ./ sum(p,dims=1)
    amp = Zygote.@ignore  updateAmplitudeAndBackground(profile,self.data,self.weight)
    Zygote.@ignore self.amplitude .= amp[:]
    return (sum(abs2,@. self.weight * (self.data - amp[1] - $(reshape(amp[2:end],1,:)) * profile)))
 end
