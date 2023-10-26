

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
    # @inbounds for o in 1:self.order
    #     λpo = (( λ - self.λ0)/self.λ0 )^(o)
    #     w += self.cλ[o + 1]  * λpo;
    #     y += self.cy[o + 1]  * λpo;
    # end

    λpo = (( λ - self.λ0)/self.λ0 ).^(1:self.order)
    w = self.cλ[1] +sum(self.cλ[2:end]  .* λpo)
    y = self.cy[1] +sum(self.cy[2:end] .* λpo)
    
    return (w,(y - x)^2)
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
