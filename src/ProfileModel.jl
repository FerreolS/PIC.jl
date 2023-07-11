

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

function UpdateProfileModel(self::ProfileModel, C::Matrix{Float64})
    @assert length(C)==(self.order+1)  "coefficients size does not match the order"
    self.cλ = C[1,:];
    self.cy = C[2,:];
    return self
end

# GaussianModel2!(ret::AbstractArray{T},fwhm, x::AbstractArray{T}) where (T<:Real)
function getProfile(pmodel::ProfileModel, dist,pixλ)  
    p = Zygote.Buffer(dist);
    #p = similar(dist)
    (w,y) = pmodel(λ)
    for (index,(d, λ)) in enumerate(zip(dist,pixλ)) 
        p[index] = GaussianModel2(w,d^2)
    end
    cp = copy(p)
    return cp./sum(cp,dims=1)
end
function getProfile!(p,pmodel::ProfileModel, dist,pixλ)  
    for (index,(d, λ)) in enumerate(zip(dist,pixλ)) 
        p[index] = GaussianModel2(pmodel(λ),d^2)
    end
    #p .= p./sum(p,dims=1)
    #return p
end
