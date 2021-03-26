
using Zygote
using TwoDimensional


# Dispersion model
mutable struct DispModel
    λ0::Float64   # reference wavelength
    order::Int32  # order of the polynomial
    cx::Array{Float64,1} # coefficients of the polynomial along the x axis
    cy::Array{Float64,1} # coefficients of the polynomial along the y axis
end

# Model of a lenslet
struct LensletModel
    bbox::BoundingBox{Int}  # Boundingbox of influence of the lenslet on the detector
    dmodel::DispModel       # dispersion model of the lenslet
    λlaser::Array{Float64,1}# wavelength of the laser 
end

"""
    lmod = LensletModel(λ0,λlaser,bbox);
    Lenslet model constructor
    * λ0 is the reference wavelength
    * λlaser is a 1D array of laser wavelengths
    * bbox is the bounding box of the lenslet on the detector 
"""
function LensletModel(λ0::Float64, λlaser::Array{Float64,1}, bbox::BoundingBox{Int})
    local order = length(λlaser) - 1; # the order of the polynomial is the number of laser -1
    cx = zeros(Float64, order + 1); # coefficients of the polynomial along the x axis
    cy = zeros(Float64, order + 1); # coefficients of the polynomial along the x axis
    LensletModel(bbox, DispModel(λ0, order, cx, cy), λlaser)
end

"""
    (x,y) = (self::DispModel)(λ) 
    compute the position (x,y)  of the wavelent λ 
    according to the dispersion law DispModel.
"""
function (self::DispModel)(λ::Float64)
    x = self.cx[1];
    y = self.cy[1];
    for o in 1:self.order
        x += self.cx[o + 1]  * (self.λ0 - λ)^(o);
        y += self.cy[o + 1]  * (self.λ0 - λ)^(o);
    end
    return (x, y)
end

"""
    D = update(D,C)
    Update the coefficients C of the DispModel D.
"""
function update(self::DispModel, C::Array{Float64,2})
    self.cx = C[1,:];
    self.cy = C[2,:];
    return self
end

"""
     GaussianModel(A,fwhm,x,y)    
    Compute the value  at position (x,y) of a 2D Gaussian of full-width at half xmaximum 
    Update the coefficients C of the DispModel D.
"""
function GaussianModel(A::Float64, fwhm::Float64, x::Float64, y::Float64)
    local fwhm2sigma = 1 / (2 * sqrt(2 * log(2.)))::Float64
    return A * exp(-(x^2 + y^2) / (2 * (fwhm * fwhm2sigma )^2));
end

"""
     GaussianModel(A,fwhm,x,y)    
    Compute the value of a 2D Gaussian of full-width at half xmaximum 
    Update the coefficients C of the DispModel D.
"""
function GaussianModel(A::Float64, fwhm::Float64, x::Array{Float64,2}, y::Array{Float64,2})
    local fwhm2sigma = 1 / (2 * sqrt(2 * log(2.)))::Float64
    return A .* exp.(.-(x.^2 + y.^2) ./ (2 * (fwhm * fwhm2sigma )^2));
end

"""
    GaussianSpotsCost(data,weight,lmodel, A,fwhm,C)    
    Compute a quadratic cost for the lenslet modelized by lmodel with parameters A fwhm and C
"""
function GaussianSpotsCost(data::Array{Float64,2}, weight::Array{Float64,2}, lmodel::LensletModel, A::Array{Float64,1}, fwhm::Array{Float64,1}, C::Array{Float64,2})
    update(lmodel.dmodel, C);
    s = 0.;
    for I in CartesianIndices(bbox)
        spotsmodel = 0;
        for (index, λ) in enumerate(lmodel.λlaser) 
            (mx, my)  = lmodel.dmodel(λ);
            spotsmodel += GaussianModel(A[index], fwhm[index], I[1] - mx, I[2] - my)  
        end    
        s += weight[I] * ( data[I] - spotsmodel)^2;
    end
    return s;
end


"""
    GaussianSpotsModel(model,lmodel, A,fwhm,C)    
    Build the model of the lenslet modelized by lmodel with parameters A fwhm and c
"""
function GaussianSpotsModel(model::Array{Float64,2}, lmodel::LensletModel, A::Array{Float64,1}, fwhm::Array{Float64,1}, C::Array{Float64,2})
    update(lmodel.dmodel, C);
    t = Zygote.Buffer(model);
    t[:] = model[:];
    for I in CartesianIndices(bbox)
        spotsmodel = 0;
        for (index, λ) in enumerate(lmodel.λlaser) 
            (mx, my)  = lmodel.dmodel(λ);
            t[I[1],I[2]] += GaussianModel(A[index], fwhm[index], I[1] - mx, I[2] - my)  
        end
    end
    return Zygote.copy(t);
end


