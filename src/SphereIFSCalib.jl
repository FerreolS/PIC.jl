module SphereIFSCalib

using Zygote, StaticArrays
using TwoDimensional


"""
    DispModel(λ0::Float64,order::Int32,cx::Array{Float64,1},cy::Array{Float64,1})

The dispersion model giving the position of a wavelength on the detector
* `λ0` is the reference wavelength
* `order` is the order of the polynomials
* `cx` is an array of coefficients of the polynomial along the x axis
* `cy` is an array of coefficients of the polynomial along the y axis
"""
mutable struct DispModel
    λ0::Float64   # reference wavelength
    order::Int  # order of the polynomial
    cx::Array{Float64,1} # coefficients of the polynomial along the x axis
    cy::Array{Float64,1} # coefficients of the polynomial along the y axis
end

"""
    (self::DispModel)(λ::Float64)

compute the position `(x,y)`  of the wavelent `λ`
according to the dispersion law `DispModel`.

### Example
```
D = DispModel(λ0, order, cx, cy);
(x,y) = D(λ)
```
"""
function (self::DispModel)(λ::Float64)
    x = self.cx[1];
    y = self.cy[1];
    for o in 1:self.order
        λpo = (( λ - self.λ0)/self.λ0 )^(o)
        x += self.cx[o + 1]  * λpo;
        y += self.cy[o + 1]  * λpo;
    end
    return (x, y)
end

"""
    UpdateDispModel(self::DispModel, C::Array{Float64,2})

Update the coefficients  of the DispModel .
* `self`: DispModel object
* `C` : array containing the polynomial coefficients.
"""
function UpdateDispModel(self::DispModel, C::Array{Float64,2})
    @assert size(C)==(2,self.order+1) "coefficients array does not have the right size"
    self.cx = C[1,:];
    self.cy = C[2,:];
    return self
end

"""
    LaserModel(nλ::Int,λlaser::Array{Float64,1},amplitude::Array{Float64,1},fwhm::Array{Float64,1})

Model of the laser illumination.
It consist on:
* `nλ` the number of laser
* `λlaser` an array of the wavelengths of the lasers
* `amplitude` an array of the amplitude of the maximum of the Gaussian spot
* `fwhm` an array of the full width at half maximum of the Gaussians
"""
mutable struct LaserModel
    nλ::Int
    λlaser::Array{Float64,1}# wavelength of the laser
    amplitude::Array{Float64,1}
    fwhm::Array{Float64,1}
end

"""
    LaserModel(λlaser::Array{Float64,1},amplitude::Array{Float64,1},fwhm::Array{Float64,1})

Constructor of the model of the laser illumination.
It consist on:
* `λlaser` an array of the wavelengths of the lasers
* `amplitude` an array of the amplitude of the maximum of each Gaussian spots
* `fwhm` an array of the full width at half maximum of the Gaussians
"""
function LaserModel(λlaser::Array{Float64,1},amplitude::Array{Float64,1},fwhm::Array{Float64,1})
    nλ = length(λlaser);
    @assert length(amplitude)==nλ "amplitude vector does not have the right size"
    @assert length(fwhm)==nλ "fwhm vector does not have the right size"
    LaserModel(nλ ,λlaser,amplitude,fwhm);
end

"""
    UpdateLaserModel(self::LaserModel,A::Array{Float64,1},fwhm::Array{Float64,1})

Update the parameters of the laser model
* `self` :  LaserModel object
* `A` : 1D  array of amplitudes of the maximum of each Gaussian spot
* `fwhm` 1D  array of the full width at half maximum of each Gaussian spot

`A` and `fwhm` must have lenth of `self.nλ`
"""
function UpdateLaserModel(self::LaserModel,A::Array{Float64,1},fwhm::Array{Float64,1})
    @assert length(A)==self.nλ "amplitude vector does not have the right size"
    @assert length(fwhm)==self.nλ "fwhm vector does not have the right size"
    self.amplitude = A;
    self.fwhm = fwhm;
end


"""
    LensletModel(bbox::BoundingBox{Int},dmodel::DispModel)

Model of a lenslet
The image of a lenslet on the detector is decribed by:
* `bbox` the boundingbox of its influence on the detector
* `dmodel` the dispersion model described by a object of type `DispModel`
"""
struct LensletModel
    bbox::BoundingBox{Int}  # Boundingbox of influence of the lenslet on the detector
    dmodel::DispModel       # dispersion model of the lenslet
 #   λlaser::Array{Float64,1}# wavelength of the laser
end


"""
    lmod = LensletModel(λ0::Float64, order::Int, bbox::BoundingBox{Int})

Lenslet model constructor
* `λ0`  : reference wavelength
* `order` : order of the polynomials
* `bbox` : bounding box of the lenslet on the detector
"""
function LensletModel(λ0::Float64, order::Int, bbox::BoundingBox{Int})
    cx = zeros(Float64, order + 1); # coefficients of the polynomial along the x axis
    cy = zeros(Float64, order + 1); # coefficients of the polynomial along the x axis
    LensletModel(bbox, DispModel(λ0, order, cx, cy))
end


"""
    lmod = LensletModel(λ0::Float64, order::Int, bbox::BoundingBox{Int},cx0::Float64,cy0::Float64)

Lenslet model constructor
* `λ0`  : reference wavelength
* `order` : order of the polynomials
* `bbox` : bounding box of the lenslet on the detector
* `cx0` :
* `cy0` :
"""
function LensletModel(λ0::Float64, order::Int, bbox::BoundingBox{Int},cx0::Float64,cy0::Float64)
    cx = zeros(Float64, order + 1); # coefficients of the polynomial along the x axis
    cy = zeros(Float64, order + 1); # coefficients of the polynomial along the x axis
    cx[1] = cx0;
    cy[1] = cy0;
    LensletModel(bbox, DispModel(λ0, order, cx, cy))
end


"""
    lmod = LensletModel(λ0::Float64, order::Int, bbox::BoundingBox{Int})

Lenslet model constructor
* `λ0`  : reference wavelength
* `order` : order of the polynomials
* `bbox` : bounding box of the lenslet on the detector
"""
function LensletModel(λ0::Float64, order::Int,cx0::Float64,cy0::Float64, widthx::Number, widthy::Number)
    cx = zeros(Float64, order + 1); # coefficients of the polynomial along the x axis
    cy = zeros(Float64, order + 1); # coefficients of the polynomial along the x axis
    cx[1] = cx0;
    cy[1] = cy0;
    bbox = round(Int,BoundingBox(xmin=cx0-widthx, ymin=cy0-widthy, xmax=cx0+widthx, ymax=cy0+widthy));
    LensletModel(bbox, DispModel(λ0, order, cx, cy))
end




"""
    GaussianModel(A,fwhm,x,y)

Compute the value  at position (x,y) of a 2D centered Gaussian
* `fwhm` : full-width at half maximum
* `A` : amplitude at (x,y) = (0,0)
* `x`, `y`: sampled postions
"""
function GaussianModel(A::Float64, fwhm::Float64, x::Float64, y::Float64)
    local fwhm2sigma = 1 / (2 * sqrt(2 * log(2.)))::Float64
    return A * exp(-(x^2 + y^2) / (2 * (fwhm * fwhm2sigma )^2));
end

"""
    GaussianModel(A::Float64, fwhm::Float64, x::AbstractArray)

Compute the value at position x 1D centered Gaussian
* `A` : amplitude at x = 0
* `fwhm` : full-width at half maximum
* `x`: array of the sampled position
"""
function GaussianModel(A::Float64, fwhm::Float64, x::AbstractArray)
    local fwhm2sigma = 1 / (2 * sqrt(2 * log(2.)))::Float64
    return A .* exp.(.-(x.^2) ./ (2 * (fwhm * fwhm2sigma )^2));
end

"""
    GaussianModel2(A::Float64, fwhm::Float64, x::AbstractArray)

Compute the value at position sqrt(r) 1D centered Gaussian
* `A` : amplitude at x = 0
* `fwhm` : full-width at half maximum
* `x`: array of the squared sampled position

Equivalent to `GaussianModel(A, fwhm, sqrt.(x))`
"""
function GaussianModel2(A::Float64, fwhm::Float64, x::AbstractArray)
    local fwhm2sigma = Float64(1) / (2 * sqrt(2 * log(2.)))
    return A .* exp.(-x ./ (2 * (fwhm * fwhm2sigma )^2));
end

function GaussianModel2(A::Float64, fwhm::Float64, x::AbstractFloat)
    local fwhm2sigma = Float64(1) / (2 * sqrt(2 * log(2.)))
    return A * exp(-x / (2 * (fwhm * fwhm2sigma )^2));
end

function GaussianModel2(fwhm::Float64, x::AbstractArray)
    local fwhm2sigma = Float64(1) / (2 * sqrt(2 * log(2.)))
    return exp.(-x ./ (2 * (fwhm * fwhm2sigma )^2));
end

function GaussianModel2(fwhm::Float64, x::AbstractFloat)
    local fwhm2sigma = Float64(1) / (2 * sqrt(2 * log(2.)))
    return exp(-x / (2 * (fwhm * fwhm2sigma )^2));
end

function GaussianModel2!(ret::AbstractArray{T},fwhm::Float64, x::AbstractArray{T}) where (T<:AbstractFloat)
        ret .= exp.(-x ./ (2 * (fwhm * Float64(1) / (2 * sqrt(2 * log(2.))) )^2));
        nothing
end

"""
    GaussianSpotsCost(data::Array{Float64,2}, weight::Array{Float64,2}, lmodel::LensletModel,  laser::LaserModel,A::Array{Float64,1}, fwhm::Array{Float64,1}, C::Array{Float64,2})

Compute a weighted quadratic cost of a lenslet model :
cost = weight .* || data - model ||^2
* `lmodel`:  model of the lenslet
* `A` : 1D array containing the amplitude  of all Gaussian spots
* `fwhm` : 1D array containing the fwhm  of all Gaussian spots
* `C` : 2D array containing the chromatic law coefficients.
"""
function GaussianSpotsCost(data::Array{Float64,2}, weight::Array{Float64,2}, lmodel::LensletModel,  laser::LaserModel,A::Array{Float64,1}, fwhm::Array{Float64,1}, C::Array{Float64,2})
    UpdateDispModel(lmodel.dmodel, C);
    UpdateLaserModel(laser,A,fwhm);
    s = 0.;
    for I in CartesianIndices(lmodel.bbox)
        spotsmodel = 0;
        for (index, λ) in enumerate(laser.λlaser)
            (mx, my)  = lmodel.dmodel(λ);
            spotsmodel += GaussianModel(laser.amplitude[index], laser.fwhm[index], I[1] - mx, I[2] - my)
        end
        s += weight[I] * ( data[I] - spotsmodel)^2;
    end
    return s;
end


"""
    GaussianSpotsModel(lmodel::LensletModel,laser::LaserModel, A::Array{Float64,1}, fwhm::Array{Float64,1}, C::Array{Float64,2})

Build the model of a lenslet
* `lmodel`:  model of the lenslet
* `A` : 1D array containing the amplitude  of all Gaussian spots
* `fwhm` : 1D array containing the fwhm  of all Gaussian spots
* `C` : 2D array containing the chromatic law coefficients.
"""
function GaussianSpotsModel(lmodel::LensletModel,laser::LaserModel, A::Array{Float64,1}, fwhm::Array{Float64,1}, C::Array{Float64,2})
    UpdateDispModel(lmodel.dmodel, C);
    UpdateLaserModel(laser,A,fwhm);
    bbox = lmodel.dmodel;
    model = zeros(Float64,round(bbox).xmax-round(bbox).xmin+1,round(bbox).ymax-round(bbox).ymin+1);
    t = Zygote.Buffer(model);
    t[:] = model[:];
    for I in CartesianIndices(bbox)
        spotsmodel = 0;
        for (index, λ) in enumerate(laser.λlaser)
            (mx, my)  = lmodel.dmodel(λ);
            t[I[1],I[2]] += GaussianModel(laser.amplitude[index], laser.fwhm[index], I[1] - mx, I[2] - my)
        end
    end
    model =  Zygote.copy(t);
end


"""
    LensletLaserImage(lmodel::LensletModel,laser::LaserModel)

Build the image of a lenslet under laser illumination
* `lmodel`: model of the lenslet
* `laser`: model of the laser illumination
"""
function LensletLaserImage(lmodel::LensletModel,laser::LaserModel)
    bbox = lmodel.bbox;
    (rx,ry) = axes(bbox) # extracting bounding box range
    spotsmodel =   zeros(Float64,size(round(bbox)));
    @inbounds for (index, λ) in enumerate(laser.λlaser)  # For all laser
        (mx, my)  = lmodel.dmodel(λ);  # center of the index-th Gaussian spot
        r = ((rx.-mx).^2) .+ ((ry.-my).^2)';
        spotsmodel = spotsmodel .+ GaussianModel2.(laser.amplitude[index], laser.fwhm[index], r)
    end
    return spotsmodel;
end

"""
    LensletLaserImageA(lmodel::LensletModel,laser::LaserModel)

Build the image of a lenslet under laser illumination
* `lmodel`: model of the lenslet
* `laser`: model of the laser illumination
"""
function LensletLaserImage!(spotsmodel::Array{Float64,3},lmodel::LensletModel,laser::LaserModel)
    bbox = lmodel.bbox;
    (rx,ry) = axes(bbox) # extracting bounding box range
  #  model = zeros(Float64,laser.nλ,size(round(bbox)...));
   # t = Zygote.Buffer(spotsmodel);
  #  spotsmodel =   zeros(Float64,size(round(bbox)));
    @inbounds for (index, λ) in enumerate(laser.λlaser)  # For all laser
        (mx, my)  = lmodel.dmodel(λ);  # center of the index-th Gaussian spot
        r = ((rx.-mx).^2) .+ ((ry.-my).^2)';
        spotsmodel[:,:,index]= GaussianModel2( laser.fwhm[index], r)
    end
    nothing
end
"""
    LensletLaserImageA(lmodel::LensletModel,laser::LaserModel)

Build the image of a lenslet under laser illumination
* `lmodel`: model of the lenslet
* `laser`: model of the laser illumination
"""
function LensletLaserImage2(lmodel::LensletModel,laser::LaserModel)
    bbox = lmodel.bbox;
    (rx,ry) = axes(bbox) # extracting bounding box range
    model = zeros(Float64,size(round(bbox))...,laser.nλ);
    t = Zygote.Buffer(model);
  #  spotsmodel =   zeros(Float64,size(round(bbox)));
    @inbounds for (index, λ) in enumerate(laser.λlaser)  # For all laser
        (mx, my)  = lmodel.dmodel(λ);  # center of the index-th Gaussian spot
        r = ((rx.-mx).^2) .+ ((ry.-my).^2)';
        t[:,:,index]= GaussianModel2.( laser.fwhm[index], r)
    end
    return copy(t)
end


"""
    LikelihoodIFS(model::LensletModel,laser::LaserModel,data::AbstractArray,weight::AbstractArray)

Build the likelihood function for a given lenslet
* `lmodel`: model of the lenslet
* `laser`: model of the laser illumination
* `data` : data
* `weight`: precision (ie inverse variance) of the data
"""
struct LikelihoodIFS{T<:AbstractFloat,N}
    nspots::Int
    model::LensletModel
    laser::LaserModel
    data::Array{T,2}
    weight::Array{T,2}
    dwd::Float64
    spots::Array{T,3}
    wspots::Array{T,3}
    b::Array{T,1}#MVector{N, Float64}
    A::Array{T,2}#MMatrix{N, N,Float64}
    # Inner constructor provided to force using outer constructors.
    function LikelihoodIFS{T}(model::LensletModel,
        laser::LaserModel,
        data::Array{T,2},
        weight::Array{T,2}) where {T<:AbstractFloat}
        N=laser.nλ;
        #@assert laser.nλ == N
        @assert laser.nλ > model.dmodel.order " the order of the law must be less than the number of laser"
        @assert size(data) == size(weight)
        dwd = sum(weight.*data.^2)
        spots = zeros(Float64,size(round(model.bbox))...,N)
        wspots = similar(spots)
        b =  zeros(Float64,N)
        A =  zeros(Float64,N,N)
        N2 = N*N
        return new{T,N}(N,model,laser,data, weight,dwd,spots,wspots,b,A)
    end
    # Inner constructor provided to force using outer constructors.
    function LikelihoodIFS{T}(model::LensletModel,
        laser::LaserModel,
        data::Array{T,2},
        weight::T) where {T<:AbstractFloat}
        N=laser.nλ;
        #@assert laser.nλ == N
        @assert laser.nλ > model.dmodel.order " the order of the law must be less than the number of laser"
        dwd = sum(weight.*data.^2)
        spots = zeros(Float64,size(round(model.bbox))...,N)
        wspots = similar(spots)
        b = @MVector zeros(Float64,N)
        A = @MMatrix zeros(Float64,N,N)
        return new{T,N}(N,model,laser,data, weight*ones(1,1),dwd,spots,wspots,b,A)
    end
end

function LikelihoodIFS(model::LensletModel,laser::LaserModel,data::AbstractArray{<:Real,2})
    T = float(eltype(data))
    LikelihoodIFS{T}(model,laser,convert(Array{T,2},data),T(1.0))
end

function LikelihoodIFS(model::LensletModel,
                        laser::LaserModel,
                        data::AbstractArray{<:Real,2},
                        weight::Union{Real,AbstractArray{<:Real,2}})
    T = float(promote_type(eltype(data),eltype(weight)))
    LikelihoodIFS{T}(model,laser,convert(Array{T,2},data), T.(weight))
end

function  (self::LikelihoodIFS)(a::Array{Float64,1},fwhm::Array{Float64,1},C::Array{Float64,2})::Float64
    UpdateDispModel(self.model.dmodel, C);
    UpdateLaserModel(self.laser,a,fwhm);
    return Float64.(sum(self.weight .* (self.data .- LensletLaserImage(self.model,self.laser)).^2))
end

"""
    (self::LikelihoodIFS)(x::Vector{Float64})
    compute the likelihood for a given lenslet for the parameters `x`

    ### Example
    ```
    laser =  LaserModel(λlaser,a0,fwhm0)
    lenslet = LensletModel(λ0,laser.nλ-1,round(bbox))
    xinit = vcat([ainit[:],fwhminit[:],cinit[:]]...)
    lkl = LikelihoodIFS(lenslet,laser,view(data,lenslet.bbox), view(weight,lenslet.bbox))
    xopt = vmlmb(lkl, xinit; verb=50)
    ```
"""
#= function  (self::LikelihoodIFS)(x::Vector{Float64})::Float64
    (a::Vector{Float64},fwhm::Vector{Float64},c::Matrix{Float64}) = (x[1:(self.laser.nλ)],x[(self.laser.nλ+1):(2*self.laser.nλ)],reshape(x[(2*self.laser.nλ+1):(4*self.laser.nλ)],2,:));
    self(a,fwhm,c)
end =#
function  (self::LikelihoodIFS)(x::Vector{Float64})::Float64
    (fwhm::Vector{Float64},c::Matrix{Float64}) = (x[1:(self.laser.nλ)],reshape(x[(self.laser.nλ+1):(3*self.laser.nλ)],2,:));
    self(fwhm,c)
end

function  (self::LikelihoodIFS)(fwhm::Array{Float64,1},C::Array{Float64,2})::Float64
    # @assert length(fwhm)== self.laser.nλ "length(fwhm) must equal to the number of lasers"
     UpdateDispModel(self.model.dmodel, C);
     bbox = self.model.bbox;
     (rx,ry) = axes(bbox) # extracting bounding box range
     m = Zygote.Buffer(self.spots);
     @inbounds for (index, λ) in enumerate(self.laser.λlaser)  # For all laser
        (mx, my)  = self.model.dmodel(λ);  # center of the index-th Gaussian spot
        r = ((rx.-mx).^2) .+ ((ry.-my).^2)';
        m[:,:,index] = GaussianModel2.( fwhm[index], r);
     end
    spots = copy(m)
    Zygote.@ignore  self.b .= updateAmplitude(self.nspots,spots,self.data,self.weight)
    sumspot =   zeros(Float64,size(round(bbox)));
    @inbounds for i =1:self.nspots
        sumspot += self.b[i] *spots[:,:,i]
    end
    return Float64.(sum(self.weight .* (self.data .-sumspot).^2))
 end

#=
function  (self::LikelihoodIFS)(fwhm::Array{Float64,1},C::Array{Float64,2})::Float64
   # @assert length(fwhm)== self.laser.nλ "length(fwhm) must equal to the number of lasers"
    UpdateDispModel(self.model.dmodel, C);
    bbox = self.model.bbox;
    (rx,ry) = axes(bbox) # extracting bounding box range
    m = Zygote.Buffer(self.spots);
    mw = Zygote.Buffer(self.wspots);
    a = Zygote.Buffer(self.A)
    b = Zygote.Buffer(self.b)
    @inbounds for (index, λ) in enumerate(self.laser.λlaser)  # For all laser
        (mx, my)  = self.model.dmodel(λ);  # center of the index-th Gaussian spot
        r = ((rx.-mx).^2) .+ ((ry.-my).^2)';
        GaussianModel2!(( m[:,:,index]), fwhm[index], r);
        mw[:,:,index] .=  m[:,:,index].* self.weight ;
        b[index] = Float64.(sum(mw[:,:,index].*self.data ));
        a[index,index] = Float64.(sum(mw[:,:,index].*m[:,:,index]));
        for i=1:index-1
            tmp = Float64.(sum(mw[:,:,index].*m[:,:,i]))
            a[i,index] = tmp
            a[index,i] = tmp;
        end
    end
    #self.A = copy(a);
  #  amplitude =  similar(self.b)
  A= copy(a)
  B= copy(b)
    amplitude =updateAmplitude(A,B)
    return Float64.(sum(amplitude .* (A * amplitude + 2 * B)) + self.dwd);
    #return Float64.(sum(self.weight .* (self.data .- LensletLaserImage(self.model,self.laser)).^2))
end =#

 Zygote.@nograd  function updateAmplitude(N::Int,spots::AbstractArray{T},data::AbstractArray{T},weight::AbstractArray{T}) where T<:AbstractFloat
    A = @MMatrix zeros(Float64,N,N)
    b = @MVector zeros(Float64,N)
    mw = similar(spots);
    @inbounds for index=1:N
        mw[:,:,index] .=  spots[:,:,index].* weight ;
        b[index] = sum(mw[:,:,index].* data );
        A[index,index] = sum(mw[:,:,index].* spots[:,:,index]);
        for i=1:index-1
            A[i,index] = A[i,index] = sum(mw[:,:,index].* spots[:,:,i])
        end
    end
    return inv(A)*b
end
end
