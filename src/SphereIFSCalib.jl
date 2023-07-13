module SphereIFSCalib

using Zygote, StaticArrays,StatsBase, LinearAlgebra
using TwoDimensional, ProgressMeter, OptimPackNextGen
include("DispModel.jl")
include("ProfileModel.jl")


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
    profile::ProfileModel   # intensity profile in the lenslet
end


"""
    lmod = LensletModel(λ0::Float64, order::Int, bbox::BoundingBox{Int})

Lenslet model constructor
* `λ0`  : reference wavelength
* `order` : order of the polynomials
* `bbox` : bounding box of the lenslet on the detector
"""
function LensletModel(λ0::Float64, disporder::Int, profileorder::Int, bbox::BoundingBox{Int})
    LensletModel(bbox, DispModel(λ0, disporder), ProfileModel(λ0, profileorder))
end


# """
#     lmod = LensletModel(λ0::Float64, order::Int, bbox::BoundingBox{Int},cx0::Float64,cy0::Float64)

# Lenslet model constructor
# * `λ0`  : reference wavelength
# * `order` : order of the polynomials
# * `bbox` : bounding box of the lenslet on the detector
# * `cx0` :
# * `cy0` :
# """
# function LensletModel(λ0::Float64, order::Int, bbox::BoundingBox{Int},cx0::Float64,cy0::Float64)
#     cx = zeros(Float64, order + 1); # coefficients of the polynomial along the x axis
#     cy = zeros(Float64, order + 1); # coefficients of the polynomial along the x axis
#     cx[1] = cx0;
#     cy[1] = cy0;
#     LensletModel(bbox, DispModel(λ0, order, cx, cy))
# end


# """
#     lmod = LensletModel(λ0::Float64, order::Int, bbox::BoundingBox{Int})

# Lenslet model constructor
# * `λ0`  : reference wavelength
# * `order` : order of the polynomials
# * `bbox` : bounding box of the lenslet on the detector
# """
# function LensletModel(λ0::Float64, order::Int,cx0::Float64,cy0::Float64, widthx::Number, widthy::Number)
#     cx = zeros(Float64, order + 1); # coefficients of the polynomial along the x axis
#     cy = zeros(Float64, order + 1); # coefficients of the polynomial along the x axis
#     cx[1] = cx0;
#     cy[1] = cy0;
#     bbox = round(Int,BoundingBox(xmin=cx0-widthx, ymin=cy0-widthy, xmax=cx0+widthx, ymax=cy0+widthy));
#     LensletModel(bbox, DispModel(λ0, order, cx, cy))
# end




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

"""
    GaussianModel2(A::Float64, fwhm::Float64, x::Real)

Compute the value at position sqrt(x) 1D centered Gaussian
* `A` : amplitude at x = 0
* `fwhm` : full-width at half maximum
* `x`: sampled position

Equivalent to `GaussianModel(A, fwhm, sqrt(x))`
"""
function GaussianModel2(A::T, fwhm::T, x::T) where (T<:Real)
    local fwhm2sigma = T(1 / (2 * sqrt(2 * log(2.))))
    return A * exp(-x / (2 * (fwhm * fwhm2sigma )^2));
end

"""
    GaussianModel2(fwhm::Float64, x::AbstractArray)

Compute the value at position sqrt(r) 1D centered Gaussian
* `fwhm` : full-width at half maximum
* `x`: array of the squared sampled position

Equivalent to `GaussianModel(1.,fwhm, sqrt.(x))`
"""
function GaussianModel2(fwhm::T, x::AbstractArray{T})where (T<:Real)
    local fwhm2sigma = T(1 / (2 * sqrt(2 * log(2.))))
    return exp.(-x ./ (2 * (fwhm * fwhm2sigma )^2));
end

"""
    GaussianModel2(fwhm::Float64, x::AbstractArray)

Compute the value at position sqrt(r) 1D centered Gaussian
* `fwhm` : full-width at half maximum
* `x`:  squared sampled position

Equivalent to `GaussianModel(1.,fwhm, sqrt(x))`
"""
function GaussianModel2(fwhm::T, x::T) where (T<:Real)
    local fwhm2sigma =T(1 / (2 * sqrt(2 * log( 2))))
    return exp(-x / (2 * (fwhm * fwhm2sigma )^2));
end

GaussianModel2(tpl::Tuple{T, T}) where (T<:Real) = GaussianModel2(tpl...)

"""
    GaussianModel2!(ret::AbstractArray{T},fwhm::Float64, x::AbstractArray)

Compute inplace the value at position sqrt(r) 1D centered Gaussian
* `ret` : output array
* `fwhm` : full-width at half maximum
* `x`:  squared sampled position

Equivalent to `GaussianModel(1.,fwhm, sqrt(x))`
"""
function GaussianModel2!(ret::AbstractArray{T},fwhm, x::AbstractArray{T}) where (T<:Real)
        @. ret = exp(-x / T(2 * (fwhm * 1 / (2 * sqrt(2 * log(2.))) )^2));
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
    LensletLaserImage!(spotsmodel::Array{Float64,3},lmodel::LensletModel,laser::LaserModel)

Build inplace the image of a lenslet under laser illumination
* `ret` : output array
* `lmodel`: model of the lenslet
* `laser`: model of the laser illumination
"""
function LensletLaserImage!(spotsmodel::Array{Float64,3},lmodel::LensletModel,laser::LaserModel)
    bbox = lmodel.bbox;
    (rx,ry) = axes(bbox) # extracting bounding box range
    @inbounds for (index, λ) in enumerate(laser.λlaser)  # For all laser
        (mx, my)  = lmodel.dmodel(λ);  # center of the index-th Gaussian spot
        r = ((rx.-mx).^2) .+ ((ry.-my).^2)';
        spotsmodel[:,:,index]= GaussianModel2( laser.fwhm[index], r)
    end
    nothing
end


"""
    LikelihoodIFS(model::LensletModel,wavelengths::AbstractArray{<:Real,1},data::AbstractArray,weight::AbstractArray)

Build the likelihood function for a given lenslet
* `lmodel`: model of the lenslet
* `laser`: wavelengths of the illumination lasers
* `data` : data
* `weight`: precision (ie inverse variance) of the data
"""
struct LikelihoodIFS{T<:Real}
    nλ::Int
    model::LensletModel
    wavelengths::Array{T,1}
    data::Array{T,2}
    weight::Array{T,2}
    spots::Array{T,3}
    amplitude::Array{T,1}#MVector{N, Float64}
    # Inner constructor provided to force using outer constructors.
    function LikelihoodIFS{T}(model::LensletModel,
        wavelengths::Array{T,1},
        data::Array{T,2},
        weight::Array{T,2}) where {T<:Real}
        nλ =length(wavelengths);
        @assert nλ > model.dmodel.order " the order of the law must be less than the number of laser"
        @assert size(data) == size(weight)
        spots = zeros(Float64,size(round(model.bbox))...,nλ)
        amplitude =  zeros(Float64,nλ)
        return new{T}(nλ,model,wavelengths,data, weight,spots,amplitude)
    end
    # Inner constructor provided to force using outer constructors.
    function LikelihoodIFS{T}(model::LensletModel,
        wavelengths::Array{T,1},
        data::Array{T,2},
        weight::T) where {T<:Real}
        nλ =length(wavelengths);
        #@assert laser.nλ == N
        @assert nλ > model.dmodel.order " the order of the law must be less than the number of laser"
        spots = zeros(Float64,size(round(model.bbox))...,nλ)
        amplitude =  zeros(Float64,nλ)
        return new{T}(nλ,model,wavelengths,data, weight*ones(1,1),spots,amplitude)
    end
end

function LikelihoodIFS(model::LensletModel,wavelengths::AbstractArray{<:Real,1},data::AbstractArray{<:Real,2})
    T = float(eltype(data))
    LikelihoodIFS{T}(model,convert(Array{T,1},wavelengths),convert(Array{T,2},data),T(1.0))
end

function LikelihoodIFS(model::LensletModel,
                        wavelengths::AbstractArray{<:Real,1},
                        data::AbstractArray{<:Real,2},
                        weight::Union{Real,AbstractArray{<:Real,2}})
    T = float(promote_type(eltype(data),eltype(weight)))
    LikelihoodIFS{T}(model,convert(Array{T,1},wavelengths),convert(Array{T,2},data), T.(weight))
end

"""
    (self::LikelihoodIFS)(x::Vector{Float64})
    compute the likelihood for a given lenslet for the parameters `x`

    ### Example
    ```
    nλ = length(λlaser)
    lenslet = LensletModel(λ0,nλ-1,round(bbox))
    xinit = vcat([fwhminit[:],cinit[:]]...)
    lkl = LikelihoodIFS(lenslet,λlaser,view(data,lenslet.bbox), view(weight,lenslet.bbox))
    xopt = vmlmb(lkl, xinit; verb=50)
    ```
"""
function  (self::LikelihoodIFS)(x::Vector{T})::Float64 where (T<:Real)
    (fwhm::Vector{T},c::Matrix{T}) = (x[1:(self.nλ)],reshape(x[(self.nλ+1):(3*self.nλ)],2,:));
    self(fwhm,c)
end

function  (self::LikelihoodIFS)(fwhm::Array{T,1},C::Array{T,2})::Float64 where (T<:Real)
    #@assert length(fwhm)== self.laser.nλ "length(fwhm) must equal to the number of lasers"
    UpdateDispModel(self.model.dmodel, C);
    bbox = self.model.bbox;
    (rx,ry) = axes(bbox) # extracting bounding box range
    m = Zygote.Buffer(self.spots);
    @inbounds for (index, λ) in enumerate(self.wavelengths)  # For all laser
        (mx, my)  = self.model.dmodel(λ);  # center of the index-th Gaussian spot
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

 SimpleGauss(x,center::Float64,fwhm::Float64) = exp(-(x-center)^2 / (2 * (fwhm * Float64(1) / (2 * sqrt(2 * log(2.))) )^2));


#= function  (self::LikelihoodIFS)(fwhm::Array{Float64,1},C::Array{Float64,2})::Float64
    # @assert length(fwhm)== self.laser.nλ "length(fwhm) must equal to the number of lasers"
    UpdateDispModel(self.model.dmodel, C);
    bbox = self.model.bbox;
    (rx,ry) = axes(bbox) # extracting bounding box range
    spotpos =  self.model.dmodel.(self.wavelengths);

   # spots  = [SimpleGauss(x,spotpos[index][1],fwhm[index]).* SimpleGauss(y, spotpos[index][2], fwhm[index]) for x in rx, y in ry, index in 1:self.nλ]
    spots  = [GaussianModel2(fwhm[index], ((x -spotpos[index][1])^2) + ((y-spotpos[index][2])^2)) for x in rx, y in ry, index in 1:self.nλ]

    Zygote.@ignore  self.amplitude .= updateAmplitude(self.nλ,spots,self.data,self.weight)
    sumspot =   zeros(Float64,size(round(bbox)));
    @inbounds for i =1:self.nλ
        sumspot += self.amplitude[i] *spots[:,:,i]
    end
    return Float64.(sum(self.weight .* (self.data .-sumspot).^2))
 end =#

#=  function  (self::LikelihoodIFS)(fwhm::Array{Float64,1},C::Array{Float64,2})::Float64
    # @assert length(fwhm)== self.laser.nλ "length(fwhm) must equal to the number of lasers"
     UpdateDispModel(self.model.dmodel, C);
     bbox = self.model.bbox;
    (rx,ry) = axes(bbox) # extracting bounding box LinRange
    r = [rx,ry];
    nx =length(rx);
    ny =length(ry);
    mx = reinterpret(reshape, Float64, self.model.dmodel.(self.wavelengths))
    p =   (broadcast((a,b,c) ->SimpleGauss.(a,b,c) ,r,mx,reshape(fwhm,(1,self.nλ))));
    spots =  reshape(hcat(broadcast((a,b) -> a*b',p[1,1:self.nλ],p[2,1:self.nλ])...),nx,ny,self.nλ)
    Zygote.@ignore  self.amplitude .= updateAmplitude(self.nλ,spots,self.data,self.weight)
    sumspot =   zeros(Float64,size(round(bbox)));
    @inbounds for i =1:self.nλ
        sumspot += self.amplitude[i] *spots[:,:,i]
    end
    return Float64.(sum(self.weight .* (self.data .-sumspot).^2))
 end =#

 """
        updateAmplitude(nλ,m,d,W)

    return the `nλ` amplitudes `a` according the the model `m`, the data and the precision `W`
    such that
    `a = argmin_a || a*m - D||^2_W`
    where
    * `nλ` : is the number of spots in the model
    * `m`:  is the model composed of `nλ` images of spots
    * `d`:  is the data
    * `W`: is the precision (inverse variance) of the data
"""
function updateAmplitude(N::Int,spots::AbstractArray{T},data::AbstractArray{T},weight::AbstractArray{T}) where T<:Real
    A = @MMatrix zeros(Float64,N,N)
    b = @MVector zeros(Float64,N)
    mw = similar(spots);
    @inbounds for index=1:N
        mw[:,:,index] .=  spots[:,:,index].* weight ;
        b[index] = sum(mw[:,:,index].* data );
        A[index,index] = sum(mw[:,:,index].* spots[:,:,index]);
        for i=1:index-1
            A[i,index] = A[index,i] = sum(mw[:,:,index].* spots[:,:,i])
        end
    end
    return inv(A)*b
end

"""
        (lenslettab, atab, fwhmtab,ctab) = fitSpectralLaw(laserData,weights,λlaser,lensletsize,cx0,cy0,cinit,fwhminit;validlenslets=true);

    fits the spectral of all lenslet identified as valid in the `valid` vector.
    * `laserData` : is the laser calibration data
    * `weight`:  is the precision (inverse variance) of the data
    * `valid`:  a boolean vector indicating valid lenslet
    * `λlaser`: is the vector of the wavelength of the lasers
    * `lensletsize` :  is a 4 Int tuple giving the size of lenslet on the detector
    * `position`: is the a priori position of the lenslets
    * `cxinit`: is the initial vector of the spectral law coefficients along x
    * `cxinit`: is the initial vector of the spectral law coefficients along y
    * `fwhminit`: is the initial vector of full width at half maximum of the spots
    * `validlenslets`: is an optionnal vector indicating the already known invalid lenslets
"""
function fitSpectralLaw(laserdata::Matrix{T},
                        weights::Matrix{T},
                        λlaser::Array{Float64,1},
                        lensletsize::NTuple{4, Int},
                        position::Matrix{Float64},
                        cxinit::Vector{Float64},
                        cyinit::Vector{Float64},
                        fwhminit::Array{Float64,1};
                        validlenslets::AbstractArray{Bool,1}=[true]
                        ) where T<:Real

    numberoflenslet = size(position)[1]
    if length(validlenslets)==1
        validlenslets = true(numberoflenslet)
    else
        numberoflenslet = min(length(validlenslets) ,numberoflenslet)
    end
    profileorder= 4
    nλ = length(λlaser)
    λ0 = mean(λlaser)# reference
    @assert length(fwhminit) == nλ

    (dxmin, dxmax,dymin,dymax) = lensletsize
    lenslettab = Array{Union{LensletModel,Missing}}(missing,numberoflenslet);
    atab = Array{Union{Float64,Missing}}(missing,nλ,numberoflenslet);
    fwhmtab = Array{Union{Float64,Missing}}(missing,nλ,numberoflenslet);
    ctab = Array{Union{Float64,Missing}}(missing,2,nλ,numberoflenslet);
    p = Progress(numberoflenslet; showspeed=true)
    Threads.@threads for i in findall(validlenslets)
        lensletbox = round(Int, BoundingBox(position[i,1]-dxmin, position[i,1]+dxmax, position[i,2]-dymin, position[i,2]+dymax));

        lenslettab[i] = LensletModel(λ0,nλ-1,profileorder, lensletbox);
        Cinit= [ [position[i,1] cxinit...]; [position[i,2] cyinit...] ];
        xinit = vcat([fwhminit[:],Cinit[:]]...);
        laserDataView = view(laserdata, lensletbox);
        weightView = view(weights,lensletbox);
        lkl = LikelihoodIFS(lenslettab[i],λlaser, laserDataView,weightView);
        cost(x::Vector{Float64}) = lkl(x);
        local xopt
        try
            xopt = vmlmb(cost, xinit; verb=false,ftol = (0.0,1e-8),maxeval=500,autodiff=true);
        catch e
            @debug showerror(stdout, e)
            @debug "Error on lenslet  $i"
            continue
        end
        (fwhmopt,copt) = (xopt[1:(nλ)],reshape(xopt[(nλ+1):(3*nλ)],2,:));
        atab[:,i] = lkl.amplitude;
        fwhmtab[:,i] = fwhmopt;
        ctab[:,:,i] = copt;
        next!(p);
    end
    ProgressMeter.finish!(p);
    return (lenslettab, atab, fwhmtab,ctab);
end

function fitSpectralLawAndProfile(laserdata::Matrix{T},
    laserweights::Matrix{T},
    lampdata::Matrix{T},
    lampweights::Matrix{T},
    λlaser::Array{Float64,1},
    lensletsize::NTuple{4, Int},
    position::Matrix{Float64},
    cxinit::Vector{Float64},
    cyinit::Vector{Float64},
    fwhminit::Array{Float64,1},
    wavelengthrange::AbstractArray{Float64,1};
    validlenslets::AbstractArray{Bool,1}=[true],
    profileorder::Int = 3
    ) where T<:Real

    numberoflenslet = size(position)[1]
    if length(validlenslets)==1
        validlenslets = true(numberoflenslet)
    else
        numberoflenslet = min(length(validlenslets) ,numberoflenslet)
    end

    nλ = length(λlaser)
    λ0 = mean(λlaser)# reference
    @assert length(fwhminit) == nλ

    (dxmin, dxmax,dymin,dymax) = lensletsize
    lenslettab = Array{Union{LensletModel,Missing}}(missing,numberoflenslet);
    laserAmplitude = Array{Union{Float64,Missing}}(missing,nλ,numberoflenslet);
    laserfwhm = Array{Union{Float64,Missing}}(missing,nλ,numberoflenslet);
    laserdist = Array{Union{Float64,Missing}}(missing,2048,2048);
    λMap =  Array{Union{Float64,Missing}}(missing,2048,2048);
    p = Progress(numberoflenslet; showspeed=true)
    Threads.@threads for i in findall(validlenslets)
        lensletbox = round(Int, BoundingBox(position[i,1]-dxmin, position[i,1]+dxmax, position[i,2]-dymin, position[i,2]+dymax));

        lenslettab[i] = LensletModel(λ0,nλ-1, profileorder,lensletbox);

        # Fit spectral law

        Cinit= [ [position[i,1] cxinit...]; [position[i,2] cyinit...] ];
        xinit = vcat([fwhminit[:],Cinit[:]]...);
        laserDataView = view(laserdata, lensletbox);
        laserWeightView = view(laserweights,lensletbox);
        spectrallkl = LikelihoodIFS(lenslettab[i],λlaser, laserDataView,laserWeightView);
        cost(x::Vector{Float64}) = spectrallkl(x);
        local xopt
        try
            xopt = vmlmb(cost, xinit; verb=false,ftol = (0.0,1e-8),maxeval=500,autodiff=true);
        catch e
            @debug showerror(stdout, e)
            @debug "Error on lenslet  $i"
            continue
        end
        fwhm = xopt[1:nλ]
        laserAmplitude[:,i] = spectrallkl.amplitude;
        laserfwhm[:,i] = fwhm
        (dist, pixλ) = distanceMap(wavelengthrange,lenslettab[i]);
        view(laserdist,lensletbox) .= dist;
        view(λMap,lensletbox) .= pixλ;

        # Fit profile
        
        lampDataView = view(lampdata, lensletbox);
        lampWeightView = view(lampweights,lensletbox);
        profilecoefs = zeros(2,profileorder)
        profilecoefs[1,1] = maximum(fwhm)
        profilecoefs[2,1:nλ] .= lenslettab[i].dmodel.cx
        pmodel  =  ProfileModel(λ0,profilecoefs)
        profilelkl = LikelihoodProfile{Float64}(pmodel,lampDataView,lampWeightView,pixλ,lensletbox)

        try
            vmlmb!(profilelkl, profilecoefs; verb=false,ftol = (0.0,1e-8),maxeval=500,autodiff=true);
        catch e
            @debug showerror(stdout, e)
            @debug "Error on lenslet  $i"
            continue
        end
        lenslettab[i] = LensletModel(lensletbox,lenslettab[i].dmodel, ProfileModel(λ0, profilecoefs))
        next!(p);
    end
    ProgressMeter.finish!(p);
    return (lenslettab, laserAmplitude, laserfwhm,laserdist, λMap);
end


function distanceMap(wavelengthrange::AbstractArray{Float64,1},
                    lenslet::LensletModel
                    )
    bbox = lenslet.bbox;
    dist = ones(Float64,size(round(bbox))).*1000;
    pixλ = ones(Float64,size(round(bbox)));
    (ax,ay) = axes(bbox)
    previous_index = 0;
    for I in CartesianIndices(dist)
        previous_index = max(1,previous_index-5);
        for  (index,λ) in enumerate(wavelengthrange[previous_index:end])
            (mx, my)  = lenslet.dmodel(λ)
            rx = ax[I[1]]-mx;
            ry = ay[I[2]]-my;
            r = sign(rx) * sqrt(rx^2 + ry^2);
            if abs(r) < abs(dist[I[1],I[2]])
                dist[I[1],I[2]] = r;
                pixλ[I[1],I[2]] = λ;
            else
                previous_index = previous_index + index-1;
                break
            end
        end
    end
    return (dist,pixλ)
end

function updateAmplitude(profile,data::Matrix{T},weight::Matrix{T}) where T<:AbstractFloat
    A = similar(data)
    b = similar(data)

    @. b = profile * data * weight
    @. A = profile^2 * weight
    A = sum(A,dims=1)
    b = sum(b,dims=1)
    zA = (A .== T(0)).||(b.<=T(0))
    if any(zA)
        A[zA] .=1
        b[zA] .=0
    end
    
    return b ./ A
end

function updateAmplitudeAndBackground(profile,data::Matrix{T},weight::Matrix{T}) where T<:AbstractFloat
    
    c = @. profile *  weight
    b = @. profile * data * weight
    a = @. profile^2 * weight
    a = sum(a,dims=1)[:]
    b = sum(b,dims=1)[:]
    c = sum(c,dims=1)[:]
    za = (a .== T(0)).||(b.<=T(0))
    if any(za)
        a[za] .=1
        b[za] .=0
    end
    

    N = length(a)
    A = Matrix{T}(undef,N+1,N+1)
    A[1,1] = sum(weight)
    A[1,2:end] .= A[2:end,1] .= c[:]
    A[2:end,2:end] .= diagm(a)

    b =  vcat(sum(data .* weight),b[:])

    return  inv(A)*b
end


end
