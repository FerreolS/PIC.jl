"""
    GaussianModel2(fwhm::Float64, x::AbstractArray)

Compute the value at lenslets_coords sqrt(r) 1D centered Gaussian
* `fwhm` : full-width at half maximum
* `x`:  squared sampled lenslets_coords

Equivalent to `GaussianModel(1.,fwhm, sqrt(x))`
"""
function GaussianModel2(fwhm::T, x::T) ::T where {T<:Real}
    fwhm2sigma = 1 / (2 * sqrt(2 * log(2)))
    exp(-x / (2 * (fwhm * fwhm2sigma)^2))
end

function GaussianModel2(t::NTuple{2,T}) ::T where {T<:Real} return GaussianModel2(t[1], t[2]) end

"""
    Spectral_LKL(model::LensletModel,wavelengths::AbstractArray{<:Real,1},data::AbstractArray,weight::AbstractArray)

Build the likelihood function for a given lenslet
* `lmodel`: model of the lenslet
* `laser`: wavelengths of the illumination lasers
* `data` : data
* `weight`: precision (ie inverse variance) of the data
"""
struct Spectral_LKL{T<:Real,MD<:AbstractMatrix{T},MW<:AbstractMatrix{T}}
    nλ::Int
    lenslet_model::LensletModel
    lasers_λs::Vector{T}
    data::MD
    weights::MW
    spots::Array{T,3}
    amplitude::Vector{T}
    function Spectral_LKL{T,MD,MW}(
        nλ, lenslet_model, lasers_λs, data, weights, spots, amplitude
    ) where {T,MD,MW}
        length(lasers_λs) == nλ     || throw(ArgumentError)
        size(data) == size(weights) || throw(ArgumentError)
        size(spots,3) == nλ         || throw(ArgumentError)
        length(amplitude) == nλ     || throw(ArgumentError)
        nλ > lenslet_model.dmodel.order || throw(ArgumentError)
        size(spots)[1:2] == size(lenslet_model.bbox) || throw(ArgumentError)
        new{T,MD,MW}(nλ, lenslet_model, lasers_λs, data, weights, spots, amplitude)
    end
end

function Spectral_LKL(
    lenslet_model::LensletModel, lasers_λs::Vector{T}, data::MD, weights::MW
) where {T<:Real,MD<:AbstractMatrix{T},MW<:AbstractMatrix{T}}
    nλ = length(lasers_λs)
    spots = zeros(T, size(lenslet_model.bbox)..., nλ)
    amplitude = zeros(T, nλ)
    Spectral_LKL{T,MD,MW}(nλ, lenslet_model, lasers_λs, data, weights, spots, amplitude)
end

"""
    (self::Spectral_LKL)(x::Vector{Float64})
compute the likelihood for a given lenslet for the parameters `x`
"""
function  (self::Spectral_LKL)(x::Vector{T}) ::T where {T<:Real}
    (fwhm::Vector{T},c::Matrix{T}) = (x[1:(self.nλ)],reshape(x[(self.nλ+1):(3*self.nλ)],2,:));
    self(fwhm,c)
end

    UpdateDispModel(self.lenslet_model.dmodel, C);
    
    bbox = self.lenslet_model.bbox
    
    (rx,ry) = axes(bbox) # extracting bounding box range
    
    m = Zygote.Buffer(self.spots)
    @inbounds for (index,λ) in enumerate(self.lasers_λs)  # For all laser
        (mx, my)  = self.lenslet_model.dmodel(λ);  # center of the index-th Gaussian spot
        r = ((rx.-mx).^2) .+ ((ry.-my).^2)';
        m[:,:,index] = GaussianModel2.(fwhm[index], r);
    end
    spots = copy(m)
    Zygote.@ignore self.amplitude .= updateAmplitude(self.nλ, spots, self.data, self.weights)
    sumspot = zeros(T, size(bbox))
    @inbounds for i =1:self.nλ
        sumspot += self.amplitude[i] *spots[:,:,i]
    end
    return sum(self.weights .* (self.data .-sumspot).^2)
end

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
function updateAmplitude(
    N::Int, spots::AbstractArray{T,3},
    data::AbstractMatrix{T}, weight::AbstractMatrix{T}
) ::Vector{T} where {T<:Real}
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

function fitSpectralLawAndProfile(
    lasers_data::Matrix{T},
    lasers_weights::Matrix{T},
    lamp_data::Matrix{T},
    lamp_weights::Matrix{T},
    lasers_λs::Vector{Float64},
    λref::Float64,
    lenslet_size::NTuple{4,Int},
    lenslets_coords::Matrix{Float64},
    cxinit::Vector{Float64},
    cyinit::Vector{Float64},
    fwhminit::Vector{Float64},
    λrange::AbstractVector{Float64};
    valid_lenslets::AbstractVector{Bool}=trues(size(lenslets_coords,1)),
    profile_order::Int = 2,
    smalltest ::Bool = false
) where {T<:Real}

    nlens = size(lenslets_coords,1)
    nλ = length(lasers_λs)

    length(fwhminit) == nλ || throw(ArgumentError)

    (dxmin, dxmax, dymin, dymax) = lenslet_size
    nrows_lamp_amplitudes = (1 + dymin + dymax) + 1 # nrows(lenslet box) + 1 additional cell
    
    lenslets_models = Vector{LensletModel}(undef, nlens)
    lasers_amplitudes = Matrix{Float64}(undef, nλ, nlens)
    lamp_amplitudes = Matrix{Float64}(undef, nrows_lamp_amplitudes, nlens)
    lasers_fwhms = Matrix{Float64}(undef, nλ, nlens)
    lasers_dists = Matrix{Float64}(undef, 2048, 2048)
    λmap =  Matrix{Float64}(undef, 2048, 2048)
    p = Progress(nlens; showspeed=true)
    
    indices = findall(valid_lenslets)
    indices = smalltest ? rand(MersenneTwister(1234), indices, 300) : indices
    
    Threads.@threads for i in indices

        lenslet_box = round(Int, BoundingBox(
            lenslets_coords[i,1]-dxmin, lenslets_coords[i,1]+dxmax,
            lenslets_coords[i,2]-dymin, lenslets_coords[i,2]+dymax),
            RoundNearestTiesUp)

        lenslets_models[i] = LensletModel(lenslet_box, λref, nλ-1, profile_order);

        # Fit spectral law

        xinit = [fwhminit...; lenslets_coords[i,:]...; cxinit[1]; cyinit[1]; cxinit[2]; cyinit[2]]
        lasers_data_view = view(lasers_data, lenslet_box);
        lasers_weights_view = view(lasers_weights,lenslet_box);
        spectral_lkl = Spectral_LKL(
            lenslets_models[i], lasers_λs, lasers_data_view, lasers_weights_view)
        try
            vmlmb!(spectral_lkl, xinit; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug showerror(stdout, e)
            @debug "Error on lenslet $i"
            continue
        end
        lasers_amplitudes[:,i] = spectral_lkl.amplitude
        lasers_fwhms[:,i] .= view(xinit, 1:nλ)
        (dist, pixλ) = distanceMap(λrange, lenslets_models[i])
        lasers_dists[lenslet_box] .= dist
        λmap[lenslet_box] .= pixλ

        # Fit profile
        
        lamp_data_view = view(lamp_data, lenslet_box)
        lamp_weights_view = view(lamp_weights, lenslet_box)
        profile_coeffs = zeros(Float64, 2, profile_order+1)
        profile_coeffs[1,1:3] .= [2.3, 2.5, 2.9] # maximum(fwhm)
        profile_coeffs[2,1] = lenslets_models[i].dmodel.cx[1]
        profile_model = ProfileModel(λref, profile_coeffs)
        profile_lkl = Profile_LKL(
            profile_model, lamp_data_view, lamp_weights_view, pixλ, lenslet_box)
        try
            vmlmb!(profile_lkl, profile_coeffs
                   ; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug showerror(stdout, e)
            @debug "Error on lenslet $i"
            continue
        end
        profile_model = ProfileModel(λref, profile_coeffs)
        lenslets_models[i] = LensletModel(lenslet_box, lenslets_models[i].dmodel, profile_model)

        profile = @. GaussianModel2(profile_model(pixλ,($(axes(lenslet_box,1)))))
        profile ./= sum(profile; dims=1)
        lamp_amplitudes[:,i] .= updateAmplitudeAndBackground(
            profile, lamp_data_view, lamp_weights_view)
            
        next!(p)
    end
    ProgressMeter.finish!(p)
    
    (lenslets_models, lasers_amplitudes, lamp_amplitudes, lasers_fwhms, lasers_dists, λmap)
end

function distanceMap(
    λrange::AbstractVector{Float64}, lenslet::LensletModel
) ::NTuple{2,Matrix{Float64}}

    dist = fill(1000e0,  size(lenslet.bbox))
    pixλ = ones(Float64, size(lenslet.bbox))
    
    (ax, ay) = axes(lenslet.bbox)
    
    previous_index = 0
    for I in CartesianIndices(dist)
        previous_index = max(1, previous_index-5)
        for (index,λ) in enumerate(λrange[previous_index:end])
            (mx, my) = lenslet.dmodel(λ)
            rx = ax[I[1]] - mx
            ry = ay[I[2]] - my
            r = sign(rx) * sqrt(rx^2 + ry^2)
            if abs(r) < abs(dist[I])
                dist[I] = r;
                pixλ[I] = λ;
            else
                previous_index += (index - 1)
                break
            end
        end
    end

    (dist, pixλ)
end

function updateAmplitude(
    nλ::Int, data::Matrix{T}, weights::Matrix{T}
) ::Vector{T} where {T<:Real}

    A = similar(data)
    b = similar(data)

    @. b = nλ * data * weights
    @. A = nλ^2 * weights
    A = sum(A; dims=1)
    b = sum(b; dims=1)
    zA = (A .== T(0)).||(b.<=T(0))
    if any(zA)
        A[zA] .= 1
        b[zA] .= 0
    end
    
    return b ./ A
end

function updateAmplitudeAndBackground(profile,data::MA,weight::MB) where {T<:AbstractFloat,MA<:AbstractMatrix{T},MB<:AbstractMatrix{T}}
    
    c = @. profile *  weight
    b = @. profile * data * weight
    a = @. profile^2 * weight
    a = sum(a,dims=1)[:]
    b = sum(b,dims=1)[:]
    c = sum(c,dims=1)[:]
    za = (a .== T(0)).||(b.<=T(0))
    if any(za)
        a[za] .=T(1)
        b[za] .=T(0)
        c[za] .=T(0)
    end
    

    N = length(a)
    A = Matrix{T}(undef,N+1,N+1)
    A[1,1] = sum(weight)
    A[1,2:end] .= A[2:end,1] .= c[:]
    A[2:end,2:end] .= diagm(a)

    b =  vcat(sum(data .* weight),b[:])

    return  inv(A)*b
end
