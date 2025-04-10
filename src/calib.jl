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

        bbox = round(Int, BoundingBox(
            lenslets_coords[i,1]-dxmin, lenslets_coords[i,1]+dxmax,
            lenslets_coords[i,2]-dymin, lenslets_coords[i,2]+dymax),
            RoundNearestTiesUp)

        lenslets_models[i] = LensletModel(bbox, λref, nλ-1, profile_order);

        # Fit Dispersion

        fitvars = [
            fwhminit...; lenslets_coords[i,:]...; cxinit[1]; cyinit[1]; cxinit[2]; cyinit[2]]
        lasers_data_view = view(lasers_data, bbox);
        lasers_weights_view = view(lasers_weights,bbox);
        disp_lkl = Disp_LKL(lenslets_models[i].bbox, lenslets_models[i].disp_model,
                            lasers_λs, lasers_data_view, lasers_weights_view)
        try
            vmlmb!(disp_lkl, fitvars; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug showerror(stdout, e)
            @debug "Error on lenslet $i"
            continue
        end
        lasers_amplitudes[:,i] = disp_lkl.amplitude
        lasers_fwhms[:,i] .= view(fitvars, 1:nλ)
        (dist, pixλ) = distanceMap(λrange, lenslets_models[i])
        lasers_dists[bbox] .= dist
        λmap[bbox] .= pixλ

        # Fit profile
        
        lamp_data_view = view(lamp_data, bbox)
        lamp_weights_view = view(lamp_weights, bbox)
        fitvars = [2.3; 2.5; 2.9; lenslets_models[i].disp_model.cx[1]; 0; 0]
        profile_lkl = Profile_LKL(bbox, lenslets_models[i].profile_model,
                                  lamp_data_view, lamp_weights_view, pixλ)
        try
            vmlmb!(profile_lkl, fitvars
                   ; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug showerror(stdout, e)
            @debug "Error on lenslet $i"
            continue
        end
        profile_model = ProfileModel(λref, fitvars)
        lenslets_models[i] = LensletModel(bbox, lenslets_models[i].disp_model, profile_model)

        profile = @. GaussianModel2(profile_model(pixλ,($(axes(bbox,1)))))
        profile ./= sum(profile; dims=1)
        lamp_amplitudes[:,i] .= updateAmplitudeAndBackground!(
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
            (mx, my) = lenslet.disp_model(λ)
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

