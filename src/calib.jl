const NLENS = 18908

const LASERS_λS = [ 987.72e-9, 1123.71e-9, 1309.37e-9, 1545.10e-9 ]
const LASERS_FWHMS_INIT = [2.3, 2.4 , 2.7]
const λRANGE = LinRange(850e-9, 1600e-9, 10000) # coarse wavelength range of the instrument

const DISP_ORDER = 2

const LENS_DX_LOWER = 2
const LENS_DX_UPPER = 2
const LENS_DY_LOWER = 21
const LENS_DY_UPPER = 18

const DISPERSION_CXY0S_INIT_PATH = joinpath(dirname(pathof(PIC)), "dispersion_cxy0s.txt")
const DISPERSION_CXY0S = readdlm(DISPERSION_CXY0S_INIT_PATH, Float64)
const DISPERSION_CX1_MEDIAN =  -0.6001811340726275
const DISPERSION_CX2_MEDIAN =  -0.3187688427580339
const DISPERSION_CY1_MEDIAN =  89.9795748752424
const DISPERSION_CY2_MEDIAN = -52.635157560302524

const PROFILE_CλS = [2.3; 2.5; 2.9]

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
    lasers_data    ::Matrix{<:Real},
    lasers_weights ::Matrix{<:Real},
    lamp_data      ::Matrix{<:Real},
    lamp_weights   ::Matrix{<:Real},
    ; nλ ::Int,
      lasers_λs ::Vector{Float64},
      lasers_fwhms_init ::Vector{Float64},
      λrange ::AbstractVector{Float64},
      λref ::Float64 = mean(lasers_λs),
      nlens ::Int = NLENS,
      disp_cxy0s ::Matrix{Float64} = DISPERSION_CXY0S,
      lens_dx_lower ::Int = LENS_DX_LOWER,
      lens_dx_upper ::Int = LENS_DX_UPPER,
      lens_dy_lower ::Int = LENS_DY_LOWER,
      lens_dy_upper ::Int = LENS_DY_UPPER,
      profile_order ::Int = 2,
      profile_cλs ::Vector{Float64} = PROFILE_CλS,
      valid_lenslets ::BitVector = trues(nlens)
)
    size(lasers_data) == size(lasers_weights) == (2048,2048) || throw(ArgumentError)
    size(lamp_data)   == size(lamp_weights)   == (2048,2048) || throw(ArgumentError)
    nλ ≥ 2                                                   || throw(ArgumentError)
    size(lasers_λs) == size(lasers_fwhms_init) == (nλ,)      || throw(ArgumentError)
    nlens ≥ 1                                                || throw(ArgumentError)
    size(disp_cxy0s) == (nlens,2)                            || throw(ArgumentError)
    lens_dx_lower ≥ 0                                        || throw(ArgumentError)
    lens_dx_upper ≥ 0                                        || throw(ArgumentError)
    lens_dy_lower ≥ 0                                        || throw(ArgumentError)
    lens_dy_upper ≥ 0                                        || throw(ArgumentError)
    profile_order ≥ 1                                        || throw(ArgumentError)
    length(profile_cλs) == profile_order + 1                 || throw(ArgumentError)
    size(valid_lenslets) == (nlens,)                         || throw(ArgumentError)

    lens_width  = lens_dx_lower + 1 + lens_dx_upper
    lens_height = lens_dy_lower + 1 + lens_dy_upper

    nrows_lamp_amplitudes = lens_height + 1 # 1 additional cell
    
    lenslets_models = Vector{LensletModel}(undef, nlens)
    lasers_amplitudes = Matrix{Float64}(undef, nλ, nlens)
    lamp_amplitudes = Matrix{Float64}(undef, nrows_lamp_amplitudes, nlens)
    lasers_fwhms = Matrix{Float64}(undef, nλ, nlens)
    lasers_dists = Matrix{Float64}(undef, 2048, 2048)
    λmap =  Matrix{Float64}(undef, 2048, 2048)

    p = Progress(nlens; showspeed=true)

    Threads.@threads for i in findall(valid_lenslets)

        bbox = round(Int, BoundingBox(
            (disp_cxy0s[i,1] - lens_dx_lower), (disp_cxy0s[i,1] + lens_dx_upper),
            (disp_cxy0s[i,2] - lens_dy_lower), (disp_cxy0s[i,2] + lens_dy_upper)),
            RoundNearestTiesUp) # rounding mode to preserve bbox size

        if size(bbox) != (lens_width,lens_height)
            @error "bbox size $(size(bbox)) should be $((lens_width,lens_height))"
            continue
        end

        ((bbox.xmin ≥ 1) & (bbox.xmax ≤ 2048) & (bbox.ymin ≥ 1) & (bbox.ymax ≤ 2048)) || continue

        lenslets_models[i] = LensletModel(bbox, λref, DISP_ORDER, profile_order);

        # Fit Dispersion

        fitvars = [ lasers_fwhms_init...                 ;
                    disp_cxy0s[i,1:2]...                 ;
                    DISPERSION_CX1_MEDIAN * (λref*1e6)   ;
                    DISPERSION_CY1_MEDIAN * (λref*1e6)   ;
                    DISPERSION_CX2_MEDIAN * (λref*1e6)^2 ;
                    DISPERSION_CY2_MEDIAN * (λref*1e6)^2 ]

        lasers_data_view = view(lasers_data, bbox)
        lasers_weights_view = view(lasers_weights, bbox)
        disp_lkl = Dispersion_LKL(bbox, lenslets_models[i].disp_model,
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
        fitvars = [profile_cλs... ; lenslets_models[i].disp_model.cx[1]; 0; 0]
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

