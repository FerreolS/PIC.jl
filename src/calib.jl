const NLENS = 18908

const LASERS_λS = [ 987.72e-9, 1123.71e-9, 1309.37e-9, 1545.10e-9 ]
const LASERS_FWHMS_INIT = [2.3, 2.4 , 2.7]
const λRANGE = LinRange(850e-9, 1600e-9, 10000) # coarse wavelength range of the instrument

const DISP_ORDER = 2

const LENS_DX_LOWER = 2
const LENS_DX_UPPER = 2
const LENS_DY_LOWER = 21
const LENS_DY_UPPER = 18

const DISPERSION_CXY0S_INIT_PATH = joinpath(dirname(pathof(PIC)), "dispersion_cxy0s_init.txt")
const DISPERSION_CXY0S = readdlm(DISPERSION_CXY0S_INIT_PATH, Float64)
const DISPERSION_CX1_MEDIAN =  -0.6001811340726275
const DISPERSION_CX2_MEDIAN =  -0.3187688427580339
const DISPERSION_CY1_MEDIAN =  89.9795748752424
const DISPERSION_CY2_MEDIAN = -52.635157560302524

const PROFILE_CλS_INIT = [2.3; 2.5; 2.9]

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
    lasers_data    ::AbstractMatrix{<:Real},
    lasers_weights ::AbstractMatrix{<:Real},
    lamp_data      ::AbstractMatrix{<:Real},
    lamp_weights   ::AbstractMatrix{<:Real},
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
      profile_cλs_init ::Vector{Float64} = PROFILE_CλS_INIT,
      valid_lenslets ::AbstractVector{Bool} = trues(nlens)
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
    length(profile_cλs_init) == profile_order + 1            || throw(ArgumentError)
    size(valid_lenslets) == (nlens,)                         || throw(ArgumentError)

    lens_width  = lens_dx_lower + 1 + lens_dx_upper
    lens_height = lens_dy_lower + 1 + lens_dy_upper

    nrows_lamp_amplitudes = lens_height + 1 # 1 additional cell
    
    lenslets_models = Vector{LensletModel}(undef, nlens)
    lasers_fwhms      = fill(NaN64, nλ, nlens)
    lasers_amplitudes = fill(NaN64, nλ, nlens)
    lasers_dists      = fill(NaN64, 2048, 2048)
    λmap              = fill(NaN64, 2048, 2048)
    lamp_amplitudes   = fill(NaN64, nrows_lamp_amplitudes, nlens)

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

        lens_lasers_data = view(lasers_data, bbox)
        lens_lasers_weights = view(lasers_weights, bbox)

        disp_lkl = Dispersion_LKL(bbox, lenslets_models[i].disp_model,
                                  lasers_λs, lens_lasers_data, lens_lasers_weights)

        disp_cxs_init = [ disp_cxy0s[i,1] ;
                          DISPERSION_CX1_MEDIAN * (λref*1e6) ;
                          DISPERSION_CX2_MEDIAN * (λref*1e6)^2 ]

        disp_cys_init = [ disp_cxy0s[i,2] ;
                          DISPERSION_CY1_MEDIAN * (λref*1e6) ;
                          DISPERSION_CY2_MEDIAN * (λref*1e6)^2 ]

        fitvars = encode_disp_lkl_fitvars(lasers_fwhms_init, disp_cxs_init, disp_cys_init)

        try
            vmlmb!(disp_lkl, fitvars; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug "Error on lenslet $i" exception=(e,catch_backtrace())
            continue
        end
        # last step of vmlmb is not necessary the chosen step
        # so we call again, to mutate fields to the chosen step values
        disp_lkl(fitvars)

        (fit_fwhm, fit_cxs, fit_cys) = decode_disp_lkl_fitvars(nλ, DISP_ORDER, fitvars)
        lasers_fwhms[:,i] .= fit_fwhm
        
        lasers_amplitudes[:,i] = disp_lkl.amplitude
        
        compute_lasers_dists_and_λmap!(λrange, lenslets_models[i], lasers_dists, λmap)

        # Fit profile
        
        lens_lamp_data = view(lamp_data, bbox)
        lens_lamp_weights = view(lamp_weights, bbox)
        profile_cxs_init = [ lenslets_models[i].disp_model.cxs[1] ; 0 ; 0 ]
        profile_lkl = Profile_LKL(bbox, lenslets_models[i].profile_model,
                                  lens_lamp_data, lens_lamp_weights, view(λmap,bbox))
        fitvars = encode_profile_lkl_fitvars(profile_cλs_init, profile_cxs_init)
        try
            vmlmb!(profile_lkl, fitvars
                   ; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug "Error on lenslet $i" exception=(e,catch_backtrace())
            continue
        end
        (fit_cλs, fit_cxs) = decode_profile_lkl_fitvars(fitvars)
        updateProfileModel!(lenslets_models[i].profile_model, fit_cλs, fit_cxs)

        lenslet_λmap = view(λmap, bbox)
        lenslet_rx = axes(bbox, 1)
        profile = @. GaussianModel2(lenslets_models[i].profile_model(lenslet_λmap, lenslet_rx))
        profile ./= sum(profile; dims=1)
        lamp_amplitudes[:,i] .= updateAmplitudeAndBackground!(
            profile, lens_lamp_data, lens_lamp_weights)
            
        next!(p)
    end
    ProgressMeter.finish!(p)
    
    (lenslets_models, lasers_fwhms, lasers_amplitudes, lasers_dists, λmap, lamp_amplitudes)
end

function compute_lasers_dists_and_λmap!(
    λrange::AbstractVector{Float64}, lenslet::LensletModel, lasers_dists, λmap
) ::Nothing

    previous_index = 0
    for I in CartesianIndices(lenslet.bbox)
        previous_index = max(1, previous_index-5)
        for (index,λ) in enumerate(λrange[previous_index:end])
            (gaussian_cx, gaussian_cy) = lenslet.disp_model(λ)
            dist_to_gaussian_cx = I[1] - gaussian_cx
            dist_to_gaussian_cy = I[2] - gaussian_cy
            r = sign(dist_to_gaussian_cx) * sqrt(dist_to_gaussian_cx^2 + dist_to_gaussian_cy^2)
            if isnan(lasers_dists[I]) || abs(r) < abs(lasers_dists[I])
                lasers_dists[I] = r;
                λmap[I] = λ;
            else
                break
            end
            previous_index += 1
        end
    end
    
    nothing
end

