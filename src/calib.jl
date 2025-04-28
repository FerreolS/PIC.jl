const NLENS = 18908

const LASERS_λS = [ 987.72e-9, 1123.71e-9, 1309.37e-9, 1545.10e-9 ]
const LASERS_FWHMS_INIT = [2.3, 2.4 , 2.7]
const λRANGE = LinRange(850e-9, 1600e-9, 10000) # coarse wavelength range of the instrument

const LASERS_ORDER = 2

const LENS_DX_LOWER = 2
const LENS_DX_UPPER = 2
const LENS_DY_LOWER = 21
const LENS_DY_UPPER = 18

const LASERS_CXY0S_INIT_PATH = joinpath(dirname(pathof(PIC)), "lasers_cxy0s_init.txt")
const LASERS_CXY0S = readdlm(LASERS_CXY0S_INIT_PATH, Float64)
const LASERS_CX1_MEDIAN =  -0.6001811340726275
const LASERS_CX2_MEDIAN =  -0.3187688427580339
const LASERS_CY1_MEDIAN =  89.9795748752424
const LASERS_CY2_MEDIAN = -52.635157560302524

const PROFILE_CλS_INIT = [2.3; 2.5; 2.9]

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
      lasers_cxy0s ::Matrix{Float64} = LASERS_CXY0S,
      lens_dx_lower ::Int = LENS_DX_LOWER,
      lens_dx_upper ::Int = LENS_DX_UPPER,
      lens_dy_lower ::Int = LENS_DY_LOWER,
      lens_dy_upper ::Int = LENS_DY_UPPER,
      lasers_order ::Int = LASERS_ORDER,
      lamp_order ::Int = 2,
      lamp_cλs_init ::Vector{Float64} = PROFILE_CλS_INIT,
      valid_lenslets ::AbstractVector{Bool} = trues(nlens)
)
    size(lasers_data) == size(lasers_weights) == (2048,2048) || throw(ArgumentError)
    size(lamp_data)   == size(lamp_weights)   == (2048,2048) || throw(ArgumentError)
    nλ ≥ 2                                                   || throw(ArgumentError)
    size(lasers_λs) == size(lasers_fwhms_init) == (nλ,)      || throw(ArgumentError)
    nlens ≥ 1                                                || throw(ArgumentError)
    size(lasers_cxy0s) == (nlens,2)                          || throw(ArgumentError)
    lens_dx_lower ≥ 0                                        || throw(ArgumentError)
    lens_dx_upper ≥ 0                                        || throw(ArgumentError)
    lens_dy_lower ≥ 0                                        || throw(ArgumentError)
    lens_dy_upper ≥ 0                                        || throw(ArgumentError)
    lamp_order ≥ 1                                        || throw(ArgumentError)
    length(lamp_cλs_init) == lamp_order + 1            || throw(ArgumentError)
    size(valid_lenslets) == (nlens,)                         || throw(ArgumentError)

    bbox_width  = lens_dx_lower + 1 + lens_dx_upper
    bbox_height = lens_dy_lower + 1 + lens_dy_upper

    lenslets_models = Vector{LensletModel}(undef, nlens)
    
    lasers_cxs          = fill(NaN64, lasers_order+1, nlens)
    lasers_cys          = fill(NaN64, lasers_order+1, nlens)
    lasers_fwhms        = fill(NaN64, nλ, nlens)
    lasers_amplitudes   = fill(NaN64, nλ, nlens)
    lasers_pixels_dists = fill(NaN64, bbox_width, bbox_height, nlens)
    lasers_pixels_λs    = fill(NaN64, bbox_width, bbox_height, nlens)
    lamp_cλs            = fill(NaN64, lamp_order+1, nlens)
    lamp_cxs            = fill(NaN64, lamp_order+1, nlens)
    lamp_backs          = fill(NaN64, nlens)
    lamp_amplitudes    = fill(NaN64, bbox_height, nlens)

    p = Progress(nlens; showspeed=true)

    assigned_lenslets = copy(valid_lenslets)

    Threads.@threads for i in findall(valid_lenslets)

        bbox = round(Int, BoundingBox(
            (lasers_cxy0s[i,1] - lens_dx_lower), (lasers_cxy0s[i,1] + lens_dx_upper),
            (lasers_cxy0s[i,2] - lens_dy_lower), (lasers_cxy0s[i,2] + lens_dy_upper)),
            RoundNearestTiesUp) # rounding mode to preserve bbox size

        size(bbox) == (bbox_width, bbox_height) || error()

        if ((bbox.xmin < 1) | (bbox.xmax > 2048) | (bbox.ymin < 1) | (bbox.ymax > 2048))
            assigned_lenslets[i] = false
            continue
        end

        lenslets_models[i] = LensletModel(bbox, nλ, λref, lamp_order);

        # Fit Lasers

        lens_lasers_data = view(lasers_data, bbox)
        lens_lasers_weights = view(lasers_weights, bbox)

        lasers_lkl = Lasers_LKL(
            nλ, lasers_order, lasers_λs, λref, bbox, lens_lasers_data, lens_lasers_weights)
        
        lasers_cxs_init = [ lasers_cxy0s[i,1] ;
                          LASERS_CX1_MEDIAN * (λref*1e6) ;
                          LASERS_CX2_MEDIAN * (λref*1e6)^2 ]

        lasers_cys_init = [ lasers_cxy0s[i,2] ;
                          LASERS_CY1_MEDIAN * (λref*1e6) ;
                          LASERS_CY2_MEDIAN * (λref*1e6)^2 ]

        fitvars = encode_lasers_lkl_fitvars(lasers_fwhms_init, lasers_cxs_init, lasers_cys_init)

        try
            vmlmb!(lasers_lkl, fitvars; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug "Error on lenslet $i" exception=(e,catch_backtrace())
            assigned_lenslets[i] = false
            continue
        end
        
        (fit_lasers_fwhms, fit_lasers_cxs, fit_lasers_cys) = decode_lasers_lkl_fitvars(
            lasers_lkl.nλ, fitvars)
            
        (cost, fit_lasers_amplitudes) = compute_lasers_cost_and_amplitudes(
            lasers_lkl, fit_lasers_cxs, fit_lasers_cys, fit_lasers_fwhms)
            
        lasers_cxs[:,i] .= fit_lasers_cxs
        lasers_cys[:,i] .= fit_lasers_cys
        lasers_fwhms[:,i] .= fit_lasers_fwhms
        lasers_amplitudes[:,i] .= fit_lasers_amplitudes

        lens_lasers_pixels_dists = view(lasers_pixels_dists,:,:,i)
        lens_lasers_pixels_λs = view(lasers_pixels_λs,:,:,i)

        compute_lasers_dists_and_λmap!(
            λrange, bbox, lasers_order, λref, fit_lasers_cxs, fit_lasers_cys,
            lens_lasers_pixels_dists, lens_lasers_pixels_λs)
            
        # Fit profile
        
        lens_lamp_data = view(lamp_data, bbox)
        lens_lamp_weights = view(lamp_weights, bbox)
        lamp_cxs_init = [ fit_lasers_cxs[1] ; 0 ; 0 ]
        lamp_lkl = Lamp_LKL(lamp_order, λref, bbox, lens_lasers_pixels_λs,
                                  lens_lamp_data, lens_lamp_weights)
        fitvars = encode_lamp_lkl_fitvars(lamp_cλs_init, lamp_cxs_init)
        try
            vmlmb!(lamp_lkl, fitvars
                   ; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug "Error on lenslet $i" exception=(e,catch_backtrace())
            assigned_lenslets[i] = false
            continue
        end
        (fit_lamp_cλs, fit_lamp_cxs) = decode_lamp_lkl_fitvars(fitvars)
        updateProfileModel!(lenslets_models[i].lamp_model, fit_lamp_cλs, fit_lamp_cxs)


        lenslet_rx = axes(bbox, 1)
        profile = @. GaussianModel2(lenslets_models[i].lamp_model(lens_lasers_pixels_λs, lenslet_rx))
        profile ./= sum(profile; dims=1)
        lamp_cλs[:,i] .= fit_lamp_cλs
        lamp_cxs[:,i] .= fit_lamp_cxs
        (fit_lamp_back, fit_lamp_amplitudes) = compute_lamp_backs_and_amps(
            profile, lens_lamp_data, lens_lamp_weights)
        lamp_backs[i] = fit_lamp_back
        lamp_amplitudes[:,i] .= fit_lamp_amplitudes
            
        next!(p)
    end
    ProgressMeter.finish!(p)
    
    (; nlens, nλ, lasers_λs, λref, lasers_order, lamp_order, lens_dx_lower, lens_dx_upper,
       lens_dy_lower, lens_dy_upper, assigned_lenslets, lenslets_models, bbox_width, bbox_height,
       lasers_cxs, lasers_cys, lasers_fwhms, lasers_amplitudes,
       lasers_pixels_dists, lasers_pixels_λs, lamp_backs, lamp_amplitudes)
end



