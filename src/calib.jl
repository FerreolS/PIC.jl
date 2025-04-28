const NLENS = 18908

const LASERS_λS = [ 987.72e-9, 1123.71e-9, 1309.37e-9, 1545.10e-9 ]
const LASERS_FWHMS_INIT = [2.3, 2.4 , 2.7]
const λRANGE = LinRange(850e-9, 1600e-9, 10000) # coarse wavelength range of the instrument

const LASERS_ORDER_DEFAULT = 2
const LAMP_ORDER_DEFAULT = 2

const BBOX_DX_LOWER_DEFAULT = 2
const BBOX_DX_UPPER_DEFAULT = 2
const BBOX_DY_LOWER_DEFAULT = 21
const BBOX_DY_UPPER_DEFAULT = 18

const LASERS_CXY0S_INIT_PATH = joinpath(dirname(pathof(PIC)), "lasers_cxy0s_init.txt")
const LASERS_CXY0S_INIT = readdlm(LASERS_CXY0S_INIT_PATH, Float64)
const LASERS_CX1_MEDIAN =  -0.6001811340726275
const LASERS_CX2_MEDIAN =  -0.3187688427580339
const LASERS_CY1_MEDIAN =  89.9795748752424
const LASERS_CY2_MEDIAN = -52.635157560302524

const LAMP_CFWHMS_INIT = [2.3; 2.5; 2.9]

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
      lasers_cxy0s_init ::Matrix{Float64} = LASERS_CXY0S_INIT,
      bbox_dx_lower ::Int = BBOX_DX_LOWER_DEFAULT,
      bbox_dx_upper ::Int = BBOX_DX_UPPER_DEFAULT,
      bbox_dy_lower ::Int = BBOX_DY_LOWER_DEFAULT,
      bbox_dy_upper ::Int = BBOX_DY_UPPER_DEFAULT,
      lasers_order ::Int = LASERS_ORDER_DEFAULT,
      lamp_order   ::Int = LAMP_ORDER_DEFAULT,
      lamp_cfwhms_init ::Vector{Float64} = LAMP_CFWHMS_INIT,
      valid_lenslets ::AbstractVector{Bool} = trues(NLENS)
)
    size(lasers_data) == size(lasers_weights) == (2048,2048) || throw(ArgumentError)
    size(lamp_data)   == size(lamp_weights)   == (2048,2048) || throw(ArgumentError)
    NLENS ≥ 1                                                || throw(ArgumentError)
    nλ ≥ 2                                                   || throw(ArgumentError)
    size(lasers_λs) == size(lasers_fwhms_init) == (nλ,)      || throw(ArgumentError)
    size(lasers_cxy0s_init) == (NLENS,2)                     || throw(ArgumentError)
    bbox_dx_lower ≥ 0                                        || throw(ArgumentError)
    bbox_dx_upper ≥ 0                                        || throw(ArgumentError)
    bbox_dy_lower ≥ 0                                        || throw(ArgumentError)
    bbox_dy_upper ≥ 0                                        || throw(ArgumentError)
    lamp_order ≥ 1                                           || throw(ArgumentError)
    length(lamp_cfwhms_init) == lamp_order + 1               || throw(ArgumentError)
    size(valid_lenslets) == (NLENS,)                         || throw(ArgumentError)

    bbox_width  = bbox_dx_lower + 1 + bbox_dx_upper
    bbox_height = bbox_dy_lower + 1 + bbox_dy_upper

    bboxes              = fill(BoundingBox{Int}(-1, -1, -1, -1), NLENS)
    lasers_cxs          = fill(NaN64, lasers_order+1, NLENS)
    lasers_cys          = fill(NaN64, lasers_order+1, NLENS)
    lasers_fwhms        = fill(NaN64, nλ, NLENS)
    lasers_amplitudes   = fill(NaN64, nλ, NLENS)
    lasers_pixels_dists = fill(NaN64, bbox_width, bbox_height, NLENS)
    lasers_pixels_λs    = fill(NaN64, bbox_width, bbox_height, NLENS)
    lamp_cfwhms         = fill(NaN64, lamp_order+1, NLENS)
    lamp_cxs            = fill(NaN64, lamp_order+1, NLENS)
    lamp_backs          = fill(NaN64, NLENS)
    lamp_amplitudes     = fill(NaN64, bbox_height, NLENS)

    p = Progress(NLENS; showspeed=true)

    assigned_lenslets = copy(valid_lenslets)

    Threads.@threads for i in findall(valid_lenslets)

        bbox = round(Int, BoundingBox(
            (lasers_cxy0s_init[i,1] - bbox_dx_lower), (lasers_cxy0s_init[i,1] + bbox_dx_upper),
            (lasers_cxy0s_init[i,2] - bbox_dy_lower), (lasers_cxy0s_init[i,2] + bbox_dy_upper)),
            RoundNearestTiesUp) # rounding mode to preserve bbox size

        bboxes[i] = bbox

        size(bbox) == (bbox_width, bbox_height) || error()

        if ((bbox.xmin < 1) | (bbox.xmax > 2048) | (bbox.ymin < 1) | (bbox.ymax > 2048))
            assigned_lenslets[i] = false
            continue
        end

        # Fit Lasers

        lens_lasers_data = view(lasers_data, bbox)
        lens_lasers_weights = view(lasers_weights, bbox)

        lasers_lkl = Lasers_LKL(
            nλ, lasers_order, lasers_λs, λref, bbox, lens_lasers_data, lens_lasers_weights)
        
        lasers_cxs_init = [ lasers_cxy0s_init[i,1] ;
                          LASERS_CX1_MEDIAN * (λref*1e6) ;
                          LASERS_CX2_MEDIAN * (λref*1e6)^2 ]

        lasers_cys_init = [ lasers_cxy0s_init[i,2] ;
                          LASERS_CY1_MEDIAN * (λref*1e6) ;
                          LASERS_CY2_MEDIAN * (λref*1e6)^2 ]

        vmlmbvars = encode_lasers_lkl_vmlmbvars(lasers_fwhms_init, lasers_cxs_init, lasers_cys_init)

        try
            vmlmb!(lasers_lkl, vmlmbvars; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug "Error on lenslet $i" exception=(e,catch_backtrace())
            assigned_lenslets[i] = false
            continue
        end
        
        (fit_lasers_fwhms, fit_lasers_cxs, fit_lasers_cys) = decode_lasers_lkl_vmlmbvars(
            lasers_lkl.nλ, vmlmbvars)
            
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
        vmlmbvars = encode_lamp_lkl_vmlmbvars(lamp_cfwhms_init, lamp_cxs_init)
        try
            vmlmb!(lamp_lkl, vmlmbvars
                   ; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
        catch e
            @debug "Error on lenslet $i" exception=(e,catch_backtrace())
            assigned_lenslets[i] = false
            continue
        end
        (fit_lamp_cfwhms, fit_lamp_cxs) = decode_lamp_lkl_vmlmbvars(vmlmbvars)

        (cost, fit_lamp_back, fit_lamp_amplitudes) = compute_lamp_cost_and_back_and_amplitudes(
            lamp_lkl, fit_lamp_cfwhms, fit_lamp_cxs)

        lamp_cfwhms[:,i] .= fit_lamp_cfwhms
        lamp_cxs[:,i] .= fit_lamp_cxs
        lamp_backs[i] = fit_lamp_back
        lamp_amplitudes[:,i] .= fit_lamp_amplitudes
            
        next!(p)
    end
    ProgressMeter.finish!(p)
    
    (; nλ, lasers_λs, λref, lasers_order, lamp_order, bbox_dx_lower, bbox_dx_upper,
       bbox_dy_lower, bbox_dy_upper, assigned_lenslets, bbox_width, bbox_height, bboxes,
       lasers_cxs, lasers_cys, lasers_fwhms, lasers_amplitudes,
       lasers_pixels_dists, lasers_pixels_λs, lamp_cfwhms, lamp_cxs, lamp_backs, lamp_amplitudes)
end
