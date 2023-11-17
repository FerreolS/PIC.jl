const ANTHONY_YJ_CX_PATH  = "resources/anthony_cx.txt"
const ANTHONY_YJ_CY_PATH  = "resources/anthony_cy.txt"
const ANTHONY_YJH_CX_PATH = "resources/coef_pol_x.txt"
const ANTHONY_YJH_CY_PATH = "resources/coef_pol_y.txt"

function get_anthony_cxy(λ0, ifsmode)
    
    src_path = pathof(PIC) |> dirname
    
    cx_path = (ifsmode == :YJ)  ? joinpath(src_path, ANTHONY_YJ_CX_PATH)  :
              (ifsmode == :YJH) ? joinpath(src_path, ANTHONY_YJH_CX_PATH) :
              error()
  
    cy_path = (ifsmode == :YJ)  ? joinpath(src_path, ANTHONY_YJ_CY_PATH)  :
              (ifsmode == :YJH) ? joinpath(src_path, ANTHONY_YJH_CY_PATH) :
              error()
    
    order = (ifsmode == :YJ) ? 2 : (ifsmode == :YJH) ? 3 : error()
    
    cx = readdlm(cx_path, header = false)
    cx0 = cx[:,1] .+ 1025
    mcxs = map(1:order) do i; median(cx[:,i+1]) * (λ0 * 1e6)^i end

    cy = readdlm(cy_path, header = false)
    cy0 = cy[:,1] .+ 1025
    mcys = map(1:order) do i; median(cy[:,i+1]) * (λ0 * 1e6)^i end
    
    (cx0, mcxs, cy0, mcys)
end

function compute_first_valid_lenses_map(cx0, cy0, dxmin, dxmax, dymin, dymax)
    return (
            ((cx0 .- dxmin) .≥ 1)
        .&  ((cx0 .+ dxmax) .≤ 2048)
        .&  ((cy0 .- dymin) .≥ 1)
        .&  ((cy0 .+ dymax) .≤ 2048)
    )
end