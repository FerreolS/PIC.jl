const ANTHONY_CX_PATH = "resources/anthony_cx.txt"
const ANTHONY_CY_PATH = "resources/anthony_cy.txt"

function get_anthony_cxy(λ0)
    
    src_path = pathof(PIC) |> dirname
    
    cx_path = joinpath(src_path, ANTHONY_CX_PATH)
    cy_path = joinpath(src_path, ANTHONY_CY_PATH)
    
    cx = readdlm(cx_path, header = false)
    cx0 = cx[:,1] .+ 1025
    mcx1 = median(cx[:,2]) *  λ0 * 1e6
    mcx2 = median(cx[:,3]) * (λ0 * 1e6)^2

    cy = readdlm(cy_path, header = false)
    cy0 = cy[:,1] .+ 1025
    mcy1 = median(cy[:,2]) *  λ0 * 1e6
    mcy2 = median(cy[:,3]) * (λ0 * 1e6)^2
    
    (cx0, mcx1, mcx2, cy0, mcy1, mcy2)
end

function compute_first_valid_lenses_map(cx0, cy0, dxmin, dxmax, dymin, dymax)
    return (
            ((cx0 .- dxmin) .≥ 1)
        .&  ((cx0 .+ dxmax) .≤ 2048)
        .&  ((cy0 .- dymin) .≥ 1)
        .&  ((cy0 .+ dymax) .≤ 2048)
    )
end