"""
    FitResult(fields...)

Structure to store the results of a call to [`fit_wavelamps_specpos`](@ref)

Some array may contain NaNs in particular the ones for which the lenslet has not been computed.

# Fields
- `lenses_boxes ::Vector{BoundingBox{Int}}`: box for each lenslet
- `λ0 ::Float64`: reference wavelength for Wave,Lamp and Spec,Pos
- `wavelamp_order ::Int`: order of the polynomes for Wave,Lamp
- `wavelamps_fits_cx ::Matrix{Float64}`: final coefficients for the x position of the center of
  the gaussian for each laser for each lenslet. Size is thus (number of lasers, number of lenslets)
- `wavelamps_fits_cy ::Matrix{Float64}`: same as `wavelamps_fits_cx` but for y position.
- `wavelamps_fits_fwhm ::Matrix{Float64}`: final fwhm for each laser for each lenslet. Size is thus
  (number of lasers, number of lenslets).
- `wavelamps_fits_background ::Vector{Float64}`: final background for each lenslet.
- `wavelamps_fits_amps ::Matrix{Float64}`: final amplitude for each laser for each lenslet. Size is
  thus (number of lasers, number of lenslets).
- `wavelamps_centers_dists ::Matrix{Float64}`: on the detector matrix, for each lenslet box, for
  each pixel of the box, the distance to the closest gaussian center. Size is thus (2048,2048).
- `wavelamps_λvals ::Matrix{Float64}`: on the detector matrix, for each lenslet box, for
  each pixel of the box, the wavelength. Size is thus (2048,2048).
- `specpos_order ::Int`: order of the polynomes for Spec,Pos.
- `specpos_fits_cx ::Matrix{Float64}`: final coefficients for the x position of the center of
  the gaussians for each lenslet Size is thus (`specpos_order`, number of lenslets).
- `specpos_fits_cλ ::Matrix{Float64}`: same as `specpos_fits_cx` but for y position.
- `specpos_fits_background ::Vector{Float64}`: final background for each lenslet.
- `sepcpos_fits_amps ::Matrix{Float64}`: final amplitude for each line of the box of each lenslet.
  Size is thus (height of the boxes, number of lenslets).
- `computed_lenses_map ::Vector{Bool}`: if `computed_lenses_map[i]` is false, then the lenslet `i`
  has not been computed, or not fully computed, because it encountered an error. The values of
  the other fields for this lenslet may thus contain NaNs.
"""
struct FitResult
    lenses_boxes ::Vector{BoundingBox{Int}}
    λ0       ::Float64

    wavelamp_order      ::Int
    wavelamps_fits_cx   ::Matrix{Float64}
    wavelamps_fits_cy   ::Matrix{Float64}
    wavelamps_fits_fwhm ::Matrix{Float64}
    wavelamps_fits_background ::Vector{Float64}
    wavelamps_fits_amps ::Matrix{Float64}

    wavelamps_centers_dists ::Matrix{Float64}
    wavelamps_λvals         ::Matrix{Float64}

    specpos_order           ::Int
    specpos_fits_cx         ::Matrix{Float64}
    specpos_fits_cλ         ::Matrix{Float64}
    specpos_fits_background ::Vector{Float64}
    sepcpos_fits_amps       ::Matrix{Float64}

    computed_lenses_map ::Vector{Bool}
end

"""
    FitResult(nb_lenses, nλ, λ0, wavelamp_order, specpos_order, box_height)

Alternative constructor that initiates arrays with correct size, filled with NaNs.
"""
function FitResult(nb_lenses, nλ, λ0, wavelamp_order, specpos_order, box_height)
    FitResult(
        Vector{BoundingBox{Int}}(undef, nb_lenses), # lenses_boxes
        λ0,
        wavelamp_order,
        fill(NaN64, wavelamp_order + 1, nb_lenses), # wavelamps_fits_cx
        fill(NaN64, wavelamp_order + 1, nb_lenses), # wavelamps_fits_cy
        fill(NaN64, nλ, nb_lenses), # wavelamps_fits_fwhm
        fill(NaN64, nb_lenses), # wavelamps_fits_background
        fill(NaN64, nλ, nb_lenses), # wavelamps_fits_amp
        fill(NaN64, 2048, 2048), # wavelamps_centers_dists
        fill(NaN64, 2048, 2048), # wavelamps_λvals
        specpos_order,
        fill(NaN64, specpos_order + 1, nb_lenses), # specpos_fits_cx
        fill(NaN64, specpos_order + 1, nb_lenses), # specpos_fits_cλ
        fill(NaN64, nb_lenses), # specpos_fits_background
        fill(NaN64, box_height, nb_lenses), # sepcpos_fits_amps
        falses(nb_lenses) # computed_lenses_map
    )
end
