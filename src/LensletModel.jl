"""
    LensletModel(bbox::BoundingBox{Int},disp_model::DispModel)

Model of a lenslet
The image of a lenslet on the detector is decribed by:
* `bbox` the boundingbox of its influence on the detector
* `disp_model` the dispersion model described by a object of type `DispModel`
"""
struct LensletModel
    bbox::BoundingBox{Int}  # Boundingbox of influence of the lenslet on the detector
    disp_model::DispModel       # dispersion model of the lenslet
    profile::ProfileModel   # intensity profile in the lenslet
    function LensletModel(bbox, disp_model, profile)
        disp_model.order  ≥ 0      || throw(ArgumentError)
        profile.order ≥ 0      || throw(ArgumentError)
        disp_model.λ0 ≈ profile.λ0 || throw(ArgumentError)
        new(bbox, disp_model, profile)
    end
end


"""
    LensletModel(::BoundingBox{Int}, λ0::Float64, disp_order::Int, profile_order::Int)

Lenslet model constructor
* `bbox` : bounding box of the lenslet on the detector
* `λ0`  : reference wavelength
* `disp_order` : DispModel order of the polynomials
* `prof_order` : ProfileModel order of the polynomials
"""
function LensletModel(bbox::BoundingBox{Int}, λ0::Float64, disp_order::Int, profile_order::Int)
    LensletModel(bbox, DispModel(λ0, disp_order), ProfileModel(λ0, profile_order))
end
