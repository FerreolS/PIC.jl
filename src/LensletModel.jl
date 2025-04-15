"""
    LensletModel(bbox::BoundingBox{Int},disp_model::DispModel)

Model of a lenslet
The image of a lenslet on the detector is decribed by:
* `bbox` the boundingbox of its influence on the detector
* `disp_model` the dispersion model described by a object of type `DispModel`
"""
struct LensletModel
    bbox::BoundingBox{Int}
    disp_model::DispersionModel
    profile_model::ProfileModel
    function LensletModel(bbox, disp_model, profile_model)
        disp_model.order  ≥ 0                || throw(ArgumentError)
        profile_model.order ≥ 0              || throw(ArgumentError)
        disp_model.λref ≈ profile_model.λref || throw(ArgumentError)
        new(bbox, disp_model, profile_model)
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
function LensletModel(bbox::BoundingBox{Int}, λref::Float64, disp_order::Int, profile_order::Int)
    LensletModel(bbox, DispersionModel(λref, disp_order), ProfileModel(λref, profile_order))
end



