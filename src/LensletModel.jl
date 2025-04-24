"""
    LensletModel(bbox::BoundingBox{Int},lasers_model::LasersModel)

Model of a lenslet
The image of a lenslet on the detector is decribed by:
* `bbox` the boundingbox of its influence on the detector
* `lasers_model` the lasers model described by a object of type `LasersModel`
"""
struct LensletModel
    bbox::BoundingBox{Int}
    lamp_model::ProfileModel
    function LensletModel(bbox, lamp_model)
        lamp_model.order ≥ 1              || throw(ArgumentError)
        new(bbox, lamp_model)
    end
end


"""
    LensletModel(::BoundingBox{Int}, λ0::Float64, lasers_order::Int, lamp_order::Int)

Lenslet model constructor
* `bbox` : bounding box of the lenslet on the detector
* `λ0`  : reference wavelength
* `lasers_order` : LasersModel order of the polynomials
* `prof_order` : ProfileModel order of the polynomials
"""
function LensletModel(
    bbox::BoundingBox{Int}, nλ::Int, λref::Float64, lamp_order::Int
)
    LensletModel(bbox, ProfileModel(λref, lamp_order))
end



