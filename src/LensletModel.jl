"""
    LensletModel(bbox::BoundingBox{Int},dmodel::DispModel)

Model of a lenslet
The image of a lenslet on the detector is decribed by:
* `bbox` the boundingbox of its influence on the detector
* `dmodel` the dispersion model described by a object of type `DispModel`
"""
struct LensletModel
    bbox::BoundingBox{Int}  # Boundingbox of influence of the lenslet on the detector
    dmodel::DispModel       # dispersion model of the lenslet
    profile::ProfileModel   # intensity profile in the lenslet
end


"""
    lmod = LensletModel(λ0::Float64, order::Int, bbox::BoundingBox{Int})

Lenslet model constructor
* `λ0`  : reference wavelength
* `order` : order of the polynomials
* `bbox` : bounding box of the lenslet on the detector
"""
function LensletModel(λ0::Float64, disporder::Int, profileorder::Int, bbox::BoundingBox{Int})
    LensletModel(bbox, DispModel(λ0, disporder), ProfileModel(λ0, profileorder))
end
