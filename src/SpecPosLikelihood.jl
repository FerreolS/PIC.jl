struct SpecPosLikelihood{
        Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64}, Mm<:AbstractMatrix{Float64} }
    λ0    ::Float64   # reference wavelength
    order ::Int       # order of the polynomial
    box     ::BoundingBox{Int}
    data    ::Md
    weights ::Mw
    λMap    ::Mm
    
    last_amps ::Vector{Float64}
    
    function SpecPosLikelihood{Md,Mw,Mm}(
            λ0, order, box, data, weights, λMap, last_amps) where {Md,Mw,Mm}
            
        size(data) == size(box) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `box`"))
        size(data) == size(weights) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `weights`"))
        size(data) == size(λMap) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `λMap`"))
        new{Md,Mw,Mm}(λ0, order, box, data, weights, λMap, last_amps)
    end
end

function SpecPosLikelihood(
    λ0::Float64, order::Int, box::BoundingBox{Int}, data::Md, weights::Mw, λMap::Mm
) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64}, Mm<:AbstractMatrix{Float64} }

    last_amps = Vector{Float64}(undef, size(box,2) + 1)
    SpecPosLikelihood{Md,Mw,Mm}(λ0, order, box, data, weights, λMap, last_amps)
end

function decode_SpecPosLikelihood_input(M::Matrix{Float64}) ::NTuple{2,Vector{Float64}}
    cx = M[:,1]
    cλ = M[:,2]
    (cx, cλ)
end

function (self::SpecPosLikelihood)(M::Matrix{Float64})::Float64
    (cx, cλ) = decode_SpecPosLikelihood_input(M)
    
    xs = polynome_with_reference(self.λ0, cx).(self.λMap)
    fwhms = polynome_with_reference(self.λ0, cλ).(self.λMap)
    
    relative_xs = (xs .- axes(self.box,1)) .^ 2
    
    model = GaussianModel2.(relative_xs, fwhms)
    
    model = model ./ sum(model; dims=1)  # normalize at 1 along x axis
    
    Zygote.@ignore begin
        self.last_amps .= updateAmplitudeAndBackground(model, self.data, self.weights)
    end
    
    resids = self.weights .* (self.data .- self.last_amps[1] .- (model .* self.last_amps[2:end]'))
    
    return sum(abs2, resids)
end

function updateAmplitudeAndBackground(profile,data::MA,weight::MB) where {T<:AbstractFloat,MA<:AbstractMatrix{T},MB<:AbstractMatrix{T}}

    c = @. profile *  weight
    b = @. profile * data * weight
    a = @. profile^2 * weight
    a = sum(a,dims=1)[:]
    b = sum(b,dims=1)[:]
    c = sum(c,dims=1)[:]
    za = (a .== T(0)).||(b.<=T(0))
    if any(za)
        a[za] .=T(1)
        b[za] .=T(0)
        c[za] .=T(0)
    end


    N = length(a)
    A = Matrix{T}(undef,N+1,N+1)
    A[1,1] = sum(weight)
    A[1,2:end] .= A[2:end,1] .= c[:]
    A[2:end,2:end] .= diagm(a)

    b =  vcat(sum(data .* weight),b[:])

    return  inv(A)*b
end

