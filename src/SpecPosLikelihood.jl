struct SpecPosLikelihood{
        Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64}, Mm<:AbstractMatrix{Float64} }
    λ0    ::Float64   # reference wavelength
    order ::Int       # order of the polynomial
    box     ::BoundingBox{Int}
    data    ::Md
    weights ::Mw
    λMap    ::Mm
    
    last_background ::Ref{Float64}
    last_amps ::Vector{Float64}

    function SpecPosLikelihood{Md,Mw,Mm}(
            λ0, order, box, data, weights, λMap, last_background, last_amps) where {Md,Mw,Mm}
            
        size(data) == size(box) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `box`"))
        size(data) == size(weights) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `weights`"))
        size(data) == size(λMap) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `λMap`"))
        size(data,2) == length(last_amps) || throw(DimensionMismatch(
            "size of `data` incompatible with size of `last_amps`"))
        new{Md,Mw,Mm}(λ0, order, box, data, weights, λMap, last_background, last_amps)
    end
end

function SpecPosLikelihood(
    λ0::Float64, order::Int, box::BoundingBox{Int}, data::Md, weights::Mw, λMap::Mm
) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64}, Mm<:AbstractMatrix{Float64} }

    last_background = Ref(NaN64)
    last_amps = Vector{Float64}(undef, size(box,2))
    SpecPosLikelihood{Md,Mw,Mm}(λ0, order, box, data, weights, λMap, last_background, last_amps)
end

function (self::SpecPosLikelihood)(M::Matrix{Float64})::Float64

    cx = M[:,1]
    cλ = M[:,2]
    
    @inbounds begin
        xs = polynome_with_reference(self.λ0, cx).(self.λMap)
        fwhms = polynome_with_reference(self.λ0, cλ).(self.λMap)
            
        relative_xs = (xs .- axes(self.box,1)) .^ 2
        
        model = Gaussian.(relative_xs, fwhms)
        
        model = model ./ sum(model; dims=1)  # normalize at 1 along x axis
    end
    
    Zygote.@ignore begin
        (background, amps) = computeAmplitudeAndBackgroundSpecPos(model, self.data, self.weights)
        self.last_background.x = background
        self.last_amps .= amps
    end
    
    @inbounds begin
        model_amped = (model .* self.last_amps') .+ self.last_background.x
        residuals   = (self.data .- model_amped) .* self.weights
        cost        = sum(abs2, residuals)
    end
    
    cost
end

# note that the gaussians are disjoint: each gaussian fills a line on the box, and has
# absolutely no data on other lines. Note that if we had used 2D gaussians,
# the matrix A would not be invertible.
function computeAmplitudeAndBackgroundSpecPos(
    profile ::AbstractMatrix{Float64},
    data    ::AbstractMatrix{Float64},
    weights ::AbstractMatrix{Float64}
) ::Tuple{Float64, Vector{Float64}}

    a = sum(profile.^2 .* weights         ; dims=1)
    c = sum(profile    .* weights         ; dims=1)
    b = sum(profile    .* weights .* data ; dims=1)
    
    a = reshape(a, :)
    b = reshape(b, :)
    c = reshape(c, :)
    
    # replacing empty data with dummy values
    za =  iszero.(a) .| (b .≤ 0e0)
    if any(za)
        a[za] .= 1e0
        b[za] .= 0e0
        c[za] .= 0e0
    end

    N = length(a)
    A = Matrix{Float64}(undef, N+1, N+1)
    A[1,1] = sum(weights)
    A[1,2:end] .= A[2:end,1] .= c
    A[2:end,2:end] .= diagm(a)

    b = [ sum(data .* weights); b ]

    v =  inv(A)*b
    background = v[1]
    amps       = v[2:end]
    
    (background, amps)
end

