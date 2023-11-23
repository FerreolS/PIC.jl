struct WaveLampLikelihood{ Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    λ0 ::Float64
    λs ::Vector{Float64}
    box ::BoundingBox{Int}
    data ::Md
    weights ::Mw

    last_amp ::Vector{Float64}
    
    function WaveLampLikelihood{Md,Mw}(
        λ0, λs, box, data, weights, last_amp) where {Md,Mw}
    
        size(box) == size(data) || throw(DimensionMismatch(
            "size of `box` must match the size of `data`"))
        size(data) == size(weights) || throw(DimensionMismatch(
            "size of `data` must match the size of `weights`"))
        length(last_amp) == length(λs) || throw(DimensionMismatch(
            "length of `last_amp` must be `length(λs)`"))
        
        new{Md,Mw}(λ0, λs, box, data, weights, last_amp)
    end
end

function WaveLampLikelihood(
    λ0::Float64, λs::Vector{Float64}, box::BoundingBox{Int}, data::Md, weights::Mw
) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    
    last_amp = similar(λs)
    WaveLampLikelihood{Md,Mw}(λ0, λs, box, data, weights, last_amp)
end

function (self::WaveLampLikelihood)(M::Matrix{Float64}) ::Float64

    fwhms = M[:,1]
    cx    = M[:,2]
    cy    = M[:,3]
        
    nλ = length(self.λs)
    (xs, ys) = axes(self.box)
    
    @inbounds begin
        # spot centers
        spots_ctrs_x = polynome_with_reference(self.λ0, cx).(self.λs)
        spots_ctrs_y = polynome_with_reference(self.λ0, cy).(self.λs)
        
        # for each laser, a matrix containing the gaussian spot
        list_spots = map(1:nλ) do i
            ctr_x = spots_ctrs_x[i]
            ctr_y = spots_ctrs_y[i]
            radiusmatrix = @. (xs - ctr_x)^2 + ((ys - ctr_y)^2)'
            Gaussian.(radiusmatrix, fwhms[i])
        end
    end
    
    # amplitude of gaussian spots
    Zygote.@ignore begin
        self.last_amp .= computeAmplitudeWaveLamp(list_spots, self.data, self.weights)
        any(isnan, self.last_amp) && error("NaN amplitude: \"$(self.last_amp)\"")
    end

    @inbounds begin
        list_spots_amped = map(i -> list_spots[i] .* self.last_amp[i], 1:nλ)

        sumspots = reduce(.+, list_spots_amped)

        cost = sum(self.weights .* (self.data .- sumspots).^2)
    end
    
    cost
end

function computeAmplitudeWaveLamp(
    spots::Vector{Matrix{Float64}}, data::AbstractArray{Float64}, weights::AbstractArray{Float64})
    
    N = length(spots)
    A = @MMatrix zeros(Float64,N,N)
    b = @MVector zeros(Float64,N)
    mw = Matrix{Float64}(undef, size(data))
    
    @inbounds for i in 1:N
        mw .= spots[i] .* weights
        b[i] = sum(mw .* data)
        A[i,i] = sum(mw .* spots[i])
        for j in 1:i-1
            A[j,i] = A[i,j] = sum(mw .* spots[j])
        end
    end
    amps = inv(A)*b
    amps
end


