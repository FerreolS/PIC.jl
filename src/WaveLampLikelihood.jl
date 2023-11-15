struct WaveLampLikelihood{ Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    λ0 ::Float64
    order ::Int
    λs ::Vector{Float64}
    box ::BoundingBox{Int}
    data ::Md
    weights ::Mw

    last_amp ::Vector{Float64}
    
    function WaveLampLikelihood{Md,Mw}(
        λ0, order, λs, box, data, weights, last_amp) where {Md,Mw}
    
        length(λs) > order || throw(DimensionMismatch(
            "`order` must be strictly less than `length(λs)`"))
        size(box) == size(data) || throw(DimensionMismatch(
            "size of `box` must match the size of `data`"))
        size(data) == size(weights) || throw(DimensionMismatch(
            "size of `data` must match the size of `weights`"))
        length(last_amp) == length(λs) || throw(DimensionMismatch(
            "length of `last_amp` must be `length(λs)`"))
        
        new{Md,Mw}(λ0, order, λs, box, data, weights, last_amp)
    end
end

function WaveLampLikelihood(
    λ0::Float64, order::Int, λs::Vector{Float64}, box::BoundingBox{Int}, data::Md, weights::Mw
) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    
    last_amp = similar(λs)
    WaveLampLikelihood{Md,Mw}(λ0, order, λs, box, data, weights, last_amp)
end

function decode_WaveLampLikelihood_input(
    nλ::Int, order::Int, V::Vector{Float64}) ::NTuple{3,Vector{Float64}}
    
    start_fwhms = 1
    end_fwhms   = nλ
    
    start_cx    = end_fwhms + 1
    end_cx      = start_cx + order
    
    start_cy    = end_cx   + 1
    end_cy      = start_cy + order
    
    fwhms = V[start_fwhms : end_fwhms]
    cx    = V[start_cx    : end_cx]
    cy    = V[start_cy    : end_cy]
    
    (fwhms, cx, cy)
end

function (self::WaveLampLikelihood)(V::Vector{Float64}) ::Float64

    nλ = length(self.λs)

    (fwhms, cx, cy) = decode_WaveLampLikelihood_input(nλ, self.order, V)
        
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
            GaussianModel2.(radiusmatrix, fwhms[i])
        end
    end
    
    # amplitude of gaussian spots
    Zygote.@ignore begin
        self.last_amp .= updateAmplitudeWaveLamp(list_spots, self.data, self.weights)
        any(isnan, self.last_amp) && error("NaN amplitude: \"$(self.last_amp)\"")
    end

    @inbounds begin
        list_spots_amped = map(i -> list_spots[i] .* self.last_amp[i], 1:nλ)

        sumspots = reduce(.+, list_spots_amped)

        cost = sum(self.weights .* (self.data .- sumspots).^2)
    end
    
    cost
end

 """
        updateAmplitude(nλ,m,d,W)

    return the `nλ` amplitudes `a` according the the model `m`, the data and the precision `W`
    such that
    `a = argmin_a || a*m - D||^2_W`
    where
    * `nλ` : is the number of spots in the model
    * `m`:  is the model composed of `nλ` images of spots
    * `d`:  is the data
    * `W`: is the precision (inverse variance) of the data
"""
function updateAmplitudeWaveLamp(
    spots::Vector{Matrix{Float64}}, data::AbstractArray{Float64}, weights::AbstractArray{Float64})
    
    N = length(spots)
    A = @MMatrix zeros(Float64,N,N)
    b = @MVector zeros(Float64,N)
    mw = Array{Float64}(undef, size(data)..., N);
    @inbounds for index=1:N
        mw[:,:,index] .=  spots[index] .* weights ;
        b[index] = sum(mw[:,:,index].* data );
        A[index,index] = sum(mw[:,:,index].* spots[index]);
        for i=1:index-1
            A[i,index] = A[index,i] = sum(mw[:,:,index].* spots[i])
        end
    end
    return inv(A)*b
end


