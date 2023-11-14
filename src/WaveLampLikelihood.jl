struct WaveLampLikelihood{ Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    order ::Int
    nλ ::Int
    poweredλs ::Matrix{Float64}
    box ::BoundingBox{Int}
    data ::Md
    weights ::Mw

    last_amp ::Vector{Float64}
    
    function WaveLampLikelihood{Md,Mw}(
        order, nλ, poweredλs, box, data, weights, last_amp) where {Md,Mw}
    
        nλ > order || throw(DimensionMismatch("`order` must be strictly less than `nλ`"))
        size(poweredλs) == (order, nλ) || throw(DimensionMismatch(
            "size of `poweredλs` must be `(order, nλ)`"))
        size(box) == size(data) || throw(DimensionMismatch(
            "size of `box` must match the size of `data`"))
        size(data) == size(weights) || throw(DimensionMismatch(
            "size of `data` must match the size of `weights`"))
        length(last_amp) == nλ || throw(DimensionMismatch(
            "length of `last_amp` must be `nλ`"))
        
        new{Md,Mw}(order, nλ, poweredλs, box, data, weights, last_amp)
    end
end

function WaveLampLikelihood(
    order::Int, nλ::Int, poweredλs::Matrix{Float64}, box::BoundingBox{Int}, data::Md, weights::Mw
) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    
    last_amp  = Vector{Float64}(undef, nλ)
    WaveLampLikelihood{Md,Mw}(order, nλ, poweredλs, box, data, weights, last_amp)
end

function decode_WaveLampLikelihood_input(
    nλ::Int, order::Int, V::Vector{Float64}) ::NTuple{3,Vector{Float64}}
    
    start_fwhm = 1
    end_fwhm   = nλ
    
    start_cx   = end_fwhm + 1
    end_cx     = start_cx + order
    
    start_cy   = end_cx   + 1
    end_cy     = start_cy + order
    
    fwhm = V[start_fwhm : end_fwhm]
    cx   = V[start_cx   : end_cx]
    cy   = V[start_cy   : end_cy]
    
    (fwhm, cx, cy)
end

function (self::WaveLampLikelihood)(V::Vector{Float64}) ::Float64

    (fwhm, cx, cy) = decode_WaveLampLikelihood_input(self.nλ, self.order, V)
        
    (xs, ys) = axes(self.box)
    
    # a matrix for each laser
    function build_matrix(i) # we avoid using "do" syntax for Zygote
        # spot centers
        spot_ctr_x = compute_polynome_aux(cx, self.poweredλs[:,i])
        spot_ctr_y = compute_polynome_aux(cy, self.poweredλs[:,i])
        
        # gaussian spot
        radiusmatrix = (xs .- spot_ctr_x).^2 .+ ((ys .- spot_ctr_y).^2)'
        GaussianModel2.(radiusmatrix, fwhm[i])
    end
    list_spots = map(build_matrix, 1:self.nλ)
    
    # amplitude of gaussian spots
    Zygote.@ignore begin
        self.last_amp .= updateAmplitudeWaveLamp(list_spots, self.data, self.weights)
        any(isnan, self.last_amp) && error("NaN amplitude: \"$(self.last_amp)\"")
    end

    list_spots_amped = map(i -> list_spots[i] .* self.last_amp[i], 1:self.nλ)

    sumspots = reduce(.+, list_spots_amped)

    sum(self.weights .* (self.data .- sumspots).^2)
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


