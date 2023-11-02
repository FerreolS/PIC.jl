struct WaveLampLikelihood{ Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    order ::Int
    nλ ::Int
    λpowers ::Matrix{Float64}
    box ::BoundingBox{Int}
    data ::Md
    weights ::Mw

    last_fwhm ::Vector{Float64}
    last_cx   ::Vector{Float64}
    last_cy   ::Vector{Float64}
    last_amp  ::Vector{Float64}
    
    function WaveLampLikelihood{Md,Mw}(
        order, nλ, λpowers, box, data, weights, last_fwhm, last_cx, last_cy, last_amp
    ) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    
        nλ > order || throw(ArgumentError("`order` must be strictly less than `nλ`"))
        size(λpowers) == (order, nλ) || throw(DimensionMismatch(
            "size of `λpowers` must be `(order, nλ)`"))
        size(box) == size(data) || throw(DimensionMismatch(
            "size of `box` must match the size of `data`"))
        size(data) == size(weights) || throw(DimensionMismatch(
            "size of `data` must match the size of `weights`"))
        length(last_fwhm) == nλ || throw(DimensionMismatch(
            "length of `last_fwhm` must be `nλ`"))
        length(last_cx) == order + 1 || throw(DimensionMismatch(
            "length of `last_cx` must be `order + 1`"))
        length(last_cy) == order + 1 || throw(DimensionMismatch(
            "length of `last_cy` must be `order + 1`"))
        length(last_amp) == nλ || throw(DimensionMismatch(
            "length of `last_amp` must be `nλ`"))
        
        new{Md,Mw}(order, nλ, λpowers, box, data, weights, last_fwhm, last_cx, last_cy, last_amp)
    end
end

function WaveLampLikelihood(
    λpowers ::Matrix{Float64}, box::BoundingBox{Int}, data::Md, weights::Mw
) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    
    (order, nλ) = size(λpowers)
    last_fwhm = Vector{Float64}(undef, nλ)
    last_cx   = Vector{Float64}(undef, order + 1)
    last_cy   = Vector{Float64}(undef, order + 1)
    last_amp  = Vector{Float64}(undef, nλ)
    
    WaveLampLikelihood{Md,Mw}(
        order, nλ, λpowers, box, data, weights, last_fwhm, last_cx, last_cy, last_amp)
end

function (self::WaveLampLikelihood)(V::Vector{Float64}) ::Float64
    
    startfwhm = 1
    endfwhm   = self.nλ
    startcx   = endfwhm + 1
    endcx     = startcx + self.order
    startcy   = endcx   + 1
    endcy     = startcy + self.order
    
    fwhm = V[startfwhm : endfwhm]
    cx   = V[startcx   : endcx]
    cy   = V[startcy   : endcy]
    
    self(fwhm, cx, cy)
end


function (self::WaveLampLikelihood)(
    fwhm::Vector{Float64}, cx::Vector{Float64}, cy::Vector{Float64}
) ::Float64
 
    Zygote.@ignore begin
        self.last_fwhm .= fwhm
        self.last_cx   .= cx
        self.last_cy   .= cy
        # we reset `last_amp`, because if the fit fails, we will be left with the
        # `last_cx` from the current run but the `last_amp` from the previous run.
        self.last_amp  .= NaN
    end
        
    (xs, ys) = axes(self.box)
    zygeachspot = Zygote.Buffer(Array{Float64}(undef, size(self.data)..., self.nλ))
    for i in 1:self.nλ
    
        # center of the laser spot
        spot_ctr_x = cx[1] + sum(cx[2:end] .* self.λpowers[:,i])
        spot_ctr_y = cy[1] + sum(cy[2:end] .* self.λpowers[:,i])
        
        # gaussian spot
        radiusmatrix = (xs .- spot_ctr_x).^2 .+ ((ys .- spot_ctr_y).^2)'
        zygeachspot[:,:,i] = GaussianModel2.(radiusmatrix, fwhm[i])
    end
    eachspot = copy(zygeachspot)
    
    # amplitude of gaussian spots
    Zygote.@ignore begin
        spot_amps = updateAmplitude(self.nλ, eachspot, self.data, self.weights)
        any(isnan, spot_amps) && error("NaN amplitude: \"$spot_amps\"")
        self.last_amp .= spot_amps
    end
    
    sumspots = zeros(Float64, size(self.data))
    for i in 1:self.nλ
        sumspots += eachspot[:,:,i] .* self.last_amp[i]
    end
    
    sum(self.weights .* (self.data .- sumspots).^2)
end
