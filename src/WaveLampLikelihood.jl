struct WaveLampLikelihood{ Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    order ::Int
    nλ ::Int
    poweredλs ::Matrix{Float64}
    box ::BoundingBox{Int}
    data ::Md
    weights ::Mw

    last_fwhm ::Vector{Float64}
    last_cx   ::Vector{Float64}
    last_cy   ::Vector{Float64}
    last_amp  ::Vector{Float64}
    
    function WaveLampLikelihood{Md,Mw}(
        order, nλ, poweredλs, box, data, weights, last_fwhm, last_cx, last_cy, last_amp
    ) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    
        nλ > order || throw(ArgumentError("`order` must be strictly less than `nλ`"))
        size(poweredλs) == (order, nλ) || throw(DimensionMismatch(
            "size of `poweredλs` must be `(order, nλ)`"))
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
        
        new{Md,Mw}(order, nλ, poweredλs, box, data, weights, last_fwhm, last_cx, last_cy, last_amp)
    end
end

function WaveLampLikelihood(
    poweredλs ::Matrix{Float64}, box::BoundingBox{Int}, data::Md, weights::Mw
) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    
    (order, nλ) = size(poweredλs)
    last_fwhm = Vector{Float64}(undef, nλ)
    last_cx   = Vector{Float64}(undef, order + 1)
    last_cy   = Vector{Float64}(undef, order + 1)
    last_amp  = Vector{Float64}(undef, nλ)
    
    WaveLampLikelihood{Md,Mw}(
        order, nλ, poweredλs, box, data, weights, last_fwhm, last_cx, last_cy, last_amp)
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
    
    spots_ctrs = get_lasers_centers(cx, cy, self.poweredλs)
    
    # a matrix for each laser
    list_spots = map(1:self.nλ) do i
    
        # spot centers
        spot_ctr_x = compute_polynome(cx, self.poweredλs[:,i])
        spot_ctr_y = compute_polynome(cy, self.poweredλs[:,i])
        
        # gaussian spot
        radiusmatrix = (xs .- spot_ctr_x).^2 .+ ((ys .- spot_ctr_y).^2)'
        GaussianModel2.(radiusmatrix, fwhm[i])
    end
    
    # amplitude of gaussian spots
    Zygote.@ignore begin
        eachspot = cat(list_spots...; dims=3)
        spot_amps = updateAmplitude(self.nλ, eachspot, self.data, self.weights)
        any(isnan, spot_amps) && error("NaN amplitude: \"$spot_amps\"")
        self.last_amp .= spot_amps
    end

    list_spots_amped = map(1:self.nλ) do i ; list_spots[i] .* self.last_amp[i] end

    sumspots = reduce(.+, list_spots_amped)

    sum(self.weights .* (self.data .- sumspots).^2)
end



