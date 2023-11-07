struct WaveLampModel{ Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    order ::Int
    nλ ::Int
    poweredλs ::Matrix{Float64}
    box ::BoundingBox{Int}
    data ::Md
    weights ::Mw
    
    function WaveLampModel{Md,Mw}(order, nλ, poweredλs, box, data, weights) where {Md,Mw}
    
        nλ > order || throw(ArgumentError("`order` must be strictly less than `nλ`"))
        size(poweredλs) == (order, nλ) || throw(DimensionMismatch(
            "size of `poweredλs` must be `(order, nλ)`"))
        size(box) == size(data) || throw(DimensionMismatch(
            "size of `box` must match the size of `data`"))
        size(data) == size(weights) || throw(DimensionMismatch(
            "size of `data` must match the size of `weights`"))
        
        new{Md,Mw}(order, nλ, poweredλs, box, data, weights)
    end
end

function WaveLampModel(
    order::Int, nλ::Int, poweredλs::Matrix{Float64}, box::BoundingBox{Int}, data::Md, weights::Mw
) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    
    WaveLampModel{Md,Mw}(order, nλ, poweredλs, box, data, weights)
end

function (self::WaveLampModel)(
    fwhm::Vector{Float64}, cx::Vector{Float64}, cy::Vector{Float64}
) ::Tuple{Matrix{Float64},Vector{Float64}}
        
    (xs, ys) = axes(self.box)
    
    # a matrix for each laser
    function build_matrix(i) # we avoid using "do" syntax for Zygote
        # spot centers
        @show cx
        @show self.poweredλs[:,i]
        spot_ctr_x = #=1014.6521027190289 =#compute_polynome_aux(cx, self.poweredλs[:,i])
        spot_ctr_y = #=1084.2108857705884 =#compute_polynome_aux(cy, self.poweredλs[:,i])
        
        # gaussian spot
        radiusmatrix = (xs .- spot_ctr_x).^2 .+ ((ys .- spot_ctr_y).^2)'
        GaussianModel2.(radiusmatrix, fwhm[i])
    end
    list_spots = map(build_matrix, 1:self.nλ)
    
    # amplitude of gaussian spots
    spots_amps = fill(1e0, self.nλ)
    Zygote.@ignore begin
        eachspot = cat(list_spots...; dims=3)
        spots_amps = updateAmplitude(self.nλ, eachspot, self.data, self.weights)
    end

    list_spots_amped =
        if any(isnan, spots_amps)
            @warn "NaN amplitude"
            list_spots
        else
            map(i -> list_spots[i] .* spots_amps[i], 1:self.nλ)
        end

    model = reduce(.+, list_spots_amped)

    (model, spots_amps)
end

#function wll_input_encode(
#    order::Int, nλ::Int, fwhm::Vector{Float64},
#    cx0::Float64, restcx::Vector{Float64}, cy0::Float64, restcy::Vector{Float64}
#) ::Matrix{Float64}
#    
#    M = fill(0e0, (nλ, 3))
#    
#    M[1:nλ, 1] .= fwhm
#    
#    M[1,         2]  = cx0
#    M[2:order+1, 2] .= restcx
#    
#    M[1,         3]  = cy0
#    M[2:order+1, 3] .= restcy
#    
#    @show M
#    
#    M
#end
#
#function wll_input_decode(order::Int, nλ::Int, M::Matrix{Float64}) ::NTuple{3,Vector{Float64}}
#
#    fwhm = M[1:nλ, 1]
#    cx   = M[1:order+1, 2]
#    cy   = M[1:order+1, 3]
#    
#    (fwhm, cx, cy)
#end

function wll_input_encode(
    fwhm::Vector{Float64},
    cx0::Float64, restcx::Vector{Float64},
    cy0::Float64, restcy::Vector{Float64}
) ::Vector{Float64}
    
    [ fwhm; cx0; restcx; cy0; restcy ]
end

function wll_input_decode(order::Int, nλ::Int, V::Vector{Float64}) ::NTuple{3,Vector{Float64}}
    
    fwhm  = V[1:nλ]
    
    start = nλ + 1
    cx    = V[start : start + order]
    
    start2 = start + order + 1
    cy    = V[start2 : start2 + order]
    
    (fwhm, cx, cy)
end

function WaveLampLikelihood(self::WaveLampModel)(V::Vector{Float64}) ::Tuple{Float64,Vector{Float64}}
    (fwhm, cx, cy) = wll_input_decode(self.model.order, self.nλ, V)
    WaveLampLikelihood(self)(fwhm, cx, cy)
end

function WaveLampLikelihood(self::WaveLampModel)(
    fwhm::Vector{Float64}, cx::Vector{Float64}, cy::Vector{Float64}
) ::Tuple{Float64,Vector{Float64}}
    
    model = self(fwhm, cx, cy)

    
    cost = sum(self.weights .* (self.data .- model).^2)
    
    cost
end


