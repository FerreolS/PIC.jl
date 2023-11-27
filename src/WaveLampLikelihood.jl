"""
    WaveLampLikelihood{Md,Mw}(λ0, λs, box, data, weights, last_background, last_amp)

Structure to store the constant data used in the fitting of the parameters for WaveLamp.

Find information about fitted parameters at [`(self::WaveLampLikelihood)`](@ref).

There is two parameters that are not fitted but computed from the others: background and amplitude.
Since they are not parameters, they are not returned to the user. That is why they are stored in
this structure. Only the values from the last run are stored. Find additional information about
these two values at [`compute_background_amplitudes_wavelamp`](@ref).

# Fields
- `λ0 ::Float64`: the reference wavelength. Usually the mean of `λs`.
- `λs ::Vector{Float64}`: the wavelength for each laser spot. In the same unity as `λ0`.
- `box ::BoundingBox{Int}`: integral pixel coordinates for the box of the lens on the detector.
- `data ::Md`: the observed data for the lens, size must be the same as `box`.
- `weights ::Mw`: the weights for `data`, size must be the same as `data`.
- `last_background ::Ref{Float64}`: stores the computed background from the last run.
- `last_amp ::Vector{Float64}`: stores the computed amplitude for each laser, from the last run.
"""
struct WaveLampLikelihood{ Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    λ0 ::Float64
    λs ::Vector{Float64}
    box ::BoundingBox{Int}
    data ::Md
    weights ::Mw
    
    last_background ::Ref{Float64}
    last_amp ::Vector{Float64}
    # Note that it was simpler to make Zygote work when these two values are stored here,
    # it seems to have less trouble to find them.
    
    function WaveLampLikelihood{Md,Mw}(
        λ0, λs, box, data, weights, last_background, last_amp) where {Md,Mw}
    
        size(box) == size(data) || throw(DimensionMismatch(
            "size of `box` must match the size of `data`"))
        size(data) == size(weights) || throw(DimensionMismatch(
            "size of `data` must match the size of `weights`"))
        length(last_amp) == length(λs) || throw(DimensionMismatch(
            "length of `last_amp` must be `length(λs)`"))
        
        new{Md,Mw}(λ0, λs, box, data, weights, last_background, last_amp)
    end
end

"""
    WaveLampLikelihood(λ0, λs, box, data, weights)
Version of the constructor with initial NaN values for background and amplitude, and with
automatic types parameters.
"""
function WaveLampLikelihood(
    λ0::Float64, λs::Vector{Float64}, box::BoundingBox{Int}, data::Md, weights::Mw
) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64} }
    
    last_background = Ref(NaN64)
    last_amp = fill(NaN64, size(λs))
    WaveLampLikelihood{Md,Mw}(λ0, λs, box, data, weights, last_background, last_amp)
end

"""
    (self::WaveLampLikelihood)(M::Matrix{Float64}) -> Float64

Compute a model from parameters `M`, and compare it to the `data`, and return the cost (WLS).

# Parameters (stored in a matrix of three columns)
- `fwhms`: FWHM for each laser. length must equal `length(self.λs)`.
- `cx`: coefficient for the polynome of the x-position of the centers of the gaussian spots of the
  lasers. The polynome order must be `length(self.λs) - 1`, therefore the number of coefficients
  must be `length(self.λs)`. The first coefficient is in absolute coordinates on the detector.
- `cy`: same as `cx` for y-position.

Background and amplitude are not fitted, find more information at
[`compute_background_amplitudes_wavelamp`](@ref).
"""
function (self::WaveLampLikelihood)(M::Matrix{Float64}) ::Float64

    fwhms = M[:,1]
    cx    = M[:,2]
    cy    = M[:,3]
        
    nλ = length(self.λs)
    size(M,1) == nλ || throw(DimensionMismatch("wrong size of parameters matrix"))

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
        (background, amps) = compute_background_amplitudes_wavelamp(
            list_spots, self.data, self.weights)
        self.last_background.x = background
        self.last_amp .= amps
        any(isnan, self.last_amp) && error("NaN amplitude: \"$(self.last_amp)\"")
    end

    @inbounds begin
        list_spots_amped = map(1:nλ) do i
                               list_spots[i] .* self.last_amp[i]
                           end

        model = reduce(.+, list_spots_amped) .+ self.last_background.x

        cost = sum(self.weights .* (self.data .- model).^2)
    end
    
    cost
end

"""
    compute_background_amplitudes_wavelamp(models, data, weights) -> Tuple{Float64,Vector{Float64}}

From a vector of models, a data matrix and a weights matrix, compute the "best" background and
amplitudes.

These two values are not fitted because they can be computed from the other parameters. For a set
of parameters, we can find a background and an amplitude such that WLS derivative is zero for the
amplitude and the background.

Proof:
- Recall that the WLS is the sum of `(data - model)² * weights` (with pixel-wise operations).
- Let us write it in a matrix-vector way.
- Let `d` be the data vector `Nx1` (flatten the matrix if needed).
- Let `C` be the number of gaussians.
- Let `gi` be the ith gaussian model vector Nx1, with amplitude `1` and background `0`.
- Let `G` be a matrix `Nx(C+1)`, with the first column filled with `1`, and second column filled
  with `g1`, third column filled with `g2`, etc.
- Let `p` be a vector `(C+1)x1 = [β, α1, α2...]`, with `β` the variable for background and `αi` the
  variable for the amplitude of the ith gaussian.
- Let `W` be a matrix `NxN`, a diagonal matrix containing the weights.
- Then the WLS can be written:  `(d - Gp)ᵀ W (d - Gp)`
- Unfolding this expression gives:  `dᵀWd - 2pᵀGᵀWd + pᵀGᵀWGp`
- Deriving this WLS expression by `p`, gives:  `-2GᵀWd + 2GᵀWGp`
- When the WLS derived by `p` is zero, it implies that:   `GᵀWGp = GᵀWd`
- We isolate `p`:   `p = (GᵀWG)⁻¹ GᵀWd   = [β, α1, α2...]`
"""
function compute_background_amplitudes_wavelamp(
    models  ::Vector{Matrix{Float64}},
    data    ::AbstractMatrix{Float64},
    weights ::AbstractMatrix{Float64}
) ::Tuple{Float64,Vector{Float64}}
    
    M = [ ones(Float64, size(data)), models... ]
    length_M = length(M)

    A = @MMatrix zeros(Float64, length_M, length_M)
    b = @MVector zeros(Float64, length_M)
    
    Mi_weights = Matrix{Float64}(undef, size(data))
    
    @inbounds for i in 1:length_M
        Mi_weights .= M[i] .* weights
        for j in 1:i
            A[i,j] = A[j,i] = sum(Mi_weights .* M[j])
        end
        b[i] = sum(Mi_weights .* data)
    end
    
    res = inv(A) * b
    (background, amps) = res[1], res[2:end]
end
