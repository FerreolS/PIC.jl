"""
    SpecPosLikelihood{Md,Mw,Mm}(λ0, order, box, data, weights, λMap, last_background, last_amps)

Structure to store the constant data used in the fitting of the parameters for SpecPos.

The model is a frame with a 1D gaussian on each line.

Find information about fitted parameters at [`compute_SpecPosLikelihood`](@ref).

There is two parameters that are not fitted but computed from the others: background and amplitudes.
Since they are not parameters, they are not returned to the user. That is why they are stored in
this structure. Only the values from the last run are stored. Find additional information about
these two values at [`compute_background_amplitudes_specpos`](@ref).

# Fields
- `λ0 ::Float64`: the reference wavelength. Usually the mean of `λs`.
- `order ::Float64`: the order of the polynomes for center x of the gaussians and fwhm of
  the gaussians.
- `box ::BoundingBox{Int}`: integral pixel coordinates for the box of the lens on the detector.
- `data ::Md`: the observed data for the lens, size must be the same as `box`.
- `weights ::Mw`: the weights for `data`, size must be the same as `data`.
- `λMap ::Mm`: the wavelength attributed for each pixel of the box.
- `last_background ::Ref{Float64}`: stores the computed background from the last run.
- `last_amps ::Vector{Float64}`: stores the computed amplitude for each line of the box,
  from the last run.

# Type parameters
`Md`, `Mw`, and `Mm` are just here to permit to give "views" instead of plain arrays.
"""
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

"""
    SpecPosLikelihood(λ0, order, box, data, weights, λMap)

Alternative constructor with NaNs as initial `last_background` and `last_amps`.
"""
function SpecPosLikelihood(
    λ0::Float64, order::Int, box::BoundingBox{Int}, data::Md, weights::Mw, λMap::Mm
) where { Md<:AbstractMatrix{Float64}, Mw<:AbstractMatrix{Float64}, Mm<:AbstractMatrix{Float64} }

    last_background = Ref(NaN64)
    last_amps = fill(NaN64, size(box,2))
    SpecPosLikelihood{Md,Mw,Mm}(λ0, order, box, data, weights, λMap, last_background, last_amps)
end

"""
    (self::SpecPosLikelihood)(M::Matrix{Float64}) -> Float64

Alias for [`compute_SpecPosLikelihood`](@ref)
"""
function (self::SpecPosLikelihood)(M::Matrix{Float64}) ::Float64
    compute_SpecPosLikelihood(self, M)
end

"""
    compute_SpecPosLikelihood(self::SpecPosLikelihood, M::Matrix{Float64}) -> Float64

Compute a model from parameters `M`, and compare it to the `data`, and return the cost (WLS).

# Parameters (stored in a matrix of two columns)
- `cx`: coefficients for the polynome of the x position of the centers of the gaussians. The
  polynome order is `self.order`, therefore the number of coefficients
  must be `self.order + 1`. The first coefficient is in absolute coordinates on the detector.
- `cλ: coefficients for the polynome of the fwhm of the gaussians. The polynome order is
  `self.order`, therefore the number of coefficients must be `self.order + 1`. The first
  coefficient is in floating-point pixel length.

Background and amplitudes are not fitted, find more information at
[`compute_background_amplitudes_specpos`](@ref).
"""
function compute_SpecPosLikelihood(self::SpecPosLikelihood, M::Matrix{Float64}) ::Float64

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
        (background, amps) = compute_background_amplitudes_specpos(model, self.data, self.weights)
        self.last_background.x = background
        self.last_amps .= amps
    end
    
    @inbounds begin
        model_amped = (model .* self.last_amps') .+ self.last_background.x
        cost        = sum(self.weights .* (self.data .- model_amped).^2)
    end
    
    cost
end

"""
    compute_background_amplitudes_wavelamp(models, data, weights) -> Tuple{Float64,Vector{Float64}}

From a model, a data matrix and a weights matrix, compute the "best" background and
amplitudes.

These two parameters are not fitted because they can be computed from the other parameters. For a
set of parameters, we can find a background and an amplitudes such that WLS derivative is zero for
the amplitudes and the background.

The proof is similar to the one of [`compute_background_amplitudes_wavelamp`](@ref), with the
difference that the gaussians here are 1D, in the sense that they only have data on their line, and
are supposed to print zero data on other lines. This allows us to put zeros in the `A` matrix.
"""
function compute_background_amplitudes_specpos(
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
    
    # substitute zero lines by dummy values
    za =  iszero.(a) .| (b .≤ 0e0)
    if any(za)
        a[za] .= 1e0
        b[za] .= 0e0
        c[za] .= 0e0
    end
    
    # note that the gaussians are disjoint: each gaussian fills a line on the box, and has
    # absolutely no data on other lines.

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

