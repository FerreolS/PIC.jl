# Polynomes

function compute_polynome_aux(coefs::Vector{Float64}, powered_xs::Vector{Float64}) ::Float64
    coefs[1] + sum(coefs[2:end] .* powered_xs)
end

function compute_polynome(coefs::Vector{Float64}, x::Float64) ::Float64
    powered_xs = x .^ (1 : length(coefs)-1)
    compute_polynome_aux(coefs, powered_xs)
end

function polynome(coefs::Vector{Float64}) ::Function
    (x::Float64) -> compute_polynome(coefs, x)
end

# Polynomes with reference for x

function apply_reference(ref::Float64, x::Float64) ::Float64
    (x - ref) / ref
end

function compute_polynome_with_reference(ref::Float64,coefs::Vector{Float64}, x::Float64) ::Float64
    refed_x = apply_reference(ref, x)
    compute_polynome(coefs, refed_x)
end

function polynome_with_reference(ref::Float64, coefs::Vector{Float64}) ::Function
    (x::Float64) -> compute_polynome_with_reference(ref, coefs, x)
end
