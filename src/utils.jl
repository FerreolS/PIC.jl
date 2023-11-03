function polynome(coefs::Vector{Float64}) ::Function
    Base.Fix1(compute_polynome, coefs)
end

function compute_polynome(coefs::Vector{Float64}, x::Float64) ::Float64
    powered_xs = x .^ (1:length(coefs)-1)
    polynome(coefs, powered_xs)
end

function compute_polynome(coefs::Vector{Float64}, powered_xs::Vector{Float64}) ::Float64
    coefs[1] + sum(coefs[2:end] .* powered_xs)
end
