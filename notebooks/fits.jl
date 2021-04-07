### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ d50f191a-1828-4db6-a41e-e8e5344ddc64
begin
	import Pkg
	Pkg.activate(mktempdir())

	using FITSIO
end

# ╔═╡ f6f7df42-97c4-11eb-2460-c55540f9416f
f = FITS("/home/user/stage/HR_4796-HD_95086/Calibration_wave_spec/IFS_res_spec.fits")

# ╔═╡ Cell order:
# ╠═d50f191a-1828-4db6-a41e-e8e5344ddc64
# ╠═f6f7df42-97c4-11eb-2460-c55540f9416f
