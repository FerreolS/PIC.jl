### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ f6fe4721-1592-4c7f-8c5f-5e7b26c5e919
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
		Pkg.PackageSpec(name="ImageIO", version="0.5"),
		Pkg.PackageSpec(name="ImageShow", version="0.2"),
		Pkg.PackageSpec(name="FileIO", version="1.6"),
		Pkg.PackageSpec(name="PNGFiles", version="0.3.6"),
	])

	using ImageShow, FileIO

end

# ╔═╡ 90852aee-c4f9-4d44-a98b-b7f2188e398d
f = open("/home/user/stage/PIC.jl/GR.png")

# ╔═╡ efd7e18c-67e2-4829-94ee-ef9a9ba0730b
image = load("GR.png")

# ╔═╡ Cell order:
# ╠═f6fe4721-1592-4c7f-8c5f-5e7b26c5e919
# ╠═90852aee-c4f9-4d44-a98b-b7f2188e398d
# ╠═efd7e18c-67e2-4829-94ee-ef9a9ba0730b
