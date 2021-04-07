### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 32f9655e-97c4-11eb-0a2b-6f0d53de22bc
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
		Pkg.PackageSpec(name="ImageIO", version="0.5"),
		Pkg.PackageSpec(name="ImageShow", version="0.2"),
		Pkg.PackageSpec(name="FileIO", version="1.6"),
		Pkg.PackageSpec(name="PNGFiles", version="0.3.6"),
		Pkg.PackageSpec(name="Colors", version="0.12"),
		Pkg.PackageSpec(name="ColorVectorSpace", version="0.8"),
		Pkg.PackageSpec(name="PlutoUI", version="0.7"), 
		Pkg.PackageSpec(name="HypertextLiteral", version="0.5")
	])

	using Colors, ColorVectorSpace, ImageShow, FileIO
	using PlutoUI
	using HypertextLiteral
end

# ╔═╡ c3b7ae52-d494-4d03-b0c1-417efc05af84
begin
	url = "https://user-images.githubusercontent.com/6933510/107239146-dcc3fd00-6a28-11eb-8c7b-41aaf6618935.png" 
	philip_filename = download(url)
	philip = load(philip_filename)
end

# ╔═╡ Cell order:
# ╠═32f9655e-97c4-11eb-0a2b-6f0d53de22bc
# ╠═c3b7ae52-d494-4d03-b0c1-417efc05af84
