### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 30a915a4-96f2-11eb-2b39-97540117c119
import Pkg

# ╔═╡ 76b7a1a0-96f2-11eb-23ef-494525c84322
Pkg.add("Distributions")

# ╔═╡ 827b0414-96f2-11eb-224c-b1f0ff602585
Pkg.add("Gadfly")

# ╔═╡ 76b81732-96f2-11eb-07c3-8798f3eb727c
using Distributions

# ╔═╡ 5b2d5600-96f2-11eb-321f-45a062c1a0d2
using Gadfly

# ╔═╡ fee70066-96f2-11eb-188c-9b9705955c8a
x_norm = rand(Normal(100, 20), 1000)

# ╔═╡ 3d51f072-96f3-11eb-3529-d7934b808942
plot(x = x_norm, Geom.density)

# ╔═╡ 91d5cba8-96f3-11eb-0e83-83ecd353bd05
plot(x = x_norm, Geom.histogram(bincount = 15))

# ╔═╡ 9f80701a-96f4-11eb-28f4-c7f51260aadf
mean(x_norm)

# ╔═╡ b22825aa-96f4-11eb-2bea-296925cdf0ee
median(x_norm)

# ╔═╡ bb376b2e-96f4-11eb-0c7d-b988fab32d90
var(x_norm)		#Variance

# ╔═╡ c3b65e9a-96f4-11eb-1913-c5dc161d052d
std(x_norm)		#Standard deviation

# ╔═╡ d5357854-96f4-11eb-0fd4-ffc512291d55
quantile(x_norm, [0.5, 0.95])		#0.5 = median

# ╔═╡ 0a690888-96f5-11eb-1ee0-d1d8d4f87574
fit(Normal, (x_norm))

# ╔═╡ Cell order:
# ╠═30a915a4-96f2-11eb-2b39-97540117c119
# ╠═76b7a1a0-96f2-11eb-23ef-494525c84322
# ╠═76b81732-96f2-11eb-07c3-8798f3eb727c
# ╠═827b0414-96f2-11eb-224c-b1f0ff602585
# ╠═5b2d5600-96f2-11eb-321f-45a062c1a0d2
# ╠═fee70066-96f2-11eb-188c-9b9705955c8a
# ╠═3d51f072-96f3-11eb-3529-d7934b808942
# ╠═91d5cba8-96f3-11eb-0e83-83ecd353bd05
# ╠═9f80701a-96f4-11eb-28f4-c7f51260aadf
# ╠═b22825aa-96f4-11eb-2bea-296925cdf0ee
# ╠═bb376b2e-96f4-11eb-0c7d-b988fab32d90
# ╠═c3b65e9a-96f4-11eb-1913-c5dc161d052d
# ╠═d5357854-96f4-11eb-0fd4-ffc512291d55
# ╠═0a690888-96f5-11eb-1ee0-d1d8d4f87574
