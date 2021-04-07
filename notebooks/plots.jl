### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 8f78c854-96e3-11eb-3cab-c387fdcf1a27
begin
	import Pkg
	Pkg.activate(mktempdir())
end

# ╔═╡ 07afb510-96b9-11eb-27b2-7b5cb6520572
begin
	Pkg.add("Plots")
	using Plots
end

# ╔═╡ f8aaeaca-96e6-11eb-226c-dfdb147a99dc
begin
	
	Pkg.add("StatsPlots")
	using StatsPlots # Required for the DataFrame user recipe
end

# ╔═╡ b9063658-96bf-11eb-3e76-6d6c57aac488
begin
	using Distributions
	plot(Normal(3, 5), lw = 3)
end

# ╔═╡ 59a02e6c-96bf-11eb-12b7-575a704131e6
begin
	x4 = 1:10; 
	y4 = rand(10, 4)
	plot(x4, y4, layout = (4, 1))
end

# ╔═╡ 820366d0-96bf-11eb-2a83-5d4accef7a7e
begin
	x = 1:10; 
	y = rand(10, 2); 
	p1 = plot(x, y, title = "Lines") # Make a line plot
	p2 = scatter(x, y) # Make a scatter plot
	p3 = plot(x, y, xlabel = "This one is labelled", lw = 1)
	p4 = histogram(x, y) # Four histograms each with 10 points? Why not!
	plot(p1, p2, p3, p4, layout = (2, 2), legend = false)
end

# ╔═╡ 433455be-96b4-11eb-3ec3-5d095cafac85
begin
	
	x1 = 1:10; 
	y1 = rand(10, 2);  # These are the plotting data
	#plot(x, y)
	plot(x, y, title = "Two Lines", label = ["Line 1" "Line 2"], lw = 2)
	xlabel!("My x label")
	ylabel!("My y label")
	
end


# ╔═╡ 9dd1e744-96e5-11eb-08fe-419f712777f2
begin
	# This plots into the web browser via Plotly
	x2 = 1:10; 
	y2 = rand(10, 2); 
	plotly() # Set the backend to Plotly
	plot(x, y, title = "This is Plotted using Plotly")
		
end

# ╔═╡ 415686f4-96be-11eb-0ea1-3b3b4231422f
begin
	x3 = 1:10; 
	y3 = rand(10, 2);
	gr() # Set the backend to GR
		# This plots using GR
		plot(x, y, title = "This is Plotted using GR")
end

# ╔═╡ 7bad1cf2-96e9-11eb-0b8d-21afa54bf89b
plot([sin], 0:0.1:4π) 

# ╔═╡ e9f7da18-96e8-11eb-0da7-959529c1c618
plot([sin,cos], 0:0.1:2π)                     # 2 series, sin.(x) and cos.(x)

# ╔═╡ Cell order:
# ╠═8f78c854-96e3-11eb-3cab-c387fdcf1a27
# ╠═07afb510-96b9-11eb-27b2-7b5cb6520572
# ╠═433455be-96b4-11eb-3ec3-5d095cafac85
# ╠═9dd1e744-96e5-11eb-08fe-419f712777f2
# ╠═415686f4-96be-11eb-0ea1-3b3b4231422f
# ╠═59a02e6c-96bf-11eb-12b7-575a704131e6
# ╠═820366d0-96bf-11eb-2a83-5d4accef7a7e
# ╠═f8aaeaca-96e6-11eb-226c-dfdb147a99dc
# ╠═b9063658-96bf-11eb-3e76-6d6c57aac488
# ╠═7bad1cf2-96e9-11eb-0b8d-21afa54bf89b
# ╠═e9f7da18-96e8-11eb-0da7-959529c1c618
