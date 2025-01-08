### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> chapter = 1
#> section = 2
#> order = 1
#> title = "Plot Correlation"
#> tags = ["analysis", "track_data", "Pluto", "PlutoUI"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ bbb110f4-eb0c-45e8-8456-85a951ecbfee
begin
     using Pkg
     Pkg.activate("../../pluto-deployment-environment")
     Pkg.instantiate()
 end

# ╔═╡ a6999fdb-d20a-4a13-b7e9-eee118eb98fa
using TidierFiles, TidierData, TidierPlots

# ╔═╡ 07f923d3-3995-4915-a7f3-c236e6a6d0ea
md"# Plot Correlation"

# ╔═╡ 8eecaeac-b0d9-4587-9a9f-61db52c31ff8
md"## Regression Plot"

# ╔═╡ 902a1f74-8d10-43c7-a1eb-d090438f0a27
ads = read_csv("../data/advertising.csv")

# ╔═╡ 52d98253-dd93-4af8-b6d3-d9d2d5b9b1e2
p1 = ggplot(ads, @aes(x = TV, y = sales)) +
     geom_point(color=:red) +
     geom_smooth(method=:lm, linewidth=3, alpha=0.5)

# ╔═╡ efb16650-04b3-47ea-8f1f-3a2436751adc
p2 = ggplot(ads, @aes(x = radio, y = sales)) +
     geom_point(color=:red) +
     geom_smooth(method=:lm, linewidth=3, alpha=0.5)

# ╔═╡ 3b604dd8-d398-4bab-91d1-72bd70a965dc
p3 = ggplot(ads, @aes(x = newspaper, y = sales)) +
     geom_point(color=:red) +
     geom_smooth(method=:lm, linewidth=3, alpha=0.5)

# ╔═╡ 636af518-6bdb-4c6c-ac03-d0971dfaa479
p13 = (p1 + p2 + p3 + plot_layout(ncol=3, nrow=1, widths=[5, 5, 5], heights=[1]))

# ╔═╡ Cell order:
# ╟─07f923d3-3995-4915-a7f3-c236e6a6d0ea
# ╠═bbb110f4-eb0c-45e8-8456-85a951ecbfee
# ╠═a6999fdb-d20a-4a13-b7e9-eee118eb98fa
# ╟─8eecaeac-b0d9-4587-9a9f-61db52c31ff8
# ╠═902a1f74-8d10-43c7-a1eb-d090438f0a27
# ╠═52d98253-dd93-4af8-b6d3-d9d2d5b9b1e2
# ╠═efb16650-04b3-47ea-8f1f-3a2436751adc
# ╠═3b604dd8-d398-4bab-91d1-72bd70a965dc
# ╠═636af518-6bdb-4c6c-ac03-d0971dfaa479
