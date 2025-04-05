### A Pluto.jl notebook ###
# v0.20.5

#> [frontmatter]
#> chapter = 1
#> section = 1
#> order = 1
#> title = "Data Manipulation"
#> tags = ["analysis", "track_data", "Pluto", "PlutoUI"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ bbb110f4-eb0c-45e8-8456-85a951ecbfee
begin
    using Pkg
    if isdir("../../pluto-deployment-environment")
        Pkg.activate("../../pluto-deployment-environment")
        Pkg.instantiate()
    else
        println(pwd())
    end
end

# ╔═╡ a6999fdb-d20a-4a13-b7e9-eee118eb98fa
using Tidier

# ╔═╡ 9d9047ab-150b-4e6f-a836-bb99625150bb
using RDatasets

# ╔═╡ 328eddd6-88c1-48ac-acdf-ef0dd90f0e81
using PlutoUI

# ╔═╡ 06e1b0b6-885a-4edf-bc55-3ffffeaee764
md"# Data Manipulation"

# ╔═╡ f54f185c-9d0a-4166-8003-c4eb19bd9c63
TableOfContents()

# ╔═╡ 1cc9ec2c-3cd7-43cc-8087-64932687f40e
md"## Glimpse"

# ╔═╡ 902a1f74-8d10-43c7-a1eb-d090438f0a27
movies = dataset("ggplot2", "movies")

# ╔═╡ 36ae5c5c-074f-4ede-9e65-795dfa0adb5e
@glimpse(movies)

# ╔═╡ Cell order:
# ╟─06e1b0b6-885a-4edf-bc55-3ffffeaee764
# ╠═bbb110f4-eb0c-45e8-8456-85a951ecbfee
# ╠═a6999fdb-d20a-4a13-b7e9-eee118eb98fa
# ╠═9d9047ab-150b-4e6f-a836-bb99625150bb
# ╠═328eddd6-88c1-48ac-acdf-ef0dd90f0e81
# ╠═f54f185c-9d0a-4166-8003-c4eb19bd9c63
# ╟─1cc9ec2c-3cd7-43cc-8087-64932687f40e
# ╠═902a1f74-8d10-43c7-a1eb-d090438f0a27
# ╠═36ae5c5c-074f-4ede-9e65-795dfa0adb5e
