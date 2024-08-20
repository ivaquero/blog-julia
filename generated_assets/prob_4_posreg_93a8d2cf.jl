### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> chapter = 3
#> section = 4
#> order = 4
#> title = "Poisson Regression"
#> tags = ["prob", "track_prob", "Pluto", "PlutoUI"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ f636175b-3889-4b73-bc2f-355cd8d0b4d1
begin
    using Pkg
    Pkg.activate("../../pluto-deployment-environment")
    Pkg.instantiate()
end

# ╔═╡ 525bc515-010c-4bbd-bd29-33e6b02bad50
using Turing, CairoMakie

# ╔═╡ 11088891-2852-4bf6-a187-98d87c56a37c
using TidierFiles, TidierData

# ╔═╡ bf89f22d-4083-4f91-86ee-c076b73b070a
using GLM

# ╔═╡ 1a299c5c-5862-41b3-ae4e-5f9ed69a3be5
using PlutoUI

# ╔═╡ ef844473-6575-4c5d-a8d2-8af433c2ff20
using PlutoUI: Slider

# ╔═╡ 6487462e-e13c-4e11-a4d0-ae772c69ce9c
md"# Bayesian Possion Regression"

# ╔═╡ ee5b43ea-66ac-43cb-b0f3-c81451241e73
TableOfContents()

# ╔═╡ e86a5399-150f-418b-9175-5b1a457f8470
md"## Pollen"

# ╔═╡ d92adcb8-4669-4da4-9942-cf68fd62ed59
pollen = read_csv("../data/pollen_meds.csv")

# ╔═╡ 6cadafed-5256-4093-a49a-21659efca463
md"### Split Data"

# ╔═╡ 9d0b036c-7281-456f-ae67-69d205641b0f
noPollen_noMeds = @chain pollen begin
    @filter(pollen == 0 && meds == 0)
end

# ╔═╡ a0ef0260-775c-4c8b-b1d2-e5c29ed3c408
noPollen_Meds = @chain pollen begin
    @filter(pollen == 0 && meds == 1)
end

# ╔═╡ 23a8037a-da6d-4666-aa84-3f088276c6b4
Pollen_noMeds = @chain pollen begin
    @filter(pollen == 1 && meds == 0)
end

# ╔═╡ 91b78fe2-9ae1-4ea8-8e66-a96e6865554a
Pollen_Meds = @chain pollen begin
    @filter(pollen == 1 && meds == 1)
end

# ╔═╡ 7f438b59-19c2-4a1f-a0f5-869d05801404
md"### Visualize Data"

# ╔═╡ a8c309e9-c37a-4cce-8756-fd123c803599


# ╔═╡ 657dc724-775e-4598-9751-b8ec39435a29
titles_p = [
    "No Pollen, No Meds",
    "No Pollen, Meds",
    "Pollen, No Meds",
    "Pollen, Meds"
]

# ╔═╡ 4c6378a0-94e0-4b03-84b9-3fb6355ec3fe
data_p = [
    noPollen_noMeds.count,
    noPollen_Meds.count,
    Pollen_noMeds.count,
    Pollen_Meds.count
]

# ╔═╡ b8a74d5f-3d15-4090-97a9-274b7b7010be
data_p[1][2]

# ╔═╡ ed28511f-19b2-482f-89bb-bb7e03f93cc5
begin
    figps = Figure()

    axs = [Axis(figps[row, col]) for row in 1:2, col in 1:2]

    for (ind, ax) in enumerate(axs)
        hist!(ax, data_p[ind])
        ax.title = titles_p[ind]
    end

    figps
end

# ╔═╡ b7cfc059-fe7a-40f8-8f90-e65a1d0d14a4
md"### Features & Labels"

# ╔═╡ 336bb024-af3e-451b-af7c-85b43cdfd759
features = Matrix(
    @chain pollen begin
        @select(pollen, meds)
    end
)

# ╔═╡ 1cc40211-bb07-4945-8932-660aa1254ddf
labels = pollen.count

# ╔═╡ 7dbeb7e9-28c8-4f35-9c1b-44c0789626ba
md"## Concepts"

# ╔═╡ 28db6864-4360-45f7-91d4-293cc2142d84
md"### Multiple Linear Regression"

# ╔═╡ 5a74531d-d86d-4037-b41f-05ae08cd7a75
md"$y_i = β_0 + β_1x_{i1} + β_2x_{i2} + … + β_px_{ip} + ϵ_i$"

# ╔═╡ f922acc8-919b-480f-8fcf-e1731db4d1f6
md"where, for $i = 1, ..., n$ observations"

# ╔═╡ 999a3a8e-b8e4-4790-ac5d-8f167a5f51c5
md"### Log Transformation"

# ╔═╡ 03536bd1-a2ce-467a-ae12-bbcb4cd89e03
md"Log transformation is a method of transforming datasets from nonlinear to linear and vice versa."

# ╔═╡ d0e941f1-da51-46f9-8d4a-39b4e408a649
md"Log transformations are often recommended for skewed data to help spead out the data."

# ╔═╡ ec810a4a-c2a0-4407-94bc-b5c4e4a12088
md"Examples include Power, Exponential, Square Root, Logarithmic and Hyperbolic."

# ╔═╡ 6f91afe1-8d26-4fbb-8be9-a3ae202e8cff
md"### Data for Plotting"

# ╔═╡ 7beb7f5b-2b4c-4520-b27e-eb12dea9352a
md"intercept $(β_0)$: $(@bind intercept Slider(-1.0:0.1:1.0, 0.0, true))"

# ╔═╡ 51debf88-cd90-421e-ad69-bc678e3369c5
md"coefficient 1 $(β_1)$: $(@bind coeff1 Slider(-0.25:0.01:0.25, 0.0, true))"

# ╔═╡ 3155971b-fbb0-4ccd-b606-72c43f6ecb3c
md"coefficient 2 $(β_2)$: $(@bind coeff2 Slider(-0.25:0.01:0.25, 0.0, true))"

# ╔═╡ a7f9029c-8996-4c62-a82d-5560fc22de22
line(x) = intercept + coeff1 * x + coeff2 * x

# ╔═╡ 1e3c2965-1f20-40a9-aec7-487cde65e57e
log_transform(x) = exp(line(x))

# ╔═╡ e4484b47-c19a-4af1-8197-7fa63c8001c9
xs = -10.0:0.1:10.0

# ╔═╡ 332546c4-be8c-45d2-8643-b52dc85b6ae9
ys1 = [line(x) for x in xs]

# ╔═╡ c7b9398d-e048-4b92-aa1a-617a2ae54190
ys2 = [log_transform(x) for x in xs]

# ╔═╡ 4aa59993-c4f0-4fdd-930f-5e0dccd1239c
begin
    figl = Figure()

    ax = Axis(figl[1, 1], title="Line")

    lines!(ax, xs, ys1, linewidth=3)
    vlines!(ax, [0], color=:black)

    axl = Axis(figl[1, 2], title="Log Transform")

    lines!(axl, xs, ys2, linewidth=3)
    vlines!(axl, [0], color=:black)

    # xlims!(-10, 10)
    # ylims!(-5, 5)

    figl
end

# ╔═╡ 4f1520c4-e249-449e-a82b-cfc40b66522f
md"## Define Model"

# ╔═╡ 0838bf77-3d27-4c9a-856d-3be38e9b2338
@model function pollen_model(X, y)
    # prior
    intercept ~ Normal(0, 1)
    pollen_x ~ Normal(0, 1)
    meds_x ~ Normal(0, 1)
    # likelihood
    for i in eachindex(y)
        line = intercept + pollen_x * X[i, 1] + meds_x * X[i, 2]
        log_transform = exp(line)
        y[i] ~ Poisson(log_transform)
    end
end

# ╔═╡ 82d96ca2-3831-44ed-9046-6736a32dc597
md"### Infer Posterior Probability"

# ╔═╡ aca0cad6-8c0b-4045-ae9c-e39f40da57a7
pollen_m = pollen_model(features, labels)

# ╔═╡ c320b899-18bf-4e27-bb01-40899aba6259
sampler = NUTS()

# ╔═╡ 419a2c4a-cda3-4f48-8a2a-9417eb70a561
samples = 1_000

# ╔═╡ 913f596f-3d7b-4810-b9d9-65cef48bc813
chain_p = sample(pollen_m, sampler, samples)

# ╔═╡ efb97d39-7166-4ef6-bafd-1597bb58ee57
begin
    figp = Figure()

    titles = ["intercept", "pollen_x", "meds_x"]

    for (ind, title) in enumerate(titles)
        ax = Axis(
            figp[ind, 1],
            title="$title",
            xlabel="Iteration",
            ylabel="Sample Values"
        )

        lines!(
            ax,
            1:samples,
            chain_p[title][:, 1],
            color=(:lightblue, 0.8),
        )

        ax2 = Axis(
            figp[ind, 2],
            title="$title",
            xlabel="Iteration",
            ylabel="Sample Values"
        )

        density!(
            ax2,
            chain_p[title][:, 1],
            color=(:lightblue, 0.8)
        )
    end

    figp
end

# ╔═╡ 3cb8d428-e3c9-4011-b0d1-ecb0e3c17797
md"### Make Predictions"

# ╔═╡ cdb7dd3e-3b26-4704-8d18-4afb641e6d97
newX = [
    0 0 # noPollen_noMeds
    0 1 # noPollen_Meds
    1 0 # Pollen_noMeds
    1 1 # Pollen_Meds
]

# ╔═╡ d2df4e22-7953-4541-9292-71c7a1d9e3e2
missing_counts = fill(missing, 4)

# ╔═╡ cca37a3e-55c5-4e37-9034-4a94abac8bae
predict_model = pollen_model(newX, missing_counts)

# ╔═╡ bfec1a7e-e4c4-4671-a6bf-3684acbb6abd
predictions = predict(predict_model, chain_p)

# ╔═╡ 56ff6e03-4ead-474b-9efa-b7a07b35ea64
begin
    figpp = Figure(size=(800, 800))

    for ind in 1:4
        ax = Axis(
            figpp[ind, 1],
            title="y[$ind]",
            xlabel="Iteration",
            ylabel="Sample Values"
        )

        lines!(
            ax,
            1:samples,
            predictions["y[$ind]"][:, 1],
            color=(:lightblue, 0.8)
        )

        ax2 = Axis(
            figpp[ind, 2],
            title="y[$ind]",
            xlabel="Iteration",
            ylabel="Sample Values"
        )

        hist!(
            ax2,
            predictions["y[$ind]"][:, 1],
            color=(:lightblue, 0.8)
        )
    end

    figpp
end

# ╔═╡ 3b588b4d-1bd4-42b2-b9f8-c37fd72a72e9
begin
    figx = Figure()

    for ind in 1:4
        ax = Axis(figx[ind, 1])
        hist!(ax, data_p[ind])
        ax.title = titles_p[ind]

        ax2 = Axis(figx[ind, 2])
        hist!(ax2, predictions["y[$ind]"][:, 1])
        ax2.title = titles_p[ind]
    end

    figx
end

# ╔═╡ Cell order:
# ╟─6487462e-e13c-4e11-a4d0-ae772c69ce9c
# ╠═f636175b-3889-4b73-bc2f-355cd8d0b4d1
# ╠═525bc515-010c-4bbd-bd29-33e6b02bad50
# ╠═11088891-2852-4bf6-a187-98d87c56a37c
# ╠═bf89f22d-4083-4f91-86ee-c076b73b070a
# ╠═1a299c5c-5862-41b3-ae4e-5f9ed69a3be5
# ╠═ef844473-6575-4c5d-a8d2-8af433c2ff20
# ╠═ee5b43ea-66ac-43cb-b0f3-c81451241e73
# ╟─e86a5399-150f-418b-9175-5b1a457f8470
# ╠═d92adcb8-4669-4da4-9942-cf68fd62ed59
# ╟─6cadafed-5256-4093-a49a-21659efca463
# ╠═9d0b036c-7281-456f-ae67-69d205641b0f
# ╠═a0ef0260-775c-4c8b-b1d2-e5c29ed3c408
# ╠═23a8037a-da6d-4666-aa84-3f088276c6b4
# ╠═91b78fe2-9ae1-4ea8-8e66-a96e6865554a
# ╟─7f438b59-19c2-4a1f-a0f5-869d05801404
# ╠═a8c309e9-c37a-4cce-8756-fd123c803599
# ╠═657dc724-775e-4598-9751-b8ec39435a29
# ╠═4c6378a0-94e0-4b03-84b9-3fb6355ec3fe
# ╠═b8a74d5f-3d15-4090-97a9-274b7b7010be
# ╠═ed28511f-19b2-482f-89bb-bb7e03f93cc5
# ╟─b7cfc059-fe7a-40f8-8f90-e65a1d0d14a4
# ╠═336bb024-af3e-451b-af7c-85b43cdfd759
# ╠═1cc40211-bb07-4945-8932-660aa1254ddf
# ╟─7dbeb7e9-28c8-4f35-9c1b-44c0789626ba
# ╟─28db6864-4360-45f7-91d4-293cc2142d84
# ╟─5a74531d-d86d-4037-b41f-05ae08cd7a75
# ╟─f922acc8-919b-480f-8fcf-e1731db4d1f6
# ╟─999a3a8e-b8e4-4790-ac5d-8f167a5f51c5
# ╟─03536bd1-a2ce-467a-ae12-bbcb4cd89e03
# ╟─d0e941f1-da51-46f9-8d4a-39b4e408a649
# ╟─ec810a4a-c2a0-4407-94bc-b5c4e4a12088
# ╟─6f91afe1-8d26-4fbb-8be9-a3ae202e8cff
# ╟─7beb7f5b-2b4c-4520-b27e-eb12dea9352a
# ╟─51debf88-cd90-421e-ad69-bc678e3369c5
# ╟─3155971b-fbb0-4ccd-b606-72c43f6ecb3c
# ╠═a7f9029c-8996-4c62-a82d-5560fc22de22
# ╠═1e3c2965-1f20-40a9-aec7-487cde65e57e
# ╠═e4484b47-c19a-4af1-8197-7fa63c8001c9
# ╠═332546c4-be8c-45d2-8643-b52dc85b6ae9
# ╠═c7b9398d-e048-4b92-aa1a-617a2ae54190
# ╠═4aa59993-c4f0-4fdd-930f-5e0dccd1239c
# ╟─4f1520c4-e249-449e-a82b-cfc40b66522f
# ╠═0838bf77-3d27-4c9a-856d-3be38e9b2338
# ╟─82d96ca2-3831-44ed-9046-6736a32dc597
# ╠═aca0cad6-8c0b-4045-ae9c-e39f40da57a7
# ╠═c320b899-18bf-4e27-bb01-40899aba6259
# ╠═419a2c4a-cda3-4f48-8a2a-9417eb70a561
# ╠═913f596f-3d7b-4810-b9d9-65cef48bc813
# ╠═efb97d39-7166-4ef6-bafd-1597bb58ee57
# ╟─3cb8d428-e3c9-4011-b0d1-ecb0e3c17797
# ╠═cdb7dd3e-3b26-4704-8d18-4afb641e6d97
# ╠═d2df4e22-7953-4541-9292-71c7a1d9e3e2
# ╠═cca37a3e-55c5-4e37-9034-4a94abac8bae
# ╠═bfec1a7e-e4c4-4671-a6bf-3684acbb6abd
# ╠═56ff6e03-4ead-474b-9efa-b7a07b35ea64
# ╠═3b588b4d-1bd4-42b2-b9f8-c37fd72a72e9
