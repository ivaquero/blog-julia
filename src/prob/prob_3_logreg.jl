### A Pluto.jl notebook ###
# v0.20.5

#> [frontmatter]
#> chapter = 3
#> section = 3
#> order = 3
#> title = "Logistic Regression"
#> tags = ["prob", "track_prob", "Pluto", "PlutoUI"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try
            Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value
        catch
            b -> missing
        end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d8134e3a-ca12-419d-81eb-d64a7e0e5969
begin
    using Pkg
    if isdir("../../pluto-deployment-environment")
        Pkg.activate("../../pluto-deployment-environment")
        Pkg.instantiate()
    else
        println(pwd())
    end
end

# ╔═╡ 9e1f8648-f895-4171-972c-1aa89d167637
using Turing, CairoMakie

# ╔═╡ 281e569f-faeb-4ce6-9294-683a53352778
using Tidier

# ╔═╡ 9ea5838a-59c8-483e-9a1c-b0f3ce9702c5
using GLM

# ╔═╡ cb48c841-1290-43ee-bd88-1ccaa69e18f4
using PlutoUI

# ╔═╡ 65488b10-bc49-4fb9-969f-d05efc44173c
using PlutoUI: Slider

# ╔═╡ 85fdd2e9-3aa2-4550-888b-3bba321987db
md"# Bayesian Logistic Regression"

# ╔═╡ 8a468b60-0628-47fa-a9f1-c7c1d3aa35f6
TableOfContents()

# ╔═╡ 18702173-4999-4205-b978-f6dce495685d
md"## Equation for a Line"

# ╔═╡ 08f4a549-1def-4d36-a972-ccc94e2df653
md"$y = mx + b$"

# ╔═╡ 95973762-9f95-43ea-8ffb-908fdde59916
md"where $m$ is the slope and $b$ is the y-intercept."

# ╔═╡ 258266f2-b801-4d32-84c5-8125c652b3f4
md"### Sigmoid Function"

# ╔═╡ fd30dc5b-d9fd-45c3-af38-461e914d54f6
md"$S(x) = \frac{1}{1 + e^{-x}}$"

# ╔═╡ 2b15594b-1c55-4e76-b323-8aa5df68f6f9
md"slope: $(@bind slope Slider(-10.0:1.0:10.0, 1.0, true))"

# ╔═╡ 1d76d6a7-44ab-4597-9c31-4111b23e2701
md"intercept: $(@bind intercept Slider(-5.0:1.0:5.0, 0.0, true))"

# ╔═╡ 54012055-634f-4a4c-8bd7-e95d9041ff59
y(x) = intercept + slope * x

# ╔═╡ 3e27ee40-e4f0-4389-9eb8-141c14491803
S(x) = 1 / (1 + exp(-y(x)))

# ╔═╡ 2a4b185e-d7c5-4637-9486-a949e5eec4fd
xs = -10:10

# ╔═╡ e9dcec74-6f74-45d0-839a-59e1497817d2
begin
    figs = Figure()

    axl = Axis(figs[1, 1], title="Line")

    lines!(axl, xs, y.(xs), linewidth=2)
    hlines!(axl, [0], color=:black)
    vlines!(axl, [0], color=:black)

    xlims!(axl, -10, 10)
    ylims!(axl, -10, 10)

    axs = Axis(figs[1, 2], title="Sigmoid Curve")

    lines!(axs, xs, S.(xs), linewidth=2)
    hlines!(axs, [0.5], color=:black, linestyle=:dash)
    vlines!(axs, [0], color=:black)

    xlims!(axs, -10, 10)
    ylims!(axs, 0, 1)

    figs
end

# ╔═╡ 1a40908a-bb1b-474e-842e-c2eb7eda2beb
md"## Wolfspider"

# ╔═╡ 76c92bce-dd07-42d9-98ca-5343127f2623
wolfspider = read_csv("../data/wolfspider.csv")

# ╔═╡ 01a83359-e819-42d1-9159-a325a2a30a7a
spider_x = wolfspider.feature

# ╔═╡ c2228232-e615-49ca-b6bd-b0847df02599
spider_y = wolfspider.class

# ╔═╡ e6393d1f-6703-47d4-a3a7-c6561ce6f102
spider_matrix = [spider_x spider_y]

# ╔═╡ a07aa69e-0013-4aa3-a9a4-11869265bdf4
begin
    figw = Figure()

    axw = Axis(figw[1, 1],
        title="Wolf Spider Present (1) or Absent (0)",
        xlabel="Median Grain Size of Sand (mm)",
        ylabel="Probability of Presence")

    sp_data = scatter!(axw,
        spider_x,
        spider_y,
        color=:red,
        markersize=15)

    xlims!(0, 1.25)
    figw
end

# ╔═╡ e2c8b634-6b2a-4c79-a52e-5c9c43ddf93f
md"### Define Model"

# ╔═╡ d34a1268-a3e7-43a1-ac29-3dd657e99467
@model function spider_model(grain_size, presence)
    # prior
    intercept ~ Uniform(-5, 0)
    slope ~ Uniform(0, 10)
    # likelihood
    line = intercept .+ slope .* grain_size
    p = 1 ./ (1 .+ exp.(-line))
    for i in eachindex(p)
        presence[i] ~ Bernoulli(p[i])
    end
end

# ╔═╡ b64bb28a-7b38-4f5c-a73f-01309a2f99d4
md"### Infer Posterior Probability"

# ╔═╡ 3cec32f6-4d85-474f-b893-9bd2f07fefd5
spider_m = spider_model(spider_x, spider_y)

# ╔═╡ 326be2ff-a0a9-476b-887a-ac5741e8b403
sampler_w = NUTS()

# ╔═╡ e6024dce-f9b0-4a59-a3d0-6490ede22468
samples = 1_000

# ╔═╡ 3bd969ee-c8c0-4db2-a2f7-f69ff132a354
chain_w = sample(spider_m, sampler_w, samples)

# ╔═╡ a2e1e3b4-da77-437a-a897-b9b3a2dc6372
md"### Visualize Results"

# ╔═╡ 69804c2d-b6f1-4e2a-af1b-a3c8c93ff627
begin
    figws = Figure()

    titles = ["intercept", "slope"]

    for ind in 1:2
        ax = Axis(figws[ind, 1],
            title=titles[ind],
            xlabel="Iteration",
            ylabel="Sample Values")

        lines!(ax, 1:samples, chain_w[titles[ind]][:, 1])

        ax2 = Axis(figws[ind, 2],
            title=titles[ind],
            xlabel="Iteration",
            ylabel="Sample Values")

        density!(ax2,
            chain_w[titles[ind]][:, 1],
            color=(:lightblue, 0.5),
            strokecolor=(:blue, 0.8),
            strokewidth=1.5)
    end

    figws
end

# ╔═╡ fbc232cf-07fa-4f77-8d7b-365d6e1791f2
xws = 0.0:0.01:1.2

# ╔═╡ 0816d71c-19a5-4dc4-a905-233638b7d480
for i in 1:samples
    intercept_w = chain_w[i, 1, 1]
    slope_w = chain_w[i, 2, 1]
    line(x) = intercept_w + slope_w * x
    pw(x) = 1 / (1 + exp(-line(x)))

    lines!(axw,
        xws,
        pw,
        color=(:blue, 0.03))
end

# ╔═╡ 63a102a1-3ba0-42c8-b60a-bf2d35bad83f
figw

# ╔═╡ 6f5feef8-d047-4f8c-b513-17fbff7ca0a8
md"### Predictions"

# ╔═╡ e586f88e-dcee-464e-8a79-7be214b42b89
newXw = [0.25, 0.5, 0.75, 1.0]

# ╔═╡ 0c18dced-2a3a-4e02-862b-069f2f983c0f
p_spider = fill(missing, length(newXw))

# ╔═╡ b95d26f7-eb65-467a-9d42-da819a4cc209
predictions_w = predict(spider_model(newXw, p_spider), chain_w)

# ╔═╡ 05051423-b852-4385-8d87-226f2f9506a1
begin
    figp = Figure()

    for ind in 1:4
        ax = Axis(figp[ind, 1],
            title="presence[$ind]",
            xlabel="Iteration",
            ylabel="Sample Values")

        lines!(ax,
            1:samples,
            predictions_w["presence[$ind]"][:, 1],
            color=(:lightblue, 0.8))

        ax2 = Axis(figp[ind, 2],
            title="presence[$ind]",
            xlabel="Iteration",
            ylabel="Sample Values")

        hist!(ax2,
            predictions_w["presence[$ind]"][:, 1],
            color=(:lightblue, 0.8))
    end

    figp
end

# ╔═╡ Cell order:
# ╟─85fdd2e9-3aa2-4550-888b-3bba321987db
# ╠═d8134e3a-ca12-419d-81eb-d64a7e0e5969
# ╠═9e1f8648-f895-4171-972c-1aa89d167637
# ╠═281e569f-faeb-4ce6-9294-683a53352778
# ╠═9ea5838a-59c8-483e-9a1c-b0f3ce9702c5
# ╠═cb48c841-1290-43ee-bd88-1ccaa69e18f4
# ╠═65488b10-bc49-4fb9-969f-d05efc44173c
# ╠═8a468b60-0628-47fa-a9f1-c7c1d3aa35f6
# ╟─18702173-4999-4205-b978-f6dce495685d
# ╟─08f4a549-1def-4d36-a972-ccc94e2df653
# ╟─95973762-9f95-43ea-8ffb-908fdde59916
# ╟─258266f2-b801-4d32-84c5-8125c652b3f4
# ╟─fd30dc5b-d9fd-45c3-af38-461e914d54f6
# ╟─2b15594b-1c55-4e76-b323-8aa5df68f6f9
# ╟─1d76d6a7-44ab-4597-9c31-4111b23e2701
# ╠═54012055-634f-4a4c-8bd7-e95d9041ff59
# ╠═3e27ee40-e4f0-4389-9eb8-141c14491803
# ╠═2a4b185e-d7c5-4637-9486-a949e5eec4fd
# ╠═e9dcec74-6f74-45d0-839a-59e1497817d2
# ╟─1a40908a-bb1b-474e-842e-c2eb7eda2beb
# ╠═76c92bce-dd07-42d9-98ca-5343127f2623
# ╠═01a83359-e819-42d1-9159-a325a2a30a7a
# ╠═c2228232-e615-49ca-b6bd-b0847df02599
# ╠═e6393d1f-6703-47d4-a3a7-c6561ce6f102
# ╠═a07aa69e-0013-4aa3-a9a4-11869265bdf4
# ╟─e2c8b634-6b2a-4c79-a52e-5c9c43ddf93f
# ╠═d34a1268-a3e7-43a1-ac29-3dd657e99467
# ╟─b64bb28a-7b38-4f5c-a73f-01309a2f99d4
# ╠═3cec32f6-4d85-474f-b893-9bd2f07fefd5
# ╠═326be2ff-a0a9-476b-887a-ac5741e8b403
# ╠═e6024dce-f9b0-4a59-a3d0-6490ede22468
# ╠═3bd969ee-c8c0-4db2-a2f7-f69ff132a354
# ╟─a2e1e3b4-da77-437a-a897-b9b3a2dc6372
# ╠═69804c2d-b6f1-4e2a-af1b-a3c8c93ff627
# ╠═fbc232cf-07fa-4f77-8d7b-365d6e1791f2
# ╠═0816d71c-19a5-4dc4-a905-233638b7d480
# ╠═63a102a1-3ba0-42c8-b60a-bf2d35bad83f
# ╟─6f5feef8-d047-4f8c-b513-17fbff7ca0a8
# ╠═e586f88e-dcee-464e-8a79-7be214b42b89
# ╠═0c18dced-2a3a-4e02-862b-069f2f983c0f
# ╠═b95d26f7-eb65-467a-9d42-da819a4cc209
# ╠═05051423-b852-4385-8d87-226f2f9506a1
