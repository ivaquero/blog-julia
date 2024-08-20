### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> chapter = 3
#> section = 5
#> order = 5
#> title = "Time Series"
#> tags = ["prob", "track_prob", "Pluto", "PlutoUI"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ fa222d31-8b6a-4df0-966b-10512c21b44c
using PlutoUI

# ╔═╡ d1c423c2-9663-4f10-90bc-a5d5c3b91eb8
using Turing, CairoMakie

# ╔═╡ 87f0a3bb-c6bd-4ece-b015-33cbcf5bd925
begin
    using Random
    Random.seed!(2)
end

# ╔═╡ de69211f-4b0c-4a32-af04-7e3378991006
md"# Bayesian Time Series"

# ╔═╡ 00bd0968-a0de-4aad-8aa3-bf5d955682c6
begin
    using Pkg
    if isdir("../../pluto-deployment-environment")
        Pkg.activate("../../pluto-deployment-environment")
        Pkg.instantiate()
    else
        println(pwd())
    end
end

# ╔═╡ 52a8f7b1-ae66-4a34-97ff-f3b0e42afd75
TableOfContents()

# ╔═╡ d8cf8725-525f-4de5-aa8e-f3bc7f36edbe
md"## Data"

# ╔═╡ 2bd6aa5b-3156-4dcb-9f32-90983645db8f
md"### Generate Data"

# ╔═╡ 78b663ea-371e-4909-a703-68fe71dbcd3b
true_ϕ_1 = -0.2

# ╔═╡ c567a213-cd16-4bba-9318-a407b3e9675a
true_ϕ_2 = 0.5

# ╔═╡ 81482ac9-a0f8-47c4-b042-7831f7ab3067
true_σ = 0.1

# ╔═╡ eab6dc87-d78a-42ba-8361-c762511cb07c
time = 50

# ╔═╡ 9806629f-a436-4e08-978b-ffc061aba6b9
X = Vector{Float64}(undef, time + 2)

# ╔═╡ 265cc300-ce80-4720-a3ec-54eb820f7874
X[1] = rand(Normal(0, true_σ))

# ╔═╡ 0a325e3d-e793-4e14-82a3-879cb0d9b151
X[2] = rand(Normal(0, true_σ))

# ╔═╡ 8d6cf4e6-5cc1-47bc-9818-ca12de1d6536
for t in 3:(time+2)
    noise = rand(Normal(0, true_σ))
    X[t] = true_ϕ_1 * X[t-1] +
           true_ϕ_2 * X[t-2] +
           noise
end

# ╔═╡ cff71227-8d53-4cbe-9fd8-04a9dbd064d9
X_data = X[3:end]

# ╔═╡ 82c5e028-6a18-4a08-8b72-23d3cdf7e19a
md"### Visualize Data"

# ╔═╡ c2edcd2f-a1a1-4517-af71-d7e0c6d511e8
begin
    fig = Figure()
    ax = Axis(fig[1, 1],
        title="Bayesian Autoregressive AR(2) Model",
        xlabel="t",
        ylabel="X_t")

    lines!(ax, X_data)

    # xlims!(0, 60)
    # ylims!(-0.6, 0.6)

    fig
end

# ╔═╡ 05f1bf71-e3a4-4790-b12c-0687a041fb5b
md"## Model"

# ╔═╡ 10cb6613-f02f-41b4-b0b8-fae77f224b63
md"### Define Model"

# ╔═╡ 6050c8aa-7642-4d00-bc2a-f95507c90ce8
@model function ts_model(time, X)
    # prior
    ϕ_1 ~ Normal(0, 1)
    ϕ_2 ~ Normal(0, 1)
    σ ~ Exponential(1)
    # likelihood
    X[1] ~ Normal(0, σ)
    X[2] ~ Normal(0, σ)
    for t in 3:(time+2)
        μ = ϕ_1 * X[t-1] + ϕ_2 * X[t-2]
        X[t] ~ Normal(μ, σ)
    end
end

# ╔═╡ d183b8b7-6a1c-4a06-9d63-af412c90bac7
md"### Infer Posterior Probability"

# ╔═╡ 5212eade-572e-417e-be69-72c16032dd3b
ts_m = ts_model(time, X)

# ╔═╡ 155326dd-8eda-4f53-a1c3-c4fa4dc34ca0
sampler = NUTS()

# ╔═╡ 389871b9-5c36-43e1-8941-f1fab46f3a7f
samples = 1_000

# ╔═╡ 45ca1b61-c986-46f5-9b13-7b87c60e3be0
chain_p = sample(ts_m, sampler, samples)

# ╔═╡ 03916377-9e04-42da-aa7d-752bcb236958
md"### Visualize Results"

# ╔═╡ b3c606d6-1dab-408d-a2aa-ef31740cd6e4
begin
    figp = Figure()

    titles = ["ϕ_1", "ϕ_2", "σ"]

    for (ind, title) in enumerate(titles)
        ax = Axis(figp[ind, 1],
            title="$title",
            xlabel="Iteration",
            ylabel="Sample Values")

        lines!(ax,
            1:samples,
            chain_p[title][:, 1],
            color=(:lightblue, 0.8))

        ax2 = Axis(figp[ind, 2],
            title="$title",
            xlabel="Iteration",
            ylabel="Sample Values")

        density!(ax2,
            chain_p[title][:, 1],
            color=(:lightblue, 0.8))
    end

    figp
end

# ╔═╡ 6abefc79-0142-498e-bbb5-d686ca13084b
md"## Make Predictions"

# ╔═╡ 43d2c00b-ef85-4516-a77e-8c713889f0c7
time_fcst = 10

# ╔═╡ 2fb53ce4-6731-4b72-bf32-e81ec3fdb949
X_fcst = Matrix{Float64}(undef, time_fcst + 2, samples)

# ╔═╡ d48ec7bc-9d35-414e-8ec2-c3c813ebd4c2
X_fcst[1, :] .= X_data[time-1]

# ╔═╡ 5adc594f-cd6c-4437-ae24-a64a80ce94e6
X_fcst[2, :] .= X_data[time]

# ╔═╡ 0c0457d7-f917-4d1b-86c0-8f137d7b3b9b
for col in 1:samples
    ϕ_1_fcst = rand(chain_p[:, 1, 1])
    ϕ_2_fcst = rand(chain_p[:, 2, 1])
    error_fcst = rand(chain_p[:, 3, 1])
    noise_fcst = rand(Normal(0, error_fcst))
    for row in 3:(time_fcst+2)
        X_fcst[row, col] = ϕ_1_fcst * X_fcst[row-1, col] +
                           ϕ_2_fcst * X_fcst[row-2, col] +
                           noise_fcst
    end
end

# ╔═╡ 910f50c6-e3a1-415c-a396-08b2a90810d0
md"### Visualize Predictions"

# ╔═╡ c5f00a0c-8636-4fd3-85e6-35e4d63cf20d
ts_fcst = time:(time+time_fcst)

# ╔═╡ 8e9e7723-22a8-47d6-88ac-17e5adb1fe9e
X_fcst_mean = [mean(X_fcst[i, 1:samples]) for i in 2:(time_fcst+2)]

# ╔═╡ 85edbc6a-ebce-4135-9e92-980aed9a8f22
begin
    for i in 1:samples
        lines!(ax,
            ts_fcst,
            X_fcst[2:end, i],
            linewidth=1,
            color=(:yellow, 0.2))
    end

    lines!(ax,
        ts_fcst,
        X_fcst_mean,
        linewidth=2,
        color=:red,
        linestyle=:dot)

    fig
end


# ╔═╡ Cell order:
# ╟─de69211f-4b0c-4a32-af04-7e3378991006
# ╠═00bd0968-a0de-4aad-8aa3-bf5d955682c6
# ╠═fa222d31-8b6a-4df0-966b-10512c21b44c
# ╠═d1c423c2-9663-4f10-90bc-a5d5c3b91eb8
# ╠═52a8f7b1-ae66-4a34-97ff-f3b0e42afd75
# ╟─d8cf8725-525f-4de5-aa8e-f3bc7f36edbe
# ╟─2bd6aa5b-3156-4dcb-9f32-90983645db8f
# ╠═87f0a3bb-c6bd-4ece-b015-33cbcf5bd925
# ╠═78b663ea-371e-4909-a703-68fe71dbcd3b
# ╠═c567a213-cd16-4bba-9318-a407b3e9675a
# ╠═81482ac9-a0f8-47c4-b042-7831f7ab3067
# ╠═eab6dc87-d78a-42ba-8361-c762511cb07c
# ╠═9806629f-a436-4e08-978b-ffc061aba6b9
# ╠═265cc300-ce80-4720-a3ec-54eb820f7874
# ╠═0a325e3d-e793-4e14-82a3-879cb0d9b151
# ╠═8d6cf4e6-5cc1-47bc-9818-ca12de1d6536
# ╠═cff71227-8d53-4cbe-9fd8-04a9dbd064d9
# ╟─82c5e028-6a18-4a08-8b72-23d3cdf7e19a
# ╠═c2edcd2f-a1a1-4517-af71-d7e0c6d511e8
# ╟─05f1bf71-e3a4-4790-b12c-0687a041fb5b
# ╟─10cb6613-f02f-41b4-b0b8-fae77f224b63
# ╠═6050c8aa-7642-4d00-bc2a-f95507c90ce8
# ╟─d183b8b7-6a1c-4a06-9d63-af412c90bac7
# ╠═5212eade-572e-417e-be69-72c16032dd3b
# ╠═155326dd-8eda-4f53-a1c3-c4fa4dc34ca0
# ╠═389871b9-5c36-43e1-8941-f1fab46f3a7f
# ╠═45ca1b61-c986-46f5-9b13-7b87c60e3be0
# ╟─03916377-9e04-42da-aa7d-752bcb236958
# ╠═b3c606d6-1dab-408d-a2aa-ef31740cd6e4
# ╟─6abefc79-0142-498e-bbb5-d686ca13084b
# ╠═43d2c00b-ef85-4516-a77e-8c713889f0c7
# ╠═2fb53ce4-6731-4b72-bf32-e81ec3fdb949
# ╠═d48ec7bc-9d35-414e-8ec2-c3c813ebd4c2
# ╠═5adc594f-cd6c-4437-ae24-a64a80ce94e6
# ╠═0c0457d7-f917-4d1b-86c0-8f137d7b3b9b
# ╟─910f50c6-e3a1-415c-a396-08b2a90810d0
# ╠═c5f00a0c-8636-4fd3-85e6-35e4d63cf20d
# ╠═8e9e7723-22a8-47d6-88ac-17e5adb1fe9e
# ╠═85edbc6a-ebce-4135-9e92-980aed9a8f22
