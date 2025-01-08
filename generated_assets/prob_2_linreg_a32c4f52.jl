### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> chapter = 3
#> section = 2
#> order = 2
#> title = "Linear Regression"
#> tags = ["prob", "track_prob", "Pluto", "PlutoUI"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ 92af5ce3-e6b7-47b4-afc0-7f68dabc66ad
begin
    using Pkg
    if isdir("../../pluto-deployment-environment")
        Pkg.activate("../../pluto-deployment-environment")
        Pkg.instantiate()
    else
        println(pwd())
    end
end

# ╔═╡ 2afc5904-5084-4d49-b2a7-9d52a58e2ba2
using Turing, CairoMakie

# ╔═╡ 6ab3eb1a-419f-4e85-84bc-f46d9a2478b4
using GLM

# ╔═╡ 3e7e62d5-38b9-4b5e-b150-38040342e167
using TidierFiles, TidierData

# ╔═╡ 5fd892e2-fd7d-41ed-bc8c-feec39dbd978
using PlutoUI

# ╔═╡ 83691d17-9db7-4f0c-ab8b-d83bcd83e0b1
using PlutoUI: Slider

# ╔═╡ 4f03a1c9-c588-43d6-871b-7dcbadecfa64
md"# Bayesian Linear Regression"

# ╔═╡ dd3e802b-d290-475b-a76f-109eda3a3339
TableOfContents()

# ╔═╡ 5a3d59c3-daf2-43ec-ab1b-edf97d035e2e
md"## Percent of Water"

# ╔═╡ 1a0591d5-7790-47cd-9854-eced324fbd21
md"### Traditional Approach"

# ╔═╡ e8729da3-efac-4857-8086-ebb1f3ffa9c9
md"inputs"

# ╔═╡ bbb6f71c-6a3f-453e-997c-08bbab68e5e4
n = 9

# ╔═╡ 38c74fa3-e378-472e-8a62-2dac5a41b4f9
p = 0.71

# ╔═╡ 45baa341-ac9b-42bc-8a9f-eecc5ac08323
md"model"

# ╔═╡ 0ec8193a-b66e-4ac6-96d9-7e9023d92644
f(n, p) = Int(round(n * p))

# ╔═╡ e983b173-3d6c-4f84-9b4d-f73e98b7cf60
md"output"

# ╔═╡ 7af3b6dd-8727-4e5a-8831-5cca7f08f31b
w = f(n, p)

# ╔═╡ c3edb027-9f82-48bb-8085-396504bef134
md"try to back into the value of p"

# ╔═╡ db5e419c-f362-4047-8748-d7c683888ac9
w / n

# ╔═╡ de9a56fe-0944-44e3-9407-e4b420b7fd13
p1 = 0.62

# ╔═╡ e9b5ba7e-8f57-4192-a14d-16caa8497640
w1 = f(n, p1)

# ╔═╡ b968fe7e-6f21-4859-8950-a78ad864011a
p2 = 0.72

# ╔═╡ a6041f4b-5af0-4591-907f-56fcb92c9463
w2 = f(n, p2)

# ╔═╡ c96b65d7-858a-48d0-805d-a1dcdae4f115
md"### Bayesian Approach"

# ╔═╡ 61f1b930-bd08-4ee0-941b-74d0bf37e692
md"observe data"

# ╔═╡ 724565dc-0044-44ea-b942-98e05b7dcb5e
tosses = 9

# ╔═╡ eca719f7-39a0-4368-b055-4489cbe70c4d
water = 6

# ╔═╡ dfc9d864-b38a-4075-90bf-de61be1384dc
md"define model"

# ╔═╡ 3e44c77d-aff2-4b8e-9dcd-51550cd2d574
@model function globe_toss(tosses, water)
    # prior
    percent_water ~ Beta(1, 1)
    # likelihood
    water ~ Binomial(tosses, percent_water)
end

# ╔═╡ 70882d33-ff15-4c43-b5ee-c9ed13bb9cda
md"infer posterior distribution"

# ╔═╡ 923971b2-ec4a-4152-b9f5-7e65143010e6
globe_m = globe_toss(tosses, water)

# ╔═╡ 8d5c22ef-6e64-4d8d-92ab-22c5d236ec10
sampler = NUTS()

# ╔═╡ ac9244c5-2f02-4942-b5e9-587e0e18ad80
samples = 1_000

# ╔═╡ 09730556-fe22-4c7f-863f-b3a871855d5a
chain = sample(globe_m, sampler, samples)

# ╔═╡ 6192f67d-e505-4db4-8e1f-c0f64be28f54
md"visualize results"

# ╔═╡ 2bf748bd-763a-4d58-8ea0-24444419ad4f
begin
    figc = Figure()

    axc = Axis(figc[1, 1],
        title="Trace Plot of $samples samples",
        xlabel="Sample",
        ylabel="Percent of Water")

    lines!(axc, chain[:percent_water][:, 1])
    xlims!(-20, 1020)
    ylims!(0.2, 1.0)

    figc
end

# ╔═╡ 1fed525e-3828-4052-a58b-2d4d9a642568
begin
    figd = Figure()

    axd = Axis(figd[1, 1],
        title="Posterior Distribution (approx)",
        xlabel="Percent of Water",
        ylabel="Density")

    density!(axd,
        chain[:percent_water][:, 1],
        color=(:lightblue, 0.5),
        strokecolor=(:darkblue, 0.9),
        strokewidth=3)
    xlims!(0, 1)
    ylims!(0, 3)

    figd
end

# ╔═╡ 56340950-73ab-487c-b7cd-baefd6955ab0
md"## Height vs. Weight"

# ╔═╡ 8294f363-2193-4f49-af23-0abd153594b0
howell1 = read_delim("../data/Howell1.csv"; delim=";")

# ╔═╡ 669f277e-37be-471b-8a85-a148cb7b67f0
howell1_adult = @chain howell1 begin
    @filter(age .>= 18)
end

# ╔═╡ 210249e4-cc9a-4c95-a9ee-48b40568f1f8
md"### Non-Bayesian Approach (GLM.jl)"

# ╔═╡ 4ee3977e-f0ca-4bf3-94ee-c2c330806661
md"define model"

# ╔═╡ 7f400d10-5251-415d-aaff-5483efa82ae6
ols = lm(@formula(weight ~ height), howell1_adult)

# ╔═╡ 6bbbd538-8b84-40d3-9136-87fa71726804
ols_intercept = coef(ols)[1]

# ╔═╡ 256e91b6-d157-4ae9-81a8-5657f393b83e
ols_slope = coef(ols)[2]

# ╔═╡ b5f35ff9-3531-4b10-93c5-58a912297656
md"make predictions"

# ╔═╡ 11964a1c-6284-4983-8f67-4d691c5d9b65
ols_newX = DataFrame(height=[140, 160, 175])

# ╔═╡ 00d1d6d3-c699-4b8e-8d20-6ce30cad4dd1
ols_predictions = predict(ols, ols_newX)

# ╔═╡ c3528e4e-b5de-4183-9dc7-6b91020290b4
md"visualize results"

# ╔═╡ d574ad83-9a99-4027-88de-40179f9495d4
begin
    fig1 = Figure()

    ax1 = Axis(fig1[1, 1],
        title="Adult Weight-Height Association (GLM)",
        xlabel="Height (cm)",
        ylabel="Weight (kg)")

    scatter!(ax1, howell1_adult.height, howell1_adult.weight)
    lines!(ax1, 1:samples, ols_intercept .+ ols_slope .* (1:samples))

    xlims!(100, 200)
    ylims!(0, 80)

    fig1
end

# ╔═╡ 81502815-d8d2-414d-bf51-52ee11820b4a
md"### Bayesian Approach (Turing.jl)"

# ╔═╡ 0a69b477-2b1b-49f1-9a17-def741beb71b
md"define model"

# ╔═╡ 6f1ecff4-7a28-434d-b8d5-923eef215b61
@model function my_model(weight, height)
    # prior
    intercept ~ Normal(0, 10)
    slope ~ Uniform(0, 1)
    error ~ Uniform(0, 10)
    # likelihood
    avg_weight = intercept .+ slope .* height
    weight ~ MvNormal(avg_weight, error)
end

# ╔═╡ 9db79764-30f7-4c0f-af9e-9596e1ef9405
model_r = my_model(howell1_adult.weight, howell1_adult.height)

# ╔═╡ 4705c138-2aff-463a-9f60-3c4d219d551f
sampler_r = NUTS()

# ╔═╡ 5f2c1f99-314b-4c18-bffc-aab774e024e9
chain_r = sample(model_r, sampler_r, samples)

# ╔═╡ a5fc300d-3816-4d53-8144-3cd0f5c0f42f
newX = [140, 160, 175]

# ╔═╡ 1c125abb-e404-402b-b175-c7c8d2115ac2
predictions = predict(my_model(missing, newX), chain_r)

# ╔═╡ d2ba032e-7387-4cb5-b239-4a2174f7b92d
begin
    fig2 = Figure()

    ax2 = Axis(fig2[1, 1],
        title="Adult Weight-Height Association (Turing)",
        xlabel="Height (cm)",
        ylabel="Weight (kg)")

    scatter!(ax2, howell1_adult.height, howell1_adult.weight, color=(:red, 0.6))

    for i in 1:samples
        intercept = chain_r[i, 1, 1]
        slope = chain_r[i, 2, 1]
        error = chain_r[i, 3, 1]
        lines!(ax2, 1:samples, intercept .+ slope .* (1:samples), color=(:orange, 0.5), linewidth=0.8)
    end

    xlims!(100, 200)
    ylims!(0, 80)

    fig2
end

# ╔═╡ 7c0c06e2-f164-467d-bbba-04b5a4f94db4
begin
    fig3 = Figure()

    titles = ["Intercept", "Slope", "Error"]

    for ind in 1:3
        ax = Axis(fig3[ind, 1],
            title=titles[ind],
            xlabel="Iteration",
            ylabel="Sample Values")

        lines!(ax,
            1:samples,
            predictions["weight[$ind]"][:, 1],
            color=(:blue, 0.8))

        ax2 = Axis(fig3[ind, 2],
            title=titles[ind],
            xlabel="Iteration",
            ylabel="Sample Values")

        density!(ax2,
            predictions["weight[$ind]"][:, 1],
            color=(:lightblue, 0.5),
            strokecolor=(:blue, 0.8),
            strokewidth=1.5)

        # xlims!(10, 90)
        # ylims!(0, 0.1)
    end

    fig3
end

# ╔═╡ Cell order:
# ╟─4f03a1c9-c588-43d6-871b-7dcbadecfa64
# ╠═92af5ce3-e6b7-47b4-afc0-7f68dabc66ad
# ╠═2afc5904-5084-4d49-b2a7-9d52a58e2ba2
# ╠═6ab3eb1a-419f-4e85-84bc-f46d9a2478b4
# ╠═3e7e62d5-38b9-4b5e-b150-38040342e167
# ╠═5fd892e2-fd7d-41ed-bc8c-feec39dbd978
# ╠═83691d17-9db7-4f0c-ab8b-d83bcd83e0b1
# ╠═dd3e802b-d290-475b-a76f-109eda3a3339
# ╟─5a3d59c3-daf2-43ec-ab1b-edf97d035e2e
# ╟─1a0591d5-7790-47cd-9854-eced324fbd21
# ╟─e8729da3-efac-4857-8086-ebb1f3ffa9c9
# ╠═bbb6f71c-6a3f-453e-997c-08bbab68e5e4
# ╠═38c74fa3-e378-472e-8a62-2dac5a41b4f9
# ╟─45baa341-ac9b-42bc-8a9f-eecc5ac08323
# ╠═0ec8193a-b66e-4ac6-96d9-7e9023d92644
# ╟─e983b173-3d6c-4f84-9b4d-f73e98b7cf60
# ╠═7af3b6dd-8727-4e5a-8831-5cca7f08f31b
# ╟─c3edb027-9f82-48bb-8085-396504bef134
# ╠═db5e419c-f362-4047-8748-d7c683888ac9
# ╠═de9a56fe-0944-44e3-9407-e4b420b7fd13
# ╠═e9b5ba7e-8f57-4192-a14d-16caa8497640
# ╠═b968fe7e-6f21-4859-8950-a78ad864011a
# ╠═a6041f4b-5af0-4591-907f-56fcb92c9463
# ╟─c96b65d7-858a-48d0-805d-a1dcdae4f115
# ╟─61f1b930-bd08-4ee0-941b-74d0bf37e692
# ╠═724565dc-0044-44ea-b942-98e05b7dcb5e
# ╠═eca719f7-39a0-4368-b055-4489cbe70c4d
# ╟─dfc9d864-b38a-4075-90bf-de61be1384dc
# ╠═3e44c77d-aff2-4b8e-9dcd-51550cd2d574
# ╟─70882d33-ff15-4c43-b5ee-c9ed13bb9cda
# ╠═923971b2-ec4a-4152-b9f5-7e65143010e6
# ╠═8d5c22ef-6e64-4d8d-92ab-22c5d236ec10
# ╠═ac9244c5-2f02-4942-b5e9-587e0e18ad80
# ╠═09730556-fe22-4c7f-863f-b3a871855d5a
# ╟─6192f67d-e505-4db4-8e1f-c0f64be28f54
# ╠═2bf748bd-763a-4d58-8ea0-24444419ad4f
# ╠═1fed525e-3828-4052-a58b-2d4d9a642568
# ╟─56340950-73ab-487c-b7cd-baefd6955ab0
# ╠═8294f363-2193-4f49-af23-0abd153594b0
# ╠═669f277e-37be-471b-8a85-a148cb7b67f0
# ╟─210249e4-cc9a-4c95-a9ee-48b40568f1f8
# ╟─4ee3977e-f0ca-4bf3-94ee-c2c330806661
# ╠═7f400d10-5251-415d-aaff-5483efa82ae6
# ╠═6bbbd538-8b84-40d3-9136-87fa71726804
# ╠═256e91b6-d157-4ae9-81a8-5657f393b83e
# ╟─b5f35ff9-3531-4b10-93c5-58a912297656
# ╠═11964a1c-6284-4983-8f67-4d691c5d9b65
# ╠═00d1d6d3-c699-4b8e-8d20-6ce30cad4dd1
# ╟─c3528e4e-b5de-4183-9dc7-6b91020290b4
# ╠═d574ad83-9a99-4027-88de-40179f9495d4
# ╟─81502815-d8d2-414d-bf51-52ee11820b4a
# ╟─0a69b477-2b1b-49f1-9a17-def741beb71b
# ╠═6f1ecff4-7a28-434d-b8d5-923eef215b61
# ╠═9db79764-30f7-4c0f-af9e-9596e1ef9405
# ╠═4705c138-2aff-463a-9f60-3c4d219d551f
# ╠═5f2c1f99-314b-4c18-bffc-aab774e024e9
# ╠═a5fc300d-3816-4d53-8144-3cd0f5c0f42f
# ╠═1c125abb-e404-402b-b175-c7c8d2115ac2
# ╠═d2ba032e-7387-4cb5-b239-4a2174f7b92d
# ╠═7c0c06e2-f164-467d-bbba-04b5a4f94db4
