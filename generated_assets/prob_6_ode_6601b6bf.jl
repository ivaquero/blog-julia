### A Pluto.jl notebook ###
# v0.20.13

#> [frontmatter]
#> chapter = 3
#> section = 6
#> order = 6
#> title = "Probabilistic ODE"
#> tags = ["prob", "track_prob", "Pluto", "PlutoUI"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 45758ef7-a24a-4e6a-8eba-b121abd70a17
begin
    using Pkg
    if isdir("../../pluto-deployment-environment")
        Pkg.activate("../../pluto-deployment-environment")
        Pkg.instantiate()
    else
        println(pwd())
    end
end

# ╔═╡ e78e9fce-4fd2-4dc8-9fd6-80c906926e45
using DifferentialEquations

# ╔═╡ 89cb583f-4e9b-4106-a282-1342a5e8d3bc
using Turing, CairoMakie

# ╔═╡ 93e15e01-1c21-424c-bfce-f0559593b7a4
using Tidier

# ╔═╡ fedfadde-ba7d-4d2e-96d7-cc513c196ebb
using PlutoUI

# ╔═╡ 590bc4e9-c263-4f46-b154-29957e3de1f5
using PlutoUI: Button, Slider

# ╔═╡ 325183d1-af22-448a-b0c6-bb72aef8b7d1
md"# Bayesian ODE"

# ╔═╡ 3ded2235-cba5-438e-a917-11ca74d271e8
TableOfContents()

# ╔═╡ f7bb2296-058c-494a-aeb2-113458a1a48b
md"## Describe Problem"

# ╔═╡ a43ed072-d89b-403f-9167-9148a6c25e4d
md"### Lotka-Volterra Equations"

# ╔═╡ 8bee0c36-afb2-4f67-b2de-b45807e0c967
md"$\begin{align}
\frac{dx}{dt} &= α x - β xy \\
\frac{dy}{dt} &= δ xy - γ y
\end{align}$"

# ╔═╡ 21481a9b-9a85-4d67-95eb-1506689ed934
md"where
-  $x$: the number of prey
-  $y$ is the number of some predator
-  $\frac{dx}{dt}$ and $\frac{dy}{dt}$: the instantaneous growth rates of the two populations
-  $t$: time"

# ╔═╡ 7220ef08-0a7e-4584-9d9f-b6ad94d719dd
md" $α, β, δ, γ$: positive real parameters describing the interaction of the two species,
-  $α$: the prey population growth rate and $β$ is the prey population decline rate.
-  $δ$ is the predator population growth rate and $γ$ is the predator population decline rate."

# ╔═╡ 458c1078-8d09-47af-887c-6bcc987e54b7
md"### DifferentialEquations.jl"

# ╔═╡ 81169a91-2cfe-46c5-9bd1-190b223289c3
md"$f(du, u, p, t) =
\begin{cases}
du[1] = p[1] u[1] - p[2] u[1] u[2] \\
du[2] = p[3] u[1] u[2] - p[4] u[2]
\end{cases}$"

# ╔═╡ 3f2ab2be-e841-4065-9212-58b5612efd90
md"Define Julia Function"

# ╔═╡ 730a4d87-58a6-4a8f-a911-1543a10c7f3e
function lotka_volterra(du, u, p, t)
    x, y = u
    α, β, δ, γ = p
    du[1] = dx = α * x - β * x * y
    du[2] = dy = δ * x * y - γ * y
end

# ╔═╡ 9ed08a32-f22c-477b-972d-d8c3ca3eaed7
md"### Assign Variables"

# ╔═╡ 72b75940-0d51-42ef-87fa-6c630734a7a0
md"- u: Populations
- x: Populations of Rabbits
- y: Population of Foxes
- p: Parameters
- α: Rabbit Population Growth Rate
- β: Rabbit Decline Rate
- δ: Fox Population Growth Rate
- γ: Fox Population Decline Rate
- t: Time
"

# ╔═╡ 811348cc-b477-47e1-a3d8-044cb61cc976
t_begin = 0

# ╔═╡ 3984f6c7-4a79-4ce7-b71e-250df72fda6b
md"## Visualize Model"

# ╔═╡ 05ee5445-251a-4582-914f-55128dc52a90
md"### Interactive Model"

# ╔═╡ c0843171-614d-4648-8337-e17ee4565717
@bind reset Button("Reset")

# ╔═╡ 6345b63f-76e1-4b1a-b8d1-ba11203cc573
begin
    reset
    md"Time (t): $(@bind t_end Slider(0.0:0.1:10.0; default=10.0, show_value=true))"
end

# ╔═╡ 3b15ae68-f579-40fb-9dd0-e30c320b7a3b
tspan = (t_begin, t_end)

# ╔═╡ 409f7e27-4f4c-4d5d-bdfd-352cab2436b5
begin
    reset
    md"""
    Initial Rabbit Population:
    	$(@bind x_begin Scrubbable(1:5, default = 1)) |
    Initial Fox Population:
    	$(@bind y_begin Scrubbable(1:5, default = 1))
    """
end

# ╔═╡ 9d54532d-75ad-40ef-a273-c422cdafdaba
u_begin = [x_begin, y_begin]

# ╔═╡ c9799f34-eee6-4bbb-8d24-90eb3404368e
begin
    reset
    md"""
    - Rabbit Growth (α): $(@bind α Scrubbable(0.5:0.01:3.0, default = 2.0, format = ".2f")) || Rabbit Decline (β): $(@bind β Scrubbable(0.5:0.01:3.0, default = 1.2, format = ".2f"))
    - Fox Growth (δ): $(@bind δ Scrubbable(0.5:0.01:3.0, default = 0.9, format = ".2f")) || Fox Decline (γ): $(@bind γ Scrubbable(0.5:0.01:3.0, default = 2.9, format = ".2f"))
    """
end

# ╔═╡ 02f55723-0802-4428-8652-5ed888124cdf
p = [α, β, δ, γ]

# ╔═╡ 8c6fc5a9-d8d7-44b5-9d58-d465d9e9a689
prob = ODEProblem(lotka_volterra, u_begin, tspan, p)

# ╔═╡ f321f07a-5434-4263-8ebc-454d1adcf64c
sol = solve(prob, Tsit5())

# ╔═╡ 7291c67e-b132-4380-84d3-becc95897500
rabbit = sol(t_end)[1]

# ╔═╡ 33e7933c-3ad4-4b67-a138-0f80d1af8f15
fox = sol(t_end)[2]

# ╔═╡ 4ea723ed-da5f-4e59-8fdf-1c85a5927dd0
md"### Plotting Model"

# ╔═╡ 2312e65f-28fb-44ce-8f23-781712e9b43f
md"""
Populations at Time t = $(t_end) |
Rabbits: $(round(Int, rabbit)) |
Foxes: $(round(Int, fox))
"""

# ╔═╡ 7aaeaf6f-39b3-4724-9f3c-ea080cb12e14
begin
    fig = Figure()

    ax = Axis(fig[1, 1],
        title="Lotka-Volterra Equations",
        xlabel="Time in Years",
        ylabel="Population")

    labels = ["Rabbits", "Foxes"]

    for (col, label) in enumerate(labels)
        lines!(ax, sol.t, sol[col, :], linewidth=2, label=label)
    end

    axislegend()

    fig
end

# ╔═╡ e31d2dbe-be69-4767-89dc-df4b6fc43100
md"""
Populations at Year $(t_end) |
Rabbits: $(round(Int, rabbit)) |
Foxes: $(round(Int, fox))
"""

# ╔═╡ 875253cb-c94b-4da1-8ef4-884b51f617ef
begin
    fig2 = Figure()

    ax2 = Axis(fig2[1, 1],
        title="Lotka-Volterra Equations (Phase Space)",
        xlabel="Rabbits Population",
        ylabel="Foxes Population")

    lines!(ax2, sol[1, :], sol[2, :], linewidth=2)
    scatter!(ax2, rabbit, fox, color=:red, markersize=15)

    fig2
end

# ╔═╡ 0f788ab0-1b80-4b9e-ad51-fe0b595bae37
md"## Bayesian Model"

# ╔═╡ a2717d8e-cfdf-4528-b013-6bb7bbb9be01
baboons = read_csv("../data/baboons_cheetahs.csv", col_names=false)

# ╔═╡ 742bf77a-53a9-4848-bb0e-b4ba18983e81
data = baboons |> Matrix

# ╔═╡ 623877ad-7c1a-4013-a602-9e4af23b530f
ts = 0:0.1:10

# ╔═╡ 1569a1da-ca65-43a0-9e0e-c5772688620d
begin
    fig3 = Figure()
    ax3 = Axis(fig3[1, 1],
        title="Bayesian Differential Equations",
        xlabel="t",
        ylabel="Population")

    for ind in 1:2
        lines!(ax3,
            ts,
            data[ind, :],
            linewidth=2,
            linestyle=:dot,
            label=labels[ind])
    end

    axislegend()

    # xlims!(0.0, 10.0)
    # ylims!(0.0, 10.0)
    # ax3.xticks(collect(0:10))
    # ax3.yticks(collect(0:10))

    fig3
end

# ╔═╡ ff8a878d-f446-46ba-9583-bdd29670213b
md"### Assign Variables"

# ╔═╡ fab5b721-565b-4d2d-97f2-a432d3b98f2f
md"
- mx_begin: initial population of rabbits
- my_begin: initial population of foxes
- mu_begin: initial population
"

# ╔═╡ b4365890-063c-447a-8416-e7bac60de041
mx_begin = data[1, 1]

# ╔═╡ 6634f7dc-a324-4fce-bde1-dabe18c95d18
my_begin = data[2, 1]

# ╔═╡ c8352fe9-bf4a-411e-b4ea-8a543bd1dcd8
mu_begin = [mx_begin, my_begin]

# ╔═╡ fe1f9435-08e2-45e6-a818-cd89ce594368
md"### Define Model"

# ╔═╡ c8d3bdbf-0f5b-48c6-9c4d-3e3ac0635d27
@model function lv_model(data, prob)
    # prior
    alpha ~ truncated(Normal(1.5, 1); lower=0, upper=5)
    beta ~ truncated(Normal(1.5, 1); lower=0, upper=5)
    delta ~ truncated(Normal(1.5, 1); lower=0, upper=5)
    gamma ~ truncated(Normal(1.5, 1); lower=0, upper=5)
    σ ~ Uniform(0, 1)
    # likelihood
    p = [alpha, beta, delta, gamma]
    sol = solve(prob, Tsit5(); p=p, saveat=0.1)
    for i in eachindex(sol)
        data[:, i] ~ MvNormal(sol[i], σ)
    end
end

# ╔═╡ 2cd7623d-0547-4617-9f55-9c4002f9ba0f
md"### Infer Posterior Probability"

# ╔═╡ 386d019a-f6be-47f8-8240-f7b5bf1b28b6
model_lv = lv_model(data, prob)

# ╔═╡ d5c945ff-443e-47de-b1bd-4f5465249b86
sampler = NUTS()

# ╔═╡ 5525d13a-8226-47e3-821c-3468a1752462
samples = 1_000

# ╔═╡ e4819bdc-7bbb-48a8-9c89-ae98c95ccb7b
chain_lv = sample(model_lv, sampler, samples)

# ╔═╡ a262b356-9749-4f76-ad2e-7eae6d53791d
begin
    fig4 = Figure(size=(800, 800))

    params = ["alpha", "beta", "delta", "gamma", "σ"]

    for (ind, param) in enumerate(params)
        ax = Axis(fig4[ind, 1],
            title="$param",
            xlabel="Iteration",
            ylabel="Sample Values")

        lines!(ax,
            1:samples,
            chain_lv[param][:, 1],
            color=(:lightblue, 0.8))

        ax2 = Axis(fig4[ind, 2],
            title="$param",
            xlabel="Iteration",
            ylabel="Sample Values")

        density!(ax2,
            chain_lv[param][:, 1],
            color=(:lightblue, 0.8))
    end

    fig4
end

# ╔═╡ 5fde2385-905d-4c16-bff5-8b2e3f5aad3e
md"### Retrodiction"

# ╔═╡ b82db2ce-6303-4cb5-b861-ac5ed8e12aa2
md"sample from posterior distributions"

# ╔═╡ 75a46c17-d19b-49fe-b09b-016c8a31b2de
posterior_samples = sample(chain_lv[[:alpha, :beta, :delta, :gamma]], 300;
    replace=false)

# ╔═╡ a5f5ff09-8415-4231-915e-747c7a113edb
md"visualize posterior samples"

# ╔═╡ 9ef96240-9140-4ca7-85a9-4a5b3bef7d81
begin
    fig5 = Figure(size=(800, 800))

    params2 = ["alpha", "beta", "delta", "gamma"]

    for (ind, param) in enumerate(params2)
        ax = Axis(fig5[ind, 1],
            title="$param",
            xlabel="Iteration",
            ylabel="Sample Values")

        lines!(ax,
            1:300,
            posterior_samples[param][:, 1],
            color=(:lightblue, 0.8))

        ax2 = Axis(fig5[ind, 2],
            title="$param",
            xlabel="Iteration",
            ylabel="Sample Values")

        density!(ax2,
            posterior_samples[param][:, 1],
            color=(:lightblue, 0.8))
    end

    fig5
end

# ╔═╡ e3a6e693-a668-4953-ab56-8725e6425c89
begin
    fig6 = Figure()

    ax6 = Axis(fig6[1, 1],
        title="Lotka-Volterra Equations (Phase Space)",
        xlabel="Rabbits Population",
        ylabel="Foxes Population")

    for ind in eachrow(Array(posterior_samples))
        sol_samples = solve(prob, Tsit5(); p=ind, saveat=0.1)
        lines!(ax6, sol_samples;
            label="",
            linewidth=0.1,
            colormap=(:Spectral, 0.5))
    end

    fig6
end

# ╔═╡ Cell order:
# ╟─325183d1-af22-448a-b0c6-bb72aef8b7d1
# ╠═45758ef7-a24a-4e6a-8eba-b121abd70a17
# ╠═e78e9fce-4fd2-4dc8-9fd6-80c906926e45
# ╠═89cb583f-4e9b-4106-a282-1342a5e8d3bc
# ╠═93e15e01-1c21-424c-bfce-f0559593b7a4
# ╠═fedfadde-ba7d-4d2e-96d7-cc513c196ebb
# ╠═590bc4e9-c263-4f46-b154-29957e3de1f5
# ╠═3ded2235-cba5-438e-a917-11ca74d271e8
# ╟─f7bb2296-058c-494a-aeb2-113458a1a48b
# ╟─a43ed072-d89b-403f-9167-9148a6c25e4d
# ╟─8bee0c36-afb2-4f67-b2de-b45807e0c967
# ╟─21481a9b-9a85-4d67-95eb-1506689ed934
# ╟─7220ef08-0a7e-4584-9d9f-b6ad94d719dd
# ╟─458c1078-8d09-47af-887c-6bcc987e54b7
# ╟─81169a91-2cfe-46c5-9bd1-190b223289c3
# ╟─3f2ab2be-e841-4065-9212-58b5612efd90
# ╠═730a4d87-58a6-4a8f-a911-1543a10c7f3e
# ╟─9ed08a32-f22c-477b-972d-d8c3ca3eaed7
# ╟─72b75940-0d51-42ef-87fa-6c630734a7a0
# ╠═811348cc-b477-47e1-a3d8-044cb61cc976
# ╟─3984f6c7-4a79-4ce7-b71e-250df72fda6b
# ╟─05ee5445-251a-4582-914f-55128dc52a90
# ╠═c0843171-614d-4648-8337-e17ee4565717
# ╟─6345b63f-76e1-4b1a-b8d1-ba11203cc573
# ╠═3b15ae68-f579-40fb-9dd0-e30c320b7a3b
# ╟─409f7e27-4f4c-4d5d-bdfd-352cab2436b5
# ╠═9d54532d-75ad-40ef-a273-c422cdafdaba
# ╟─c9799f34-eee6-4bbb-8d24-90eb3404368e
# ╠═02f55723-0802-4428-8652-5ed888124cdf
# ╠═8c6fc5a9-d8d7-44b5-9d58-d465d9e9a689
# ╠═f321f07a-5434-4263-8ebc-454d1adcf64c
# ╠═7291c67e-b132-4380-84d3-becc95897500
# ╠═33e7933c-3ad4-4b67-a138-0f80d1af8f15
# ╟─4ea723ed-da5f-4e59-8fdf-1c85a5927dd0
# ╟─2312e65f-28fb-44ce-8f23-781712e9b43f
# ╠═7aaeaf6f-39b3-4724-9f3c-ea080cb12e14
# ╟─e31d2dbe-be69-4767-89dc-df4b6fc43100
# ╠═875253cb-c94b-4da1-8ef4-884b51f617ef
# ╟─0f788ab0-1b80-4b9e-ad51-fe0b595bae37
# ╠═a2717d8e-cfdf-4528-b013-6bb7bbb9be01
# ╠═742bf77a-53a9-4848-bb0e-b4ba18983e81
# ╠═623877ad-7c1a-4013-a602-9e4af23b530f
# ╠═1569a1da-ca65-43a0-9e0e-c5772688620d
# ╟─ff8a878d-f446-46ba-9583-bdd29670213b
# ╟─fab5b721-565b-4d2d-97f2-a432d3b98f2f
# ╠═b4365890-063c-447a-8416-e7bac60de041
# ╠═6634f7dc-a324-4fce-bde1-dabe18c95d18
# ╠═c8352fe9-bf4a-411e-b4ea-8a543bd1dcd8
# ╟─fe1f9435-08e2-45e6-a818-cd89ce594368
# ╠═c8d3bdbf-0f5b-48c6-9c4d-3e3ac0635d27
# ╟─2cd7623d-0547-4617-9f55-9c4002f9ba0f
# ╠═386d019a-f6be-47f8-8240-f7b5bf1b28b6
# ╠═d5c945ff-443e-47de-b1bd-4f5465249b86
# ╠═5525d13a-8226-47e3-821c-3468a1752462
# ╠═e4819bdc-7bbb-48a8-9c89-ae98c95ccb7b
# ╠═a262b356-9749-4f76-ad2e-7eae6d53791d
# ╟─5fde2385-905d-4c16-bff5-8b2e3f5aad3e
# ╟─b82db2ce-6303-4cb5-b861-ac5ed8e12aa2
# ╠═75a46c17-d19b-49fe-b09b-016c8a31b2de
# ╟─a5f5ff09-8415-4231-915e-747c7a113edb
# ╠═9ef96240-9140-4ca7-85a9-4a5b3bef7d81
# ╠═e3a6e693-a668-4953-ab56-8725e6425c89
