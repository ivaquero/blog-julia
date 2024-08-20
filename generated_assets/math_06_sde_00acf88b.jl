### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> chapter = 2
#> section = 6
#> order = 6
#> title = "Stochastic Differential Equation"
#> tags = ["math", "track_comp", "Pluto", "PlutoUI"]
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

# ╔═╡ d3fde8a5-8cfe-49d4-af5d-2b5bd340f528
begin
    using Pkg
    Pkg.activate("../../pluto-deployment-environment")
    Pkg.instantiate()
end

# ╔═╡ 408a51cf-e35c-47ce-9dd6-c636bc938683
using DifferentialEquations

# ╔═╡ c4a89b4a-9b2e-417e-88f5-309482f2dcc1
using CairoMakie

# ╔═╡ f0d620bd-2606-4d52-9483-1ead2737f5aa
using PlutoUI

# ╔═╡ dbe0670d-1bb1-45f4-8673-af9b805cfec9
using PlutoUI: Button, Slider

# ╔═╡ 549c7352-6c31-11ed-2ca7-f577b5b16505
md"# Stochastic Differential Equation"

# ╔═╡ 1f3dcdb5-cbc8-4daa-8a9a-a098e7e3b936
md"Sources:"

# ╔═╡ f3cb1131-debd-49af-a931-26fa170bfc7c
md"[Wikipedia Article on Stochastic Differential Equation](https://en.wikipedia.org/wiki/Stochastic_differential_equation)"

# ╔═╡ e914b53c-a778-49fa-8e15-5aa8598b07a6
md"[DifferentialEquations.jl Documentation](https://diffeq.sciml.ai/stable/)"

# ╔═╡ 1a286802-8d7b-412f-bc0a-f177c0ccc4cf
TableOfContents()

# ╔═╡ 6b682619-f397-4fa2-8f3a-76b17651441c
md"## SDE Version of Verhulst Equation"

# ╔═╡ e22cfc96-3ae3-4402-88f3-6b38d0481eb9
md"### Define Julia Functions"

# ╔═╡ 340389f2-d8a6-42e5-992c-0aaead399c13
md"#### ODE"

# ╔═╡ f20c1489-a6b8-43d3-850b-3ac8f4b4c45c
function verhulst(uv, p, t)
    N = uv
    r, K = p
    r * N * (1 - (N / K))
end

# ╔═╡ 650d2a91-2da2-412f-b374-d27deb4e0095
md"#### Noise"

# ╔═╡ 0689f410-30a5-4cda-8512-2e43bbc9e807
md"### Assign Variables"

# ╔═╡ f60c92ed-bcc9-4827-bef6-260fe040a8bd
md"#### uv: Population"

# ╔═╡ 7a88e186-91a2-48e8-88bb-61ee0e48f61a
md"#### pv: Parameters"

# ╔═╡ 46ee364d-9e0d-4af6-8888-816db1ea0509
md"##### $r$: Growth Rate per Month"

# ╔═╡ 3d7b9676-8265-40aa-ae52-5eeb4efb85be
md"##### $K$: Maximum Population"

# ╔═╡ 1267372f-d3f2-46fe-b292-f2695c2a42c8
md"#### tv: Time in Months"

# ╔═╡ 05898844-4c13-4b6d-bd43-eee4a8550ba1
tv_begin = 0.0

# ╔═╡ 8b4618ff-b910-43e3-8757-dc144152c48e
md"Note: Time Step is required for SDE."

# ╔═╡ 8ba08b34-80d3-4191-b75c-9666f15e52b5
tv_step = 0.01

# ╔═╡ 9efab18e-3901-47cf-9a4e-e29adfa70799
md"### Define Problems"

# ╔═╡ fda74c49-5c88-4e56-8976-62e3d90be8da
md"#### Define ODE Problem"

# ╔═╡ e710b57e-866f-4b02-8d96-dce68c1ad8af
md"#### Define SDE Problem"

# ╔═╡ a10a6e86-3698-4c18-887c-23c3a48d1505
md"### Solve Problems"

# ╔═╡ e2b94f37-2c40-40da-8155-589fac67cd7a
md"#### Solve ODE Problem"

# ╔═╡ 18df1ba0-994c-47c5-a87d-5c0a1a3c1106
md"#### Solve SDE Problem"

# ╔═╡ 6a502139-535e-4b5b-9592-56b12c5bcd66
md"### Plot Solutions"

# ╔═╡ 0878c2ca-2403-485a-a0ca-c54b20054a7d
md"#### Control Panel"

# ╔═╡ 6bbfe9be-0fd2-4bdb-a8b9-8fcf8884c9ef
md"Maximum Population ($K$): $(@bind K Slider(100.0:200.0, 200.0, true))"

# ╔═╡ 329103af-22b2-479c-9f27-7f4e142ec7d3
md"Months: $(@bind tv_end Slider(0.0:72.0, 72.0, true))"

# ╔═╡ d4be6a98-54c3-45d6-830a-fc32780ea838
tvspan = (tv_begin, tv_end)

# ╔═╡ 259c1360-3e41-42e9-992d-31424b29108e
md"Beginning Population: $(@bind N_begin Slider(1.0:100.0, 1.0, true))"

# ╔═╡ ecf5f3ec-72d3-496a-9941-6c359f8b0de3
@bind new_solv Button("New Solution")

# ╔═╡ 34f7cd32-0456-4748-844c-8572e3724f78
md"Growth Rate ($r$): $(@bind r Slider(0.1:0.01:1.0, 0.14, true))"

# ╔═╡ d0c5e98c-23fe-4282-8d86-0b326bea853e
pv = [r, K]

# ╔═╡ 6824788b-80eb-4f45-bd8e-f3b33a358bd0
probv_ODE = ODEProblem(verhulst, N_begin, tvspan, pv)

# ╔═╡ 54209c09-e84e-48af-9f68-a5805e4c3f1d
solv_ODE = solve(probv_ODE)

# ╔═╡ f2f2ba06-b074-43bd-9a28-2729cb1a3c0c
md"Volatility (beta): $(@bind βv Slider(0.0:0.01:0.1, 0.0, true))"

# ╔═╡ fc5437c8-16d9-45ec-a39a-a6c4f862e081
gv(uv, p, t) = βv * uv

# ╔═╡ 414a92c5-cf92-4bd3-aff5-5a3df7c2b4cf
probv_SDE = SDEProblem(verhulst, gv, N_begin, tvspan, pv)

# ╔═╡ f2219c3f-5457-4358-b74f-0ab079bc2dd1
begin
    new_solv
    solv_SDE = solve(probv_SDE, EM(), dt=tv_step)
end

# ╔═╡ 4afe23e4-df7a-413c-a12e-44f837ecbc96
begin
    figv = Figure()

    axv = Axis(figv[1, 1],
        title="SDE Version of Verhulst Equation",
        xlabel="Time in Months",
        ylabel="Rabbit Population")

    lines!(axv, solv_ODE.t, solv_ODE.u, linewidth=2, label="ODE")
    lines!(axv, solv_SDE.t, solv_SDE.u, linewidth=2, label="SDE")
    axislegend(position=:rb)
    figv
end

# ╔═╡ 09eee192-3df3-47e3-8ac0-127b1c8b5d80
md"Ending Population: $(round(Int, solv_SDE(tv_end)))"

# ╔═╡ Cell order:
# ╟─549c7352-6c31-11ed-2ca7-f577b5b16505
# ╟─1f3dcdb5-cbc8-4daa-8a9a-a098e7e3b936
# ╟─f3cb1131-debd-49af-a931-26fa170bfc7c
# ╟─e914b53c-a778-49fa-8e15-5aa8598b07a6
# ╠═d3fde8a5-8cfe-49d4-af5d-2b5bd340f528
# ╠═408a51cf-e35c-47ce-9dd6-c636bc938683
# ╠═c4a89b4a-9b2e-417e-88f5-309482f2dcc1
# ╠═f0d620bd-2606-4d52-9483-1ead2737f5aa
# ╠═dbe0670d-1bb1-45f4-8673-af9b805cfec9
# ╠═1a286802-8d7b-412f-bc0a-f177c0ccc4cf
# ╟─6b682619-f397-4fa2-8f3a-76b17651441c
# ╟─e22cfc96-3ae3-4402-88f3-6b38d0481eb9
# ╟─340389f2-d8a6-42e5-992c-0aaead399c13
# ╠═f20c1489-a6b8-43d3-850b-3ac8f4b4c45c
# ╟─650d2a91-2da2-412f-b374-d27deb4e0095
# ╠═fc5437c8-16d9-45ec-a39a-a6c4f862e081
# ╟─0689f410-30a5-4cda-8512-2e43bbc9e807
# ╟─f60c92ed-bcc9-4827-bef6-260fe040a8bd
# ╟─7a88e186-91a2-48e8-88bb-61ee0e48f61a
# ╟─46ee364d-9e0d-4af6-8888-816db1ea0509
# ╟─3d7b9676-8265-40aa-ae52-5eeb4efb85be
# ╠═d0c5e98c-23fe-4282-8d86-0b326bea853e
# ╟─1267372f-d3f2-46fe-b292-f2695c2a42c8
# ╠═05898844-4c13-4b6d-bd43-eee4a8550ba1
# ╠═d4be6a98-54c3-45d6-830a-fc32780ea838
# ╟─8b4618ff-b910-43e3-8757-dc144152c48e
# ╠═8ba08b34-80d3-4191-b75c-9666f15e52b5
# ╟─9efab18e-3901-47cf-9a4e-e29adfa70799
# ╟─fda74c49-5c88-4e56-8976-62e3d90be8da
# ╠═6824788b-80eb-4f45-bd8e-f3b33a358bd0
# ╟─e710b57e-866f-4b02-8d96-dce68c1ad8af
# ╠═414a92c5-cf92-4bd3-aff5-5a3df7c2b4cf
# ╟─a10a6e86-3698-4c18-887c-23c3a48d1505
# ╟─e2b94f37-2c40-40da-8155-589fac67cd7a
# ╠═54209c09-e84e-48af-9f68-a5805e4c3f1d
# ╟─18df1ba0-994c-47c5-a87d-5c0a1a3c1106
# ╠═f2219c3f-5457-4358-b74f-0ab079bc2dd1
# ╟─6a502139-535e-4b5b-9592-56b12c5bcd66
# ╠═4afe23e4-df7a-413c-a12e-44f837ecbc96
# ╟─0878c2ca-2403-485a-a0ca-c54b20054a7d
# ╟─6bbfe9be-0fd2-4bdb-a8b9-8fcf8884c9ef
# ╟─329103af-22b2-479c-9f27-7f4e142ec7d3
# ╟─259c1360-3e41-42e9-992d-31424b29108e
# ╠═ecf5f3ec-72d3-496a-9941-6c359f8b0de3
# ╟─34f7cd32-0456-4748-844c-8572e3724f78
# ╟─f2f2ba06-b074-43bd-9a28-2729cb1a3c0c
# ╟─09eee192-3df3-47e3-8ac0-127b1c8b5d80
