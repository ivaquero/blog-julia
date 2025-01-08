### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> chapter = 2
#> section = 2
#> order = 2
#> title = "First Order ODE"
#> tags = ["math", "track_comp", "Pluto", "PlutoUI"]
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

# ╔═╡ 42f2dc7e-b3db-4c45-8ba6-3a1470e01edc
begin
    using Pkg
    if isdir("../../pluto-deployment-environment")
        Pkg.activate("../../pluto-deployment-environment")
        Pkg.instantiate()
    else
        println(pwd())
    end
end

# ╔═╡ 17c89fe7-4d64-4b3f-9d5e-5f57f932b3eb
using DifferentialEquations

# ╔═╡ 4ea98aee-b137-43ee-abb0-c6e843f7f1ef
using PlutoUI

# ╔═╡ 813e4fc2-f47c-4c89-ba74-8a1bdc0c3b5f
using PlutoUI: Slider

# ╔═╡ c1ec8f90-e684-4f4e-aa09-63c7401ad5f5
using CairoMakie

# ╔═╡ 519a539d-e6e6-4dd3-9e62-7b7a23cb4697
md"# Ordinary Differential Equation - 1st Order"

# ╔═╡ bb18a2be-70aa-468b-8302-5aa3d853f491
md"## Add Packages"

# ╔═╡ e02a76c8-60b6-4a0d-9fa8-95ca65e39625
TableOfContents()

# ╔═╡ a1ba8b4e-404d-11ed-28bd-93e136588e06
md"## Continuous Compound Interest"

# ╔═╡ 8aa21f28-3ac2-40f9-91c1-d1664883dbc0
md"### Describe Problem"

# ╔═╡ 4a4a9d78-85c3-4ffc-9f5d-8867b65c7d28
md"#### ODE in DifferentialEquations.jl"

# ╔═╡ d81f7702-032f-4a9b-b93b-8e98f066797f
md"### Define Function"

# ╔═╡ 247ee8fe-d88f-40de-9353-0e3f7acfd1a4
md"#### Continuous Compound Interest in Math"

# ╔═╡ 5918ce14-8161-4c72-972b-336c7f2952a6
md"$\frac{du}{dt} = pu$"

# ╔═╡ cd7af692-c403-40bc-a9b1-afedcf084efe
md"#### Continuous Compound Interest in Julia"

# ╔═╡ 9b8c6636-91ab-4dbc-a9c8-beb6100c9878
f(u, p, t) = p * u

# ╔═╡ c0b079e3-ff90-40e4-a46d-b87e17883f84
md"### Assign Variables"

# ╔═╡ 3639bbff-9572-4b3c-902f-5a3ad11376a2
md"#### u: Bank Account Balance"

# ╔═╡ e1be3ed7-7325-4a69-ae2c-72feaea13e57
md"#### p: Annual Interest Rate"

# ╔═╡ f20019d3-c4fb-405e-966b-6c9a16aab9a0
md"#### t: Time in Years"

# ╔═╡ b0e36de7-07dd-4c03-a814-32f2b81d8625
t_begin = 0.0

# ╔═╡ 0f52dcf4-22c1-4e86-946f-e742bae934ad
md"### Define ODE Problem"

# ╔═╡ 6c61c8bf-b06b-483c-850a-137f7296861d
md"### Solve Problem"

# ╔═╡ cd06f8bb-7c29-4924-9244-62dacd44bc2d
md"### Plot Solution"

# ╔═╡ f7124c0e-b647-42bf-a3b0-4859ff457f1f
md"Years: $(@bind t_end Slider(1.0:50.0, 1.0, true))"

# ╔═╡ 503e6643-069a-440b-b8f2-21c6a9f623a2
tspan = (t_begin, t_end)

# ╔═╡ 242a0d36-7b52-4f34-a084-5c7c2d85550d
md"Interest Rate: $(@bind p Slider(0.01:0.01:1.0, 1.0, true))"

# ╔═╡ 2f8255fc-59ae-4c49-bd77-5702c743cfa2
md"Beginning Balance: $(@bind u_begin Slider(1.0:1000.0, 1.0, true))"

# ╔═╡ dae87a87-e3e2-480e-8ca3-04e2d9b11a96
prob = ODEProblem(f, u_begin, tspan, p)

# ╔═╡ d33379e9-ddf5-4cb7-8208-20d4513c38c8
sol = solve(prob)

# ╔═╡ 36a2f7d4-2ee2-43d9-a466-4014d233b79c
md"Ending Balance: $(round(Int, sol(t_end)))"

# ╔═╡ 15ffebdf-7835-4170-8757-73a254058b7d
begin
    fig = Figure()

    ax = Axis(fig[1, 1],
        title="Continuous Compound Interest",
        xlabel="Time in Years",
        ylabel="Account Balance in Dollars")

    lines!(ax, sol, linewidth=2)

    # ax.xlims=0.0:50.0
    # ax.ylims=0.0:1_000_000.0
    fig
end

# ╔═╡ 4e6ccd72-7dd5-40a3-8d63-825c8afc78f6
sol(1)

# ╔═╡ 84fb94d7-6244-4349-b998-b5fdf34c50e1
md"### Appendix: Euler's Number"

# ╔═╡ d37b6161-36fa-46e5-85c0-77c81be10bd0
EULER = ℯ

# ╔═╡ 5ce0e4ae-a3e4-492a-95c9-6c927d39a765
e_est = sol(1)

# ╔═╡ 79bb60af-0d72-4ac1-96a5-7cf45b67984d
EULER - e_est

# ╔═╡ 44ba2f50-f27e-446f-b9a0-f3c990de3531
md"## Verhulst Equation"

# ╔═╡ 467203ea-3ff4-4c61-b2b6-0861580d3cf2
md"Reference: [Wikipedia Article](https://en.wikipedia.org/wiki/Pierre_Fran%C3%A7ois_Verhulst)"

# ╔═╡ f905649d-2306-4584-8f6f-bb90ed00e7f6
md"### Define Function"

# ╔═╡ a998e20e-a1d6-4df3-a572-715f8fc31e1b
md"#### Verhulst Equation in Math"

# ╔═╡ 64eb3362-0286-4fb0-9d6a-2ba1823f5f94
md"Version of equation popularized by Raymond Pearl and Lowell Reed:"

# ╔═╡ 3639f9c2-8c01-462b-8a54-de6fd29428fb
md"$\frac{dN}{dt} = rN \left( 1 - \frac{N}{K} \right)$"

# ╔═╡ 5f62077b-b7d2-4b18-91b2-5e3ef9df5cfa
md"where $N$ is the population at any time $t$,"

# ╔═╡ 6d8429fc-f619-409e-be6a-42a5376bdb22
md"and $r$ is the growth rate of the population,"

# ╔═╡ eb1b5533-3081-466f-80af-534756a2459a
md"and $K$ is the maximum population that the environment can support."

# ╔═╡ a61fed92-0b7c-47c4-908a-a181ae54dfde
md"#### Verhulst Equation in Julia"

# ╔═╡ 3a3a222c-7146-448a-a772-bcacf47e2458
md"$f(u, p, t) = p[1]u \left( 1 - \frac{u}{p[2]} \right)$"

# ╔═╡ 4ca3be54-5127-4737-b993-e9153984a8e8
# f(u, p, t) = p[1] * u * (1 - (u / p[2]))

# ╔═╡ 2af33e75-8779-474d-a3e2-c53da7e338f6
function verhulst(u, p, t)
    N = u
    r, K = p
    r * N * (1 - (N / K))
end

# ╔═╡ c8863b1e-94d1-44a6-bf04-c2f238ad3b02
md"#### Assign Variables"

# ╔═╡ cb80d2dd-ffc4-48ba-abaf-22e0121db6c4
md"#### u: Population"

# ╔═╡ 82ea56f7-8bdf-4dad-b465-4e085af15de3
md"#### p: Parameters"

# ╔═╡ 948bc4db-3063-424f-9aff-eb0a4dac66a9
md"#### $r$: Growth Rate per Month"

# ╔═╡ f22f48e8-c580-4299-93cc-8084320a2e26
md"#### $K$: Maximum Population"

# ╔═╡ 5247452f-c73b-497a-9d28-420e48d66a41
md"#### t: Time in Months"

# ╔═╡ 84adc789-646f-48f6-9e04-006063728d23
tv_begin = 0.0

# ╔═╡ 45bd6b85-000b-4e5d-a03f-cceeff350662
md"### Define ODE Problem"

# ╔═╡ 8d5116e9-8944-4d73-a560-7d3a4d4c4493
md"### Solve Problem"

# ╔═╡ e3cbef76-759a-463a-a363-66931f10ab09
md"### Plot Solution"

# ╔═╡ 1b2e68fc-cbb1-4da2-8ccf-f6d4f0c01cb1
md"Maximum Population ($K$): $(@bind K Slider(100.0:200.0, 200.0, true))"

# ╔═╡ 18e16cea-3474-48f2-8871-4a25700ed0da
md"Months: $(@bind tv_end Slider(1.0:72.0, 72.0, true))"

# ╔═╡ fe38ae5a-4985-491e-8eef-b1689dcc8f9c
tvspan = (tv_begin, tv_end)

# ╔═╡ 71070eab-2de5-4899-b3cb-8a460ba13afe
md"Growth Rate ($r$): $(@bind r Slider(0.01:0.01:1.0, 0.14, true))"

# ╔═╡ af35dad9-dc2f-4c89-b52e-79c188fc777d
pv = [r, K]

# ╔═╡ e8280117-1592-4399-a335-23165329bf4d
md"Beginning Population: $(@bind N_begin Slider(1.0:20.0, 1.0, true))"

# ╔═╡ eeef9e2d-9d8e-4e54-92fe-2f002c03acd3
probv = ODEProblem(verhulst, N_begin, tvspan, pv)

# ╔═╡ 2ac53aa0-7eee-4a4e-8b8d-088952753b99
solv = solve(probv)

# ╔═╡ 4ffca09c-f841-47bd-b0c9-61d64b7f47d2
md"Ending Population: $(round(Int, sol(tv_end)))"

# ╔═╡ 3e6ce89d-d076-4cb7-a7b3-2614412c7db1
begin
    figv = Figure()

    axv = Axis(figv[1, 1],
        title="Verhulst Equation",
        xlabel="Time in Months",
        ylabel="Rabbit Population")

    lines!(axv, solv, linewidth=2)

    # ax.xlims=0.0:50.0
    # ax.ylims=0.0:1
    figv
end

# ╔═╡ Cell order:
# ╟─519a539d-e6e6-4dd3-9e62-7b7a23cb4697
# ╟─bb18a2be-70aa-468b-8302-5aa3d853f491
# ╠═42f2dc7e-b3db-4c45-8ba6-3a1470e01edc
# ╠═17c89fe7-4d64-4b3f-9d5e-5f57f932b3eb
# ╠═4ea98aee-b137-43ee-abb0-c6e843f7f1ef
# ╠═813e4fc2-f47c-4c89-ba74-8a1bdc0c3b5f
# ╠═c1ec8f90-e684-4f4e-aa09-63c7401ad5f5
# ╠═e02a76c8-60b6-4a0d-9fa8-95ca65e39625
# ╟─a1ba8b4e-404d-11ed-28bd-93e136588e06
# ╟─8aa21f28-3ac2-40f9-91c1-d1664883dbc0
# ╟─4a4a9d78-85c3-4ffc-9f5d-8867b65c7d28
# ╟─d81f7702-032f-4a9b-b93b-8e98f066797f
# ╟─247ee8fe-d88f-40de-9353-0e3f7acfd1a4
# ╟─5918ce14-8161-4c72-972b-336c7f2952a6
# ╟─cd7af692-c403-40bc-a9b1-afedcf084efe
# ╠═9b8c6636-91ab-4dbc-a9c8-beb6100c9878
# ╟─c0b079e3-ff90-40e4-a46d-b87e17883f84
# ╟─3639bbff-9572-4b3c-902f-5a3ad11376a2
# ╟─e1be3ed7-7325-4a69-ae2c-72feaea13e57
# ╟─f20019d3-c4fb-405e-966b-6c9a16aab9a0
# ╠═b0e36de7-07dd-4c03-a814-32f2b81d8625
# ╠═503e6643-069a-440b-b8f2-21c6a9f623a2
# ╟─0f52dcf4-22c1-4e86-946f-e742bae934ad
# ╠═dae87a87-e3e2-480e-8ca3-04e2d9b11a96
# ╟─6c61c8bf-b06b-483c-850a-137f7296861d
# ╠═d33379e9-ddf5-4cb7-8208-20d4513c38c8
# ╟─cd06f8bb-7c29-4924-9244-62dacd44bc2d
# ╟─f7124c0e-b647-42bf-a3b0-4859ff457f1f
# ╟─242a0d36-7b52-4f34-a084-5c7c2d85550d
# ╟─2f8255fc-59ae-4c49-bd77-5702c743cfa2
# ╟─36a2f7d4-2ee2-43d9-a466-4014d233b79c
# ╠═15ffebdf-7835-4170-8757-73a254058b7d
# ╠═4e6ccd72-7dd5-40a3-8d63-825c8afc78f6
# ╟─84fb94d7-6244-4349-b998-b5fdf34c50e1
# ╠═d37b6161-36fa-46e5-85c0-77c81be10bd0
# ╠═5ce0e4ae-a3e4-492a-95c9-6c927d39a765
# ╠═79bb60af-0d72-4ac1-96a5-7cf45b67984d
# ╟─44ba2f50-f27e-446f-b9a0-f3c990de3531
# ╟─467203ea-3ff4-4c61-b2b6-0861580d3cf2
# ╟─f905649d-2306-4584-8f6f-bb90ed00e7f6
# ╟─a998e20e-a1d6-4df3-a572-715f8fc31e1b
# ╟─64eb3362-0286-4fb0-9d6a-2ba1823f5f94
# ╟─3639f9c2-8c01-462b-8a54-de6fd29428fb
# ╟─5f62077b-b7d2-4b18-91b2-5e3ef9df5cfa
# ╟─6d8429fc-f619-409e-be6a-42a5376bdb22
# ╟─eb1b5533-3081-466f-80af-534756a2459a
# ╟─a61fed92-0b7c-47c4-908a-a181ae54dfde
# ╟─3a3a222c-7146-448a-a772-bcacf47e2458
# ╠═4ca3be54-5127-4737-b993-e9153984a8e8
# ╠═2af33e75-8779-474d-a3e2-c53da7e338f6
# ╟─c8863b1e-94d1-44a6-bf04-c2f238ad3b02
# ╟─cb80d2dd-ffc4-48ba-abaf-22e0121db6c4
# ╟─82ea56f7-8bdf-4dad-b465-4e085af15de3
# ╟─948bc4db-3063-424f-9aff-eb0a4dac66a9
# ╟─f22f48e8-c580-4299-93cc-8084320a2e26
# ╠═af35dad9-dc2f-4c89-b52e-79c188fc777d
# ╟─5247452f-c73b-497a-9d28-420e48d66a41
# ╠═84adc789-646f-48f6-9e04-006063728d23
# ╠═fe38ae5a-4985-491e-8eef-b1689dcc8f9c
# ╟─45bd6b85-000b-4e5d-a03f-cceeff350662
# ╠═eeef9e2d-9d8e-4e54-92fe-2f002c03acd3
# ╟─8d5116e9-8944-4d73-a560-7d3a4d4c4493
# ╠═2ac53aa0-7eee-4a4e-8b8d-088952753b99
# ╟─e3cbef76-759a-463a-a363-66931f10ab09
# ╟─1b2e68fc-cbb1-4da2-8ccf-f6d4f0c01cb1
# ╟─18e16cea-3474-48f2-8871-4a25700ed0da
# ╟─71070eab-2de5-4899-b3cb-8a460ba13afe
# ╟─e8280117-1592-4399-a335-23165329bf4d
# ╟─4ffca09c-f841-47bd-b0c9-61d64b7f47d2
# ╠═3e6ce89d-d076-4cb7-a7b3-2614412c7db1
