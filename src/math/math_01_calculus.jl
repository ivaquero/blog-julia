### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> chapter = 2
#> section = 1
#> order = 1
#> title = "Calculus Basics"
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

# ╔═╡ 680cfeeb-3159-465d-9df2-f632edbe8ceb
begin
    using Pkg
    if isdir("../../pluto-deployment-environment")
        Pkg.activate("../../pluto-deployment-environment")
        Pkg.instantiate()
    else
        println(pwd())
    end
end

# ╔═╡ 34e60bb1-7e43-459a-ba02-08b63ea092d1
using ForwardDiff, ReverseDiff, QuadGK

# ╔═╡ e0440410-e7fc-4dcb-ac92-28f6c6c109fb
using CairoMakie

# ╔═╡ a99851d5-4eb3-4f0c-8c47-8d6b933cc91e
using PlutoUI

# ╔═╡ e63a9c4a-bed7-4ca0-84fd-85139270d4a1
using PlutoUI: Slider

# ╔═╡ c4dd6dd0-8896-4294-b64c-c6be3d30f30c
md"# Calculus"

# ╔═╡ 69768169-ef40-4e30-80dd-28dfb7a19d37
TableOfContents()

# ╔═╡ 5436ca50-ae35-11ed-3a4f-3930783feae1
md"## Forward Mode Automatic Differentiation"

# ╔═╡ 50cedd6d-90d1-4b18-b20a-9447711081ec
md"### Define Function"

# ╔═╡ 1ae44b2a-acdd-488d-97d8-c675e9351a89
f(x) = x^3 - x

# ╔═╡ 596746e8-23c9-479d-acfc-c7510850e7c9
md"### Find Derivative"

# ╔═╡ 4c4f8984-4beb-4cc3-b8d2-51114dfbf4bc
md"$\frac{dy}{dx} = f'(x) = \dot{y}$"

# ╔═╡ de59c452-8f1a-4469-9917-3b1f8d1fcc4b
md"$f'(x) = \lim_{h → 0} \frac{f(x + h) - f(x)}{h}$"

# ╔═╡ 58eb7e22-be78-491a-b3fa-2754704ce616
md"where $x$ is a fixed value of the function $f(x)$ and $h$ is a number close to zero."

# ╔═╡ 648c01de-c97e-434f-8ed7-846ebaa515bf
md"### Define Tangent Line Segment"

# ╔═╡ 549b91d9-e9ed-4dde-8b69-815c4d5c9c0a
md"$y = mx + b$"

# ╔═╡ cc35bdef-5828-455a-926e-5d9503d1c02c
md"where $m$ is the slope and $b$ is the y-intercept."

# ╔═╡ 04905701-d33a-4619-97f9-6d6b8a288447
x_range = -2.0:0.01:2.0

# ╔═╡ de0c7162-d323-45bc-8316-cc04f98b2c1f
md"### Control Panel"

# ╔═╡ 63259faa-6e5d-4181-bf62-192ee73129d5
md"x: $(@bind x Slider(x_range, 0.0, true))"

# ╔═╡ b63e80e3-ef66-4535-9965-0e54d9cb9e86
x_left = x - 1

# ╔═╡ 3eccc181-17a9-476f-af81-bbc4cbf6d583
x_right = x + 1

# ╔═╡ 6c0dc657-8759-4a54-b273-b8264706391d
md"## Reverse Mode Automatic Differentiation"

# ╔═╡ 95fc32a9-c2f1-4071-8dd6-1057cfab5061
md"Code Source: [ReverseDiff.jl Gradient Example in GitHub](https://github.com/JuliaDiff/ReverseDiff.jl/blob/master/examples/gradient.jl)"

# ╔═╡ 33c7b54d-f3d4-4b4d-a29c-bda4650dc7cf
md"### Define Function"

# ╔═╡ fd9eb183-2e41-4f74-8e4a-e70a12c64c36
f(A, B) = sum(A' * B + A * B')

# ╔═╡ d004b65c-c401-4571-99a5-79e2c5bc242a
y = f(x)

# ╔═╡ baba864f-155c-4714-b913-353593a23f1d
dydx = slope = ForwardDiff.derivative(f, x)

# ╔═╡ 76fb0c38-7124-4ef0-b700-eed7cb0097ca
bias = y - slope * x

# ╔═╡ ca5a9d16-df2f-4b36-aebd-dd3180380e30
y_left = slope * x_left + bias

# ╔═╡ 9357fac5-70ab-421c-8b1a-8c884ffcfb08
y_right = slope * x_right + bias

# ╔═╡ 9bb9483d-2822-44b6-aafc-243e294fb4c1
md"""
x: $(x) | f(x): $(y) | f'(x): $(dydx)
"""

# ╔═╡ 5ff7c523-755a-492e-aa8a-029d58be7534
begin
    figf = Figure()

    axf = Axis(figf[1, 1],
        title="Relationship between x, f(x) and f'(x)",
        xlabel="x",
        ylabel="f(x)")

    lines!(axf, x_range, f.(x_range))
    lines!(axf, [x_left, x_right], [y_left, y_right], linewidth=2, color=:red)
    scatter!(x, y, markersize=20, strokecolor=:black, strokewidth=2, color=:lightgreen)

    figf
end

# ╔═╡ 66ab1dba-19ed-4b89-9458-3b5872109ad5
md"### Generate Inputs by using Random Numbers"

# ╔═╡ 8fe27ce6-13fc-4e8b-8f07-65c992f8f3d9
A, B = rand(100, 100), rand(100, 100)

# ╔═╡ 88d6515f-3880-4755-aa60-d2adb475a6d5
c = f(A, B)

# ╔═╡ f22d072f-7c0a-478b-9e42-e3b62a0099ea
inputs = (A, B)

# ╔═╡ 6d6255d5-de3c-4228-a885-f305cf98a615
md"### Find Gradient"

# ╔═╡ 22f98339-fd40-40a5-9cf3-5181fdb6eb97
md"""
Note: Conceptually, a Gradient is the same as a Derivative, but the term Gradient is typically used for functions with several inputs and a single output.
"""

# ╔═╡ 43d33e42-7a6a-4d76-965a-931ae428d2dc
ReverseDiff.gradient(f, inputs)

# ╔═╡ d3a971f3-4346-43be-9440-35f927945154
md"Note: See GitHub repository for more optimal examples."

# ╔═╡ 075c34ea-f38b-424a-955e-40523a605039
md"## Numerical Integration"

# ╔═╡ 180266dc-3236-45c6-9c81-5b32a5f8699b
md"### Define Function"

# ╔═╡ a3b22cff-831d-46b5-981b-5c076f0f9216
f2(x) = sqrt(x)

# ╔═╡ 26ce383b-61ac-45f6-ae38-ed85e1d35394
md"### Assign Limits of Integration"

# ╔═╡ cbd222f3-5c67-4bbf-a697-50eae35d8721
md"$[a, b]$"

# ╔═╡ 38c99667-b150-4400-ab62-7a1d4093e3e9
md"### Evaluate Integrals"

# ╔═╡ 53a01aeb-ab29-4dc3-b1cb-67e2cb048445
md"### Area Under a Curve"

# ╔═╡ c784d200-2167-47f3-9218-3ddcfc28b164
md"$A = \int_{a}^{b} f(x) \text{d}x$"

# ╔═╡ 43b8cae1-31ab-47cc-b82c-7e8d7aad843d
md"""
"The integral from a to b of f of x, dx."
"""

# ╔═╡ 333348ea-c28a-4bbb-b520-9eb0c7afb331
md"### Control Panel"

# ╔═╡ 19cd4a95-810e-4949-8989-8ee49672389a
md"""
a: $(@bind a Slider(0.0:0.01:4.0, 0.0, true)) |
b: $(@bind b Slider(0.0:0.01:4.0, 4.0, true))
"""

# ╔═╡ 18128faa-f189-4d38-8e7c-d006b70274bc
xs = a < b ? (a:0.01:b) : (b:0.01:a)

# ╔═╡ 5e8a4df3-5333-4632-85ff-db37d3299ae0
ya = round(f2(a), digits=2)

# ╔═╡ 960b3129-c171-453c-a0fd-0ae9eb1cd9ff
yb = round(f2(b), digits=2)

# ╔═╡ b7a108c5-1285-45fb-8e49-f7f30ed9b5bd
dyda = round(ForwardDiff.derivative(f2, a), digits=2)

# ╔═╡ 46853a59-61f6-47f1-b6a9-e56b95a43f1b
dydb = round(ForwardDiff.derivative(f, b), digits=2)

# ╔═╡ c9c26468-66c0-4ab5-ae3e-a01ba940669a
area, area_error = quadgk(f2, a, b)

# ╔═╡ 2725cbf6-7002-4701-897f-6debdbc393d2
A2 = round(area, digits=2)

# ╔═╡ f5753df4-7186-4b29-b825-e496d8ae3461
begin
    figi = Figure()

    axi = Axis(figi[1, 1],
        title="Numerical Integration",
        xlabel="x",
        ylabel="f(x)")

    lines!(axi, xs, f2.(xs), linewidth=10, color=:red)
    band!(axi, xs, 0, f2.(xs), color=:lightblue)
    hlines!([0], linewidth=2, color=:green)
    xlims!(0, 4)
    ylims!(0, 2)

    figi
end

# ╔═╡ 32ac36a6-435c-42af-a42b-18d5c955a8d2
md"### Volume of a Solid of Revolution (around x-axis)"

# ╔═╡ fd31b4b4-c66c-4c18-a4f5-d52c58e4a977
md"$V = \int_{a}^{b} π [f(x)]^2 \text{d}x$"

# ╔═╡ c0181955-4b63-4161-b794-f20aafd35aa1
md"Arc Length"

# ╔═╡ ad3daca1-2443-42f9-9add-d655c6cb762d
md"$s = \int_{a}^{b} \sqrt{1 + \left( \frac{dy}{dx} \right)^2} \text{d}x$"

# ╔═╡ 0194dfa5-0cdc-4934-a9af-49d0b23ab167
md"Area of a Surface of Revolution (around x-axis)"

# ╔═╡ 56a5d53c-3e6d-43e7-9e87-693fe4434f23
md"$S = \int_{a}^{b} 2 π f(x) \sqrt{1 + \left( \frac{dy}{dx} \right)^2} \text{d}x$"

# ╔═╡ d7e3e0ca-c4aa-48b0-97f9-a3a0539528d3
volume, volume_error = quadgk(x -> π * f2(x)^2, a, b)

# ╔═╡ 058ccc42-6bd0-4d89-96f7-ec9ff99ce215
V = round(volume, digits=2)

# ╔═╡ 48fcc02a-b14c-45a8-a855-8a3ed385525b
md"""
a: $(a) | b: $(b) | f(a): $(ya) | f(b): $(yb) | f'(a): $(dyda) | f'(b): $(dydb) \
Area: $(A2) | Volume (x-axis): $(V)
"""

# ╔═╡ 756c841c-9fb0-4cd6-a21b-e154558c0749
function plotcircles!(ax, a, b)
    for x in LinRange(a, b, 20)
        θ = LinRange(0, 2π, 360)
        lines!(ax,
            fill(x, 360),
            f2(x) * cos.(θ),
            f2(x) * sin.(θ),
            linewidth=2,
            color=:dodgerblue)
    end
end

# ╔═╡ 38d5a7d3-d60d-4cdb-be14-37f7d1bc0636
begin
    figi3 = Figure()

    axi3 = Axis3(figi3[1, 1],
        title="Revolution Around X-Axis",
        xlabel="x",
        ylabel="y",
        zlabel="z",
        aspect=:data,
        azimuth=-0.3π)

    plotcircles!(axi3, a, b)
    lines!(axi3, [0, 4], [0, 0], [0, 0], color=:black)
    xlims!(0, 4),
    ylims!(-2, 2),
    zlims!(-2, 2)

    figi3
end

# ╔═╡ 9c7c27c7-6808-42a6-8555-8673898213d0
md"### Arc Length"

# ╔═╡ ec4a7031-b1ef-4502-9757-ed6bc7ef006e
md"$s = ∫_{a}^{b}\sqrt{1 + \left(\frac{dy}{dx}\right)^2} dx$"

# ╔═╡ 4398d534-bdc9-431f-b28a-40ca27a2d6af
arc_length, arc_error = quadgk(x -> sqrt(1 + ForwardDiff.derivative(f2, x)^2),
    a, b)

# ╔═╡ f6fc8569-3a55-489d-8f02-cfbbc6210d78
s = round(arc_length, digits=2)

# ╔═╡ b15f1c01-6d26-4e45-8ff0-61a41778aa81
md"### Area of a Surface of Revolution (around x-axis)"

# ╔═╡ 872aff75-1f6d-4cd6-a452-dceb49ee1889
md"$S = ∫_{a}^{b} 2πf(x) \sqrt{1 + \left(\frac{dy}{dx}\right)^2}dx$"

# ╔═╡ efbdee04-f922-432e-a374-10044b77bab1
surface_area, surface_error = quadgk(x -> 2 * pi * f2(x) * sqrt(1 + ForwardDiff.derivative(f2, x)^2),
    a, b)

# ╔═╡ e5fb18a0-01ef-4584-8c4c-dcd7eb7c9ee2
S = round(surface_area, digits=2)

# ╔═╡ 65250eb6-92a8-4f8c-80b6-bc3ad8b2cde9
md"""
a: $(a) | b: $(b) | f(a): $(ya) | f(b): $(yb) | f'(a): $(dyda) | f'(b): $(dydb) \
Area: $(A2) | Volume (x-axis): $(V) | Arc Length: $(s) | Surface Area (x-axis): $(S)
"""

# ╔═╡ Cell order:
# ╟─c4dd6dd0-8896-4294-b64c-c6be3d30f30c
# ╠═680cfeeb-3159-465d-9df2-f632edbe8ceb
# ╠═34e60bb1-7e43-459a-ba02-08b63ea092d1
# ╠═e0440410-e7fc-4dcb-ac92-28f6c6c109fb
# ╠═a99851d5-4eb3-4f0c-8c47-8d6b933cc91e
# ╠═e63a9c4a-bed7-4ca0-84fd-85139270d4a1
# ╠═69768169-ef40-4e30-80dd-28dfb7a19d37
# ╟─5436ca50-ae35-11ed-3a4f-3930783feae1
# ╟─50cedd6d-90d1-4b18-b20a-9447711081ec
# ╠═1ae44b2a-acdd-488d-97d8-c675e9351a89
# ╠═d004b65c-c401-4571-99a5-79e2c5bc242a
# ╟─596746e8-23c9-479d-acfc-c7510850e7c9
# ╟─4c4f8984-4beb-4cc3-b8d2-51114dfbf4bc
# ╟─de59c452-8f1a-4469-9917-3b1f8d1fcc4b
# ╟─58eb7e22-be78-491a-b3fa-2754704ce616
# ╠═baba864f-155c-4714-b913-353593a23f1d
# ╟─648c01de-c97e-434f-8ed7-846ebaa515bf
# ╟─549b91d9-e9ed-4dde-8b69-815c4d5c9c0a
# ╟─cc35bdef-5828-455a-926e-5d9503d1c02c
# ╠═76fb0c38-7124-4ef0-b700-eed7cb0097ca
# ╠═b63e80e3-ef66-4535-9965-0e54d9cb9e86
# ╠═3eccc181-17a9-476f-af81-bbc4cbf6d583
# ╠═ca5a9d16-df2f-4b36-aebd-dd3180380e30
# ╠═9357fac5-70ab-421c-8b1a-8c884ffcfb08
# ╠═04905701-d33a-4619-97f9-6d6b8a288447
# ╟─de0c7162-d323-45bc-8316-cc04f98b2c1f
# ╟─63259faa-6e5d-4181-bf62-192ee73129d5
# ╟─9bb9483d-2822-44b6-aafc-243e294fb4c1
# ╠═5ff7c523-755a-492e-aa8a-029d58be7534
# ╟─6c0dc657-8759-4a54-b273-b8264706391d
# ╟─95fc32a9-c2f1-4071-8dd6-1057cfab5061
# ╟─33c7b54d-f3d4-4b4d-a29c-bda4650dc7cf
# ╠═fd9eb183-2e41-4f74-8e4a-e70a12c64c36
# ╟─66ab1dba-19ed-4b89-9458-3b5872109ad5
# ╠═8fe27ce6-13fc-4e8b-8f07-65c992f8f3d9
# ╠═88d6515f-3880-4755-aa60-d2adb475a6d5
# ╠═f22d072f-7c0a-478b-9e42-e3b62a0099ea
# ╟─6d6255d5-de3c-4228-a885-f305cf98a615
# ╟─22f98339-fd40-40a5-9cf3-5181fdb6eb97
# ╠═43d33e42-7a6a-4d76-965a-931ae428d2dc
# ╟─d3a971f3-4346-43be-9440-35f927945154
# ╟─075c34ea-f38b-424a-955e-40523a605039
# ╟─180266dc-3236-45c6-9c81-5b32a5f8699b
# ╠═a3b22cff-831d-46b5-981b-5c076f0f9216
# ╟─26ce383b-61ac-45f6-ae38-ed85e1d35394
# ╟─cbd222f3-5c67-4bbf-a697-50eae35d8721
# ╟─38c99667-b150-4400-ab62-7a1d4093e3e9
# ╟─53a01aeb-ab29-4dc3-b1cb-67e2cb048445
# ╟─c784d200-2167-47f3-9218-3ddcfc28b164
# ╟─43b8cae1-31ab-47cc-b82c-7e8d7aad843d
# ╟─333348ea-c28a-4bbb-b520-9eb0c7afb331
# ╠═19cd4a95-810e-4949-8989-8ee49672389a
# ╠═18128faa-f189-4d38-8e7c-d006b70274bc
# ╠═5e8a4df3-5333-4632-85ff-db37d3299ae0
# ╠═960b3129-c171-453c-a0fd-0ae9eb1cd9ff
# ╠═b7a108c5-1285-45fb-8e49-f7f30ed9b5bd
# ╠═46853a59-61f6-47f1-b6a9-e56b95a43f1b
# ╠═c9c26468-66c0-4ab5-ae3e-a01ba940669a
# ╠═2725cbf6-7002-4701-897f-6debdbc393d2
# ╟─48fcc02a-b14c-45a8-a855-8a3ed385525b
# ╠═f5753df4-7186-4b29-b825-e496d8ae3461
# ╟─32ac36a6-435c-42af-a42b-18d5c955a8d2
# ╟─fd31b4b4-c66c-4c18-a4f5-d52c58e4a977
# ╟─c0181955-4b63-4161-b794-f20aafd35aa1
# ╟─ad3daca1-2443-42f9-9add-d655c6cb762d
# ╟─0194dfa5-0cdc-4934-a9af-49d0b23ab167
# ╟─56a5d53c-3e6d-43e7-9e87-693fe4434f23
# ╠═d7e3e0ca-c4aa-48b0-97f9-a3a0539528d3
# ╠═058ccc42-6bd0-4d89-96f7-ec9ff99ce215
# ╠═756c841c-9fb0-4cd6-a21b-e154558c0749
# ╠═38d5a7d3-d60d-4cdb-be14-37f7d1bc0636
# ╟─9c7c27c7-6808-42a6-8555-8673898213d0
# ╟─ec4a7031-b1ef-4502-9757-ed6bc7ef006e
# ╠═4398d534-bdc9-431f-b28a-40ca27a2d6af
# ╠═f6fc8569-3a55-489d-8f02-cfbbc6210d78
# ╟─b15f1c01-6d26-4e45-8ff0-61a41778aa81
# ╟─872aff75-1f6d-4cd6-a452-dceb49ee1889
# ╠═efbdee04-f922-432e-a374-10044b77bab1
# ╠═e5fb18a0-01ef-4584-8c4c-dcd7eb7c9ee2
# ╟─65250eb6-92a8-4f8c-80b6-bc3ad8b2cde9
