### A Pluto.jl notebook ###
# v0.20.3

#> [frontmatter]
#> chapter = 2
#> section = 5
#> order = 5
#> title = "Attractors"
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

# ╔═╡ 42059750-855e-49e4-87e3-0c2238aff2af
begin
    using Pkg
    if isdir("../../pluto-deployment-environment")
        Pkg.activate("../../pluto-deployment-environment")
        Pkg.instantiate()
    else
        println(pwd())
    end
end

# ╔═╡ 0a619df8-1d09-4fc4-bfb3-ea29d07f6e22
using DifferentialEquations

# ╔═╡ d0768431-0b3a-47ca-bbe7-9721a6630ddf
using CairoMakie

# ╔═╡ 15a108e1-8f63-464a-8060-5a466e9e7ccb
using PlutoUI

# ╔═╡ 3c700cb2-3194-4ed1-ba03-def708a04867
using PlutoUI: Button, Slider

# ╔═╡ 8e784e49-e0b4-424d-8fee-1cf29fcfb67f
md"# Attractors"

# ╔═╡ 9026f467-8f1c-4aaa-8b42-c3b2843b4949
TableOfContents()

# ╔═╡ 28b77353-7fdd-4587-8455-6aa050e95b76
md"## Lorenz System"

# ╔═╡ ba349aaa-7cf7-4604-be8c-208f8b716e7a
md"Source: [Wikipedia Article on Lorenz System](https://en.wikipedia.org/wiki/Lorenz_system)"

# ╔═╡ 0223e4d3-76d5-4fad-999b-de922a8d1872
md"[Lorenz, Edward Norton (1963). \"Deterministic Nonperiodic Flow\". Journal of the Atmospheric Sciences.](https://journals.ametsoc.org/view/journals/atsc/20/2/1520-0469_1963_020_0130_dnf_2_0_co_2.xml)"

# ╔═╡ 182dbe40-297a-49dc-8d09-ab454f7c10e7
md"[Harris, William (2014). \"How Chaos Theory Works\". Updated 2020. HowStuffWorks.com.](https://science.howstuffworks.com/math-concepts/chaos-theory4.htm)"

# ╔═╡ 6d2c7897-6f30-4de2-982e-d19b0cd1b264
md"### Describe Problem"

# ╔═╡ 48c52331-f17d-4bd2-b10c-a02941ac9582
md"#### System of First Order Ordinary Differential Equations"

# ╔═╡ b04c7bee-1804-4f2b-ab4a-22f658ab5665
md"##### Lorenz System"

# ╔═╡ 8dfc4c1d-3f7f-4d2b-84c6-544dce385651
md"$\begin{align}
\frac{dx}{dt} &= σ(y - x) \\
\frac{dy}{dt} &= x(ρ - z) - y \\
\frac{dz}{dt} &= xy - β z
\end{align}$"

# ╔═╡ d317a1ad-3ee0-424b-95a8-8a6cf283521a
md"\"The equations relate the properties of a two-dimensional fluid layer uniformly warmed from below and cooled from above. In particular, the equations describe the rate of change of three quantities with respect to time.\" -- Wikipedia"

# ╔═╡ 74e56cb1-dd49-4417-9807-c4b31252427d
md"\"In these equations $x$ is proportional to the intensity of the convective motion, while $y$ is proportional to the temperature difference between the ascending and descending currents...The variable $z$ is proportional to the distortion of the vertical temperature profile from linearity...\" -- Lorenz (1963)"

# ╔═╡ 53415544-c11c-4648-8d21-9c9f0290f782
md"\"...where $σ$ represents the ratio of fluid viscosity to thermal conductivity, $ρ$ represents the difference in temperature between the top and bottom of the system and $β$ is the ratio of box width to box height.\" -- Harris (2014, 2020)"

# ╔═╡ bcfd1d61-bd62-4df4-b050-d36c836b18dc
md"##### DifferentialEquations.jl"

# ╔═╡ 635797b6-0113-4398-849a-c2db628b0168
md"$f(du, u, p, t) =
\begin{cases}
du[1] = p[1] (u[2] - u[1]) \\
du[2] = u[1] (p[2] - u[3]) - u[2] \\
du[3] = u[1] u[2] - p[3] u[3]
\end{cases}$"

# ╔═╡ 6488b2d3-0ef0-4701-864e-6a19da585dda
md"### Define Julia Function"

# ╔═╡ ca34cd21-d354-4d24-acbf-7bddd393f3b7
function lorenz(du, u, p, t)
    x, y, z = u
    σ, ρ, β = p
    du[1] = dx = σ * (y - x)
    du[2] = dy = x * (ρ - z) - y
    du[3] = dz = x * y - β * z
end

# ╔═╡ c9229d01-43c9-4643-bd4d-811a65acfade
md"### Assign Variables"

# ╔═╡ 56075c58-3b32-4bc6-a6b5-e48a266d34db
md"#### u: Unknown Functions x, y and z"

# ╔═╡ ae40b3b9-6be5-42bc-a817-ac50273cf86e
md"#### p: Parameters σ, ρ and β"

# ╔═╡ 7c5b8443-a8e5-44e6-a368-d1e350119447
md"#### t: Time"

# ╔═╡ d5835fe1-bccc-49e0-911d-cb503ceaf69f
lt_begin = 0.0

# ╔═╡ cd7be3cd-26b8-4ae4-bd5f-d0b69a452d0e
lt_end = 100.0

# ╔═╡ 46cd7fd2-5a5b-4574-aa86-fe3b593d80d9
ltspan = (lt_begin, lt_end)

# ╔═╡ 8b0c2b96-04f9-4dde-86db-d2f890d968b9
md"### Define System of First Order ODEs Problem"

# ╔═╡ 723f31cf-fde6-410f-bd4b-bbf2b630c69d
md"### Solve Problem"

# ╔═╡ 329e160c-848c-405e-8c82-554a1c595c93
md"### Prep for Plotting"

# ╔═╡ d205c38b-eed2-46de-8c11-465ccc3e9db3
md"#### Break Out $t$"

# ╔═╡ 6d9ff372-9819-4e84-9cf2-93b103744835
md"#### Break Out $u$"

# ╔═╡ 4751657d-6344-4e5d-93ad-49ca358b5753
md"### Plot Solution"

# ╔═╡ 81fcbc6b-46dd-47a0-b0eb-26ebd66498e4
md"#### $u$ with Respect to Time Plots"

# ╔═╡ 38d1a9d7-232f-4a8d-8e78-81d42ab63db2
md"#### 2-Dimensional Phase Space Plots"

# ╔═╡ 69d24f15-4542-472f-a4fd-ff205b7e12a7
md"#### 3-Dimensional Phase Space Plot"

# ╔═╡ c7d6af18-cc74-43ce-998e-2909048d47c0
md"##### Control Panel"

# ╔═╡ bef53fbe-1740-4def-bc11-eab65befe283
@bind reset Button("Reset")

# ╔═╡ a52939bc-9448-4003-9333-ffaf87fef207
begin
    reset
    md"""
    Initial ||
    x: $(@bind lx_begin NumberField(0:1, default = 0)) ||
    y: $(@bind ly_begin NumberField(0:1, default = 1)) ||
    z: $(@bind lz_begin NumberField(0:1, default = 0))
    """
end

# ╔═╡ 3667c4fa-6237-4264-a7dd-11fb76043f8e
lu_begin = [lx_begin, ly_begin, lz_begin]

# ╔═╡ a71d0c7f-d449-4dc3-b3ff-da62d77c43d7
begin
    reset
    md"""
    Parameters ||
    σ: $(@bind lσ Scrubbable(0.0:0.1:20.0, default = 10.0)) |
    ρ: $(@bind lρ Scrubbable(0.0:0.1:50.0, default = 28.0)) |
    β: $(@bind lβ Scrubbable(0.00:0.01:5.00, default = 8/3, format = ".2f"))
    """
end

# ╔═╡ 2c7370ff-0a38-4087-9f63-33b7ab29cc92
lp = [lσ, lρ, lβ]

# ╔═╡ 34a26341-53ed-4d7e-9726-dc1ad573e8b6
lprob = ODEProblem(lorenz, lu_begin, ltspan, lp)

# ╔═╡ 79cda05d-d676-4c51-b0c3-ba7e278b76ef
soll = solve(lprob, reltol=1e-8, abstol=1e-8)

# ╔═╡ 1563f1c1-4160-498d-a38d-a155133670e2
soll.t

# ╔═╡ d5816687-2c6d-4f67-923f-4fbda1f36ba6
n = length(soll.t)

# ╔═╡ 57c12cb6-2e0a-48e0-a9e7-f957a60e6cb5
begin
    reset
    md"Index: $(@bind idx Slider(1:n, n, true))"
end

# ╔═╡ f813b1f3-10eb-4606-8707-aeacff6d734c
t = @view soll.t[1:idx]

# ╔═╡ 8da2a737-8f2e-4b79-a021-32a7537fe5a1
time = t[idx]

# ╔═╡ 1772ba28-c51f-4645-a43d-3840036a63bf
soll.u

# ╔═╡ e599666b-7d11-4010-b6d1-9a13d9efd443
u = @view soll.u[1:idx]

# ╔═╡ 801d509d-e372-472c-8a63-08c88f24e299
soll_u_matrix = reduce(hcat, u)'

# ╔═╡ bf4e99c7-9509-4fe5-a717-2db23bfca5fb
x = soll_u_matrix[:, 1]

# ╔═╡ 1d2ae2b7-7706-4662-95bb-7c2daaa54776
convection = x[idx]

# ╔═╡ b0c7e81e-fcd0-4a4b-9c09-f71020f8824a
begin
    figx = Figure()

    axx = Axis(figx[1, 1],
        title="Lorenz System (Convection)",
        xlabel="t",
        ylabel="Convection (x)")

    lines!(axx, t, x, linewidth=2)

    figx
end

# ╔═╡ 4af0fb0c-5505-4ab6-8215-4af4870c090b
y = soll_u_matrix[:, 2]

# ╔═╡ 69e1c979-9f8b-438e-a7f9-5b180f5d03fa
horizontal = y[idx]

# ╔═╡ 962ebd44-e2dc-44aa-b170-0037864f9e2b
begin
    figy = Figure()

    axy = Axis(figy[1, 1],
        title="Lorenz System (Horizontal)",
        xlabel="t",
        ylabel="Horizontal (y)")

    lines!(axy, t, y, linewidth=2)

    figy
end

# ╔═╡ 8bdc8356-a276-47e9-a0cf-a2285f21ed36
begin
    figxy = Figure()

    axxy = Axis(figxy[1, 1],
        title="Lorenz Attractor (Phase Space x y)",
        xlabel="Convection (x)",
        ylabel="Horizontal (y)")

    lines!(axxy, x, y, linewidth=2)

    figxy
end

# ╔═╡ 9eb53234-6bc4-4757-b8ca-33abc5bebab1
z = soll_u_matrix[:, 3]

# ╔═╡ ba4cd136-3b0f-449f-82f1-39bbae270ff1
vertical = z[idx]

# ╔═╡ 27496fe5-84f5-4e98-b253-4bc9774c7247
begin
    figz = Figure()

    axz = Axis(figz[1, 1],
        title="Lorenz System (Vertical)",
        xlabel="t",
        ylabel="Vertical (z)")

    lines!(axz, t, z, linewidth=2)

    figz
end

# ╔═╡ a33809a9-35ac-40f8-ad52-6fb70600c206
begin
    figxz = Figure()

    axxz = Axis(figxz[1, 1],
        title="Lorenz Attractor (Phase Space x z)",
        xlabel="Convection (x)",
        ylabel="Vertical (z)")

    lines!(axxz, x, z, linewidth=2)

    figxz
end

# ╔═╡ 1be6adca-42f9-4c4a-a684-624022127265
begin
    figyz = Figure()

    axyz = Axis(figyz[1, 1],
        title="Lorenz Attractor (Phase Space y z)",
        xlabel="Horizontal (y)",
        ylabel="Vertical (z)")

    lines!(axyz, y, z, linewidth=2)

    figyz
end

# ╔═╡ 991ada5d-3df6-4571-b000-8af04a5f3629
begin
    figxyz = Figure()

    axxyz = Axis3(figxyz[1, 1],
        title="Lorenz Attractor (Phase Space x y z)\n\n",
        xlabel="Convection (x)",
        ylabel="Horizontal (y)",
        zlabel="Vertical (z)",
        azimuth=-0.2π)

    lines!(axxyz, x, y, z, linewidth=2)
    scatter!(axxyz, convection, horizontal, vertical, color=:red, markersize=2)
    figxyz
end

# ╔═╡ 23ffe24e-95dd-496a-b0ba-fd2420543dc9
begin
    figl = Figure()

    axl = Axis(figl[1, 1],
        title="Lorenz System",
        xlabel="t",
        ylabel="u")

    labels2 = ["Convection (x)" "Horizontal (y)" "Vertical (z)"]

    for (col, label) in enumerate(labels2)
        lines!(axl, t, soll_u_matrix[:, col], linewidth=2, label=label)
    end

    figl
end

# ╔═╡ 0a80db8d-7d77-4fb8-8f1c-a1820939d027
md"""
Time (t): $(round(time; digits = 2)) |
Convection (x): $(round(convection; digits = 2)) |
Horizontal (y): $(round(horizontal; digits = 2)) |
Vertical (z): $(round(vertical; digits = 2))
"""

# ╔═╡ Cell order:
# ╟─8e784e49-e0b4-424d-8fee-1cf29fcfb67f
# ╠═42059750-855e-49e4-87e3-0c2238aff2af
# ╠═0a619df8-1d09-4fc4-bfb3-ea29d07f6e22
# ╠═d0768431-0b3a-47ca-bbe7-9721a6630ddf
# ╠═15a108e1-8f63-464a-8060-5a466e9e7ccb
# ╠═3c700cb2-3194-4ed1-ba03-def708a04867
# ╠═9026f467-8f1c-4aaa-8b42-c3b2843b4949
# ╟─28b77353-7fdd-4587-8455-6aa050e95b76
# ╟─ba349aaa-7cf7-4604-be8c-208f8b716e7a
# ╟─0223e4d3-76d5-4fad-999b-de922a8d1872
# ╟─182dbe40-297a-49dc-8d09-ab454f7c10e7
# ╟─6d2c7897-6f30-4de2-982e-d19b0cd1b264
# ╟─48c52331-f17d-4bd2-b10c-a02941ac9582
# ╟─b04c7bee-1804-4f2b-ab4a-22f658ab5665
# ╟─8dfc4c1d-3f7f-4d2b-84c6-544dce385651
# ╟─d317a1ad-3ee0-424b-95a8-8a6cf283521a
# ╟─74e56cb1-dd49-4417-9807-c4b31252427d
# ╟─53415544-c11c-4648-8d21-9c9f0290f782
# ╟─bcfd1d61-bd62-4df4-b050-d36c836b18dc
# ╟─635797b6-0113-4398-849a-c2db628b0168
# ╟─6488b2d3-0ef0-4701-864e-6a19da585dda
# ╠═ca34cd21-d354-4d24-acbf-7bddd393f3b7
# ╟─c9229d01-43c9-4643-bd4d-811a65acfade
# ╟─56075c58-3b32-4bc6-a6b5-e48a266d34db
# ╠═3667c4fa-6237-4264-a7dd-11fb76043f8e
# ╟─ae40b3b9-6be5-42bc-a817-ac50273cf86e
# ╠═2c7370ff-0a38-4087-9f63-33b7ab29cc92
# ╟─7c5b8443-a8e5-44e6-a368-d1e350119447
# ╠═d5835fe1-bccc-49e0-911d-cb503ceaf69f
# ╠═cd7be3cd-26b8-4ae4-bd5f-d0b69a452d0e
# ╠═46cd7fd2-5a5b-4574-aa86-fe3b593d80d9
# ╟─8b0c2b96-04f9-4dde-86db-d2f890d968b9
# ╠═34a26341-53ed-4d7e-9726-dc1ad573e8b6
# ╟─723f31cf-fde6-410f-bd4b-bbf2b630c69d
# ╠═79cda05d-d676-4c51-b0c3-ba7e278b76ef
# ╟─329e160c-848c-405e-8c82-554a1c595c93
# ╟─d205c38b-eed2-46de-8c11-465ccc3e9db3
# ╠═1563f1c1-4160-498d-a38d-a155133670e2
# ╠═d5816687-2c6d-4f67-923f-4fbda1f36ba6
# ╠═f813b1f3-10eb-4606-8707-aeacff6d734c
# ╠═8da2a737-8f2e-4b79-a021-32a7537fe5a1
# ╟─6d9ff372-9819-4e84-9cf2-93b103744835
# ╠═1772ba28-c51f-4645-a43d-3840036a63bf
# ╠═e599666b-7d11-4010-b6d1-9a13d9efd443
# ╠═801d509d-e372-472c-8a63-08c88f24e299
# ╠═bf4e99c7-9509-4fe5-a717-2db23bfca5fb
# ╠═4af0fb0c-5505-4ab6-8215-4af4870c090b
# ╠═9eb53234-6bc4-4757-b8ca-33abc5bebab1
# ╠═1d2ae2b7-7706-4662-95bb-7c2daaa54776
# ╠═69e1c979-9f8b-438e-a7f9-5b180f5d03fa
# ╠═ba4cd136-3b0f-449f-82f1-39bbae270ff1
# ╟─4751657d-6344-4e5d-93ad-49ca358b5753
# ╟─81fcbc6b-46dd-47a0-b0eb-26ebd66498e4
# ╠═23ffe24e-95dd-496a-b0ba-fd2420543dc9
# ╠═b0c7e81e-fcd0-4a4b-9c09-f71020f8824a
# ╠═962ebd44-e2dc-44aa-b170-0037864f9e2b
# ╠═27496fe5-84f5-4e98-b253-4bc9774c7247
# ╠═38d1a9d7-232f-4a8d-8e78-81d42ab63db2
# ╠═8bdc8356-a276-47e9-a0cf-a2285f21ed36
# ╠═a33809a9-35ac-40f8-ad52-6fb70600c206
# ╠═1be6adca-42f9-4c4a-a684-624022127265
# ╟─69d24f15-4542-472f-a4fd-ff205b7e12a7
# ╠═991ada5d-3df6-4571-b000-8af04a5f3629
# ╟─c7d6af18-cc74-43ce-998e-2909048d47c0
# ╟─bef53fbe-1740-4def-bc11-eab65befe283
# ╟─57c12cb6-2e0a-48e0-a9e7-f957a60e6cb5
# ╟─a52939bc-9448-4003-9333-ffaf87fef207
# ╟─a71d0c7f-d449-4dc3-b3ff-da62d77c43d7
# ╟─0a80db8d-7d77-4fb8-8f1c-a1820939d027
