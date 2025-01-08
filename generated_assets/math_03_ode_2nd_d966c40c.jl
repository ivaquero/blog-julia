### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> chapter = 2
#> section = 3
#> order = 3
#> title = "Second Order ODE"
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

# ╔═╡ aaaa0553-0b07-4fcf-a249-d48fd3a623cf
begin
    using Pkg
    if isdir("../../pluto-deployment-environment")
        Pkg.activate("../../pluto-deployment-environment")
        Pkg.instantiate()
    else
        println(pwd())
    end
end

# ╔═╡ 85838d71-5f65-4e17-8c00-87cf7a1b58ee
using DifferentialEquations

# ╔═╡ a15ac9be-089a-47cd-b076-9f3b32963831
using CairoMakie

# ╔═╡ 08cc349e-aed9-46ff-b5a3-a750682cf56d
using PlutoUI

# ╔═╡ 436151b7-e035-457b-8aa3-903ff582ccd3
using PlutoUI: Button, Slider

# ╔═╡ c19e6b15-bdad-4af9-a507-5aaf94cf29c4
md"# Ordinary Differential Equation - 2nd Order"

# ╔═╡ c64bf6da-952d-45a1-85d5-241f8961d090
md"## Add Packages"

# ╔═╡ 6da01cda-e68b-4a32-9812-7f1ee8c72d47
TableOfContents()

# ╔═╡ 93baba10-4b48-11ed-15e3-57b798e16e05
md"## Projectile Motion"

# ╔═╡ 893d028f-542e-4efd-9be8-9fe4f81757e9
md"### Describe Problem"

# ╔═╡ 38b81041-f663-4428-bb05-d85e3027a80b
md"#### 1st Order Ordinary Differential Equation"

# ╔═╡ 3c5c3bcf-4c78-490f-a6ae-ea5557138ca3
md"#### Relationship between Position, Velocity and Acceleration"

# ╔═╡ 0baae2d0-0199-4edf-b284-31755f90a5f3
md"$u = position$"

# ╔═╡ 105266a1-9ec6-47c0-a279-5be812bafb15
md"First Order Derivative of Position with respect to Time:"

# ╔═╡ 87026eaf-167a-47c7-b895-056154ad49ed
md"$\frac{du}{dt} = velocity$"

# ╔═╡ 9e667814-7977-47bf-95d1-bd27a887693e
md"First Order Derivative of Velocity with respect to Time:"

# ╔═╡ 3ff45669-557a-4f8d-b52b-81d8dd2b9a03
md"$\frac{dv}{dt} = acceleration$"

# ╔═╡ b1b768d1-e421-4d0d-8a89-ef4564d8fa52
md"Example: Meters per Second per Second:"

# ╔═╡ f81b56aa-ab26-4c74-bcb4-06715d27cb62
md"$\frac{m}{s^2}$"

# ╔═╡ 94dbdf02-efca-45e2-abc6-dc9a20659b17
md"#### 2nd Order Ordinary Differential Equation"

# ╔═╡ 3e09e75b-2701-43c9-84e0-1ed393ad9ab1
md"Second Derivative of Position with respect to Time:"

# ╔═╡ 4c0a2697-0b91-4faa-a1d6-9d8eba804c4e
md"$\frac{d^2u}{dt^2} = acceleration$"

# ╔═╡ 61137347-8c11-48d9-b44f-c578cc25881b
md"#### DifferentialEquations.jl"

# ╔═╡ cfed74e5-5035-45c8-8062-978bb25e263c
md"$\frac{d^2u}{dt^2} = f(du, u, p, t)$"

# ╔═╡ 3c68f18b-8344-4879-a173-eab036675c79
md"where $u$ is the unknown function and $du$ is the first derivative of $u$ with respect to $t$."

# ╔═╡ 0e76740a-be26-47d0-a7ab-c2688ddf18bf
md"##### Projectile Motion Problem"

# ╔═╡ adeedcef-b65a-46d6-82d0-4347b2a0843b
md"$f(du, u, p, t) = [0, g]$"

# ╔═╡ bd61576c-bca7-49b1-b30e-353e2af43eb0
md"where $g$ is the frictionless, free fall gravitational acceleration."

# ╔═╡ 00509e7f-c194-4e26-ae41-d80e5f9337e0
md"### Define Function"

# ╔═╡ e047d82a-97af-4d12-a5fe-0bd2878f8c0a
md"#### Gravitational Accelerations in $\frac{m}{s^2}$"

# ╔═╡ 2c60a984-9360-4cf6-b2fa-6de06068921c
md"Source: [Wikipedia Article](https://en.wikipedia.org/wiki/Gravitational_acceleration)"

# ╔═╡ 23e6b417-ff94-47da-897c-30405f4d43e1
const g_EARTH = -9.8067

# ╔═╡ 54307131-bc75-4bd2-8dc6-fbefa826e3a9
const g_MOON = -1.625

# ╔═╡ c7611782-80ba-4ace-847b-ae97f2df8fa7
const g_MARS = -3.728

# ╔═╡ 2a6d3c67-ce18-4763-94ac-1ffdb52185d3
md"### Assign Variables"

# ╔═╡ 090bd835-8321-4562-8ead-3b836ec2187c
md"#### du: Velocity in Meters per Second"

# ╔═╡ 7be99dc3-a454-4ee2-8781-17bb798ed3ef
md"$\cos θ = \frac{\text{adjacent}}{\text{hypotenuse}} = \frac{vx}{v}$"

# ╔═╡ 535c3260-0f85-461b-851c-e67aa08e5d7b
md"$\sin θ = \frac{\text{opposite}}{\text{hypotenuse}} = \frac{vy}{v}$"

# ╔═╡ a301fb45-1d35-4db5-8522-aadce039f554
md"#### u: Position in Meters"

# ╔═╡ 558df6c3-f439-4e80-9dc4-8fa14248da41
x_begin = 0.0

# ╔═╡ 107e3db6-5ebb-41db-b74f-41d5b6adb03b
y_begin = 0.0

# ╔═╡ 190f00a9-8f36-4f3c-9592-af8b2dbc672f
u_begin = [x_begin, y_begin]

# ╔═╡ d6b66898-cdd4-4a70-93bf-b74bfd51f87e
md"#### t: Time in Seconds"

# ╔═╡ 122443c1-6751-4717-ac35-9262ca2e5583
t_begin = 0.0

# ╔═╡ 9c204420-7d1c-457e-bfbd-4465a2e1ea66
md"### Define ODE Problem"

# ╔═╡ 864d0a33-1ebf-4855-a686-bd0dcdc5e7bf
md"### Solve Problem"

# ╔═╡ 3661663f-8f66-4e13-a990-69f129c432bc
md"### Plot Solution"

# ╔═╡ 91d831c2-125f-49b4-916e-46bb9385b9ec
@bind select_g Select(["Earth", "Moon", "Mars"])

# ╔═╡ 53b7e346-9b05-4470-9b11-02b10c98a0d8
if select_g == "Moon"
    g = g_MOON
elseif select_g == "Mars"
    g = g_MARS
else
    g = g_EARTH
end

# ╔═╡ 45365e48-fba2-4d3a-8d01-4f4f3637d0ac
f(du, u, p, t) = [0, g]

# ╔═╡ 4ce7a287-79e1-44e5-8457-cd7c2c1873e0
md"Seconds: $(@bind t_end Slider(0.0:0.1:20.0, 20.0, true))"

# ╔═╡ 5b3db5f4-47e8-485d-a157-bf8d1d2bb8f8
tspan = (t_begin, t_end)

# ╔═╡ b78d2f18-5826-4ac4-9b1e-2a9f28e8b342
md"Initial Velocity (m/s): $(@bind v_begin Slider(50.0:150.0, 100.0, true))"

# ╔═╡ 77e4774b-e677-4483-8a36-9ca71f14abad
md"Angle in Degrees: $(@bind θ Slider(0.0:90.0, 45.0, true))"

# ╔═╡ 4889b9c3-9e88-4b2a-8537-ddb415a39fe6
vx_begin = cosd(θ) * v_begin

# ╔═╡ 71225132-d38b-4293-affa-f9ccc8b4d995
vy_begin = sind(θ) * v_begin

# ╔═╡ 0ac7ef2e-f6c8-4842-b159-950f2b51124d
du_begin = [vx_begin, vy_begin]

# ╔═╡ 93548b0e-ff09-482b-ac4a-3ce7ef8c1c43
prob = SecondOrderODEProblem(f, du_begin, u_begin, tspan)

# ╔═╡ fec0a303-a585-4459-9390-774c62ea777e
sol = solve(prob)

# ╔═╡ e504738e-e820-412d-8358-795b3485ce55
sol(20)[1]

# ╔═╡ 79428bd2-64dc-4a2c-b13b-293033c02e32
sol(20)[4]

# ╔═╡ be1c2e2a-815c-4f47-885f-8d5b385fe24a
md"""
x Position (m): $(round(sol(t_end)[3]; digits = 2)) |
y Position (m): $(round(sol(t_end)[4]; digits = 2)) |
x Velocity (m/s): $(round(sol(t_end)[1]; digits = 2)) |
y Velocity (m/s): $(round(sol(t_end)[2]; digits = 2))
"""

# ╔═╡ 7a472404-d82d-41fe-b2d8-6812065a11cf
begin
    fig = Figure()

    ax = Axis(fig[1, 1],
        title="Projectile Motion (All Values)",
        xlabel="Time in Seconds",
        ylabel="Velocity (m/s) | Position (m)")

    labels = ["x Velocity", "y Velocity", "x Position", "y Position"]

    for (col, label) in enumerate(labels)
        lines!(ax, sol.t, sol[col, :], linewidth=2, label=label)
    end

    axislegend()
    fig
end

# ╔═╡ c75c69ff-221a-4aca-a583-c1511f6cb265
begin
    fig2 = Figure()

    ax2 = Axis(fig2[1, 1],
        title="Projectile Motion (Trajectory)",
        xlabel="x Position (m)",
        ylabel="y Position (m)")

    lines!(ax2, sol[3, :], sol[4, :], linewidth=2)

    fig2
end

# ╔═╡ c9482118-4c14-4a27-a7e4-ed3d2b51e0c5
md"## Simple Harmonic Motion"

# ╔═╡ fd310adf-b672-48e2-9418-9d84481f2a06
md"Source: [Wikipedia Article](https://en.wikipedia.org/wiki/Simple_harmonic_motion)"

# ╔═╡ 395cd47b-531c-4521-824d-0399e40c08b5
md"### Describe Problem"

# ╔═╡ db7f8b0d-8abd-46bb-ad36-1304428d4756
md"#### Second Order Ordinary Differential Equation"

# ╔═╡ 90780034-995d-4abb-992f-cf5e408ce7b1
md"From Hooke's Law:"

# ╔═╡ 62c348ed-c0fc-45a5-ae2b-7578f92d54c5
md"$F_s = -kx$"

# ╔═╡ 350c25c6-3ea2-4191-ad89-a4412e02195c
md"where $F_s$ is the restoring force of the spring, $k$ is a positive number, characteristic of the spring and $x$ is the displacement from the equilibrium (or mean) position."

# ╔═╡ 381d7225-c51b-4666-aced-a2de9df42841
md"From Newton's Second Law of Motion:"

# ╔═╡ 5dce232a-eded-4f2d-b033-e73130142ef4
md"$F = ma$"

# ╔═╡ 58c69620-5430-4e74-88dc-86a81afa1e0c
md"where $m$ is the mass of the object and $a$ is the acceleration."

# ╔═╡ d36463ca-117c-4c5e-bf8d-1fc028fed64d
md"From the last tutorial:"

# ╔═╡ 68f87501-5fdc-475c-86b1-7f591b8e13b2
md"$a = \frac{d^2x}{dt^2}$"

# ╔═╡ 9a16de72-6899-4c1d-9f18-0e095c9e1d84
md"where the acceleration, $a$, is the second derivative of position, $x$, with respect to time."

# ╔═╡ 0ec2381b-77f8-450a-818d-851ee7380354
md"Combine equations:"

# ╔═╡ f310cc54-5e8b-4731-9772-70f5166c9310
md"$F = m \frac{d^2x}{dt^2} = -kx$"

# ╔═╡ 922869e4-ba94-45ad-a1dd-61a719cf9e6b
md"$\frac{d^2x}{dt^2} = - \frac{k}{m} x$"

# ╔═╡ edaf6150-cb0e-4db3-bdac-c06c9016aed3
md"\"where $m$ is the inertial mass of the oscillating body, $x$ is the displacement from the equilibrium (or mean) position and $k$ is a constant (the spring constant for a mass on a spring).\" -- Wikipedia"

# ╔═╡ 23ad3a7a-3230-4b1a-8130-675f9c5538a6
md"#### DifferentialEquations.jl"

# ╔═╡ c5ed2565-234f-4b94-b05c-b53d80b07c47
md"$\frac{d^2u}{dt^2} = f(du, u, p, t)$"

# ╔═╡ 27216dd3-e8ba-4ffc-9014-e20624598af7
md"#### Simple Harmonic Motion Problem"

# ╔═╡ 041059c9-dcfb-4b89-ab69-9bf965c2e012
md"$f(du, u, p, t) = - \frac{p[1]}{p[2]} u$"

# ╔═╡ 1b12edd0-2015-49c9-b3f9-841919ef38c3
md"### Define Julia Function"

# ╔═╡ 12684fb2-b778-42a8-923a-616691b03f93
function harmonic(du, u, p, t)
    dx = du
    x = u
    k, m = p
    -(k / m) * x
end

# ╔═╡ 9f793759-296a-4d30-85d4-9a7a910136f6
md"### Assign Variables"

# ╔═╡ d147328a-f41d-40ec-9786-29b85bf770e4
md"#### du: Velocity in Meters per Second"

# ╔═╡ 2ef4d36c-1d0f-4378-8646-9937121f9a2c
dx_begin = 0.0

# ╔═╡ 62681afb-e3b5-4d39-8ef2-117171bbab09
md"#### u: Position in Meters"

# ╔═╡ eefd97bb-f73f-4793-afa1-196c26c01430
md"#### p: Parameters"

# ╔═╡ 5c91e7ba-049e-4434-b5eb-b50fa3ad3d14
md"#### k: Spring Constant in Newtons per Meter"

# ╔═╡ 8ad8b31f-5e5d-400e-ab95-09e90b76a089
md"#### m: Mass in Kilograms"

# ╔═╡ d694b1ed-80fa-429d-b594-8ca47e8662b4
md"#### t: Time in Seconds"

# ╔═╡ e372c42d-0487-400c-929c-48f86bb81228
md"### Define Second Order ODE Problem"

# ╔═╡ 5905065f-9b3f-4a4a-86b6-43f91b100ece
md"### Solve Problem"

# ╔═╡ c8315485-5835-436b-bb72-aa3a003fb51d
md"### Plot Solution"

# ╔═╡ 6e661aed-d69c-4178-a1a4-bf6663729f42
@bind reset Button("Reset Sliders")

# ╔═╡ 2ba83fcc-ea8f-4962-9365-cac38790a9ff
begin
    reset
    md"""
    Time (seconds): $(@bind tm_end Slider(0.0:0.01:20.0, 20.0, true)) |
    Initial Position (m): $(@bind xx_begin Slider(-0.25:0.01:0.25, 0.25, true))
    """
end

# ╔═╡ ee5bdd91-16a0-462c-ad19-90fbfe9a6298
tmspan = (t_begin, tm_end)

# ╔═╡ 47a4fdc2-cd81-44c7-be6f-1c688cdf4df1
begin
    reset
    md"""
    Spring Constant (N/m): $(@bind k Slider(0.5:0.01:2.0, 1.0, true)) |
    Mass (kg): $(@bind m Slider(0.5:0.01:2.0, 1.0, true))
    """
end

# ╔═╡ 84ddfbe3-11bb-4f39-8ef1-bd7c1382611b
p = [k, m]

# ╔═╡ 979beae9-a4c4-41d3-ac87-22ab7de3c95f
probm = SecondOrderODEProblem(harmonic, dx_begin, xx_begin, tmspan, p)

# ╔═╡ 96d45b34-e92e-4684-9d9c-b47ce33961ba
solm = solve(probm)

# ╔═╡ 0577f6c0-24fa-4a24-8417-427e87d23a17
position = solm(tm_end)[2]

# ╔═╡ 6a176b28-86e4-4c6e-9120-f428bdf729db
velocity = solm(tm_end)[1]

# ╔═╡ 286a45e7-615d-43f2-b9d1-89ce9e1cf9b0
acceleration = -k / m * position

# ╔═╡ 99b4bef8-85be-4380-9680-7af0db0cb4f9
begin
    fig3 = Figure()

    ax3 = Axis(fig3[1, 1],
        title="Simple Harmonic Motion (All Values)",
        xlabel="Time (seconds)",
        ylabel="Velocity (m/s) | Position (m)")

    labels2 = ["Position", "Velocity"]

    for (col, label) in enumerate(labels2)
        lines!(ax3, solm.t, solm[col, :], linewidth=2, label=label)
    end

    fig3
end

# ╔═╡ d1df8bd7-6069-4081-8eb5-d33b1372a385
begin
    fig4 = Figure()

    ax4 = Axis(fig4[1, 1],
        title="Simple Harmonic Motion (Phase Space)",
        xlabel="Position (m)",
        ylabel="Velocity (m/s)")

    lines!(ax4, solm[2, :], solm[1, :], linewidth=2)
    scatter!(ax4, position, velocity, color=:red, markersize=10)
    xlims!(ax4, -0.5, 0.5)
    ylims!(ax4, -0.5, 0.5)

    fig4
end

# ╔═╡ 33f9126d-1fd0-4b44-88cd-7e29dc5f6116
md"### Appendix: Period (T)"

# ╔═╡ ef292aaf-2b46-48bf-83ae-d7de78b9729b
md"$T = 2 π \sqrt{\frac{m}{k}}$"

# ╔═╡ d3fcc368-cf6d-4596-bff9-53e48ee28185
period = 2 * pi * sqrt(m / k)

# ╔═╡ c992c858-455c-40b5-ad8d-f6df216b6087
md"""
Position (m): $(round(position; digits = 2)) |
Velocity (m/s): $(round(velocity; digits = 2)) |
Acceleration (m/s/s): $(round(acceleration; digits = 2)) |
Period (s): $(round(period; digits = 2))
"""

# ╔═╡ 5faad3c0-e91b-4a2e-aa8b-5712c6efe67f
md"""
Position (m): $(round(position; digits = 2)) |
Velocity (m/s): $(round(velocity; digits = 2)) |
Acceleration (m/s/s): $(round(acceleration; digits = 2)) |
Period (s): $(round(period; digits = 2))
"""

# ╔═╡ Cell order:
# ╟─c19e6b15-bdad-4af9-a507-5aaf94cf29c4
# ╟─c64bf6da-952d-45a1-85d5-241f8961d090
# ╠═aaaa0553-0b07-4fcf-a249-d48fd3a623cf
# ╠═85838d71-5f65-4e17-8c00-87cf7a1b58ee
# ╠═a15ac9be-089a-47cd-b076-9f3b32963831
# ╠═08cc349e-aed9-46ff-b5a3-a750682cf56d
# ╠═436151b7-e035-457b-8aa3-903ff582ccd3
# ╠═6da01cda-e68b-4a32-9812-7f1ee8c72d47
# ╟─93baba10-4b48-11ed-15e3-57b798e16e05
# ╟─893d028f-542e-4efd-9be8-9fe4f81757e9
# ╟─38b81041-f663-4428-bb05-d85e3027a80b
# ╟─3c5c3bcf-4c78-490f-a6ae-ea5557138ca3
# ╟─0baae2d0-0199-4edf-b284-31755f90a5f3
# ╟─105266a1-9ec6-47c0-a279-5be812bafb15
# ╟─87026eaf-167a-47c7-b895-056154ad49ed
# ╟─9e667814-7977-47bf-95d1-bd27a887693e
# ╟─3ff45669-557a-4f8d-b52b-81d8dd2b9a03
# ╟─b1b768d1-e421-4d0d-8a89-ef4564d8fa52
# ╟─f81b56aa-ab26-4c74-bcb4-06715d27cb62
# ╟─94dbdf02-efca-45e2-abc6-dc9a20659b17
# ╟─3e09e75b-2701-43c9-84e0-1ed393ad9ab1
# ╟─4c0a2697-0b91-4faa-a1d6-9d8eba804c4e
# ╟─61137347-8c11-48d9-b44f-c578cc25881b
# ╟─cfed74e5-5035-45c8-8062-978bb25e263c
# ╟─3c68f18b-8344-4879-a173-eab036675c79
# ╟─0e76740a-be26-47d0-a7ab-c2688ddf18bf
# ╟─adeedcef-b65a-46d6-82d0-4347b2a0843b
# ╟─bd61576c-bca7-49b1-b30e-353e2af43eb0
# ╟─00509e7f-c194-4e26-ae41-d80e5f9337e0
# ╠═45365e48-fba2-4d3a-8d01-4f4f3637d0ac
# ╟─e047d82a-97af-4d12-a5fe-0bd2878f8c0a
# ╟─2c60a984-9360-4cf6-b2fa-6de06068921c
# ╠═23e6b417-ff94-47da-897c-30405f4d43e1
# ╠═54307131-bc75-4bd2-8dc6-fbefa826e3a9
# ╠═c7611782-80ba-4ace-847b-ae97f2df8fa7
# ╠═53b7e346-9b05-4470-9b11-02b10c98a0d8
# ╟─2a6d3c67-ce18-4763-94ac-1ffdb52185d3
# ╟─090bd835-8321-4562-8ead-3b836ec2187c
# ╟─7be99dc3-a454-4ee2-8781-17bb798ed3ef
# ╠═4889b9c3-9e88-4b2a-8537-ddb415a39fe6
# ╟─535c3260-0f85-461b-851c-e67aa08e5d7b
# ╠═71225132-d38b-4293-affa-f9ccc8b4d995
# ╠═0ac7ef2e-f6c8-4842-b159-950f2b51124d
# ╟─a301fb45-1d35-4db5-8522-aadce039f554
# ╠═558df6c3-f439-4e80-9dc4-8fa14248da41
# ╠═107e3db6-5ebb-41db-b74f-41d5b6adb03b
# ╠═190f00a9-8f36-4f3c-9592-af8b2dbc672f
# ╟─d6b66898-cdd4-4a70-93bf-b74bfd51f87e
# ╠═122443c1-6751-4717-ac35-9262ca2e5583
# ╠═5b3db5f4-47e8-485d-a157-bf8d1d2bb8f8
# ╟─9c204420-7d1c-457e-bfbd-4465a2e1ea66
# ╠═93548b0e-ff09-482b-ac4a-3ce7ef8c1c43
# ╟─864d0a33-1ebf-4855-a686-bd0dcdc5e7bf
# ╠═fec0a303-a585-4459-9390-774c62ea777e
# ╠═e504738e-e820-412d-8358-795b3485ce55
# ╠═79428bd2-64dc-4a2c-b13b-293033c02e32
# ╟─3661663f-8f66-4e13-a990-69f129c432bc
# ╠═91d831c2-125f-49b4-916e-46bb9385b9ec
# ╟─4ce7a287-79e1-44e5-8457-cd7c2c1873e0
# ╟─b78d2f18-5826-4ac4-9b1e-2a9f28e8b342
# ╟─77e4774b-e677-4483-8a36-9ca71f14abad
# ╟─be1c2e2a-815c-4f47-885f-8d5b385fe24a
# ╠═7a472404-d82d-41fe-b2d8-6812065a11cf
# ╠═c75c69ff-221a-4aca-a583-c1511f6cb265
# ╟─c9482118-4c14-4a27-a7e4-ed3d2b51e0c5
# ╟─fd310adf-b672-48e2-9418-9d84481f2a06
# ╟─395cd47b-531c-4521-824d-0399e40c08b5
# ╟─db7f8b0d-8abd-46bb-ad36-1304428d4756
# ╟─90780034-995d-4abb-992f-cf5e408ce7b1
# ╟─62c348ed-c0fc-45a5-ae2b-7578f92d54c5
# ╟─350c25c6-3ea2-4191-ad89-a4412e02195c
# ╟─381d7225-c51b-4666-aced-a2de9df42841
# ╟─5dce232a-eded-4f2d-b033-e73130142ef4
# ╟─58c69620-5430-4e74-88dc-86a81afa1e0c
# ╟─d36463ca-117c-4c5e-bf8d-1fc028fed64d
# ╟─68f87501-5fdc-475c-86b1-7f591b8e13b2
# ╟─9a16de72-6899-4c1d-9f18-0e095c9e1d84
# ╟─0ec2381b-77f8-450a-818d-851ee7380354
# ╟─f310cc54-5e8b-4731-9772-70f5166c9310
# ╟─922869e4-ba94-45ad-a1dd-61a719cf9e6b
# ╟─edaf6150-cb0e-4db3-bdac-c06c9016aed3
# ╟─23ad3a7a-3230-4b1a-8130-675f9c5538a6
# ╟─c5ed2565-234f-4b94-b05c-b53d80b07c47
# ╟─27216dd3-e8ba-4ffc-9014-e20624598af7
# ╟─041059c9-dcfb-4b89-ab69-9bf965c2e012
# ╟─1b12edd0-2015-49c9-b3f9-841919ef38c3
# ╠═12684fb2-b778-42a8-923a-616691b03f93
# ╟─9f793759-296a-4d30-85d4-9a7a910136f6
# ╟─d147328a-f41d-40ec-9786-29b85bf770e4
# ╠═2ef4d36c-1d0f-4378-8646-9937121f9a2c
# ╟─62681afb-e3b5-4d39-8ef2-117171bbab09
# ╟─eefd97bb-f73f-4793-afa1-196c26c01430
# ╟─5c91e7ba-049e-4434-b5eb-b50fa3ad3d14
# ╟─8ad8b31f-5e5d-400e-ab95-09e90b76a089
# ╠═84ddfbe3-11bb-4f39-8ef1-bd7c1382611b
# ╟─d694b1ed-80fa-429d-b594-8ca47e8662b4
# ╠═ee5bdd91-16a0-462c-ad19-90fbfe9a6298
# ╟─e372c42d-0487-400c-929c-48f86bb81228
# ╠═979beae9-a4c4-41d3-ac87-22ab7de3c95f
# ╟─5905065f-9b3f-4a4a-86b6-43f91b100ece
# ╠═96d45b34-e92e-4684-9d9c-b47ce33961ba
# ╠═0577f6c0-24fa-4a24-8417-427e87d23a17
# ╠═6a176b28-86e4-4c6e-9120-f428bdf729db
# ╠═286a45e7-615d-43f2-b9d1-89ce9e1cf9b0
# ╟─c8315485-5835-436b-bb72-aa3a003fb51d
# ╟─c992c858-455c-40b5-ad8d-f6df216b6087
# ╟─6e661aed-d69c-4178-a1a4-bf6663729f42
# ╟─2ba83fcc-ea8f-4962-9365-cac38790a9ff
# ╟─47a4fdc2-cd81-44c7-be6f-1c688cdf4df1
# ╟─5faad3c0-e91b-4a2e-aa8b-5712c6efe67f
# ╠═99b4bef8-85be-4380-9680-7af0db0cb4f9
# ╠═d1df8bd7-6069-4081-8eb5-d33b1372a385
# ╟─33f9126d-1fd0-4b44-88cd-7e29dc5f6116
# ╟─ef292aaf-2b46-48bf-83ae-d7de78b9729b
# ╠═d3fcc368-cf6d-4596-bff9-53e48ee28185
