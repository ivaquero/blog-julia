### A Pluto.jl notebook ###
# v0.20.3

#> [frontmatter]
#> chapter = 2
#> section = 3
#> order = 3
#> title = "ODE Systems"
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

# ╔═╡ 835e13a7-badc-4bb4-92ee-61afa3972a71
begin
    using Pkg
    if isdir("../../pluto-deployment-environment")
        Pkg.activate("../../pluto-deployment-environment")
        Pkg.instantiate()
    else
        println(pwd())
    end
end

# ╔═╡ 07582e60-4d07-49b8-a6ab-e962eec8d81b
using DifferentialEquations

# ╔═╡ 1e7ae697-d85a-4c1c-b0da-e863f7dc3e74
using CairoMakie

# ╔═╡ 3b7e547f-3c75-4624-a79b-d9ea6ceb0ead
using PlutoUI

# ╔═╡ 4bc6152e-3201-4678-a876-e5577984cd57
using PlutoUI: Button, Slider

# ╔═╡ ee40337c-6e10-42a6-8aba-8a69a1ace228
using Colors

# ╔═╡ 6eb6bd53-f864-4a82-8173-f533ea7beeaa
md"# ODE Systems"

# ╔═╡ cedd2cac-46f9-4875-afe1-b1b3a936135d
md"## Add Packages"

# ╔═╡ d3ed21d2-d3c2-416d-84e4-560fdee7fb0a
TableOfContents()

# ╔═╡ 9ee9c82c-7d06-4e0c-a353-c1c2e86f71a3
md"## Double Pendulum"

# ╔═╡ ba9eb20a-8c8f-455b-9d4d-f44601b4fb44
md"Sources:"

# ╔═╡ cfb79dad-db47-4a66-9d23-cd1cb78e16c5
md"[Wikipedia Article on Double Pendulum](https://en.wikipedia.org/wiki/Double_pendulum)"

# ╔═╡ 7aa78e57-7b0d-43b2-b0a3-da4ba3a9d70c
md"[Neumann, Erik (2002). Double Pendulum Simulation at myphysicslab.com.](https://www.myphysicslab.com/pendulum/double-pendulum-en.html)"

# ╔═╡ 445b7e01-3867-4814-9237-999ec6600e39
md"[Ma, Yingbo; Rackauckas, Chris (2020). \"Differential Equations Models\".](https://nextjournal.com/sosiris-de/diffeq-models)"

# ╔═╡ 60969d34-b75d-4862-93bf-46addf6be92d
md"### Describe Problem"

# ╔═╡ 4b6da8ef-c845-44b4-b980-976671240d26
md"#### System of 1st Order ODEs"

# ╔═╡ 04b935a7-6c6d-4fe1-ba48-149de8e64046
md"Note: see Numerical Solution Equations for Double Pendulum Simulation at myphysicslab.com."

# ╔═╡ cc092311-1167-4cfa-aef7-4f63beb62b40
md"Double Pendulum Problem"

# ╔═╡ 36a0f007-6071-4224-a6ab-6cbcecf5c9ef
md"$f(du, u, p, t) =
\begin{cases}
du[1] = {\theta_1}' = \omega_1 \\
du[2] = {\theta_2}' = \omega_2 \\
du[3] = {\omega_1}' \\
du[4] = {\omega_2}'
\end{cases}$"

# ╔═╡ b0195989-346a-41fd-8f7b-a9bf87424b29
md"Define Julia Function"

# ╔═╡ d8db785d-cdf4-4b94-91f3-5bbf9496878f
function double_pendulum(du, u, p, t)
    theta1, theta2, omega1, omega2 = u

    m1, m2, L1, L2, g = p

    c = cos(theta1 - theta2)
    s = sin(theta1 - theta2)

    du[1] = omega1
    du[2] = omega2

    du[3] = (-g * (2 * m1 + m2) * sin(theta1) -
             m2 * g * sin(theta1 - 2 * theta2) -
             2 * s * m2 * (omega2^2 * L2 + omega1^2 * L1 * c)) /
            (L1 * (2 * m1 + m2 - m2 * cos(2 * theta1 - 2 * theta2)))

    du[4] = (2 * s * (omega1^2 * L1 * (m1 + m2) +
                      g * (m1 + m2) * cos(theta1) +
                      omega2^2 * L2 * m2 * c)) /
            (L2 * (2 * m1 + m2 - m2 * cos(2 * theta1 - 2 * theta2)))
end

# ╔═╡ 6211c0cd-c282-450c-816d-921da466e267
md"### Assign Variables"

# ╔═╡ f328cad6-72c8-46b4-8368-5ec386a5fc6e
md"#### u: θ1, θ2, ω1 and ω2"

# ╔═╡ 3aa7e90e-e889-4309-8814-5d0b8406b288
ω1 = 0.0

# ╔═╡ b1dd866e-044a-47e4-b2e9-e9a540dc59ae
ω2 = 0.0

# ╔═╡ b5f1925c-1eaa-479f-afd3-583eddc84007
md"#### p: Parameters m1, m2, L1, L2 and g"

# ╔═╡ 1fa233ba-873e-4066-91f0-2375651733f3
md"#### t: Time in Seconds"

# ╔═╡ 06c64392-0e4c-4afa-8181-498aab2d586a
t_begin = 0.0

# ╔═╡ 5b638093-fefa-418f-99c4-2588448da696
md"### Solve Problem"

# ╔═╡ d37a967d-64c2-46ad-9fd6-0169be8c2b50
md"### Prep for Plotting"

# ╔═╡ 550f5e88-7aaf-44a7-8ea9-243e3d3d6600
md"### Plot Solution"

# ╔═╡ 367fea2b-84e6-42da-addd-040a64345cc2
md"##### Control Panel"

# ╔═╡ 77d77b90-1726-46d2-9cfd-d85115c14290
@bind preset Button("Reset")

# ╔═╡ ecff4935-d3cd-4ae4-b92a-00c172971b28
begin
    preset
    md"Time (s): $(@bind t_end Slider(0.0:0.01:100.0, 100.0, true))"
end

# ╔═╡ c6ced053-87a0-4608-91d6-bc3ce8c5b574
tspan = (t_begin, t_end)

# ╔═╡ 0880ed90-9970-4606-b45a-d9b58cd1f739
begin
    preset
    md"""
    Initial Conditions ||
    θ1 (deg): $(@bind θ1d Scrubbable(0.0:0.1:180.0, default = 90.0)) |
    θ2 (deg): $(@bind θ2d Scrubbable(0.0:0.1:180.0, default = 120.0))
    """
end

# ╔═╡ 9128427c-ecc1-450a-85ce-fdd37c45305f
θ1 = deg2rad(θ1d)

# ╔═╡ 9f41835f-2cdf-4fda-b955-e4666d762ea1
θ2 = deg2rad(θ2d)

# ╔═╡ 536da0f6-c2c5-4987-9e3c-2124825be7a4
u_begin = [θ1, θ2, ω1, ω2]

# ╔═╡ d8fea5d8-4b66-4593-b994-b5d92f30a9ad
begin
    reset
    md"""
    Parameters ||
    m1 (kg): $(@bind m1 Scrubbable(0.1:0.1:2.0, default = 1.0)) |
    m2 (kg): $(@bind m2 Scrubbable(0.1:0.1:2.0, default = 1.5)) |
    L1 (m): $(@bind L1 Scrubbable(0.1:0.1:2.0, default = 1.5)) |
    L2 (m): $(@bind L2 Scrubbable(0.1:0.1:2.0, default = 2.0)) |
    g (m/s/s): $(@bind g Scrubbable(9.8:0.0001:9.81, default = 9.8067,
    	format = ".4f"))
    """
end

# ╔═╡ bd285b60-e43f-4754-9446-2ba293626819
pp = [m1, m2, L1, L2, g]

# ╔═╡ 078bb957-3d1d-477f-9e64-cd94f7d01e6b
prob = ODEProblem(double_pendulum, u_begin, tspan, pp)

# ╔═╡ 18471ae5-52f1-4549-86e0-ad1f00d4bcdb
sol = solve(prob, Vern7(), reltol=1e-10, abstol=1e-10)

# ╔═╡ 6cdd8d51-2c4e-4ef0-9d5b-5cd982f7aa6f
begin
    fig = Figure()

    ax = Axis(fig[1, 1],
        xlabel="Time")

    for col in range(1, 4)
        lines!(ax, sol.t, sol[col, :], linewidth=2)
    end

    fig
end

# ╔═╡ 2ebe56db-8d9a-488a-831e-7e8a5cbd5548
function polar2cart(sol; vars=(1, 2))
    idx = sol.t[1]:0.01:sol.t[end]

    p1 = map(x -> x[vars[1]], sol.(idx))
    p2 = map(y -> y[vars[2]], sol.(idx))

    x1 = L1 * sin.(p1)
    y1 = -L1 * cos.(p1)

    x2 = L2 * sin.(p2)
    y2 = -L2 * cos.(p2)

    ((x1, y1), (x1 + x2, y1 + y2))
end

# ╔═╡ d0511ed6-7cb1-4b1d-a246-bfeea8d25a9c
pendulum1, pendulum2 = polar2cart(sol)

# ╔═╡ 0cbcc92b-c817-47cb-bd5d-24244b78cc70
x1 = pendulum1[1][end]

# ╔═╡ c177c824-04f4-49e0-be7f-aec181595632
y1 = pendulum1[2][end]

# ╔═╡ 9b8949f6-acfd-4668-ae07-a51161cbf789
x2 = pendulum2[1][end]

# ╔═╡ f2d2c387-5f23-47d6-9d8b-3d56f7d5ce8a
y2 = pendulum2[2][end]

# ╔═╡ 624eb883-9ec3-4c02-8824-dfa9a0099eea
begin
    figx1y1 = Figure()

    axx1y1 = Axis(figx1y1[1, 1],
        title="Double Pendulum Phase Space x2 y2",
        xlabel="x",
        ylabel="y")

    lines!(axx1y1, pendulum1[1], pendulum1[2], linewidth=2)

    figx1y1
end

# ╔═╡ ee43de08-3ec1-4517-98b6-456b2db9d2b3
begin
    figx2y2 = Figure()

    axx2y2 = Axis(figx2y2[1, 1],
        title="Double Pendulum Phase Space x2 y2",
        xlabel="x",
        ylabel="y")

    lines!(axx2y2, pendulum2[1], pendulum2[2], linewidth=2)

    lines!(axx2y2, [0, x1, x2], [0, y1, y2], color=:black, linewidth=2)

    for (pair, color) in zip(zip([0, x1, x2], [0, y1, y2]), [:black, :green, :red])
        scatter!(axx2y2, pair[1], pair[2], color=color, markersize=20)
    end

    figx2y2
end

# ╔═╡ 4a7ff5e7-e08e-4788-9dba-00e3bd2005a2
md"""
Time: $(t_end) ||
x1: $(round(x1; digits = 2)) |
y1: $(round(y1; digits = 2)) ||
x2: $(round(x2; digits = 2)) |
y2: $(round(y2; digits = 2)) ||
"""

# ╔═╡ 9cc0df8c-422e-40f3-b727-f71ed682682b
md"## Three-Body Problem"

# ╔═╡ 604b2c88-f8fb-47cd-b07b-45833a188696
md"Source: [Wikipedia Article on Three-Body Problem](https://en.wikipedia.org/wiki/Three-body_problem)"

# ╔═╡ 88e2cba3-934c-4cc5-86ad-f0c16b6126cc
md"### Describe Problem"

# ╔═╡ 3ec52742-6da5-4e9e-a466-19eacb7e52f0
md"#### System of First Order Ordinary Differential Equations"

# ╔═╡ 90dbd390-48b9-4420-aca9-4dceca5fe89c
md"$\begin{align}

\mathbf{{\ddot{r}}_1} &=
	-G m_2 \frac{{\mathbf{r_1}} - {\mathbf{r_2}}}
		{{|{\mathbf{r_1}} - {\mathbf{r_2}}|}^3} -
	G m_3 \frac{{\mathbf{r_1}} - {\mathbf{r_3}}}
		{{|{\mathbf{r_1}} - {\mathbf{r_3}}|}^3} \\

\mathbf{{\ddot{r}}_2} &=
	-G m_3 \frac{{\mathbf{r_2}} - {\mathbf{r_3}}}
		{{|{\mathbf{r_2}} - {\mathbf{r_3}}|}^3} -
	G m_1 \frac{{\mathbf{r_2}} - {\mathbf{r_1}}}
		{{|{\mathbf{r_2}} - {\mathbf{r_1}}|}^3} \\

\mathbf{{\ddot{r}}_3} &=
	-G m_1 \frac{{\mathbf{r_3}} - {\mathbf{r_1}}}
		{{|{\mathbf{r_3}} - {\mathbf{r_1}}|}^3} -
	G m_2 \frac{{\mathbf{r_3}} - {\mathbf{r_1}}}
		{{|{\mathbf{r_3}} - {\mathbf{r_1}}|}^3} \\

\end{align}$"

# ╔═╡ 18ab2e3d-4304-45da-af0d-0ac79030fcf3
md"where the mass of the body is $m_i$,"

# ╔═╡ 1d5287ff-1c3f-4c78-9650-30fe711d3c42
md"the position of the body is $\mathbf{r_i} = (x_i, y_i, z_i)$,"

# ╔═╡ ef6d63a4-3b7e-4ae9-be75-77c26dc5429a
md"and $G$ is the gravitational constant."

# ╔═╡ a7210da8-b481-4604-94f3-4bed550e600e
md"##### DifferentialEquations.jl"

# ╔═╡ 89def556-a2b2-4123-b0ee-ae0bbcdfac2c
md"$f(du, u, p, t) =

\begin{cases}

du[1] &= vx_1 \\
du[2] &= vy_1 \\
du[3] &= vx_2 \\
du[4] &= vy_2 \\
du[5] &= vx_3 \\
du[6] &= vy_3 \\

du[7] &= {vx_1}' \\
du[8] &= {vy_1}' \\
du[9] &= {vx_2}' \\
du[10] &= {vy_2}' \\
du[11] &= {vx_3}' \\
du[12] &= {vy_3}' \\

\end{cases}$"

# ╔═╡ 795d2ec7-2ae8-477c-9e3b-5c240877f4b7
md"Define Julia Function"

# ╔═╡ fbf6dce6-3f37-4c71-99b9-9e0cb5568ddc
function three_body(du, u, p, t)

    # position, velocity
    x1, y1, x2, y2, x3, y3,
    vx1, vy1, vx2, vy2, vx3, vy3 = u

    # mass
    m1, m2, m3 = p

    # velocity
    du[1] = vx1
    du[2] = vy1

    du[3] = vx2
    du[4] = vy2

    du[5] = vx3
    du[6] = vy3

    # distance from body to body
    r12 = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    r13 = sqrt((x1 - x3)^2 + (y1 - y3)^2)
    r23 = sqrt((x2 - x3)^2 + (y2 - y3)^2)

    # acceleration (set G = 1)
    du[7] = -m2 * ((x1 - x2) / (r12^3)) - m3 * ((x1 - x3) / (r13^3))
    du[8] = -m2 * ((y1 - y2) / (r12^3)) - m3 * ((y1 - y3) / (r13^3))

    du[9] = -m3 * ((x2 - x3) / (r23^3)) - m1 * ((x2 - x1) / (r12^3))
    du[10] = -m3 * ((y2 - y3) / (r23^3)) - m1 * ((y2 - y1) / (r12^3))

    du[11] = -m1 * ((x3 - x1) / (r13^3)) - m2 * ((x3 - x2) / (r23^3))
    du[12] = -m1 * ((y3 - y1) / (r13^3)) - m2 * ((y3 - y2) / (r23^3))
end

# ╔═╡ d04f66f0-9bb5-4504-88d2-3266c0999413
md"### Assign Variables"

# ╔═╡ f01158e4-7774-4292-8c23-0bfd41c3234a
md"#### u: x1, y1, x2, y2, x3, y3, vx1, vy1, vx2, vy2, vx3, vy3"

# ╔═╡ 247a08e4-6fc5-420f-bd57-779d377d9195
md"#### Initial Positions"

# ╔═╡ 5f76da8a-a9a6-42c7-943f-e6691c10969d
md"##### Initial Conditions"

# ╔═╡ d5b7eb79-d465-4ccd-ba92-0e607ccedc25
md"#### p3: Parameters m1, m2, m3"

# ╔═╡ d9d9aa6d-ab33-48f1-be4e-8acbf527bd38
md"#### t3: Time"

# ╔═╡ eac02b2c-9975-4b75-a9d0-20a43c9c6006
t3_begin = 0.0

# ╔═╡ 65356d11-14dc-462d-b442-f4a115a44ea1
md"### Solve Problem"

# ╔═╡ f42dd494-ab86-41db-abc5-e9b24abff373
md"### Prep for Plotting"

# ╔═╡ 2476846e-f10a-4088-a17d-24b1efe54ad8
md"#### Obtain Position Values from Solution"

# ╔═╡ 2c63d96d-67a6-4404-99b6-f234b005a7ce
function values(sol; dt=0.001)
    idx = sol.t[1]:dt:sol.t[end]

    x1_values = map(value -> value[1], sol.(idx))
    y1_values = map(value -> value[2], sol.(idx))
    x2_values = map(value -> value[3], sol.(idx))
    y2_values = map(value -> value[4], sol.(idx))
    x3_values = map(value -> value[5], sol.(idx))
    y3_values = map(value -> value[6], sol.(idx))

    return ([x1_values, x2_values, x3_values],
        [y1_values, y2_values, y3_values])
end

# ╔═╡ 6fdded87-8502-4726-80ad-97c8a5021bc1
md"### Plot Solution"

# ╔═╡ ccebd5a3-192f-40d8-a27c-9a19dd46f946
julia = Colors.JULIA_LOGO_COLORS

# ╔═╡ 27243123-cc92-4726-abc7-0ca2aac127e8
md"#### Control Panel"

# ╔═╡ 8bb25845-80b2-4d5a-89b6-11199e833282
@bind reset3 Button("Reset")

# ╔═╡ 7b671c88-c649-455f-8149-7ac335b1753d
begin
    reset3
    @bind t3_end Slider(0.0:0.001:100.0, 100.0, true)
end

# ╔═╡ cd44fe85-8cda-4ca2-b97c-26147041d81c
t3span = (t3_begin, t3_end)

# ╔═╡ 32694d3e-3d29-4787-a918-1ca3fc17c272
begin
    reset3
    md"""
    Positions ||
    x1: $(@bind x1_0 Scrubbable(-2.0:0.1:2.0, default = -1.0)) |
    y1: $(@bind y1_0 Scrubbable(-2.0:0.1:2.0, default =  0.0)) |
    x2: $(@bind x2_0 Scrubbable(-2.0:0.1:2.0, default =  1.0)) |
    y2: $(@bind y2_0 Scrubbable(-2.0:0.1:2.0, default =  0.0)) |
    x3: $(@bind x3_0 Scrubbable(-2.0:0.1:2.0, default =  0.0)) |
    y3: $(@bind y3_0 Scrubbable(-2.0:0.1:2.0, default = sqrt(3))) |
    """
end

# ╔═╡ 48383a1f-9bbe-49b9-b190-0e311a09e6fc
begin
    reset3
    md"""
    Velocities ||
    vx1: $(@bind vx1_0 Scrubbable(-2.0:0.0001:2.0, default = -0.148 * cosd(60),
    	format = ".4f")) |
    vy1: $(@bind vy1_0 Scrubbable(-2.0:0.0001:2.0, default =  0.148 * sind(60),
    	format = ".4f")) |
    vx2: $(@bind vx2_0 Scrubbable(-2.0:0.0001:2.0, default = -0.148 * cosd(60),
    	format = ".4f")) |
    vy2: $(@bind vy2_0 Scrubbable(-2.0:0.0001:2.0, default = -0.148 * sind(60),
    	format = ".4f")) |
    vx3: $(@bind vx3_0 Scrubbable(-2.0:0.0001:2.0, default =  0.148,
    	format = ".4f")) |
    vy3: $(@bind vy3_0 Scrubbable(-2.0:0.0001:2.0, default =  0.0,
    	format = ".4f")) |
    """
end

# ╔═╡ 3736adf4-95d2-4116-96c4-ebb341b8b4cc
u3_begin = [
    x1_0, y1_0, x2_0, y2_0, x3_0, y3_0,
    vx1_0, vy1_0, vx2_0, vy2_0, vx3_0, vy3_0,
]

# ╔═╡ ebf34acc-7c11-42eb-887a-578d99b23345
begin
    reset3
    md"""
    Masses ||
    m1: $(@bind m31 Scrubbable(0.01:0.01:2.0, default = 0.1, format = ".2f")) |
    m2: $(@bind m32 Scrubbable(0.01:0.01:2.0, default = 0.1, format = ".2f")) |
    m3: $(@bind m33 Scrubbable(0.01:0.01:2.0, default = 0.1, format = ".2f")) |
    """
end

# ╔═╡ 2362a0e7-aeef-4324-a06b-4196ad32a805
p3 = [m31, m32, m33]

# ╔═╡ aa67af0c-3047-4412-8ee0-1a4a90845d87
prob3 = ODEProblem(three_body, u3_begin, t3span, p3)

# ╔═╡ 01836ec2-d4ff-4d97-85dd-05c719619d9c
sol3 = solve(prob3, Vern7(), reltol=1e-16, abstol=1e-16)

# ╔═╡ 9acc05f3-c010-40c4-a18f-0d210fb63585
x_values, y_values = values(sol3)

# ╔═╡ d80ca788-b226-48c3-9912-6c58f90f62c8
x1_end = x_values[1][end]

# ╔═╡ 892b2d80-14db-4a2e-910d-bd1c141f4000
y1_end = y_values[1][end]

# ╔═╡ fd1679b9-0df0-4348-9011-72b5c9899a77
x2_end = x_values[2][end]

# ╔═╡ 34b6de5e-d74e-4a9a-8d01-290c700ad35e
y2_end = y_values[2][end]

# ╔═╡ 5705d937-d7b8-4162-99f7-7b442aa21499
x3_end = x_values[3][end]

# ╔═╡ d5fe4b95-13dc-4561-b18d-302501f5850a
y3_end = y_values[3][end]

# ╔═╡ 7634a4fe-94a8-4df0-bd36-44a99c331ca3
x_values[1]

# ╔═╡ d011a472-042e-4dc1-9a03-6d2ee3affadc
begin
    fig3xy = Figure()

    ax3xy = Axis(fig3xy[1, 1],
        title="Three-Body Problem",
        xlabel="x",
        ylabel="y")

    colors = [julia.red julia.purple julia.green]

    for (ind, color) in enumerate(colors)
        lines!(ax3xy, x_values[ind], y_values[ind], linewidth=2, color=color)
    end

    for pair in zip([x1_end, x2_end, x3_end], [y1_end, y2_end, y3_end], colors)
        scatter!(ax3xy, pair[1], pair[2], markersize=25, color=pair[3])
    end

    fig3xy
end

# ╔═╡ ec990ef6-267b-4e50-9a82-feb426183a49
begin
    fig3a = Figure()

    ax3a = Axis(fig3a[1, 1],
        xlabel="Time")

    for col in range(1, 4)
        lines!(ax3a, sol3.t, sol3[col, :], linewidth=2)
    end

    fig3a
end

# ╔═╡ Cell order:
# ╟─6eb6bd53-f864-4a82-8173-f533ea7beeaa
# ╟─cedd2cac-46f9-4875-afe1-b1b3a936135d
# ╠═835e13a7-badc-4bb4-92ee-61afa3972a71
# ╠═07582e60-4d07-49b8-a6ab-e962eec8d81b
# ╠═1e7ae697-d85a-4c1c-b0da-e863f7dc3e74
# ╠═3b7e547f-3c75-4624-a79b-d9ea6ceb0ead
# ╠═4bc6152e-3201-4678-a876-e5577984cd57
# ╠═ee40337c-6e10-42a6-8aba-8a69a1ace228
# ╠═d3ed21d2-d3c2-416d-84e4-560fdee7fb0a
# ╟─9ee9c82c-7d06-4e0c-a353-c1c2e86f71a3
# ╟─ba9eb20a-8c8f-455b-9d4d-f44601b4fb44
# ╟─cfb79dad-db47-4a66-9d23-cd1cb78e16c5
# ╟─7aa78e57-7b0d-43b2-b0a3-da4ba3a9d70c
# ╟─445b7e01-3867-4814-9237-999ec6600e39
# ╟─60969d34-b75d-4862-93bf-46addf6be92d
# ╟─4b6da8ef-c845-44b4-b980-976671240d26
# ╟─04b935a7-6c6d-4fe1-ba48-149de8e64046
# ╟─cc092311-1167-4cfa-aef7-4f63beb62b40
# ╟─36a0f007-6071-4224-a6ab-6cbcecf5c9ef
# ╟─b0195989-346a-41fd-8f7b-a9bf87424b29
# ╠═d8db785d-cdf4-4b94-91f3-5bbf9496878f
# ╟─6211c0cd-c282-450c-816d-921da466e267
# ╟─f328cad6-72c8-46b4-8368-5ec386a5fc6e
# ╠═9128427c-ecc1-450a-85ce-fdd37c45305f
# ╠═9f41835f-2cdf-4fda-b955-e4666d762ea1
# ╠═3aa7e90e-e889-4309-8814-5d0b8406b288
# ╠═b1dd866e-044a-47e4-b2e9-e9a540dc59ae
# ╠═536da0f6-c2c5-4987-9e3c-2124825be7a4
# ╟─b5f1925c-1eaa-479f-afd3-583eddc84007
# ╠═bd285b60-e43f-4754-9446-2ba293626819
# ╟─1fa233ba-873e-4066-91f0-2375651733f3
# ╠═06c64392-0e4c-4afa-8181-498aab2d586a
# ╠═c6ced053-87a0-4608-91d6-bc3ce8c5b574
# ╠═078bb957-3d1d-477f-9e64-cd94f7d01e6b
# ╟─5b638093-fefa-418f-99c4-2588448da696
# ╠═18471ae5-52f1-4549-86e0-ad1f00d4bcdb
# ╟─d37a967d-64c2-46ad-9fd6-0169be8c2b50
# ╠═2ebe56db-8d9a-488a-831e-7e8a5cbd5548
# ╠═d0511ed6-7cb1-4b1d-a246-bfeea8d25a9c
# ╠═0cbcc92b-c817-47cb-bd5d-24244b78cc70
# ╠═c177c824-04f4-49e0-be7f-aec181595632
# ╠═9b8949f6-acfd-4668-ae07-a51161cbf789
# ╠═f2d2c387-5f23-47d6-9d8b-3d56f7d5ce8a
# ╟─550f5e88-7aaf-44a7-8ea9-243e3d3d6600
# ╠═6cdd8d51-2c4e-4ef0-9d5b-5cd982f7aa6f
# ╠═624eb883-9ec3-4c02-8824-dfa9a0099eea
# ╠═ee43de08-3ec1-4517-98b6-456b2db9d2b3
# ╟─367fea2b-84e6-42da-addd-040a64345cc2
# ╠═77d77b90-1726-46d2-9cfd-d85115c14290
# ╟─ecff4935-d3cd-4ae4-b92a-00c172971b28
# ╟─0880ed90-9970-4606-b45a-d9b58cd1f739
# ╟─d8fea5d8-4b66-4593-b994-b5d92f30a9ad
# ╟─4a7ff5e7-e08e-4788-9dba-00e3bd2005a2
# ╟─9cc0df8c-422e-40f3-b727-f71ed682682b
# ╟─604b2c88-f8fb-47cd-b07b-45833a188696
# ╟─88e2cba3-934c-4cc5-86ad-f0c16b6126cc
# ╟─3ec52742-6da5-4e9e-a466-19eacb7e52f0
# ╟─90dbd390-48b9-4420-aca9-4dceca5fe89c
# ╟─18ab2e3d-4304-45da-af0d-0ac79030fcf3
# ╟─1d5287ff-1c3f-4c78-9650-30fe711d3c42
# ╟─ef6d63a4-3b7e-4ae9-be75-77c26dc5429a
# ╟─a7210da8-b481-4604-94f3-4bed550e600e
# ╟─89def556-a2b2-4123-b0ee-ae0bbcdfac2c
# ╟─795d2ec7-2ae8-477c-9e3b-5c240877f4b7
# ╠═fbf6dce6-3f37-4c71-99b9-9e0cb5568ddc
# ╟─d04f66f0-9bb5-4504-88d2-3266c0999413
# ╟─f01158e4-7774-4292-8c23-0bfd41c3234a
# ╟─247a08e4-6fc5-420f-bd57-779d377d9195
# ╟─5f76da8a-a9a6-42c7-943f-e6691c10969d
# ╠═3736adf4-95d2-4116-96c4-ebb341b8b4cc
# ╟─d5b7eb79-d465-4ccd-ba92-0e607ccedc25
# ╠═2362a0e7-aeef-4324-a06b-4196ad32a805
# ╟─d9d9aa6d-ab33-48f1-be4e-8acbf527bd38
# ╠═eac02b2c-9975-4b75-a9d0-20a43c9c6006
# ╠═cd44fe85-8cda-4ca2-b97c-26147041d81c
# ╠═aa67af0c-3047-4412-8ee0-1a4a90845d87
# ╟─65356d11-14dc-462d-b442-f4a115a44ea1
# ╠═01836ec2-d4ff-4d97-85dd-05c719619d9c
# ╟─f42dd494-ab86-41db-abc5-e9b24abff373
# ╟─2476846e-f10a-4088-a17d-24b1efe54ad8
# ╠═2c63d96d-67a6-4404-99b6-f234b005a7ce
# ╠═9acc05f3-c010-40c4-a18f-0d210fb63585
# ╠═d80ca788-b226-48c3-9912-6c58f90f62c8
# ╠═892b2d80-14db-4a2e-910d-bd1c141f4000
# ╠═fd1679b9-0df0-4348-9011-72b5c9899a77
# ╠═34b6de5e-d74e-4a9a-8d01-290c700ad35e
# ╠═5705d937-d7b8-4162-99f7-7b442aa21499
# ╠═d5fe4b95-13dc-4561-b18d-302501f5850a
# ╟─6fdded87-8502-4726-80ad-97c8a5021bc1
# ╠═ec990ef6-267b-4e50-9a82-feb426183a49
# ╠═ccebd5a3-192f-40d8-a27c-9a19dd46f946
# ╠═7634a4fe-94a8-4df0-bd36-44a99c331ca3
# ╠═d011a472-042e-4dc1-9a03-6d2ee3affadc
# ╟─27243123-cc92-4726-abc7-0ca2aac127e8
# ╠═8bb25845-80b2-4d5a-89b6-11199e833282
# ╟─7b671c88-c649-455f-8149-7ac335b1753d
# ╟─32694d3e-3d29-4787-a918-1ca3fc17c272
# ╟─48383a1f-9bbe-49b9-b190-0e311a09e6fc
# ╟─ebf34acc-7c11-42eb-887a-578d99b23345
