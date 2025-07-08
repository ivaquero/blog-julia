### A Pluto.jl notebook ###
# v0.20.13

#> [frontmatter]
#> chapter = 2
#> section = 7
#> order = 7
#> title = "Partial Differential Equation"
#> tags = ["math", "track_comp", "Pluto", "PlutoUI"]
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

# ╔═╡ 9f1b8143-4882-4a55-8d0a-e2a09e1e6269
begin
    using Pkg
    if isdir("../../pluto-deployment-environment")
        Pkg.activate("../../pluto-deployment-environment")
        Pkg.instantiate()
    else
        println(pwd())
    end
end

# ╔═╡ b9bc3e57-74cd-4c9c-b5d9-726deec6c79d
using DifferentialEquations

# ╔═╡ eb9a290e-1f8c-40ef-8c22-385d4081b539
using CairoMakie

# ╔═╡ 6a74f193-a026-45c9-aeec-636af4764777
using PlutoUI

# ╔═╡ 0bc89c15-5f27-45de-b4cc-320e66703912
using PlutoUI: Button, Slider

# ╔═╡ e603881a-d6bb-4f16-ba6a-c414d45b9ccb
md"# Partial Differential Equation (PDE)"

# ╔═╡ a8f8881a-c501-41a6-8661-2051c01aa368
md"Sources:"

# ╔═╡ 6563abc2-5913-403b-b58e-306f7744171f
md"[Wikipedia Article on Heat Equation](https://en.wikipedia.org/wiki/Heat_equation)"

# ╔═╡ 1c1ff0ab-3840-4ba7-9086-8ae2e0579bf2
md"[Wikipedia Article on Partial Differential Equation](https://en.wikipedia.org/wiki/Partial_differential_equation)"

# ╔═╡ fd414020-ea77-43bc-a8c1-11587223f11d
md"[Wikipedia Article on Method of Lines](https://en.wikipedia.org/wiki/Method_of_lines)"

# ╔═╡ d9aaa8a9-4b3f-46df-8372-68018391a9cc
md"[SciML Notes on PDEs](https://book.sciml.ai/notes/14-PDEs_Convolutions_and_the_Mathematics_of_Locality/)"

# ╔═╡ fc898119-fc6d-4c32-863f-baa70e30080a
md"[3Blue1Brown YouTube Video on PDE (DE2)](https://youtu.be/ly4S0oi3Yz8?list=PLZHQObOWTQDNPOjrT6KVlfJuKtYTftqH6)"

# ╔═╡ 53dc7233-fc51-4b46-8f7b-ce81f4b8e4b6
md"## Add Packages"

# ╔═╡ f0393b21-9dbe-4017-bf40-6e574004b36e
TableOfContents()

# ╔═╡ 33426f6e-71b3-11ed-092c-910e0f4bd7f0
md"## 1-Dimensional Heat Equation"

# ╔═╡ 4c318d00-9eb5-492a-aace-4f54d9c2cbca
md"### Describe Problem"

# ╔═╡ 4f5703ec-7d60-4764-96c7-8918a22ea95a
md"#### Heat Equation"

# ╔═╡ c8f1ea66-4b51-4c6c-9cc4-b9e3002f981c
md"$\frac{∂ u}{∂ t} = α \frac{∂^2 u}{∂ x^2}$"

# ╔═╡ 94071c80-a838-4fea-ae57-29324e6013e3
md"where $u$ is Temperature, $x$ is position, $t$ is time,"

# ╔═╡ 63b827df-13e9-4c0d-8b1b-e88b4535889b
md"and $α$ is a positive coefficient called the thermal diffusivity of the medium."

# ╔═╡ 3b883662-6eec-4921-86a2-1d90246c2af2
md"$\frac{∂^2 u}{∂ x^2} = \frac{u_{i-1} - 2 u_i + u_{i+1}}{dx^2}$"

# ╔═╡ 0dd536c0-ba2c-4ee7-8e06-51e339a93aa0
md"#### DifferentialEquations.jl"

# ╔═╡ 01985ac8-fe8b-40b3-846d-3646622153cf
md"$f(du, u, p, t) =
\begin{cases}
du[1] = p[1] \left( \frac{-u[1] + u[2]}{p[2]^2} \right) \\
du[2:(\text{end} - 1)] = p[1] \left( \frac{u[i - 1] - 2 u[i] + u[i + 1]}{p[2]^2} \right) \\
du[\text{end}] = p[1] \left( \frac{u[\text{end} - 1] - u[\text{end}]}{p[2]^2} \right)
\end{cases}$"

# ╔═╡ 408c8fd9-533f-4916-89cb-02167366a2f0
md"where $du[1]$ is the left boundary, $du[\text{end}]$ is the right boundary, $du[2:(\text{end} - 1 )]$ are interior functions,"

# ╔═╡ 47348043-c8fa-4318-8f81-ba6df2d51e7a
md"and $p[1]$ is $\alpha$ and $p[2]$ is $dx$, which is a discrete step in position."

# ╔═╡ 003f6317-f1cc-43b5-b193-98be6b96895e
md"### Define Julia Function"

# ╔═╡ 2d233e27-5083-4567-88e1-58e91d82bdc4
function heat_equation(du, u, p, t)

    # Parameters
    α, dx = p

    # Left Boundary
    du[1] = α * (-u[1] + u[2]) / dx^2

    # Interior Points
    for i in 2:(length(u)-1)
        du[i] = α * (u[i-1] - 2 * u[i] + u[i+1]) / dx^2
    end

    # Right Boundary
    du[end] = α * (u[end-1] - u[end]) / dx^2
end

# ╔═╡ 36b14a8c-9acb-48d9-84d7-682f57bcc317
md"### Assign Variables"

# ╔═╡ 3e27ffce-19b6-45eb-b725-083f3b9e8c21
md"#### u: Temperature T is a function of position x and time t"

# ╔═╡ d1d68c7f-df1f-4d6e-802c-5df4b93fefa3
md"#### Position Along 1-D Wire"

# ╔═╡ eb5e2ae2-4ca8-494a-a2ee-657f228f8b03
L = 1.0

# ╔═╡ f0d62914-2213-43a2-acd3-def383c421ad
dx = 0.01

# ╔═╡ fe8f51c2-e304-458f-9ca2-fa0a6729de80
x = 0.0:dx:L

# ╔═╡ 994a11d7-0d2a-4a49-8ab5-b588d9260374
md"#### Initial Temperatures Along 1-D Wire"

# ╔═╡ db6f27d2-ccad-477f-b299-d958f90e9425
md"#### p: Parameters α and dx"

# ╔═╡ bd841703-fd38-4fbc-b76d-536c8234b50d
md"Note: dx already defined above."

# ╔═╡ b78c836f-7abb-46eb-b161-9bff3ae545c3
md"#### t: Time"

# ╔═╡ ef143c00-d28d-436b-bbba-e6ad96c7f8c4
t_begin = 0.0

# ╔═╡ 3dc6ab16-dba5-478b-b65f-9e3a3b0c1492
md"### Define System of First Order ODEs Problem"

# ╔═╡ c68506fc-0b58-4187-8884-92392ed78f0d
md"### Solve Problem"

# ╔═╡ 4550b812-6951-4207-9a0d-4b1459c8b3f5
md"### Plot Solution"

# ╔═╡ dfe84d47-a51a-46a4-9d3c-c97e776e2574
md"### Control Panel"

# ╔═╡ 7fd17ac7-4f34-49d2-af6d-7628b75bb5ff
@bind reset Button("Reset")

# ╔═╡ dd28ee23-cfc4-4b61-a5e4-f386e5974cfb
begin
    reset
    md"Time: $(@bind t_end Slider(0.0:0.001:1.0; default=0.0, show_value=true))"
end

# ╔═╡ e4696021-f124-4293-b824-ada27d2ecd5a
tspan = (t_begin, t_end)

# ╔═╡ 0559872d-4cce-4396-9ba1-a539aeeb8f79
begin
    reset
    md"""
    Temperature ||
    Left: $(@bind left_temp Slider(0.0:0.01:1.0; default=0.9, show_value=true)) |
    Right: $(@bind right_temp Slider(0.0:0.01:1.0; default=0.1, show_value=true)) ||
    """
end

# ╔═╡ 03b4477e-4cdf-45a0-a51b-4d6067d1ee8f
u_begin = [i < L / 2 ? left_temp : right_temp for i in x]

# ╔═╡ b13e8490-6b04-493a-8bca-17fcdd994484
begin
    reset
    md"α: $(@bind α Slider(0.0:0.001:2.0; default=0.5, show_value=true))"
end

# ╔═╡ c231f801-6423-4726-9324-be36a9f1a375
p = [α, dx]

# ╔═╡ 34252e38-3e9d-4038-9dcb-dfa130dd7442
prob = ODEProblem(heat_equation, u_begin, tspan, p)

# ╔═╡ 6a1e3f99-0cf5-4773-9093-c8888b65fab0
sol = solve(prob)

# ╔═╡ b15e0cca-d357-4062-9464-dddb623b4880
T = sol[end]

# ╔═╡ c9bb9cf3-e82b-494e-95b0-45c68c4e918b
Tleft = T[1]

# ╔═╡ b9eae908-52e3-4cdf-8af7-5dd7c6147aee
T25 = T[25]

# ╔═╡ ccba2904-2dcc-42e1-a90d-b6ef608aeef5
T50 = T[50]

# ╔═╡ c8a37386-3d11-473f-81c9-5ebf7b7a1757
T51 = T[51]

# ╔═╡ 6b7784fa-c42c-4d3a-bb0c-f46a3ae89b79
T76 = T[76]

# ╔═╡ e5376c51-5df6-48a8-bae1-de66d7939d93
Tright = T[end]

# ╔═╡ 4914a69e-082f-4e43-b26c-ce3cf04a66fd
md"""
Temperatures along wire ||
Tleft: $(round(Tleft; digits = 2)) |
T25: $(round(T25; digits = 2)) |
T50: $(round(T50; digits = 2)) |
T51: $(round(T51; digits = 2)) |
T76: $(round(T76; digits = 2)) |
Tright: $(round(Tright; digits = 2)) ||
"""

# ╔═╡ 09c5a78b-48b8-49fd-ad82-2d8b8ddca5f6
begin
    fig = Figure()

    ax = Axis(fig[1, 1],
        title="1-D Heat Equation",
        xlabel="Position (x)",
        ylabel="Temperature (T)")

    lines!(ax, x, T, linewidth=2)

    fig
end

# ╔═╡ Cell order:
# ╟─e603881a-d6bb-4f16-ba6a-c414d45b9ccb
# ╟─a8f8881a-c501-41a6-8661-2051c01aa368
# ╟─6563abc2-5913-403b-b58e-306f7744171f
# ╟─1c1ff0ab-3840-4ba7-9086-8ae2e0579bf2
# ╟─fd414020-ea77-43bc-a8c1-11587223f11d
# ╟─d9aaa8a9-4b3f-46df-8372-68018391a9cc
# ╟─fc898119-fc6d-4c32-863f-baa70e30080a
# ╟─53dc7233-fc51-4b46-8f7b-ce81f4b8e4b6
# ╠═9f1b8143-4882-4a55-8d0a-e2a09e1e6269
# ╠═b9bc3e57-74cd-4c9c-b5d9-726deec6c79d
# ╠═eb9a290e-1f8c-40ef-8c22-385d4081b539
# ╠═6a74f193-a026-45c9-aeec-636af4764777
# ╠═0bc89c15-5f27-45de-b4cc-320e66703912
# ╠═f0393b21-9dbe-4017-bf40-6e574004b36e
# ╟─33426f6e-71b3-11ed-092c-910e0f4bd7f0
# ╟─4c318d00-9eb5-492a-aace-4f54d9c2cbca
# ╟─4f5703ec-7d60-4764-96c7-8918a22ea95a
# ╟─c8f1ea66-4b51-4c6c-9cc4-b9e3002f981c
# ╟─94071c80-a838-4fea-ae57-29324e6013e3
# ╟─63b827df-13e9-4c0d-8b1b-e88b4535889b
# ╟─3b883662-6eec-4921-86a2-1d90246c2af2
# ╟─0dd536c0-ba2c-4ee7-8e06-51e339a93aa0
# ╟─01985ac8-fe8b-40b3-846d-3646622153cf
# ╟─408c8fd9-533f-4916-89cb-02167366a2f0
# ╟─47348043-c8fa-4318-8f81-ba6df2d51e7a
# ╟─003f6317-f1cc-43b5-b193-98be6b96895e
# ╠═2d233e27-5083-4567-88e1-58e91d82bdc4
# ╟─36b14a8c-9acb-48d9-84d7-682f57bcc317
# ╟─3e27ffce-19b6-45eb-b725-083f3b9e8c21
# ╟─d1d68c7f-df1f-4d6e-802c-5df4b93fefa3
# ╠═eb5e2ae2-4ca8-494a-a2ee-657f228f8b03
# ╠═f0d62914-2213-43a2-acd3-def383c421ad
# ╠═fe8f51c2-e304-458f-9ca2-fa0a6729de80
# ╟─994a11d7-0d2a-4a49-8ab5-b588d9260374
# ╠═03b4477e-4cdf-45a0-a51b-4d6067d1ee8f
# ╟─db6f27d2-ccad-477f-b299-d958f90e9425
# ╟─bd841703-fd38-4fbc-b76d-536c8234b50d
# ╠═c231f801-6423-4726-9324-be36a9f1a375
# ╟─b78c836f-7abb-46eb-b161-9bff3ae545c3
# ╠═ef143c00-d28d-436b-bbba-e6ad96c7f8c4
# ╠═e4696021-f124-4293-b824-ada27d2ecd5a
# ╟─3dc6ab16-dba5-478b-b65f-9e3a3b0c1492
# ╠═34252e38-3e9d-4038-9dcb-dfa130dd7442
# ╟─c68506fc-0b58-4187-8884-92392ed78f0d
# ╠═6a1e3f99-0cf5-4773-9093-c8888b65fab0
# ╠═b15e0cca-d357-4062-9464-dddb623b4880
# ╠═c9bb9cf3-e82b-494e-95b0-45c68c4e918b
# ╠═b9eae908-52e3-4cdf-8af7-5dd7c6147aee
# ╠═ccba2904-2dcc-42e1-a90d-b6ef608aeef5
# ╠═c8a37386-3d11-473f-81c9-5ebf7b7a1757
# ╠═6b7784fa-c42c-4d3a-bb0c-f46a3ae89b79
# ╠═e5376c51-5df6-48a8-bae1-de66d7939d93
# ╟─4550b812-6951-4207-9a0d-4b1459c8b3f5
# ╟─dfe84d47-a51a-46a4-9d3c-c97e776e2574
# ╠═7fd17ac7-4f34-49d2-af6d-7628b75bb5ff
# ╟─dd28ee23-cfc4-4b61-a5e4-f386e5974cfb
# ╟─0559872d-4cce-4396-9ba1-a539aeeb8f79
# ╟─b13e8490-6b04-493a-8bca-17fcdd994484
# ╟─4914a69e-082f-4e43-b26c-ce3cf04a66fd
# ╠═09c5a78b-48b8-49fd-ad82-2d8b8ddca5f6
