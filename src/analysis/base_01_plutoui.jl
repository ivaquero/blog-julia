### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> chapter = 1
#> section = 0
#> order = 0
#> title = "Pluto UI"
#> tags = ["analysis", "track_data", "Pluto", "PlutoUI"]
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

# ╔═╡ b721b6dc-c66a-400c-8f8f-253ea8acb18f
begin
    using Pkg
    Pkg.activate("../../pluto-deployment-environment")
    Pkg.instantiate()
end

# ╔═╡ d0884674-dc0e-11ed-188a-8d07ea88e62b
using PlutoUI

# ╔═╡ a87e5dbe-c309-4cf7-85bb-86e3ced6b055
using PlutoUI: Slider, Button

# ╔═╡ 4ab44cd0-72e5-47cb-a1ed-60e6985bddae
using CairoMakie

# ╔═╡ 1db55175-f4b5-4095-91b6-442bf48f4287
md"# PlutoUI Cheatsheet"

# ╔═╡ aced29d6-75ac-4c24-b497-9ea54b8c9d15
TableOfContents()

# ╔═╡ ba3f2bdf-2993-4d8b-a948-5f30987b8617
md"## Intro to PlutoUI"

# ╔═╡ 51d2100c-137e-4c2b-ad84-aa2a3beac96c
md"General tools:"

# ╔═╡ 1be81d6c-e7e4-4f13-aadb-e0264d5b698a
TableOfContents()

# ╔═╡ d7d60a95-3d40-4c04-8fd2-dfef7caa03c3
md"### with_terminal()"

# ╔═╡ 0af5bbf7-1d94-46d6-a02e-47903aedd556
with_terminal() do
    println("Hello, World!")
end

# ╔═╡ ef29ecb3-a3af-4592-a5fc-fd155e5d4e2e
md"## Interactive PlutoUI Widgets"

# ╔═╡ 8b2e6dad-327e-46c8-bbd7-ee910c8abf91
md"### Single Slider"

# ╔═╡ 75f80a8b-4164-4e46-a158-658f5bbb98e3
md"Fahrenheit to Celsius equation:"

# ╔═╡ 0ec94ad7-fe8e-400b-83f0-92ec2c07bbc2
md"$°C = \frac{5}{9} × (°F - 32)$"

# ╔═╡ 123a515c-6028-4eb9-911a-fec7f9a45bc0
@bind sliderF Slider(50.0:90.0, 72.0, true)

# ╔═╡ fc3a5370-2110-4444-842e-f13e1f18942b
sliderC = round(5 / 9 * (sliderF - 32), digits=1)

# ╔═╡ 239a38eb-e9d5-4748-af23-7255dcf5f29f
md"Temperature: $(sliderF) °F | $(sliderC) °C"

# ╔═╡ 748b8ca5-000c-4307-bf2d-b1e09b10fe6b
md"### Scrubbable"

# ╔═╡ febc3943-6723-48e7-8940-8143f8427145
md"Thermostat: $(@bind scrubF Scrubbable(50.0:90.0, default = 72.0)) °F"

# ╔═╡ a4e792ae-4251-4286-8a59-525641da53da
scrubC = round(5 / 9 * (scrubF - 32), digits=1)

# ╔═╡ 0868a272-9acb-48c1-88e5-de6c77e38cbc
md"Temperature: $(scrubF) °F | $(scrubC) °C"

# ╔═╡ 620caaeb-2474-4c99-b614-a5c012028040
md"### NumberField"

# ╔═╡ 3efcd737-7ac5-41df-9811-e5cb935dcd64
@bind numberF NumberField(50.0:90.0, default=72.0)

# ╔═╡ 46971901-b892-4973-98c9-42a102c50dd0
numberC = round(5 / 9 * (numberF - 32), digits=1)

# ╔═╡ b4df93eb-9919-4432-92b7-c84595526f83
md"Temperature: $(numberF) °F | $(numberC) °C"

# ╔═╡ db253112-b3e3-41d5-82c9-956418e724e0
md"### CheckBox"

# ╔═╡ e222646a-3b17-46dd-a715-bac895249526
md"Check the box to see the temperature in Celsius: $(@bind checkC CheckBox())"

# ╔═╡ d771d98f-6ace-4f3c-9c74-0950b7bd77d3
if checkC
    md"Temperature: $(numberC) °C"
else
    md"Temperature: $(numberF) °F"
end

# ╔═╡ fc216d91-71a1-4b7e-8da0-997490f78fcd
md"### TextField"

# ╔═╡ 4ae37cc0-664f-4744-9e2b-78240b5e029a
@bind textF TextField(default="72")

# ╔═╡ fe83c6c1-a596-4273-8e72-7cc12f2e41ac
parseF = parse(Float64, textF)

# ╔═╡ f1f3a5bc-cacc-40f8-882d-51bd39395eaf
textC = round(5 / 9 * (parseF - 32), digits=1)

# ╔═╡ b10f13fd-6fa0-4b4c-a4ba-4bb062b7dbda
md"Temperature: $(parseF) °F | $(textC) °C"

# ╔═╡ ce9fbc37-3833-40a0-a366-d9e923197bce
md"### Select"

# ╔═╡ 3ea2de2d-aad3-4ac4-9eb8-c78169e13c11
@bind selectScale Select(["Fahrenheit", "Celsius"])

# ╔═╡ ced2b4da-7278-40e5-86be-be29a4436dff
if selectScale == "Celsius"
    md"Temperature: $(numberC) °C"
else
    md"Temperature: $(numberF) °F"
end

# ╔═╡ ba8d6414-d9d1-4b1f-a4b5-0bf7ebe355b1
md"### MultiSelect"

# ╔═╡ 8361cd51-6b1c-4fc7-81b2-75d00c8ca1ce
@bind multiScales MultiSelect(["Fahrenheit", "Celsius"])

# ╔═╡ e19fbdcc-76f0-439d-9ba3-5a39afc0ca4e
if multiScales == ["Fahrenheit", "Celsius"]
    md"Temperature: $(numberF) °F | $(numberC) °C"
elseif multiScales == ["Fahrenheit"]
    md"Temperature: $(numberF) °F"
elseif multiScales == ["Celsius"]
    md"Temperature: $(numberC) °C"
else
    md"Please make selection(s)."
end

# ╔═╡ bcc5d6e1-7c75-4e7a-b23e-d4e420029bdd
md"### MultiCheckBox"

# ╔═╡ c0db1698-5c3d-494b-ad4b-4d25b9dbbd47
@bind boxScales MultiCheckBox(["Fahrenheit", "Celsius"])

# ╔═╡ 38a976d0-9378-4524-8212-165be96ae0da
if boxScales == ["Fahrenheit", "Celsius"]
    md"Temperature: $(numberF) °F | $(numberC) °C"
elseif boxScales == ["Fahrenheit"]
    md"Temperature: $(numberF) °F"
elseif boxScales == ["Celsius"]
    md"Temperature: $(numberC) °C"
else
    md"Please make selection(s)."
end

# ╔═╡ 56b8ac1c-42c9-43d2-a811-f7f94a05bd7f
md"### Button"

# ╔═╡ e73d3c85-d33f-4c45-998a-6b08a06dc783
@bind random Button("Random Data Generator")

# ╔═╡ 907ded68-4cbd-4b27-8663-15d630b4a51a
typeof(random)

# ╔═╡ b56bb17b-45fe-44cc-bd1f-e6ce4daca011
md"## Interactive Plot"

# ╔═╡ b5ec371d-dec5-4f40-9ff6-158e49a98f73
md"Equation of a straight line:"

# ╔═╡ 5e08562e-9ac3-4369-92d9-5ef5b7a4f5e0
md"$y = mx + b$"

# ╔═╡ d9d78b46-fe6e-4f93-8dcd-0c45ca2c1037
md"where $m$ is the rise-over-run slope and $b$ is the y-intercept."

# ╔═╡ 5aa696a3-df8a-484c-9dd1-c5e6076f5c13
begin
    slope = @bind m Slider(-10.0:0.1:10.0, 1.0, true)
    intercept = @bind b Slider(-10.0:0.1:10.0, 0.0, true)

    md"Parameters: \
    Slope: $(slope) \
    y-intercept: $(intercept)"
end

# ╔═╡ 8ab22640-bf43-4f6e-8a27-f54bc25e054c
x = range(-10, 10)

# ╔═╡ e9dcbe7c-3cf7-4f13-b59f-97acbe636087
scatter(x, markersize=3)

# ╔═╡ 720b426b-11f5-4016-8336-cef66aa988b8
f(x) = m .* x .+ b

# ╔═╡ 8a47b5a9-be3b-41e7-bdab-2a4ace9b4d4a
lines(f(x), linewidth=3)

# ╔═╡ Cell order:
# ╟─1db55175-f4b5-4095-91b6-442bf48f4287
# ╠═b721b6dc-c66a-400c-8f8f-253ea8acb18f
# ╠═d0884674-dc0e-11ed-188a-8d07ea88e62b
# ╠═a87e5dbe-c309-4cf7-85bb-86e3ced6b055
# ╠═4ab44cd0-72e5-47cb-a1ed-60e6985bddae
# ╠═aced29d6-75ac-4c24-b497-9ea54b8c9d15
# ╟─ba3f2bdf-2993-4d8b-a948-5f30987b8617
# ╟─51d2100c-137e-4c2b-ad84-aa2a3beac96c
# ╟─1be81d6c-e7e4-4f13-aadb-e0264d5b698a
# ╟─d7d60a95-3d40-4c04-8fd2-dfef7caa03c3
# ╠═0af5bbf7-1d94-46d6-a02e-47903aedd556
# ╟─ef29ecb3-a3af-4592-a5fc-fd155e5d4e2e
# ╟─8b2e6dad-327e-46c8-bbd7-ee910c8abf91
# ╟─75f80a8b-4164-4e46-a158-658f5bbb98e3
# ╟─0ec94ad7-fe8e-400b-83f0-92ec2c07bbc2
# ╠═123a515c-6028-4eb9-911a-fec7f9a45bc0
# ╠═fc3a5370-2110-4444-842e-f13e1f18942b
# ╟─239a38eb-e9d5-4748-af23-7255dcf5f29f
# ╟─748b8ca5-000c-4307-bf2d-b1e09b10fe6b
# ╠═febc3943-6723-48e7-8940-8143f8427145
# ╠═a4e792ae-4251-4286-8a59-525641da53da
# ╠═0868a272-9acb-48c1-88e5-de6c77e38cbc
# ╟─620caaeb-2474-4c99-b614-a5c012028040
# ╠═3efcd737-7ac5-41df-9811-e5cb935dcd64
# ╠═46971901-b892-4973-98c9-42a102c50dd0
# ╠═b4df93eb-9919-4432-92b7-c84595526f83
# ╟─db253112-b3e3-41d5-82c9-956418e724e0
# ╠═e222646a-3b17-46dd-a715-bac895249526
# ╠═d771d98f-6ace-4f3c-9c74-0950b7bd77d3
# ╟─fc216d91-71a1-4b7e-8da0-997490f78fcd
# ╠═4ae37cc0-664f-4744-9e2b-78240b5e029a
# ╠═fe83c6c1-a596-4273-8e72-7cc12f2e41ac
# ╠═f1f3a5bc-cacc-40f8-882d-51bd39395eaf
# ╠═b10f13fd-6fa0-4b4c-a4ba-4bb062b7dbda
# ╟─ce9fbc37-3833-40a0-a366-d9e923197bce
# ╠═3ea2de2d-aad3-4ac4-9eb8-c78169e13c11
# ╠═ced2b4da-7278-40e5-86be-be29a4436dff
# ╠═ba8d6414-d9d1-4b1f-a4b5-0bf7ebe355b1
# ╠═8361cd51-6b1c-4fc7-81b2-75d00c8ca1ce
# ╠═e19fbdcc-76f0-439d-9ba3-5a39afc0ca4e
# ╟─bcc5d6e1-7c75-4e7a-b23e-d4e420029bdd
# ╠═c0db1698-5c3d-494b-ad4b-4d25b9dbbd47
# ╠═38a976d0-9378-4524-8212-165be96ae0da
# ╟─56b8ac1c-42c9-43d2-a811-f7f94a05bd7f
# ╠═e73d3c85-d33f-4c45-998a-6b08a06dc783
# ╠═907ded68-4cbd-4b27-8663-15d630b4a51a
# ╠═e9dcbe7c-3cf7-4f13-b59f-97acbe636087
# ╟─b56bb17b-45fe-44cc-bd1f-e6ce4daca011
# ╟─b5ec371d-dec5-4f40-9ff6-158e49a98f73
# ╟─5e08562e-9ac3-4369-92d9-5ef5b7a4f5e0
# ╟─d9d78b46-fe6e-4f93-8dcd-0c45ca2c1037
# ╠═5aa696a3-df8a-484c-9dd1-c5e6076f5c13
# ╠═8ab22640-bf43-4f6e-8a27-f54bc25e054c
# ╠═720b426b-11f5-4016-8336-cef66aa988b8
# ╠═8a47b5a9-be3b-41e7-bdab-2a4ace9b4d4a
