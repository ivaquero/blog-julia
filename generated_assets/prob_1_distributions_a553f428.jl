### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> chapter = 3
#> section = 1
#> order = 1
#> title = "Distributions"
#> tags = ["prob", "track_prob", "Pluto", "PlutoUI"]
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

# ╔═╡ da32c4aa-e936-4f1d-b327-d3396a0d06c3
begin
    using Pkg
	if isdir("../../pluto-deployment-environment")
    	Pkg.activate("../../pluto-deployment-environment")
    	Pkg.instantiate()
	else
		println(pwd())
	end
end

# ╔═╡ aaf27acb-a9e5-4674-8d66-4bb6bdbc44ba
using Distributions

# ╔═╡ 73663d4c-52af-448c-9058-2c6e87d49e0c
using CairoMakie

# ╔═╡ fe318e8d-77e2-4ab5-b60d-c3879e7ae85e
using PlutoUI

# ╔═╡ 8e2814b8-436e-4156-bfd5-d1ba90f6fe54
using PlutoUI: Slider

# ╔═╡ 49bc4b90-c42e-11ed-25e4-bbb8e4953cf1
md"# Probability Distributions"

# ╔═╡ c6449c89-9128-49a2-bf54-6d7d5b51d862
TableOfContents()

# ╔═╡ 11c4e8aa-d35d-47ef-ae94-60e4f91c8ee1
md"## Continuous Distributions"

# ╔═╡ 3d069ab4-d11e-46b0-8245-e2cf8c3ae7bb
md"### Normal Distribution"

# ╔═╡ 0e4baa34-34d8-47ad-9b4a-34cff7863615
md"Mean (μ): $(@bind μ Slider(-3.0:0.1:3.0, 0.0, true))"

# ╔═╡ fb376e87-15a8-4363-942b-e36b83d73fdb
md"Standard Deviation (σ): $(@bind σ Slider(0.5:0.01:3.0, 1.0, true))"

# ╔═╡ d72c319f-63cf-4ef7-9c39-8107a67b1828
begin
    fign = Figure()

    axn = Axis(fign[1, 1],
        title="Normal Distribution, N($μ, $σ)",
        xlabel="x",
        ylabel="pdf(x)")

    lines!(axn, Normal(μ, σ))
    vlines!([μ])
    vlines!([μ + σ, μ - σ])
    xlims!(-5, 5)
    ylims!(0, 0.9)

    fign
end

# ╔═╡ 412c36a7-c779-4b62-9176-a0488fefac11
md"### Beta Distribution"

# ╔═╡ b0128411-7b22-45c7-8bfb-3e2205cfd9a1
md"α: $(@bind α Slider(0.1:0.1:5.0, 1.0, true))"

# ╔═╡ 6a2808f4-59a6-4182-81f2-567705493449
md"β: $(@bind β Slider(0.1:0.1:5.0, 1.0, true))"

# ╔═╡ 4ed0a3c5-5124-4d5a-8d46-8e1a721f7660
begin
    figb = Figure()

    axb = Axis(figb[1, 1],
        title="Beta Distribution, Beta($α, $β)",
        xlabel="x",
        ylabel="pdf(x)")

    lines!(axb, Beta(α, β))
    bx = 0:0.01:1
    lines!(axb,
        bx,
        zeros(101, 1)[:, 1],
        pdf.(Beta(α, β), bx),
        color=:lightblue)
    xlims!(0, 1)
    ylims!(0, 5)

    figb
end

# ╔═╡ 2bdba2e2-64de-45e6-9146-ff1913e3f529
md"### Uniform Distribution"

# ╔═╡ 6bc3bbd2-ab8e-4db1-949f-53ceaa10463a
md"a: $(@bind ua Slider(-5.:1.0:3.0, 1.0, true))"

# ╔═╡ 279df7fc-219a-461a-b359-c5b2cc9e7676
md"b: $(@bind ub Slider(4:1.0:9, 1.0, true))"

# ╔═╡ 2d52058a-e947-4f0d-a98c-0fee250d2def
begin
    figu = Figure()
    axu = Axis(figu[1, 1],
        title="Uniform Distribution, U($ua, $ub)",
        xlabel="x",
        ylabel="pdf(x)")

    ux = -5.0:9.0
    lines!(axu,
        ux,
        pdf.(Uniform(ua, ub), ux),
        linewidth=3)

    # xlims!(-5, 12)
    # ylims!(0, 1)

    figu
end

# ╔═╡ 1cc23123-8039-417b-a804-3c68065c2954
md"## Discrete Distributions"

# ╔═╡ 0daa3875-a08f-44df-8e23-0cffb4fa4b78
md"### Bernoulli Distribution"

# ╔═╡ ebbc484b-c5d7-4263-a0a3-d4e5d367f850
md"Success Rate (p): $(@bind pb Slider(0.00:0.01:1.0, 0.5, true))"

# ╔═╡ 1f76995f-7ab3-4739-b37d-47af87a14cff
begin
    figd = Figure()
    axd = Axis(figd[1, 1],
        title="Bernoulli Distribution, B($pb)",
        xlabel="x",
        ylabel="pmf(x)",
        xticks=-1:1:10)

    stem!(axd,
        Bernoulli(pb),
        linewidth=3,
        markersize=12,
        color=:white,
        stemcolor=:blue,
        strokecolor=:blue,
        strokewidth=2)
    scatter!(axd, Bernoulli(pb), markersize=8)

    xlims!(-0.2, 1.2)
    ylims!(0, 1)

    figd
end

# ╔═╡ 279c7fe4-8a46-465b-a647-73bcc36a33b1
md"### Binomial Distribution"

# ╔═╡ bc2dd470-8c03-41c1-ae9a-84f3f92882c2
md"Number of Trials (n): $(@bind n Slider(1:1:10, 1, true))"

# ╔═╡ 51a0937c-1e56-46ce-bb74-4e35b51e8bd3
md"Probability of Success (p): $(@bind p Slider(0.00:0.01:1.0, 0.5, true))"

# ╔═╡ cb546aaf-e442-4d5c-b6ac-6c0f41db54a1
begin
    figs = Figure()
    axs = Axis(figs[1, 1],
        title="Binomial Distribution, B($n, $p)",
        xlabel="x",
        ylabel="pmf(x)")

    stem!(axs,
        Binomial(n, p),
        linewidth=3,
        markersize=12,
        color=:white,
        stemcolor=:blue,
        strokecolor=:blue,
        strokewidth=2)
    scatter!(axs, Binomial(n, p), markersize=8)

    # xlims!(-0.5, 10.5)
    # ylims!(0, 1)

    figs
end

# ╔═╡ bf400397-6449-46ec-9127-7c759c48adcb
md"### Poisson Distribution"

# ╔═╡ 1462928b-3c6a-4106-9c68-ce45f8267bc1
md"Average Rate of Occurrences (λ): $(@bind λ Slider(1.0:0.1:10.0, 1.0, true))"

# ╔═╡ 507fe99c-4308-4916-8fbb-3aaf52d781f8
begin
    figpo = Figure()
    axpo = Axis(figpo[1, 1],
        title="Poisson Distribution, B($λ)",
        xlabel="x",
        ylabel="pmf(x)")

    stem!(axpo,
        Poisson(λ),
        linewidth=3,
        markersize=12,
        color=:white,
        stemcolor=:blue,
        strokecolor=:blue,
        strokewidth=2)
    scatter!(axpo, Poisson(λ), markersize=8)

    # xlims!(0, 10)
    # ylims!(0, 1)

    figpo
end

# ╔═╡ Cell order:
# ╟─49bc4b90-c42e-11ed-25e4-bbb8e4953cf1
# ╠═da32c4aa-e936-4f1d-b327-d3396a0d06c3
# ╠═aaf27acb-a9e5-4674-8d66-4bb6bdbc44ba
# ╠═73663d4c-52af-448c-9058-2c6e87d49e0c
# ╠═fe318e8d-77e2-4ab5-b60d-c3879e7ae85e
# ╠═8e2814b8-436e-4156-bfd5-d1ba90f6fe54
# ╠═c6449c89-9128-49a2-bf54-6d7d5b51d862
# ╟─11c4e8aa-d35d-47ef-ae94-60e4f91c8ee1
# ╟─3d069ab4-d11e-46b0-8245-e2cf8c3ae7bb
# ╟─0e4baa34-34d8-47ad-9b4a-34cff7863615
# ╟─fb376e87-15a8-4363-942b-e36b83d73fdb
# ╠═d72c319f-63cf-4ef7-9c39-8107a67b1828
# ╟─412c36a7-c779-4b62-9176-a0488fefac11
# ╟─b0128411-7b22-45c7-8bfb-3e2205cfd9a1
# ╟─6a2808f4-59a6-4182-81f2-567705493449
# ╠═4ed0a3c5-5124-4d5a-8d46-8e1a721f7660
# ╟─2bdba2e2-64de-45e6-9146-ff1913e3f529
# ╠═6bc3bbd2-ab8e-4db1-949f-53ceaa10463a
# ╠═279df7fc-219a-461a-b359-c5b2cc9e7676
# ╠═2d52058a-e947-4f0d-a98c-0fee250d2def
# ╟─1cc23123-8039-417b-a804-3c68065c2954
# ╟─0daa3875-a08f-44df-8e23-0cffb4fa4b78
# ╟─ebbc484b-c5d7-4263-a0a3-d4e5d367f850
# ╠═1f76995f-7ab3-4739-b37d-47af87a14cff
# ╟─279c7fe4-8a46-465b-a647-73bcc36a33b1
# ╟─bc2dd470-8c03-41c1-ae9a-84f3f92882c2
# ╟─51a0937c-1e56-46ce-bb74-4e35b51e8bd3
# ╠═cb546aaf-e442-4d5c-b6ac-6c0f41db54a1
# ╟─bf400397-6449-46ec-9127-7c759c48adcb
# ╟─1462928b-3c6a-4106-9c68-ce45f8267bc1
# ╠═507fe99c-4308-4916-8fbb-3aaf52d781f8
