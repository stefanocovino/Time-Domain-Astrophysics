### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 6a1315d1-9a6d-4ce0-b1c0-3fe22beb1ec2
begin
	using CommonMark
	using PlutoUI
end

# ╔═╡ 4d477519-c44f-434c-b7e0-8daaa5009358
md"""
**What is this?**


*This jupyter notebook is part of a collection of notebooks on various topics discussed during the Time Domain Astrophysics course delivered by Stefano Covino at the [Università dell'Insubria](https://www.uninsubria.eu/) in Como (Italy). Please direct questions and suggestions to [stefano.covino@inaf.it](mailto:stefano.covino@inaf.it).*
"""

# ╔═╡ ec67de24-d88a-46fa-ae46-c0cd7b797adc
md"""
**This is a `Pluto` notebook**
"""

# ╔═╡ 3ddd0f61-79d0-473c-8fec-a0e0c3fc72bf
TableOfContents()

# ╔═╡ 5029a214-0841-40fb-b397-4a2e1047bfb7
md"""
$(LocalResource("Pics/TimeDomainBanner.jpg"))
"""

# ╔═╡ 404060d3-23ec-400b-84cf-779e63b90293
md"""
# Baeysian view of the LS Periodogram
***

- What we want to be able to do is to detect variability and measure the period in the face of both noisy and incomplete data. Instead we'll use Fourier decomposition to get a more useful tool for actual data analysis.

- For a periodic signal we have:

$$y(t+P)=y(t),$$ where $P$ is the period.

- We can create a *phased light curve* that plots the data as function of phase:
$$\phi=\frac{t}{P} − {\rm int}\left(\frac{t}{P}\right),$$

- where ${\rm int}(x)$ returns the integer part of $x$.
"""

# ╔═╡ 72cf6fcf-2e2f-4dee-994b-ce2bfef51802
md"""
### A Single Sinusoid
***

- Let's take the case where the data are drawn from a single sinusoidal signal:

$$y(t)=A \sin(\omega t+\phi)+\epsilon$$

- and determine whether or not the data are indeed consistent with periodic variability and, if so, what is the period.

- This model is **non-linear** in the frequency term, $\omega$ and the phase, $\phi$ and therefore We rewrite the argument as $\omega(t−t_0)$ (reexpressing the phase term) and use trigonometrics identies to rewrite the model as:

$$y(t)=a \sin(\omega t)+b \cos(\omega t)$$

- where

$$A=(a^2+b^2)^{1/2} \text{ and } \phi=\tan^{−1}(b/a)$$

- The model is now linear with respect to coefficients $a$ and $b$ (and nonlinear only with respect to frequency, $\omega$).
"""

# ╔═╡ 3941f9fb-5b32-4152-8e3f-65767e63e055
md"""
- Assuming constant uncertainties on the data, we can write a likelihood function down:

$$L =\prod^N_{j=1}\frac{1}{\sqrt{2\pi}\sigma} \exp \left(\frac{−[y_j−a \sin(\omega t_j)−b \cos(\omega t_j)]^2}{2\sigma^2} \right) $$

- where $y_i$ is the measurement (e.g., the brightness of a star) taken at time $t_i$.

- With a lof of math we do not report here, and assuming uniform priors on $a, b, \omega$, and $\sigma$ (which gives nonuniform priors on $A$ and $\phi$), the posterior distribution of parameters can be simplified to:

$$p(\omega,a,b,\sigma|{t,y}) \propto \sigma^{−N} \exp \left(\frac{−NQ}{2\sigma^2} \right)$$

"""

# ╔═╡ 37e1a8cb-742a-4199-a534-34d4976e0f34
cm"""
- with

```math
Q= V - {2\over N} \left[ a \, I(\omega) + b \, R(\omega) - a\, b\, M(\omega) - {1 \over 2} a^2 \, S(\omega) - {1 \over 2} b^2 \,C(\omega)\right]
```

- and

```math
V = {1\over N} \sum_{j=1}^N y_j^2
```

```math
I(\omega) = \sum_{j=1}^N y_j   \sin(\omega t_j)
```

```math
R(\omega) = \sum_{j=1}^N y_j  \cos(\omega t_j)
```

```math
M(\omega) = \sum_{j=1}^N \sin(\omega t_j) \, \cos(\omega t_j)
```

```math
S(\omega) = \sum_{j=1}^N \sin^2(\omega t_j)
```

```math
C(\omega) = \sum_{j=1}^N  \cos^2(\omega t_j)
```

- *Note that I, R, M, S, C only depend on ``\omega`` and the data*.
"""

# ╔═╡ febb4a7e-792a-4d98-92d5-09ec53091be1
md"""
- If $N>>1$ and we have data that extends longer than the period:

$$S(\omega) \approx C(\omega) \approx N/2$ and $M(\omega) \ll N/2$$

- and

$$Q \approx V - {2\over N} \left[ a \, I(\omega) + b \, R(\omega)\right]  + {1 \over 2} (a^2 + b^2)$$
"""

# ╔═╡ c03057ef-1d71-46fa-b1af-f84347a95cda
cm"""
### The posterior for many, randomly spaced, observations
***

- If we marginalize over ``a`` and ``b`` (as we are interested in the period):

```math
p(\omega,\sigma|\{t,y\}) \propto  \sigma^{-(N-2)} \exp \left( { - N V \over 2 \sigma^2} + { P(\omega) \over \sigma^2}       \right)
```

- with

```math
P(\omega) = {1 \over N} [ I^2(\omega) + R^2(\omega)]
```

```math
V = {1\over N} \sum_{j=1}^N y_j^2
```

```math
I(\omega) = \sum_{j=1}^N y_j   \sin(\omega t_j)
```

```math
R(\omega) = \sum_{j=1}^N y_j  \cos(\omega t_j)
```
"""

# ╔═╡ fc6786d5-7d35-46e5-b904-0e66b087708a
md"""
-  we know the noise $\sigma$ then

```math
p(\omega|\{t,y\}, \sigma) \propto \exp \left( { P(\omega) \over \sigma^2} \right)
```

- and we now have the posterior for $\omega$!
"""

# ╔═╡ 0eb79796-8eb4-42c3-86ab-2ee3fbdbf21c
cm"""
## Significance of the peaks in the periodogram
***

- Let's compute the ``\chi^2`` for the LS periodogram:

```math
\chi^2(\omega) \equiv {1 \over \sigma^2} \sum_{j=1}^N [y_j-y(t_j)]^2 =
  {1 \over \sigma^2} \sum_{j=1}^N [y_j- a_0\, \sin(\omega t_j) - b_0 \, \cos(\omega t_j)]^2
```
  
- which we can simplify to:

```math
\chi^2(\omega) =  \chi_0^2 \, \left[1 - {2 \over N \, V}  \, P(\omega) \right]
```

- where, again, ``P(\omega)`` is the periodogram and ``\chi_0^2`` is the ``\chi^2`` for a model with ``y(t)``=constant:

```math
\chi_0^2 = {1 \over \sigma^2} \sum_{j=1}^N y_j^2 = {N \, V \over \sigma^2}
```
"""

# ╔═╡ 61307ab8-8148-4459-9a34-fd2e0af112a8
md"""
- We'll now renormalise the periodogram as:

$$P_{\rm LS}(\omega) = \frac{2}{N V} P(\omega),$$  

- where $0 \le P_{\rm LS}(\omega) \le 1$.

- With this renormalization, the reduction in $\chi^2(\omega)$ for the harmonic model, relative to $\chi^2$ for the pure noise model, $\chi^2_0$ is:

$${\chi^2(\omega) \over \chi^2_0}=  1 - P_{LS}(\omega).$$

- To determine if our source is variable or not, we first compute $P_{\rm LS}(\omega)$ and then model the odds ratio for our variability model vs. a no-variability model.

- If our variability model is "correct", then the peak of $P(\omega)$ gives the best $\omega$ and the $\chi^2$ at $\omega = \omega_0$ is $N$.
"""

# ╔═╡ 41b8576d-21fc-4683-a379-8a6bfb835a9a
md"""
- If the true frequency is $\omega_0$ then the maximum peak in the periodogram should have a height:

```math
P(\omega_0) = {N \over 4} (a_0^2 + b_0^2)
```

- and standard deviation:
  
```math
\sigma_P(\omega_0)  = {\sqrt{2} \over 2} \, \sigma^2.
```
"""

# ╔═╡ 232d90be-8b74-4158-a32a-a69d5122fc80
md"""
### Credits
***

This notebook contains material obtained from [https://github.com/gnarayan/ast596_2023_Spring](https://github.com/gnarayan/ast596_2023_Spring).
"""

# ╔═╡ b36fd613-95c8-44bf-876d-4eb345c26f08
cm"""
## Course Flow

<table>
  <tr>
	<td></td>
    <td>Previous lecture</td>
    <td>Next lecture</td>
  </tr>
  <tr>
	<td>notebook</td>
    <td><a href="./open?path=Lectures/Lecture-Lomb-Scargle/Lecture-Lomb-Scargle.jl">Irregular sampling</a></td>
    <td><a href="./open?path=Lectures/Lecture-Lomb-Scargle/Lecture-Lomb-Scargle.jl">Irregular sampling</a></td>
  </tr>
  <tr>
	<td>html</td>
    <td><a href="../../Lectures/Lecture-Lomb-Scargle/Lecture-Lomb-Scargle.html">Irregular sampling</a></td>
    <td><a href="../../Lectures/Lecture-Lomb-Scargle/Lecture-Lomb-Scargle.html">Irregular sampling</a></td>
  </tr>
 </table>


"""

# ╔═╡ 206474b8-0811-4785-8a71-acdcfd20b76c
md"""
**Copyright**

This notebook is provided as [Open Educational Resource](https://en.wikipedia.org/wiki/Open_educational_resources). Feel free to use the notebook for your own purposes. The text is licensed under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/), the code of the examples, unless obtained from other properly quoted sources, under the [MIT license](https://opensource.org/licenses/MIT). Please attribute the work as follows: *Stefano Covino, Time Domain Astrophysics - Lecture notes featuring computational examples, 2026*.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CommonMark = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
CommonMark = "~0.8.15"
PlutoUI = "~0.7.61"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.5"
manifest_format = "2.0"
project_hash = "5eaed3ce1c5a576233ab4f3159d6c7feab12cee4"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CommonMark]]
deps = ["Crayons", "PrecompileTools"]
git-tree-sha1 = "3faae67b8899797592335832fccf4b3c80bb04fa"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.15"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "1833212fd6f580c20d4291da9c1b4e8a655b128e"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.0.0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "7e71a55b87222942f0f9337be62e26b1f103d3e4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.61"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╟─4d477519-c44f-434c-b7e0-8daaa5009358
# ╟─ec67de24-d88a-46fa-ae46-c0cd7b797adc
# ╟─6a1315d1-9a6d-4ce0-b1c0-3fe22beb1ec2
# ╟─3ddd0f61-79d0-473c-8fec-a0e0c3fc72bf
# ╟─5029a214-0841-40fb-b397-4a2e1047bfb7
# ╟─404060d3-23ec-400b-84cf-779e63b90293
# ╟─72cf6fcf-2e2f-4dee-994b-ce2bfef51802
# ╟─3941f9fb-5b32-4152-8e3f-65767e63e055
# ╟─37e1a8cb-742a-4199-a534-34d4976e0f34
# ╟─febb4a7e-792a-4d98-92d5-09ec53091be1
# ╟─c03057ef-1d71-46fa-b1af-f84347a95cda
# ╟─fc6786d5-7d35-46e5-b904-0e66b087708a
# ╟─0eb79796-8eb4-42c3-86ab-2ee3fbdbf21c
# ╟─61307ab8-8148-4459-9a34-fd2e0af112a8
# ╟─41b8576d-21fc-4683-a379-8a6bfb835a9a
# ╟─232d90be-8b74-4158-a32a-a69d5122fc80
# ╟─b36fd613-95c8-44bf-876d-4eb345c26f08
# ╟─206474b8-0811-4785-8a71-acdcfd20b76c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
