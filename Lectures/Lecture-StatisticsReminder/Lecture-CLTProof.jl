### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 0ea989ec-dce0-4802-b117-413c6ee3e969
begin
	using CommonMark
	using PlutoUI
end

# ╔═╡ 3926f03f-e9f3-4bb3-a3b7-5f0541d65d09
md"""
**What is this?**


*This jupyter notebook is part of a collection of notebooks on various topics discussed during the Time Domain Astrophysics course delivered by Stefano Covino at the [Università dell'Insubria](https://www.uninsubria.eu/) in Como (Italy). Please direct questions and suggestions to [stefano.covino@inaf.it](mailto:stefano.covino@inaf.it).*
"""

# ╔═╡ 4968fd4d-b86b-4071-9b18-07aa8b8b0738
md"""
**This is a `pluto` notebook**
"""

# ╔═╡ a9b1d316-dbcf-464d-83ab-fff7ac01f6b1
TableOfContents()

# ╔═╡ da8b93b1-c1b1-4d2d-bd93-00b5ea97dbf6
md"""
$(LocalResource("Pics/TimeDomainBanner.jpg"))
"""

# ╔═╡ 5a7eb5d1-ec5f-48e0-ac1d-e84322d11767
md"""
# A proof of the central limit theorem
***
"""

# ╔═╡ 009d6216-81b2-4823-bcb0-ba662d8c6a73
md"""
- Let $X_1$, $X_2$,…,$X_N$ be i.i.d. random variables that form a random sample of size ’N’. Assume that we have drawn this sample from a population that has a mean $μ$ and variance $σ²$.

- Let $\bar{X}_N$ be the sample mean: $\bar{X}_N = \frac{X_1+X_2...+X_N}{N}$

- Let $\bar{Z}_N$ be the standardized sample mean: $\bar{Z}_N = \frac{\bar{X}_N-\mu}{\sigma/\sqrt{N}}$

- The Central Limit Theorem states that as N tends to infinity, $\bar{Z}_N$ *converges in distribution* to $N(0,1)$, i.e. the CDF of $\bar{Z}_N$ becomes identical to the CDF of $N(0, 1)$. 

- To prove this statement, we use the property of the Moment Generating Function (MGF) that if the MGFs of $X$ and $Y$ are identical, then so are their CDFs.
"""

# ╔═╡ 979a339d-0c70-468e-aa53-8867244016e4
md"""
### Moment Generating Function
***

- The k-th moment of a random variable $X$ is the expected value of $X$ raised to the k-th power: $\mu_k(X) = \mathbb{E}(X^k) = \sum_i x_i^k P(X=x_i)$

- The k-th moment of $X$ around some value $c$ is known as the k-th central moment of $X$: $\mu_k(X) = \mathbb{E}((X-c)^k) = \sum_i (x_i-c)^k P(X=x_i)$

- The k-th standardized moment of $X$ is the k-th central moment of $X$ divided by k-th power of the standard deviation of $X$: $\frac{\mu_k(X)}{\sigma^k} = \frac{\mathbb{E}((X-c)^k)}{\sigma^k}$

- For the record, the first 5 moments of $X$ have specific values or meanings attached to them:

    - The zeroth’s raw and central moments of X are $\mathbb{E}(X^0)$ and $\mathbb{E}[(X — c)^0]$ respectively. Both equate to 1.
    - The 1st raw moment of $X$ is $\mathbb{E}(X)$. It’s the mean of $X$.
    - The second central moment of $X$ around its mean is $\mathbb{E}[X — \mathbb{E}(X)]^2$. It’s the variance of X.
    - The third and fourth standardized moments of $X$ are $\mathbb{E}[X — \mathbb{E}(X)]^3/σ^3$, and $\mathbb{E}[X — \mathbb{E}(X)]^4/σ^4$. They are the skewness and kurtosis of $X$ respectively. 
    
- Let's now define a new random variable $tX$ where $t$ is a real number. Here’s the Taylor series expansion of $e$ to the power $tX$ evaluated at $t = 0$:

```math
 e^{tX} = \sum_{k=0}^\infty \frac{(tX)^k}{k!} = 1 + \frac{t}{1!}X + \frac{t^2}{2!}X^2 ... 
```

- If we apply the *Expectation operator* on both sides of the above equation we get from the linearity of the operator:

```math
\mathbb{E}(e^{tX}) = \sum_{k=0}^\infty \mathbb{E}(\frac{(tX)^k}{k!}) = \sum_{k=0}^\infty \frac{t^k}{k!} \mathbb{E}(X^k) 
```

- If we now write the general form of an *Exponential Generating Function* (EGF):

```math
\sum_{k=0}^N \frac{a_k}{n!} x^k = a_0 + a_1 x + \frac{a_2}{2!} x^2 + ... 
```

- We see that $\mathbb{E}(X^k)$ are the coefficients $a_k$ in the EGF. Thus we have the *Moment Generating Function* (MGF):

```math
M_X(t) = \mathbb{E}(e^{tX}) = \sum_{k=0}^\infty \frac{t^k}{k!} \mathbb{E}(X^k) = 1 + t\mathbb{E}(X) + \frac{t^2}{2!}\mathbb{E}(X^2) + ... 
```

- Now, it is easy to realize that the k-th derivative of the EGF, when evaluated at $x = 0$ gives us the k-th coefficient of the underlying sequence. So, $M_X^0(t=0) = 1, M_X^1(t=0) = \mathbb{E}(X), ..., M_X^k(t=0) = \mathbb{E}(X^k)$

- If two random variables $X$, $Y$ have identical moments (i.e. identic MGF) they must unavoidably have identical CDF.

- We can say more. If $Y = aX + b$, then $M_Y(t) = \mathbb{E}(e^{(aX+b)t}) = e^{bt}M_X(t)$.

- Besides, if $Y$ the sum of $N$ independent, identically distributed random variables, $Y = X_1 + X_2 + ...$ then $M_Y(t) = \mathbb{E}(e^{t(X_1+X_2+..)}) = \mathbb{E}(e^{tX_1})\mathbb{E}(e^{tX_2})...$, due to the independency. And, finally $M_Y(t) = [\mathbb{E}(e^{tX})]^N = [M_X(t)]^N$ since they are identically distributed.

- One last useful result related to MGF is that if $X \sim \mathcal{N}(0,1)$ then $M_X(t) = e^{\frac{t^2}{2}}$, since the given distribution has mean = 0, variance = 1, skew = 0 and kurtosis = 1.
"""

# ╔═╡ 6f051bea-173a-4fff-a4ed-df724cc33a1a
md"""
- Coming back to the proof of the CLT, we now understand that we need to prove that the MGF of $\bar{Z}_N$ converges to the MGF of $\mathcal{N}(0,1)$, i.e.:

```math
\lim_{N\to\infty} M_{\bar{Z}_N} (t) = e^{\frac{t^2}{2}}
```

- Let's define $Z_k = \frac{X_k - \mu}{\sigma}$. Then, $\bar{Z}_N = \frac{1}{\sqrt N}\sum_{k=1}^N Z_k$. We may write:

```math
M_{\bar{Z}_N} (t) = M_{\frac{1}{\sqrt N}\sum_{k=1}^N Z_k}(t) = M_{\frac{Z_1}{\sqrt(N)}+\frac{Z_2}{\sqrt(N)}+ ...}(t)
```

- The various $Z_k/sqrt(N)$ are, by construction, independent random variable, so that:

```math
M_{\bar{Z}_N} (t) = [M_{Z/\sqrt{N}} (t)]^N = [M_{Z} (t/\sqrt{N})]^N
```
"""

# ╔═╡ 0fdbce17-cf81-414a-a018-9197df77db76
md"""
- Now, let's create a Taylor series expansion of $M_Z(t/\sqrt{N})$ at $t = 0$:

```math
M_{Z} (t/\sqrt{N}) = \sum_{k=0}^\infty M_Z^k(0) \frac{(t/\sqrt{N})^k}{k!}
```

- Next, we split this expansion into two parts. The first part is a finite series of three terms corresponding to $k = 0, k = 1$, and $k = 2$. The second part is the remainder of the infinite series.

- $M^0$, $M^1$, $M^2$, etc. are the 0-th, 1st, 2nd, and so on derivatives of the MGF $M_Z(t/\sqrt{N})$ evaluated at ($t/\sqrt{N}) = 0$. We know that these derivatives of the MGF are the 0-th, 1st, 2nd, etc. moments of $Z$.

- The $0$-th moment, $M^0(0)$, is always 1. $Z$ is, by its construction, a standard normal random variable. Hence, its first moment (mean), $M^1(0) = 0$, and its second moment (variance), $M^2(0) = 1$. So:

```math
M_{Z} (t/\sqrt{N}) = 1 + 0 + \frac{(t/\sqrt{N})^2}{2!} + \sum_{k=3}^\infty M_Z^k(0) \frac{(t/\sqrt{N})^k}{k!}
```

- Now, let's summarize where we are:

```math
\lim_{N\to\infty} M_{\bar{Z}_N} (t) = \lim_{N\to\infty} \left(M_{Z} (t/\sqrt{N})\right)^N
```

```math
\lim_{N\to\infty} M_{\bar{Z}_N} (t) = \lim_{N\to\infty} \left( 1 + 0 + \frac{(t/\sqrt{N})^2}{2!} + R_2(t/\sqrt{N})\right)^N
```

- Where $R_r(x)$ is the remainder or residual from approximating a function using the Taylor polynomial of order $r$.

- Taylor's theporem states that $\lim_{x \to a} \frac{R_r(x)}{(x-a)^r} = 0$. In out case this turns out to be:

```math
\lim_{\frac{t}{\sqrt{N}} \to a} \frac{R_2(t/\sqrt{N})}{(t/\sqrt{N}-a)^2} = 0
```

- But, with $a=0$, this is equivalent to the limit:

```math
\lim_{N \to \infty} \frac{R_2(t/\sqrt{N})}{(t/\sqrt{N})^2} = 0
```

- In the above limit, the L.H.S. will tend to zero as long as $N$ tends to infinity independent of what value $t$ has as long as it’s finite. 

- Therefore, we can write:

```math
\lim_{N \to \infty} n R_2(t/\sqrt{N}) = 0
```

"""

# ╔═╡ 326de926-4b20-44a6-bf9f-3d7c7b813de3
md"""
- With little algebra, we have that:

```math
\lim_{N\to\infty} M_{\bar{Z}_N} (t) = \lim_{N\to\infty} \left( 1 + \frac{1}{N} [(t^2/2) + N R_2(t/\sqrt{N})]\right)^N
```

- And, therefore:

```math
\lim_{N\to\infty} M_{\bar{Z}_N} (t) = \lim_{N\to\infty} \left( 1 + \frac{(t^2/2)}{N} \right)^N
```

- Since $\lim_{N\to\infty} \left( 1 + \frac{x}{N} \right)^N = e^x$, we can finally write:

```math
\lim_{N\to\infty} M_{\bar{Z}_N} (t) = e^\frac{t^2}{2}
```

- And the CLT is proved!
"""

# ╔═╡ 85d1de10-a1ea-4991-8d81-65d12171ecdd
md"""
### Credits
***

This notebook contains material obtained by [https://towardsdatascience.com/a-proof-of-the-central-limit-theorem-8be40324da83](https://towardsdatascience.com/a-proof-of-the-central-limit-theorem-8be40324da83).
"""

# ╔═╡ 66d96dd5-2a8a-46c3-8897-0dbed4c69986
cm"""
## Course Flow
***

<table>
  <tr>
    <td>Previous lecture</td>
    <td>Next lecture</td>
  </tr>
  <tr>
      <td><a href="./open?path=Lectures/Lecture-StatisticsReminder/Lecture-StatisticsReminder.jl">Reminder of frequentist statistics</a></td>
    <td><a href="./open?path=Lectures/Lecture-StatisticsReminder/Lecture-StatisticsReminder.jl">Reminder of frequentist statistics</a></td>
  </tr>
 </table>


"""

# ╔═╡ a0bfbb2f-da37-4943-9791-1b4b3b534f9e
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
# ╟─3926f03f-e9f3-4bb3-a3b7-5f0541d65d09
# ╟─4968fd4d-b86b-4071-9b18-07aa8b8b0738
# ╟─0ea989ec-dce0-4802-b117-413c6ee3e969
# ╟─a9b1d316-dbcf-464d-83ab-fff7ac01f6b1
# ╟─da8b93b1-c1b1-4d2d-bd93-00b5ea97dbf6
# ╟─5a7eb5d1-ec5f-48e0-ac1d-e84322d11767
# ╟─009d6216-81b2-4823-bcb0-ba662d8c6a73
# ╟─979a339d-0c70-468e-aa53-8867244016e4
# ╟─6f051bea-173a-4fff-a4ed-df724cc33a1a
# ╟─0fdbce17-cf81-414a-a018-9197df77db76
# ╟─326de926-4b20-44a6-bf9f-3d7c7b813de3
# ╟─85d1de10-a1ea-4991-8d81-65d12171ecdd
# ╟─66d96dd5-2a8a-46c3-8897-0dbed4c69986
# ╟─a0bfbb2f-da37-4943-9791-1b4b3b534f9e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
