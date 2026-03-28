### A Pluto.jl notebook ###
# v0.20.24

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
# Mutual Information and Time-Series Analysis
***

> We see here an algorithm introduced by [Huijse et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018ApJS..236...12H/abstract), to identify light curve periods based on quadratic mutual information (MI). We invite any intersted reader to check the original publication for any further detail.

- This is a non-parametric method since it does not rely on sinusoidal models for the data. 

- Instead, a metric on the phase diagram of the light curve ``\{\phi_i, m_i\}_{i=1,\ldots,N}`` is optimized, where ``m_i``, magnitudes, and ``\phi_i``, phases, are obtained from the time instants ``t_i`` given a certain trial period ``P`` as:

```math
\phi_i = \frac{\text{mod}(t_i, P)}{P} ~ \in ~ [0, 1], 
```

- where ``\text{mod}(\cdot, \cdot)`` stands for the division remainder operator.





- We'll make an extensive use of the information theoretic concept of Mutual Information (MI). 
    - In a broad sense, MI measures the reduction of the uncertainty of a random variable (RV) given that we know a second RV. MI can also be seen as a measure of dependence although, unlike correlation, MI is able to capture non-linear dependence between RVs. 
    - More formally, MI is posed as the divergence (statistical distance) between the joint probability density function (PDF) of the RVs and the product of their marginal PDFs.

- Several definitions of MI exist in the literature, being Shannon's MI the most well known ([Gray 2023](https://ee.stanford.edu/~gray/it.pdf)). 
- Shannon's MI for continuous RVs ``X`` and ``Y`` with joint PDF ``f_{X, Y}(\cdot, \cdot)`` is defined as: 

```math
\text{MI}_S(X, Y) = D_{KL}(f_{X,Y} || f_X f_Y) = 
```
```math 
= \iint f_{X,Y} \log f_{X,Y} \,dx \,dy - \int f_{X} \log f_X  \,dx  - \int f_{Y} \log f_Y \,dy 
```

- where ``D_{KL}(\cdot || \cdot)`` is the Kullback-Leibler divergence and ``f_X (x)= \int f_{X,Y} (x,y)\,dy``, ``f_Y (y) = \int f_{X,Y} (x, y)\,dx`` are the marginal PDFs of ``X`` and ``Y``, respectively.




- We intend to avoid the estimation of the PDF by using MI definitions arising from generalized divergences. Such MI estimators have been proposed in the Information Theoretic Learning (ITL, e.g., [Principe 2000](https://link.springer.com/article/10.1023/A:1008143417156)) literature. 
- In what follows we present the derivation of two MI definitions for continuous RVs from the ITL framework. 

- Starting from the Euclidean distance between probability density functions:

```math
D_{ED}(f(x) || g(x)) = \int (f(x) - g(x))^2 \,dx,  
```

- the Euclidean distance Quadratic MI between RVs ``X`` and ``Y`` is defined as: 

```math
\text{QMI}_{ED}(X, Y) = D_{ED}(f_{X,Y} (x,y)|| f_{X}(x) f_{Y}(y)) = 
```
```math
= \iint  f_{X,Y}^2 \,dx \,dy - 2 \iint f_{X,Y} f_X f_Y \,dx \,dy +  \int f_X^2 \,dx \int f_Y ^2 \,dy 
```
```math
= V_J - 2 V_C + V_M 
```

- where ``f_{X,Y}(\cdot, \cdot)`` is the joint PDF of ``X`` and ``Y`` while ``f_X(\cdot)`` and ``f_Y(\cdot)`` are the marginal PDFs, respectively.

- The terms ``V_J``, ``V_M`` and ``V_C`` correspond to the integrals of the squared joint PDF, squared product of the marginal PDFs and product of joint PDF and marginal PDFs, respectively. 

> The ITL framework provides an estimator of these quantities that can be computed directly from data samples. This estimator is called the Information Potential (IP) of an RV and it corresponds to the expected value of its PDF. 

- Note that the expected value of a PDF is equivalent to the integral of the squared PDF.





- In ITL a strong emphasis is given to the estimation of these quantities directly from data in a non-parametric way. 

- As an example consider the ITL estimation of the so-called [Renyi](https://en.wikipedia.org/wiki/R%C3%A9nyi_entropy)'s second order generalization ``H_{2} (X)`` of Shannon's entropy of a continuous RV: 

```math
H_{2} (X) = - \log \int f_X(x)^2 \,dx, 
```

- where ``f_X(x)`` is the RV's PDF. 

- Assuming that we have ``\{x_i\}_{i=1,\dots,N}`` realizations of the RV its PDF can be computed using a kernel density estimator (KDE):

```math
f_X(x) = \frac{1}{N} \sum_{i=1}^N \text{G}_h \left( x-x_i\right) = \frac{1}{N\sqrt{2\pi}h} \sum_{i=1}^N \exp \left( \frac{\|x-x_i\|^2}{2h^2} \right), 
```

- where ``\text{G}_h(\cdot)`` is the Gaussian kernel with bandwidth ``h``.

- Using the Gaussian convolution property, i.e. the convolution of two Gaussian functions is also a Gaussian, we obtain:

```math
H_{2} (X) = - \log \frac{1}{N^2} \int  \sum_{i=1}^N \sum_{j=1}^N  \text{G}_h \left( x-x_i \right)  \text{G}_h \left(x-x_j \right) \,dx =
```
```math
= - \log \frac{1}{N^2}  \sum_{i=1}^N \sum_{j=1}^N  \text{G}_{\sqrt{2}h} \left( x_i-x_j \right) = - \log \text{IP}_X, 
```

- where ``\text{IP}_X`` is the Information Potential (IP), an estimator of the expected value of the PDF of ``X`` estimated directly from the data samples bypassing the estimation of the PDF. 




- Assuming that we have ``\{x_i, y_i\}_{i=1,\ldots,N}`` i.i.d. realizations of RVs ``X`` and ``Y`` and using the IP estimator we get:
```math
V_M = \text{IP}_X \text{IP}_Y = \left (\frac{1}{N^2} \sum_{i,j=1}^{N,N}  \text{G}_{\sqrt{2}h} \left( x_i-x_j \right) \right) \left ( \frac{1}{N^2} \sum_{i,j=1}^{N,N} \text{G}_{\sqrt{2}h} \left( y_i-y_j \right) \right), 
```

```math
V_J = \text{IP}_{X,Y} = \frac{1}{N^2} \sum_{i=1}^{N} \sum_{j=1}^{N} \text{G}_{\sqrt{2}h} \left( x_i-x_j \right)  \text{G}_{\sqrt{2}h} \left( y_i-y_j \right), 
```

```math
V_C = \text{IP}_{X\times Y} = \frac{1}{N} \sum_{i=1}^{N}  \left ( \frac{1}{N}\sum_{j=1}^{N} \text{G}_{\sqrt{2}h} \left( x_i-x_j \right) \right) \left ( \frac{1}{N} \sum_{j=1}^{N} \text{G}_{\sqrt{2}h} \left( y_i-y_j \right) \right), 
```

- where ``\text{G}_{h} \left( x \right) = \frac{1}{\sqrt{2\pi}h} \exp \left( \frac{\|x\|^2}{2h^2} \right)`` is the Gaussian kernel with bandwidth ``h``. 




- The second ITL quadratic MI that we consider is obtained by defining a divergence measure based on the Cauchy-Schwarz inequality:

```math
D_{CS}(f(x) || g(x)) = -\log \frac{\left(\int f(x)g(x) \,dx\right)^2}{\int f(x)^2 \,dx \int g(x)^2 \,dx}, 
```

- then the Cauchy-Schwarz Quadratic MI for continuous RVs ``X`` and ``Y`` becomes:

```math
\text{QMI}_{CS}(X, Y) = D_{CS}(f_{X,Y} (x,y)|| f_{X}(x) f_{Y}(y)) = 
```
```math
= \log \iint f_{X,Y}^2 \,dx \,dy -2 \log \iint  f_{X,Y} f_X f_Y \,dx \,dy + 
```
```math
+ \log \int f_X^2 \,dx \int f_Y ^2 \,dy = \log V_J - 2 \log V_C + \log V_M, 
```

"""

# ╔═╡ 72cf6fcf-2e2f-4dee-994b-ce2bfef51802
md"""
## Period Estimation by Maximizing Mutual Information
***

- A typical analysis starts by applying the epoch folding transformation for a certain trial period to the (possibly) unevenly sampled time-series to obtain the phase diagram ``\{\phi_i, m_i, \sigma_i\}_{i=1,\ldots,N}``. 
- We assume that the light curve is periodic with an unknown period. The phases ``\{\phi_i\}`` correspond to our non-parametric model of the periodicity, while ``\{m_i\}`` correspond to our noisy observations. As usual ``\{\sigma_i\}`` are the estimated errors on our observations. 
- If the light curve is periodic with period ``P_T``, then folding with this period will yield the model that best explains our observations. This can be measured by calculating the MI between phases and magnitudes, i.e. the amount of information shared by model and observation. 
- We can test several models (foldings) and find the one which maximizes MI to detect the best period. 
    - Notice that MI requires independent and identically distributed (iid) realizations of the RVs. Although light curves are time series and hence there exist serial correlations in time, these correlations are broken in the phase diagram.
    - Phase is a function of time and period, and several periods are tested per light curve. If the period is not related to the underlying periodicity of the data the phase diagram is filled uniformly and serial correlations in the joint space are broken. 
 
- A second interpretation on using MI for periodicity detection lays on MI's definition as the divergence (statistical distance) between the joint PDF and the marginal PDF of the RVs. 
    - If the light curve is folded with a wrong period, the structure in the joint PDF will be almost equal to the product of the marginal PDFs, i.e. magnitudes (or fluxes) are independent of the phases. 
    - On the other hand, if the correct period is chosen the joint PDF will present structure that is not captured by the product of the marginals. By maximizing MI we are maximizing the dependency between model and observations.



- Let's denote ``M`` and ``\Phi`` as the RVs associated to magnitude and phase, respectively. We can estimate the PDF of ``M`` given its realizations using KDE as follows:

```math
f_M(m) = \frac{1}{N} \sum_{i=1}^N \text{G}_{\sqrt{\sigma_i^2+h_m^2}}(m-m_i) = \frac{1}{N} \sum_{i=1}^N \frac{1}{\sqrt{2 \pi (\sigma_i^2 + h_m^2)}} \exp \left( - \frac{1}{2} \frac{(m - m_i )^2}{(\sigma_i^2 + h_m^2)} \right), 
```
- where each sample ``m_i`` has a bandwidth that incorporates the KDE bandwidth ``h_m`` and its given uncertainty ``\sigma_i``. 

- As ``\Phi`` is a periodic RV we need a periodic kernel to appropriately estimate its PDF. We consider a kernel arising from the Wrapped Cauchy (WC) distribution \citep{jammalamadaka2001topics} and estimate ``\Phi``'s PDF as:

```math
f_\Phi(\phi) = \frac{1}{N} \sum_{i=1}^N \text{WC}_{h_\phi}(\phi-\phi_i)  = \frac{1}{2 \pi N} \sum_{i=1}^N \frac{1 - e^{-2 h_\phi}}{1 + e^{-2 h_\phi} - 2 e^{- h_\phi} \cos(2\pi (\phi - \phi_i))}, 
```

- where ``h_\phi \in (0, \infty)`` is the scale of the Cauchy distribution. 

- For ``h_\phi \to \infty`` the WC kernel behaves like the circular uniform distribution, while for ``h_\phi \to 0`` it concentrates on its mean. 

- The joint PDF of ``\Phi`` and ``M`` is estimated as:

```math
f_{\Phi, M}(\phi, m) = \frac{1}{N} \sum_{i=1}^N \text{G}_{\sqrt{\sigma_i^2+h_m^2}}(m-m_i) \cdot \text{WC}_{h_\phi}(\phi-\phi_i), 
```

- because the multiplication of valid kernel functions is also a kernel.




- Using the Gaussian kernel for the magnitudes (or fluxes) and the WC kernel for phases we obtain:
```math
\text{IP}_M = \frac{1}{N^2} \sum_{i=1}^N \sum_{j=1}^{N}  \text{G}_{\sqrt{2h_m^2 + \sigma_i^2+ \sigma_j^2}} \left( m_i-m_j \right), 
```
```math
\text{IP}_\Phi = \frac{1}{N^2} \sum_{i=1}^N \sum_{j=1}^{N}  \text{WC}_{2 h_\phi} \left( \phi_i-\phi_j \right), 
```
```math   
\text{IP}_{\Phi,M} = \frac{1}{N^2} \sum_{i=1}^N \sum_{j=1}^{N}  \text{G}_{\sqrt{2h_m^2 + \sigma_i^2+ \sigma_j^2}} \left( m_i-m_j \right) \text{WC}_{2 h_\phi} \left( \phi_i-\phi_j \right), 
```

- and therefore:
```math
\text{IP}_{\Phi \times M} = \frac{1}{N} \sum_{i=1}^{N}  \left ( \frac{1}{N}\sum_{j=1}^{N}  \text{G}_{\sqrt{2h_m^2 + \sigma_i^2+ \sigma_j^2}} \left( m_i-m_j \right) \right) \left ( \frac{1}{N} \sum_{j=1}^{N} \text{WC}_{2 h_\phi} \left( \phi_i-\phi_j \right) \right), 
```



- Through these potentials we restate the QMI estimators as:
```math
\text{QMI}_{ED}(\Phi, M) =   \text{IP}_{\Phi,M} - 2 \text{IP}_{\Phi \times M} + \text{IP}_\Phi \text{IP}_M, 
```
```math
\text{QMI}_{CS}(\Phi, M) =   \log \text{IP}_{\Phi,M} - 2 \log \text{IP}_{\Phi \times M} + \log \text{IP}_\Phi  + \log \text{IP}_M, 
```

> If the period of a light curve is estimated by maximizing the QMI for a range of trial periods. This yields a QMI periodogram!



- An interesting problem by itself is the choice of the KDE bandwidth. In this case we have two parameters ``h_\phi`` and ``h_m``. 
    - The former is associated to the phases which are always constrained to ``[0, 2\pi]``, i.e. the dynamic range of this variable is fixed. QMI is not too sensitive to ``h_\phi`` as long as it is not extremely small or large. It was found empirically that ``h_\phi = 1`` is a good choice and we keep it constant to make comparisons between QMI values easier. 
    - The second bandwidth ``h_m`` is more difficult to set as the dynamic range of the magnitudes is not known \emph{a priori}. Following [Silverman (1986)](https://www.taylorfrancis.com/books/mono/10.1201/9781315140919/density-estimation-statistics-data-analysis-bernard-silverman), we may write:
```math
h_m = 0.9 \cdot \text{min} ( \sqrt{\text{VAR}[m]}, ~\text{IQR}[m]/1.349) \cdot N^{-1/5}, 
```
- where ``\text{VAR}[m]`` is the variance of the magnitudes, ``\text{IQR}[m]`` is the interquartile range of the magnitudes and ``N`` is the number of samples. To avoid overestimation of ``h_m`` we use the weighted versions of variance and IQR, with weights ``w_i = \sigma_i^{-2}``, ``i = 1,\ldots,N``. o frequency, $\omega$).
"""

# ╔═╡ 2d596c28-74bc-4ff8-a030-fbac18dbceb0
md"""
## Reference & Material

Material and papers related to the topics discussed in this lecture.

- [Huijse et al. (2018) - Robust Period Estimation Using Mutual Information for Multiband Light Curves in the Synoptic Survey](https://ui.adsabs.harvard.edu/abs/2018ApJS..236...12H/abstract)
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
    <td><a href="./open?path=Lectures/Lecture-NonParametricAnalysis/Lecture-NonParametricAnalysis.jl">Lecture about non-parametric analysis</a></td>
    <td><a href="./open?path=Lectures/Lecture-NonParametricAnalysis/Lecture-NonParametricAnalysis.jl">Lecture about non-parametric analysis</a></td>
  </tr>
  <tr>
	<td>html</td>
    <td><a href="Lectures/Lecture-NonParametricAnalysis/Lecture-NnParametricAnalysis.html">Lecture about non-parametric analysis</a></td>
<td><a href="Lectures/Lecture-NonParametricAnalysis/Lecture-NnParametricAnalysis.html">Lecture about non-parametric analysis</a></td>
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
# ╟─2d596c28-74bc-4ff8-a030-fbac18dbceb0
# ╟─b36fd613-95c8-44bf-876d-4eb345c26f08
# ╟─206474b8-0811-4785-8a71-acdcfd20b76c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
