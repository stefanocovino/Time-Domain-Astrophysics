### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 3c2986d8-2379-4cd5-a747-435a5180076f
begin
	using CairoMakie
	using CommonMark
	using PlutoUI
end

# ╔═╡ 72a2571d-dc92-41d6-b73b-ba9e985ca4f1
md"""
**What is this?**


*This jupyter notebook is part of a collection of notebooks on various topics discussed during the Time Domain Astrophysics course delivered by Stefano Covino at the [Università dell'Insubria](https://www.uninsubria.eu/) in Como (Italy). Please direct questions and suggestions to [stefano.covino@inaf.it](mailto:stefano.covino@inaf.it).*
"""

# ╔═╡ 12200cd2-5e1e-4d96-9284-60fb950fd70a
md"""
**This is a `Julia` notebook**
"""

# ╔═╡ a81e5916-c0e7-410d-8c6e-917303ab087a
md"""
$(LocalResource("Pics/TimeDomainBanner.jpg"))
"""

# ╔═╡ 17903372-3749-4bfc-9519-7181bc0539f1
md"""
# The [Gregory-Loredo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract) algorithm

***

- We assume to have a set of $N$ arrival times, $D = {t_j}$, with $j = 1$ to $N$, over a some observing interval of duration $T$.

- We do not put any constrain on the arrival rate, $r(t)$.


"""

# ╔═╡ 0e4c318a-511c-4646-9b33-634e5cb87c6d
cm"""
## The Likelihood function
***

- The probability of ``D`` given some rate function, ``r(t)``, the likelihoof function is computed as follows:
     - The observing interval is divided in small time intervals ``\Delta t_i``, each containing either 1 or zero evewnts.
     - Adopting the Poissonian distribution, the probability of seeing ``n`` events in an interval ``\Delta t`` about time ``t`` is:

```math
p_n(t) = \frac{[r(t) \Delta t]^n e^{-r(t) \Delta t}}{n!}
```

- We also assume the rate can be considered constant within ``\Delta t``

- If ``N`` and ``Q`` are the number of time intervals with 1 or no events we can write:

```math
p(D|r,I) = \prod_{i=1}^N p_1(t_i)  \prod_{k=1}^Q p_0(t_k)
```

- For no event we have: ``p_0(t) = e^{-r(t) \Delta t}`` and for one event we have: ``p_1(t) = r(t) \Delta t e^{-r(t) \Delta t}``, so that:

```math
p(D|r,I) = \Delta t^N \prod_{i=1}^N r(t_i) e^{-\sum_{k=1}^{N+Q} r(t_k) \Delta t} = \Delta t^N \prod_{i=1}^N r(t_i) e^{\int_T r(t) dt}
```

- where the sum of ``r(t)\Delta t`` over all the observed intervals is expressed by the integral of the rate over the intervals with a range of integration ``T=(N+Q)\Delta t``.

- ``T`` is the total duration of the observtions and includes all the intervals. These intervals do not need to be contiguous, therefore allowing to remove intervals when, e.g., the detectors were not operational.

"""

# ╔═╡ 78e65245-8981-4bdf-863a-e8497bb8e2e6
cm"""

## Constant model
***

- The simplest model for the data has a constant event rate, ``A``. We denote this one-parameter model as ``M_1``.

- Setting ``r(t) = A`` we have:

```math 
p(D | A, M_1) = \Delta t^N A^N e^{-AT}
```

- We also assume a uniform prior for ``A`` from ``0`` to an upper limit ``A_{\rm max}``: ``~ p(A | M_1) = \frac{1}{A_{\rm max}}``. The dependence on ``A_{\rm max}`` turns out to be very weak. 

"""

# ╔═╡ c7fc2cd6-6715-46d7-a9e0-c6e0bb7b13eb
md"""
## Reference & Material

Material and papers related to the topics discussed in this lecture.


- [Gregory & Loredo (1992) - "A New Method for the Detection of a Periodic Signal of Unknown Shape and Period”](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract)
"""

# ╔═╡ aeaad396-5e2a-4ec5-8b05-e160e5712f35
cm"""
## Course Flow

<table>
  <tr>
    <td>Previous lecture</td>
    <td>Next lecture</td>
  </tr>
  <tr>
    <td><a href="./open?path=Lectures/Lecture - Wavelet Analysis/Lecture-Wavelet-Analysis.jl">Lecture about wavelet analysis</a></td>
    <td><a href="./open?path=Lectures/Science Case - FRBs/Lecture-FRBs.j">Science case about FRBs</a></td>
  </tr>
 </table>

"""

# ╔═╡ 7aecb81a-a226-4133-9182-53e159f5fc43
md"""
**Copyright**

This notebook is provided as [Open Educational Resource](https://en.wikipedia.org/wiki/Open_educational_resources). Feel free to use the notebook for your own purposes. The text is licensed under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/), the code of the examples, unless obtained from other properly quoted sources, under the [MIT license](https://opensource.org/licenses/MIT). Please attribute the work as follows: *Stefano Covino, Time Domain Astrophysics - Lecture notes featuring computational examples, 2026*.
"""

# ╔═╡ Cell order:
# ╟─72a2571d-dc92-41d6-b73b-ba9e985ca4f1
# ╟─12200cd2-5e1e-4d96-9284-60fb950fd70a
# ╠═3c2986d8-2379-4cd5-a747-435a5180076f
# ╟─a81e5916-c0e7-410d-8c6e-917303ab087a
# ╟─17903372-3749-4bfc-9519-7181bc0539f1
# ╟─0e4c318a-511c-4646-9b33-634e5cb87c6d
# ╟─78e65245-8981-4bdf-863a-e8497bb8e2e6
# ╟─c7fc2cd6-6715-46d7-a9e0-c6e0bb7b13eb
# ╟─aeaad396-5e2a-4ec5-8b05-e160e5712f35
# ╟─7aecb81a-a226-4133-9182-53e159f5fc43
