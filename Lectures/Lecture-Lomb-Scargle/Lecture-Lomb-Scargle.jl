### A Pluto.jl notebook ###
# v0.20.24

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

# ╔═╡ 93f1b40e-8990-405d-bba1-9f98abd054ca
begin
	using CairoMakie
	using CommonMark
	using CSV
	using DataFrames
	using FFTW
	using LinearAlgebra
	using LombScargle
	using PlutoUI
	using Random
	using Statistics
end

# ╔═╡ ff837227-c0bc-433c-b462-9e18ff3d947b
md"""
**What is this?**


*This jupyter notebook is part of a collection of notebooks on various topics discussed during the Time Domain Astrophysics course delivered by Stefano Covino at the [Università dell'Insubria](https://www.uninsubria.eu/) in Como (Italy). Please direct questions and suggestions to [stefano.covino@inaf.it](mailto:stefano.covino@inaf.it).*
"""

# ╔═╡ 53a01399-1f0b-4fbb-9966-97124e9eff42
md"""
**This is a `Pluto` notebook**
"""

# ╔═╡ 8f9cd30e-5598-4e1b-b379-1473232ea867
TableOfContents()

# ╔═╡ f7100fe5-2cfd-448d-abcc-66a623a8f11f
# ╠═╡ show_logs = false
md"""
$(LocalResource("Pics/TimeDomainBanner.jpg"))
"""

# ╔═╡ 31dd4ceb-8afd-4b55-aacb-2d888681ad1c
# ╠═╡ show_logs = false
md"""
# Unevenly sampled data

***

- Harmonic analysis of unevenly spaced data is problematic due to the loss of information and increase in aliasing.

- A (partial) solution to this problem was proposed in astronomy and it is now widely used and known as *Lomb-Scargle periodogram*.

- In case of regular sampling it is equivalent to the Fourier periodogram.

- An example of LS periodogram compute analying optical observation of the LINEAR object [11375941](https://simbad.u-strasbg.fr/simbad/sim-id?Ident=linear11375941&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id) is shown below:

$(LocalResource("Pics/linear11375941lc.png"))
$(LocalResource("Pics/linear11375941ls.png"))
"""

# ╔═╡ e6386bbc-ea18-49a6-b64a-51b627d31565
md"""
## Many possible techniques
***

- Fourier methods:
    - are based on the Fourier transform, power spectra and correlation functions (Schuster periodogram, Lomb-Scargle periodogram, correlation-based methods, wavelet transform, etc.)

- Phase-folding methods:
    - depend on folding observations as a function of phase (String Length, Analysis of Variance, Phase Dispersion Minimization, Gregory & Loredo method, conditional entropy method, correntropy methods, etc.)

- Least-squares methods:
    - involve fitting a model to the data at each candidate frequency (Lomb-Scargle periodogram, Supersmoother approach, orthogonal polynomial fits, etc.)

- Bayesian approaches:
    - apply Bayesian probability theory to the problem (Lomb-Scargle generalisation, Gregory & Loredo, Gaussian process models, stochastic process models, etc.)

"""

# ╔═╡ 477704cc-f022-4cc5-a4a4-679387b819bd
# ╠═╡ show_logs = false
cm"""
## A reminder of the CFT
***

- Given a continuous signal ``g(t)`` the Fourier transform and its inverse are:

```math
\hat{g}(f) \equiv \int_{-\infty}^\infty g(t) e^{-2\pi i f t} dt \quad  g(t) \equiv \int_{-\infty}^\infty \hat{g}(f) e^{+2\pi i f t} df
```

- We can also define the Fourier transform operator:

```math
\mathcal{F}\{g\} = \hat{g} \quad    \mathcal{F}^{-1}\{\hat{g}\} = g
```

- ``g`` and ``ĝ`` are know as a Fourier pair: ``g \Longleftrightarrow \hat{g}``.

- The Fourier Transform (FT) is a linear operator:

```math
\mathcal{F}\{f(t) + g(t)\} = \mathcal{F}\{f(t)\} + \mathcal{F}\{g(t)\}\nonumber \quad
       \mathcal{F}\{A f(t)\} = A\mathcal{F}\{f(t)\}
```

- The FT of a sinusoid with frequency ``f_0`` is a sum of delta functions at ``\pm f_0``, where ``\delta(f)\equiv\int_{-\infty}^\infty e^{-2\pi i x f}df``.

- We can write: `` \mathcal{F}\{e^{2\pi f_0 t}\} = \delta(f - f_0)``, and:

```math
\mathcal{F}\{\cos(2\pi f_0 t)\} = \frac{1}{2}\left[\delta(f - f_0) + \delta(f + f_0)\right] \quad
      \mathcal{F}\{\sin(2\pi f_0 t)\} = \frac{1}{2i}\left[\delta(f - f_0) - \delta(f + f_0)\right]
```

- Relations that can be derived from Euler’s formula: ``e^{ix} = \cos x + i\sin x``

- A time shift imparts a phase in the FT:  ``  \mathcal{F}\{g(t - t_0)\} = \mathcal{F}\{g(t)\} e^{-2\pi i ft_0}``.

- And, as we know, the squared amplitude of the FT of a continuous signal is known as the power spectral density (PSD):  `` \mathcal{P}_g \equiv \left|\mathcal{F}\{g\}\right|^2 ``.

    - Note that if ``g`` is real-valued, it follows that ``P_g`` is an even function, {i.e.} ``\mathcal{P}_g(f) = \mathcal{P}_g(-f)``.

$(LocalResource("Pics/FTpairs.png"))
"""

# ╔═╡ f95ec17d-9c56-4897-b5b2-65c3470fc294
# ╠═╡ show_logs = false
md"""
### The convolution theorem
***

- A convolution of two functions, usually denoted by the $\ast$ symbol, is defined as follows: $[f \ast g](t) \equiv \int_{-\infty}^\infty f(\tau)g(t - \tau) d\tau$

$(LocalResource("Pics/convolution.png"))

- It can be shown that the FT of a convolution is the point-wise product of the individual FTs:

```math\mathcal{F}\{f \ast g\} = \mathcal{F}\{f\} \cdot \mathcal{F}\{g\} \qquad   \mathcal{F}\{f \cdot g\} = \mathcal{F}\{f\} \ast \mathcal{F}\{g\}
```

"""

# ╔═╡ 0f01200a-665e-4ae2-a6b9-b5c192955f3f
# ╠═╡ show_logs = false
cm"""
### Effect of the observing window on a FT
***

- If ``g_{obs}(t) = g(t)W(t)``, by the convolution theorem we have: ``\mathcal{F}\{g_{obs}\} = \mathcal{F}\{g\} \ast \mathcal{F}\{W\}``.

$(LocalResource("Pics/observingwindow.png"))

- If we assume to observe a continuous periodic signal over a limited span of time the observed signal is the pointwise product of the underlying infinite periodic signal with a rectangular window function.

- By the convolution theorem, the FT is given by the convolution of the transform of the underlying function (here a set of delta functions at the component frequencies) and the transform of the window function (here a sinc function).

- For the purely periodic signal  this convolution has the effect of replacing each delta function with a sinc function.

- Because of the inverse relationship between the width of the window and the width of its transform, it follows that a wider observing window leads to proportionally less spread in the FT of the observed function.

| $(LocalResource("Pics/densesampling.png")) | $(LocalResource("Pics/sparsesampling.png")) |
| ----------------------------------------- | ------------------------------------------- |

- The observed FT is a convolution of the true transform (here a localized Gaussian) and the window transform (here another Dirac comb). Howver, if the sampling rate is "too low", the result is that the FT of the window function has spacing narrower than the FT of the signal.

    - The observed FT suffers from aliasing of signals, such that not all frequency information can be recovered.

- This is the famous Nyquist sampling theorem. If we have a regularly-sampled function with a sampling rate of ``f_0 = 1/T``, we can only fully recover the frequency information if the signal is *band-limited* between frequencies ``\pm f_0/2``.
"""

# ╔═╡ a24629ad-7302-4aff-b7bf-fcc7635b9691
md"""
### The discrete Fourier transform (DFT)
***

- Given and infinitely long and continuous signal $g(t)$ observed on a regular grid with spacing $Δt$ we have:

```math
g_{obs} = g(t) III_{\Delta t}(t) \qquad \hat{g}_{obs}(f) = \sum_{n=-\infty}^\infty g(n\Delta t) e^{-2\pi i f n \Delta t}
```

- In the real world we never have an infinite number of observations, defining $g_n \equiv g(n\Delta t)$ we may write:

```math
\hat{g}_{obs}(f) = \sum_{n=0}^N g_n e^{-2\pi i f n \Delta t}
```

- or, with $\Delta f = 1 / (N\Delta t)$ and $\hat{g}_k \equiv \hat{g}_{obs}(k\Delta f)$:

```math
\hat{g}_k = \sum_{n=0}^N g_n e^{-2\pi i k n / N}
```

- which is the standard form of the DFT.
"""

# ╔═╡ 533e41e2-08c3-450b-83eb-ae19664a6cfe
md"""
### The Classical Periodogram
***

- Applying the definition of Fourier spectrum $P_g \equiv |F{g}|^2$ we then have:

```math
P_S(f) = \frac{1}{N}\left|\sum_{n=1}^N g_n e^{-2\pi i f t_n}\right|^2
```

- This is the periodogram, which is an “estimator” of the true power spectrum of the underlying continuum function $g(t)$. It is sometimes called the *Schuster periodogram*.

- This point should be well emphasized: the *periodogram* and the *power spectrum* are conceptually different things, i.e. the periodogram is the statistic we compute from our data, and it is an *estimator* of the power spectrum, the underlying continuous function of interest.

- As we know, it is not a consistent estimator since it suffers from intrinsic variance even in the limit of an infinite number of observations.

"""

# ╔═╡ ae7829b7-b652-4e54-9231-f7bc3c3b6093
# ╠═╡ show_logs = false
cm"""
### Nonuniform sampling
***

- In the real world, particularly in fields like Astronomy where observations are subject to influences of weather and diurnal, lunar, or seasonal cycles, the sampling rate is generally far from uniform.

- In the general non-uniform case, we measure some signal at a set of ``N`` times, denoted ``\{t_n\}``, and defining the following observing window:

```math
W_{\{t_n\}}(t) = \sum_{n=1}^{N} \delta(t - t_n)
```

- Applying this window to the true underlying signal ``g(t)``, we get:

```math
g_{obs}(t) = g(t) W_{\{t_n\}}(t) = \sum_{n=1}^{N} g(t_n)\delta(t - t_n)
```

- And, just as in the evenly-sampled case, the FT of the observed signal is a convolution of the transforms of the true signal and the window:

```math
\mathcal{F}\{g_{obs}\} = \mathcal{F}\{g\} \ast \mathcal{F}\{W_{\{t_n\}}\}
```

- Unlike in the uniform case, the window transform ``\mathcal{F}\{W_{\{t_n\}}\}`` will generally not be a straightforward sequence of delta functions, and the symmetry present in the Dirac comb is broken by the uneven sampling, leading the transform to be much more "noisy".

| $(LocalResource("Pics/sparseirregular.jpg", :height=>300)) | $(LocalResource("Pics/denseirregular.png", :height=>300)) |
| ------------------------------------------------------ | ---------------------------------------------------- |

- The irregular spacing within the observing window translates to irregular frequency peaks in its transform, causing the observed transform to be noisier. Even with very dense sampling of the function, the FT cannot be exactly recovered due to the imperfect aliasing present in the window transform.

- Eventualy, the FT of a non-uniformly spaced delta functions looks like random noise, and partly it is.
    - It reflects the “random” distribution of the sampling time.

- A denser sampling helps, but to some extent “noise” is unavoidable.

"""

# ╔═╡ 9a422815-783a-473c-8538-6ecf5c208bb0
# ╠═╡ show_logs = false
md"""
### However, (almost) no Nyquist limit!
***

- The Nyquist limit is a direct consequence of the symmetry in the Dirac comb window function that describes evenly-sampled data, and uneven sampling destroys the symmetry that underlies its definition.

- For unevenly-sampled data, a "Nyquist limit" might or might not exist, and even in cases where it does exist it tends to be far larger (and thus far less relevant) than in the evenly-sampled case.

$(LocalResource("Pics/nonyquist.jpg"))

- An example of data for which the various poorly-motivated "pseudo-Nyquist" approaches fail.
    - The upper panels show the data, a noisy sinusoid with a frequency of 100 (i.e. a period of 0.01).
    - The lower left panel shows a histogram of spacings between observations: the minimum spacing is 2.55, meaning that the signal has over 250 full cycles between the closest pair of observations.
- Nevertheless, the periodogram (lower right) clearly identifies the correct period, though it is orders of magnitude larger than pseudo-Nyquist estimates based on average or minimum sampling rate.
"""

# ╔═╡ 6fa31a81-cbbe-401b-acb0-0d44b0980363
cm"""
## The Lomb-Scargle periodogram
***

- The classical periodogram can be rewritten as:

```math
P(f) = \frac{1}{N}\left|\sum_{n=1}^N g_n e^{-2\pi i f t_n} \right|^2 = \frac{1}{N}\left[
    \left(\sum_n g_n \cos(2\pi f t_n)\right)^2
    + \left(\sum_n g_n \sin(2\pi f t_n)\right)^2
    \right]
```
    
- In principle this formula could be used for non-uniform sampling too, yet the obtained periodogram does not offer, in general, some of the useful statistical properties of the even sampling case.

- The problem was addressed by proposing a generalized expression:

```math
P(f) = \frac{A^2}{2}\left(\sum_n g_n \cos(2\pi f [t_n-\tau])\right)^2
       + \frac{B^2}{2} \left(\sum_n g_n \sin(2\pi f [t_n-\tau])\right)^2
``

- where ``A``, ``B``, and ``\tau`` are arbitrary functions of the frequency ``f`` and observing times ``\{t_i\}`` (but not the values ``\{g_n\}``).

- It is possible to define ``A``, ``B``, and ``\tau`` such as:
    - the periodogram reduces to the classical form in the case of equally spaced observations;
    - the periodogram distribution is analytically computable;
    - the periodogram is insensitive to global time shifts in the data.

- The values of $A$ and $B$ leading to these properties result in the following
form of the generalized periodogram:

```math
P_{LS}(f) =
  \frac{1}{2} \Bigg\{
  \bigg(\sum_n g_n \cos(2\pi f [t_n-\tau])\bigg)^2 \bigg/
  \sum_n \cos^2(2\pi f [t_n-\tau])  + ~ \bigg(\sum_n g_n \sin(2\pi f [t_n-\tau])\bigg)^2 \bigg/
  \sum_n \sin^2(2\pi f [t_n-\tau])  \Bigg\}
```

- where ``\tau`` is specified for each ``f`` to ensure time-shift invariance:

```math
\tau = \frac{1}{4\pi f}\tan^{-1}\Bigg(
  \frac{\sum_n \sin(4\pi f t_n)}{\sum_n \cos(4\pi f t_n)}\Bigg)
```

- This modified periodogram differs from the classical periodogram only to the extent that the denominators ``\sum_n \sin^2(2\pi f t_n)`` and ``\sum_n \cos^2(2\pi f t_n)`` differ from ``N/2``, which is the expected value of each of these quantities in the limit of complete phase sampling at each frequency.

- Thus, in many cases of interest the Lomb-Scargle periodogram only differs slightly (but not negligibly) from the classical/Schuster periodogram.

- An important result is that the LS periodogram could be obtained fitting a simple sinusoidal model to the data at each given frequency, and deriving the periodogram from the goodness of fit (``\chi^2``).

    - We will see later that this possible different derivation has important consequences.

"""

# ╔═╡ 367c08ff-cedd-4e63-b53a-3178b0ec45cb
md"""
#### Exercise about classical vs Lomb-Scargle periodograms
***

- Let's study a curve observed at 30 epochs with some noise. Although with a sparse sampling, it is a sinusoidal curve with period.
"""

# ╔═╡ ef919d43-ffb2-404f-94e4-a10b9a543e09
begin
	
	t = [48.68293613, 43.68819185, 95.71437369, 71.00603639, 48.54022079,
	       17.820341  , 89.19924435, 61.75940713, 95.14511764, 80.92812141,
	       68.32096035, 11.61911367, 76.40351048, 24.05532227, 91.73672731,
	       36.27323062, 39.10355853, 49.25845392, 93.18918459, 30.02253181,
	       29.31289813, 88.33416194, 94.6922681 , 59.80890946, 52.82566667,
	       90.27610636, 12.10326369,  4.5971768 , 68.08097096, 38.76998251]
	
	y = [-1.06749625, -0.25882508,  0.54888049,  0.50273813, -1.27428505,
	       -1.17429603,  0.79366548, -0.64418343,  0.59271595, -0.65160916,
	       -0.46201257, -0.88418838,  0.74365945, -1.13026806, -0.78967783,
	       -0.75395011,  0.72986805, -0.58586657, -0.60369567, -0.9487803 ,
	       -0.49721427,  1.21028735,  0.32603806,  0.0210877 ,  0.09798819,
	        0.78991859, -0.3939756 , -1.13497054, -1.0160708 ,  0.84690192]
	
	
	fg1 = Figure()
	
	ax1fg1 = Axis(fg1[1, 1])
	
	scatter!(t,y,)
	
	fg1
end

# ╔═╡ ce5b3710-ac14-44f5-85d8-abc9756c1f11
md"""
- Now, let's write one of the possible versions of the Fourier periodogram. This is not optimized as a FFT, yet it works for any input frequency set.
"""

# ╔═╡ b985cda2-8242-468b-b5de-2299bfa50b33
md"""
- And then compute the Lomb-Scargle periodogram and Schuster periodograms, together with the functions characterized the Lomb-Scargle formula.
"""

# ╔═╡ 2e3e0747-48c2-43d1-abd4-b15f61a402bd
begin
	lsplan = LombScargle.plan(t,y,normalization=:psd,minimum_frequency=0,maximum_frequency=0.5,samples_per_peak=20)
	
	lsprg = lombscargle(lsplan)
end;

# ╔═╡ 18958f04-3eb0-4e4e-9957-ca0c5540e5a5
function SchusterPeriodogram(t, mag, freq)
    pwr = abs.(exp.(-2im .* π .* lsprg.freq * t') * y).^2
    return pwr/length(t)
end;

# ╔═╡ ccb2f0e6-517b-4860-bac2-dac964e8d388
begin
	p_schuster = SchusterPeriodogram(t, y, lsprg.freq)
	
	tau = 1 ./ (4 .* π .* lsprg.freq[2:end]) .* atan.(sum(sin.(4 .* π .* lsprg.freq[2:end] * t'),dims=2), sum(cos.(4 .* π .* lsprg.freq[2:end] * t'),dims=2))
end;

# ╔═╡ 0b098634-7890-4551-8343-34dbba3159ba
begin
	sw = 0.
	cw = 0.
	for ep in t
	    global sw = sw .+ sin.(2 .* π .* lsprg.freq[2:end] .* (ep .- tau)).^2
	    global cw = cw .+ cos.(2 .* π .* lsprg.freq[2:end] .* (ep .- tau)).^2
	end
	sin_window = sw
	cos_window = cw
end;

# ╔═╡ 2503892b-3bcf-41b8-a636-562d2613bd59
md"""
- And, finally, let's plot the results.
"""

# ╔═╡ f652aeb7-bcf0-481c-924c-412f406b3542
begin
	fg2 = Figure(size=(1000,1000))
	
	ax1fg2 = Axis(fg2[1, 1],
	    title="Spectral Power"
	    )
	
	lines!(lsprg.freq[2:end],lsprg.power[2:end],label="Lomb-Scargle Periodogram",color=:blue)
	lines!(lsprg.freq[2:end],p_schuster[2:end],label="Classical Periodogram",color=:red)
	
	axislegend()
	
	xlims!(0,0.5)
	
	ax2fg2 = Axis(fg2[2, 1],
	    xlabel="Frequency"
	    )
	
	
	lines!(lsprg.freq[2:end],vec(cos_window),label=L"$\sum_n\ \cos^2 [2\pi f (t_n-\tau)]$",color=:red)
	lines!(lsprg.freq[2:end],vec(sin_window),label=L"$\sum_n\ \sin^2 [2\pi f (t_n-\tau)]$",color=:blue)
	hlines!(length(t)/2,linestyle=:dash,color=:black)
	
	axislegend()
	
	xlims!(0,0.5)
	ylims!(0,30)
	
	fg2
	
	
end

# ╔═╡ 637108d6-c977-4a2d-b312-e40bf375d5f1
md"""
- As anticipated, LS and Fourier periodogram are usually rather similar, and their differences rely in how much the quantities plotted in the bottom plot are different from N/2.
"""

# ╔═╡ 325a6b84-af05-4d78-94a2-3674a501c1c2
# ╠═╡ show_logs = false
md"""
## LS extensions
***

- Considering the LS formula as the result of a regular fitting procedure allows one to discussn interesting generalizations.

- In the least squares interpretation of the periodogram, a sinusoidal model is proposed at each candidate frequency $f$: $y(t;f) = A_f \sin(2 \pi f (t - \phi_f))$, where the amplitude $A_f$ and phase $\phi_f$ can vary as a function of frequency.

- These model parameters are fit to the data by constructing the $\chi^2$ statistic at each frequency: $\chi^2(f) \equiv \sum_n \big(y_n - y(t_n;f)\big)^2$ and minimizing $\chi^2(f)$ at with respect to $A_f$ and $\phi_f$.

- Denoting the minimum value as $\hat{\chi}^2(f)$ the LS periodogram can be equivalently written:

```math
P(f) = \frac{1}{2}\big[\hat{\chi}^2_0 - \hat{\chi}^2(f)\big]
```

- where $\hat{\chi}^2_0$ is the non-varying reference model.

- Basing on this view a trivial yet fundamental extension is to include in the fitting procedure the (Gaussian) errors on the data through the standard change to the $\chi^2$ expression:

```math
\chi^2(f) \equiv \sum_n \left(\frac{y_n - y_{model}(t_n;f)}{\sigma_n}\right)^2
```

- This generalization of $\chi^2(f)$ also suggests a convenient way to construct a periodogram in the presence of correlated observational noise. If we let $\Sigma$ denote the  $N\times N$ noise covariance matrix for $N$ observations:

```math
\vec{y} = [y_1, y_2,\cdots y_n]^T \qquad \vec{y}_{model} = [y_{model}(t_1),y_{model}(t_1),\cdots y_{model}(t_n)]^T
```

- and therefore:

```math
\chi^2(f) = (\vec{y}-\vec{y}_{model})^T\Sigma^{-1}(\vec{y}-\vec{y}_{model})
```

- which reduces to the previous form if noise is uncorrelated (i.e., if the off-diagonal terms of $\Sigma$ are zero).

- Another import extension involves adding an offset term to the sinusoidal model at each frequency (the "floating mean" model):

```math
y_{model}(t;f) = y_0(f) + A_f \sin(2 \pi f (t - \phi_f))
```

- This avoids the need to “zero-center” the light-curve, particularly important for incomplete (e.g. flux-limited) list curve.

$(LocalResource("Pics/floatingmean.png"))

-  The previous plot shows a comparison of the standard and floating mean periodograms for data with a frequency of 0.3 and a selection effect which removes faint observations. In this case the mean estimated from the observed data is not close to the true mean, which leads to the failure of the standard periodogram to recover the correct frequency. A floating mean model correctly recovers the true frequency of 0.3.

"""

# ╔═╡ 50a1fb43-f72d-43bd-bea7-99bac76db19e
md"""
- The interpretation of a periodogram as a regular fit procedure opens the way to even more extensions. It is possible, for instance, to model the curve with multiple Fourier components, allowing a greater flexibility:

```math
y_{model}(t;f) = A_f^{0} + \sum_{k=1}^K A_f^{(k)} \sin(2\pi k f (t - \phi_f^{(k)}))
```
"""

# ╔═╡ 7cd277fd-6b32-4ca8-9c8e-1eab55659c74
md"""
#### Exercise about the effect of uncertainties in periodogram computation (homoscedastic vs heteroskedastic case)
***

- We now want to evaluate the effect of heteroskedasticity on the computaton of a periodogram.
"""

# ╔═╡ 46aea63c-b5a8-4977-8b87-f1cc479498e0
md"""
Noise (base) level: $( @bind σ_base PlutoUI.Slider(0.:0.1:2., default=0.3) ) 
"""

# ╔═╡ 1c023f88-4e2c-494b-9c03-002591955b60
begin
	Random.seed!(42)
	
	# ─────────────────────────────────────────────────────────────────────────────
	# 1. GENERA IL SEGNALE
	# ─────────────────────────────────────────────────────────────────────────────
	
	N_total   = 300          # campioni su griglia regolare (per il periodogramma classico)
	N_sparse  = 90           # campioni irregolari (per LS)
	T_obs     = 50.0         # durata totale [unità di tempo]
	f_true    = 0.25         # frequenza del segnale [Hz]  → periodo = 4
	A_signal  = 2.0          # ampiezza del segnale
	
	# Griglia regolare (usata dal periodogramma FFT)
	t_reg   = range(0.0, T_obs; length=N_total)
	noise_reg = randn(N_total) .* 0.8
	y_reg   = A_signal .* sin.(2π .* f_true .* t_reg) .+ noise_reg
	
	# Campionamento irregolare: rimuovi cluster di punti (simula osservazione astronomica)
	idx_keep = sort(randperm(N_total)[1:N_sparse])
	t_irr    = collect(t_reg)[idx_keep]
	
	# Incertezze eterogenee: alcune misure molto più rumorose
	#σ_base  = 0.3
	σ_vals  = σ_base .* (1.0 .+ 4.0 .* rand(N_sparse).^2)   # distribuzione asimmetrica
	
	# Rumore proporzionale all'incertezza locale
	noise_irr = σ_vals .* randn(N_sparse)
	y_irr     = A_signal .* sin.(2π .* f_true .* t_irr) .+ noise_irr
end;

# ╔═╡ d5cd7643-6c82-4be1-a995-6e1f79547288
begin
	# ─────────────────────────────────────────────────────────────────────────────
	# 2. PERIODOGRAMMA CLASSICO (FFT su griglia regolare)
	# ─────────────────────────────────────────────────────────────────────────────
	
	Y     = fft(y_reg .- mean(y_reg))
	freqs_fft = fftfreq(N_total, N_total / T_obs)
	
	# Solo frequenze positive
	pos   = freqs_fft .> 0
	freqs_fft_pos = freqs_fft[pos]
	power_fft     = (2 / N_total) .* abs2.(Y[pos])
end;

# ╔═╡ ae2ac755-825b-4687-aeec-1de9a4ff24fb
begin
	# ─────────────────────────────────────────────────────────────────────────────
	# 3. PERIODOGRAMMA DI LOMB-SCARGLE (dati irregolari + incertezze)
	# ─────────────────────────────────────────────────────────────────────────────
	
	# Versione SENZA pesi (ignora le incertezze)
	plan_nw = LombScargle.plan(t_irr, y_irr; normalization=:standard,
	                            minimum_frequency=0.01, maximum_frequency=2.0,
	                            samples_per_peak=10)
	pgram_nw = lombscargle(plan_nw)
	
	# Versione CON pesi (usa le incertezze)
	plan_w = LombScargle.plan(t_irr, y_irr, σ_vals; normalization=:standard,
	                           minimum_frequency=0.01, maximum_frequency=2.0,
	                           samples_per_peak=10)
	pgram_w = lombscargle(plan_w)
	
	freqs_ls  = freqpower(pgram_nw)[1]
	power_nw  = freqpower(pgram_nw)[2]
	power_w   = freqpower(pgram_w)[2]
	
	# Livello di falso allarme al 1%
	fap_nw = LombScargle.fap(pgram_nw, 0.01)
	fap_w  = LombScargle.fap(pgram_w,  0.01)
end;

# ╔═╡ 981eb94f-2bb9-4330-892b-9f0a5f929a00
begin
	# ─────────────────────────────────────────────────────────────────────────────
	# 4. FIGURA
	# ─────────────────────────────────────────────────────────────────────────────
	
	fig = Figure(size=(1100, 860), backgroundcolor=:white)
	
	# Palette colori
	col_signal  = RGBf(0.15, 0.45, 0.75)
	col_noisy   = RGBf(0.85, 0.35, 0.15)
	col_fft     = RGBf(0.30, 0.30, 0.30)
	col_ls_nw   = RGBf(0.80, 0.50, 0.10)
	col_ls_w    = RGBf(0.10, 0.65, 0.35)
	col_fap     = RGBf(0.70, 0.10, 0.10)
	
	# ── Titolo globale ──────────────────────────────────────────────────────────
	Label(fig[0, 1:2],
	    "Classic Periodogram vs Lomb-Scargle Periodogram with Heteroskedastic Uncertainties";
	    fontsize=17, font=:bold, tellwidth=false)
	
	# ── Pannello A: serie temporale irregolare ──────────────────────────────────
	ax_ts = Axis(fig[1, 1:2];
	    title  = "A  –  Time series with irreguylar sampling",
	    xlabel = "Time",
	    ylabel = "Amplitude",
	    titlesize=13, xlabelsize=11, ylabelsize=11)
	
	# Segnale "vero" sottostante
	t_fine  = range(0.0, T_obs; length=1000)
	y_true  = A_signal .* sin.(2π .* f_true .* t_fine)
	lines!(ax_ts, t_fine, y_true; color=col_signal, linewidth=1.5,
	       linestyle=:dash, label="True signal  (f = $(f_true) Hz)")
	
	# Dati con incertezze
	errorbars!(ax_ts, t_irr, y_irr, σ_vals;
	           color=(col_noisy, 0.40), linewidth=0.8, whiskerwidth=4)
	scatter!(ax_ts, t_irr, y_irr;
	         color=col_noisy, markersize=5, label="Measurements with variable σ")
	
	axislegend(ax_ts; position=:rt, framevisible=false, labelsize=10)
	
	# ── Pannello B: distribuzione delle incertezze ──────────────────────────────
	ax_sig = Axis(fig[2, 1];
	    title  = "B  –  Distribution of σ",
	    xlabel = "σ",
	    ylabel = "Counts",
	    titlesize=13, xlabelsize=11, ylabelsize=11)
	
	hist!(ax_sig, σ_vals; bins=20, color=(col_noisy, 0.70),
	      strokecolor=:white, strokewidth=0.5)
	vlines!(ax_sig, [σ_base]; color=col_signal, linewidth=2,
	        linestyle=:dash, label="Base σ")
	axislegend(ax_sig; position=:rt, framevisible=false, labelsize=10)
	
	# ── Pannello C: periodogramma FFT ───────────────────────────────────────────
	ax_fft = Axis(fig[2, 2];
	    title  = "C  –  FFT Periodogram (regular sampling, σ = cost.)",
	    xlabel = "Frequency [Hz]",
	    ylabel = "Power",
	    titlesize=13, xlabelsize=11, ylabelsize=11)
	
	lines!(ax_fft, freqs_fft_pos, power_fft; color=col_fft, linewidth=1.2)
	vlines!(ax_fft, [f_true]; color=col_signal, linewidth=2,
	        linestyle=:dash, label="True f = $(f_true)")
	xlims!(ax_fft, 0.0, 1.5)
	axislegend(ax_fft; position=:rt, framevisible=false, labelsize=10)
	
	# ── Pannello D: Lomb-Scargle senza pesi ─────────────────────────────────────
	ax_nw = Axis(fig[3, 1];
	    title  = "D  –  Lomb-Scargle Periodogram (Uncertainties not used)",
	    xlabel = "Frequency [Hz]",
	    ylabel = "LS Power",
	    titlesize=13, xlabelsize=11, ylabelsize=11)
	
	lines!(ax_nw, freqs_ls, power_nw; color=col_ls_nw, linewidth=1.2)
	hlines!(ax_nw, [fap_nw]; color=col_fap, linewidth=1.5,
	        linestyle=:dot, label="FAP 1%")
	vlines!(ax_nw, [f_true]; color=col_signal, linewidth=2,
	        linestyle=:dash, label="True f")
	xlims!(ax_nw, 0.0, 1.5)
	axislegend(ax_nw; position=:rt, framevisible=false, labelsize=10)
	
	# ── Pannello E: Lomb-Scargle con pesi ───────────────────────────────────────
	ax_w = Axis(fig[3, 2];
	    title  = "E  –  Lomb-Scargle Periodogram (Uncertainties included)",
	    xlabel = "Frequency [Hz]",
	    ylabel = "LS Power",
	    titlesize=13, xlabelsize=11, ylabelsize=11)
	
	lines!(ax_w, freqs_ls, power_w; color=col_ls_w, linewidth=1.2)
	hlines!(ax_w, [fap_w]; color=col_fap, linewidth=1.5,
	        linestyle=:dot, label="FAP 1%")
	vlines!(ax_w, [f_true]; color=col_signal, linewidth=2,
	        linestyle=:dash, label="True f")
	xlims!(ax_w, 0.0, 1.5)
	axislegend(ax_w; position=:rt, framevisible=false, labelsize=10)
	
	# ── Pannello F: confronto diretto LS con vs senza pesi ──────────────────────
	ax_cmp = Axis(fig[4, 1:2];
	    title  = "F  – Direct comparison: LS with and without uncertainties",
	    xlabel = "Frequency [Hz]",
	    ylabel = "Normalized LS Power",
	    titlesize=13, xlabelsize=11, ylabelsize=11)
	
	# Normalizza al picco per confronto visivo
	lines!(ax_cmp, freqs_ls, power_nw ./ maximum(power_nw);
	       color=col_ls_nw, linewidth=1.5, label="LS without uncertainties")
	lines!(ax_cmp, freqs_ls, power_w  ./ maximum(power_w);
	       color=col_ls_w,  linewidth=1.5, label="LS with uncertainties")
	vlines!(ax_cmp, [f_true]; color=col_signal, linewidth=2,
	        linestyle=:dash, label="True f = $(f_true) Hz")
	xlims!(ax_cmp, 0.0, 1.5)
	axislegend(ax_cmp; position=:rt, framevisible=false, labelsize=10)
	
	# ── Note a piè di figura ────────────────────────────────────────────────────
	note = """
	Methodological notes:
	• FFT (panel C): requires a uniform time grid; ignores constructional uncertainties.
	• Unweighted LS (D): handles irregular sampling but treats all measurements as equivalent.
	• Weighted LS (E, F): assigns weight wᵢ = 1/σᵢ² → more precise measurements matter more.
	With heterogeneous uncertainties, the peak at the true frequency is sharper and the noise floor is lower.
	"""
	Label(fig[5, 1:2], note; fontsize=9, tellwidth=false,
	      halign=:left, justification=:left,
	      color=RGBf(0.35, 0.35, 0.35))
	
	rowsize!(fig.layout, 0, Auto(0.3))
	rowsize!(fig.layout, 5, Auto(0.5))
	rowgap!(fig.layout, 8)
	colgap!(fig.layout, 12)
	
	fig
end

# ╔═╡ 910a90a8-7b85-4d69-ba8d-3cabbfe89bc4
md"""
### The Window Function
***

- It is clear that a careful evaluation of the window function is crucial in interpreting the results of a LS analysis.

- The LS periodogram can offer a sufficiently reliable way to compute it:

```math
\mathcal{P}_W(f;\{t_n\}) = \left|\sum_{n=1}^{N} e^{-2\pi i f t_n}\right|^2
```

- It is, in its essence, the periodogram for data $g_n=1$ at all times $t_n$. No floating-mean model should be used in this case.

"""

# ╔═╡ 8ca6fd17-566d-4aae-961c-cc4a4e6ab1d0
md"""
#### Exercise about multiple component periodgram
***

- We analyse the light-curve of th LINEAR object [14752041](http://simbad.u-strasbg.fr/simbad/sim-id?Ident=linear+14752041&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id), an eclipsing binary star.
"""

# ╔═╡ 2b0e312c-4285-4e63-9c68-b737646869d9
begin
	lnr = DataFrame(CSV.File("LINEAR_14752041.csv"))
	first(lnr,5)
end

# ╔═╡ 171b2e3c-058e-4d82-9d0e-b5dd5f521404
begin
	fg3 = Figure()
	
	ax1fg3 = Axis(fg3[1, 1],
	    xlabel="time (MJD)",
	    ylabel="magnitude",
	    title="LINEAR object 14752041",
	    )
	
	scatter!(lnr[!,:t],lnr[!,:mag])
	errorbars!(lnr[!,:t],lnr[!,:mag],lnr[!,:magerr])
	
	ax1fg3.yreversed=true
	
	fg3
end

# ╔═╡ ec859add-e8e7-413b-b9d4-7d2b8694a46d
md"""
- This looks like a rich light-curve, althgough with a rather irregular sampling and large gaps.

- Before analysing the light-curve let's study the effect of the observing window on the periodogram.
"""

# ╔═╡ eae2ed93-31fb-417a-9621-5586cf4f9732
begin
	ls_window = lombscargle(lnr[!,:t], ones(size(lnr[!,:t])), fit_mean=false, center_data=false, minimum_frequency=0.001, maximum_frequency=10)
	
	
	fg4 = Figure()
	
	ax1fg4 = Axis(fg4[1, 1],
	    xlabel="period (day)",
	    ylabel="power",
	    title="Observing window periodogram",
	    xscale=log10
	    )
	
	lines!(1 ./ ls_window.freq,ls_window.power)
	
	fg4
end

# ╔═╡ f38c7158-5206-467b-987d-38b4aa80e935
md"""
- The window periodogram is indeed quite intersting: togetehr with a rather large peak corresponding to the $\sim 1$ year periodicity of the observations, there is power at several hours and days, reflecting the observing strategy.

- Let's move now to the data periodogram:
"""

# ╔═╡ ce90af7b-ef2c-4cdf-84cb-de76f5dd5c66
begin
	lsper = lombscargle(lnr[!,:t],lnr[!,:mag],lnr[!,:magerr],minimum_frequency=0.001,maximum_frequency=10)
	
	fg5 = Figure(size=(1000,1000))
	
	ax1fg5 = Axis(fg5[1, 1],
	    xlabel="period (day)",
	    ylabel="power",
	    title="LINEAR object 14752041 periodogram",
	    xscale=log10
	    )
	
	lines!(1 ./ lsper.freq,lsper.power)
	
	
	fg5
end

# ╔═╡ 0a342241-26b4-425e-b90f-fc555efb0b67
# ╠═╡ show_logs = false
md"""
- The periodogram shows a forest of interesting features, mainly for periods of a few hours.
    - This often happensa when the shape of variability is not close to a sinusoid.
    
- This case can be deal with a multi-term analysis.

$(LocalResource("Pics/higherorder.png"))

- 1-term and 6-term Lomb-Scargle models fit to an eclipsing binary. The standard periodogram finds an alias of the true 17.5-hour frequency, because a simple sinusoidal model cannot closely fit both the primary and secondary eclipse. A six-term Fourier model, on the other hand, does find the true period.

- In principle it is not even necessary to think to sinusoidal terms only. If a more complex signal shape is physically motivated, there are no intrinsic limitations.

- Model comparison is better carried out in a full Bayesian framework, with properly formalized prior knowledge, etc.

- Quite interestingly, it can be shown that the LS periodogram is in fact the optimal statistics for detecting a stationary sinusoidal signal in the presence of Gaussian noise.

"""

# ╔═╡ 17444051-6ee8-44af-8cfc-275bcb495d86
# ╠═╡ show_logs = false
md"""
### Bayesian LS periodogram
***

- The least squares view of the Lomb-Scargle periodogram creates a natural bridge, via maximum likelihood, to Bayesian periodic analysis. In fact, in the Bayesian view, the Lomb-Scargle periodogram is the optimal statistic for detecting a stationary sinusoidal signal in the presence of Gaussian noise.

- For the standard, simple-sinusoid model, the Bayesian periodogram is given by the posterior probability of frequency ``f`` given the data ``D`` and sinusoidal model ``M``:

```math
p(f\mid D, M) \propto e^{P_{LS}(f)}
```

- where ``P_{LS}(f)`` is the LS power. The effect of this exponentiation is to suppress side-lobes and alias peaks in relation to the largest one of the spectrum.

$(LocalResource("Pics/bayesianls.png"))

- Be aware that ``p(f\mid D,M)`` is the probability data are drawn from a sinusoidal model. Not the probability that data are periodic in general.

    - More technically, this is indeed the probability conditioned on the assumption that the data are drawn from a sinusoidal model.

- Some more detail about the derivation of the LS periodogram in a Bayesian framework can be found here ([notebook](./open?path=Lectures/Lecture-Lomb-Scargle/Lecture-BayesLS.jl), [html](../../Lectures/Lecture-Lomb-Scargle/Lecture-BayesLS.html)). 
"""

# ╔═╡ 467e6316-fa6d-4487-86cc-8ba408577fe7
# ╠═╡ show_logs = false
cm"""
### LS frequency grid
***

- At variance with the DFT, the frequency grid for LS analysis is not determined by the sampling.

- The most important warning is to pay attention not to choose a too sparse grid, with the risk to miss important features.

$(LocalResource("Pics/lsgrid.png"))

- In general, for an observation length ``T``, we have sinc-shaped peaks of width ``\sim 1/T``. A rule of the thumb advise is typically to oversample the grid by a factor ``\sim 5``.

"""

# ╔═╡ bd68480e-e2f6-4759-8333-066001c2ae19
md"""
### LS (or DFT) peak uncertainty
***

- It happens frequenyly thay period uncertainties are quoted basing on the peak width. This is somehow ill-defined, since the width of a periodogram peak depends, essentialy, on the length of the observation.
"""

# ╔═╡ ec7fa3e8-2635-4af3-b2e7-6435e96ebcb2
md"""
#### Exercize about LS (and DFT) periodogram peak uncertainty
***

- First, let's define a function to generate light-curves of different lengths and S/N.
"""

# ╔═╡ 0dafbc5d-d383-44a8-8dc8-704f8a36f7fe
function create_data(N; T=4, signal_to_noise=5, period=1.0, random_state=None)
    rng = Random.seed!(random_state)
    t = T .* rand(N)
    dy = 0.5 ./ signal_to_noise .* ones(size(t))
    y = sin.(2 * π * t ./ period) .+ dy .* randn(N)
    return t, y, dy
end;

# ╔═╡ fe1aa970-abc8-4720-9f09-5fe931100263
md"""
- And plot several light-curves with different parameters.
"""

# ╔═╡ 9d413d79-c9fb-4c30-ac02-92380510b607
# ╠═╡ show_logs = false
begin
	fg6 = Figure(size=(700,350))
	
	ax1fg6 = Axis(fg6[1, 1],
	    title="Peak scaling with number of data points (fixed S/N=10)",
	    ylabel=L"$P_{LS}$ (normalized)"
	    )
	
	ax2fg6 = Axis(fg6[2, 1],
	    ylabel=L"$P_{LS} / (S/N)^2$ (PSD-normalized)",
	    xlabel="frequency"
	    )
	
	ax3fg6 = Axis(fg6[1, 2],
	    title="Peak scaling with signal-to-noise ratio (fixed N=1000)",
	    ylabel=L"$P_{LS}$ (normalized)"
	    )
	
	ax4fg6 = Axis(fg6[2, 2],
	    ylabel=L"$P_{LS} / (S/N)^2$ (PSD-normalized)",
	    xlabel="frequency"
	    )
	
	
	SN = 10
	for N in [1000, 100, 10]
	    t, y, dy = create_data(N, signal_to_noise=SN, random_state=68345)
	    freq = LinRange(0.01, 4, 2000)
	
	    ls1 = lombscargle(t, y, dy, normalization=:standard, frequencies=freq)
	    ls2 = lombscargle(t, y, dy, normalization=:psd, frequencies=freq)
	
	    lines!(ax1fg6,freq,ls1.power,label="N="*string(N))
	    lines!(ax2fg6,freq,ls2.power ./ SN, label="N="*string(N))
	end
	
	
	N = 1000
	for SN in [10, 1, 0.1]
	    t, y, dy = create_data(N, signal_to_noise=SN, random_state=68345)
	    freq = LinRange(0.01, 4, 2000)
	
	    ls3 = lombscargle(t, y, dy, normalization=:standard, frequencies=freq)
	    ls4 = lombscargle(t, y, dy, normalization=:psd, frequencies=freq)
	
	    lines!(ax3fg6,freq,ls3.power,label="S/N="*string(SN))
	    lines!(ax4fg6,freq,ls4.power ./ SN^2, label="S/N="*string(SN))
	end
	
	axislegend(ax1fg6)
	axislegend(ax2fg6)
	axislegend(ax3fg6)
	axislegend(ax4fg6)
	
	
	
	fg6
end

# ╔═╡ 26d6e75f-dd44-4876-ad72-e5606d27798f
md"""
- The effect of the number of points $N$ and the signal-to-noise ratio $S/N$ on the expected width and height of the periodogram peak.
    - Top panels show the normalized periodogram, while bottom panels show the PSD-normalized periodogram scaled by noise variance. Perhaps surprisingly, neither the number of points nor the signal-to-noise ratio affects the peak width.

- Within the Bayesian interpretation, one can derive an approximate expression relating the uncertainty to the number of samples, $N$, and the average $S/N$ ratio, $Σ$:

```math
\sigma_f \approx f_{1/2} \sqrt{\frac{2}{N\Sigma^2}}
```

"""

# ╔═╡ 66da0ee4-4d50-47e6-ab43-3385b2a9f9fa
# ╠═╡ show_logs = false
md"""
### LS Peak Significance
***

- The typical approach to quantifying the significance of a peak is the *False Alarm Probability* (FAP), which measures the probability that a dataset with no signal would, due to coincidental alignment among the random


- Even for the LS periodogram, and pure Gaussian noise, it can be proved that the values of the unnormalized periodogram follow a $\chi^2$ distribution with two degrees of freedom.

- If $Z=P(f_0)$ is the periodogram value at a peak with frequency $f_0$, then the cumulative probability to observe a value less then $Z$ is:

```math
P_{single}(Z) = 1 - \exp(-Z)
```

- Analogously to the DFT, we are generally not interested in the distribution of one particular randomly chosen frequency, but rather the distribution of the highest peak of the periodogram.

- The problem here, at variance with the DFT case, is not simple to solve since the value at one frequency is correlated with the value at other frequencies in a way that is quite difficult to analytically express. These correlations come from the convolution with the survey window.

- One common approach is to assume that it can be modeled on some *effective number* of independent frequencies $_{Neff}$, so that the FAP can be estimated as:

```math
FAP(z) \approx 1 - \big[P_{single}(z)\big]^{N_{eff}}
```

- A very simple estimate for $_{Neff}$ can based on the arguments about the expected peak width, $δf = 1/T$. In this approximation, the number of independent peaks in a range $0 ≤ f ≤ f_{max}$ is assumed to be $N_{eff} = f_{max} T$.

- Finding a reliable way to accurately compute $_{Neff}$ analytically is still and open problem. Often the problem is solved computationally by, e.g., a bootstrap procedure.

$(LocalResource("Pics/independentfreqs.png"))

- In the previous plot a comparison of various approaches to computing the FAP for simulated observations with both structured and unstructured survey windows is shown.
"""

# ╔═╡ e30cc0a2-ce7f-4e9f-9f3f-ebc14a76e936
md"""
#### Exercize about peak periodogram significance evaluation
***

- Let's generate a short time-series with the function defined above and compute the sigificance level by means of a boostrap analysis.
"""

# ╔═╡ 2ef299bf-e198-4bde-83cd-ecd86e566683
begin
	t7, y7, dy7 = create_data(30, signal_to_noise=2, random_state=583)
	
	lsper7 = lombscargle(t7, y7, dy7, maximum_frequency=4, samples_per_peak=10)
	
	lsboot7 = LombScargle.bootstrap(1000, t7, y7, dy7, maximum_frequency=4, samples_per_peak=10)
	
	p85 = LombScargle.fapinv(lsboot7,0.85)
	p95 = LombScargle.fapinv(lsboot7,0.95)
	p99 = LombScargle.fapinv(lsboot7,0.99)
	
	
	fg7 = Figure(size=(700,350))
	
	ax1fg7 = Axis(fg7[1, 1],
	    )
	
	ax2fg7 = Axis(fg7[1, 2],
	    )
	
	scatter!(ax1fg7,t7,y7)
	errorbars!(ax1fg7,t7,y7,dy7)
	
	lines!(ax2fg7,lsper7.freq,lsper7.power)
	hlines!(ax2fg7,p85,label="85%")
	hlines!(ax2fg7,p95,label="95%")
	hlines!(ax2fg7,p99,label="99%")
	
	axislegend(ax2fg7)
	
	fg7
end

# ╔═╡ 6d2ccb03-4455-4249-9d3d-be4cb5d7a45f
md"""
## Various possible PSD normalizations
***

- When considering the periodogram from the Fourier perspective, it is useful to normalize the periodogram such that in the special case of equally spaced data it recovers the standard Fourier power spectrum.

- This is the so-called "psd" normalization, and the equivalent least-squares expression is:

```math
P(f) = \frac{1}{2}\big[\hat{\chi}^2_0 - \hat{\chi}^2(f)\big]
```

- For equally spaced data this becomes:

```math
P(f) = \frac{1}{N} \left| FFT(y_n) \right|^2
```

- Under the "psd" normalization periodogram units are “unit(y)$^2$”. And can be interpreted as squared amplitudes of the Fourier component at each frequency.

- However, if uncertainties are included in the analysis, the periodogram becomes unitless: it is essentially a measure of periodic content in signal-to-noise ratio rather than in signal itself.

- In the least-squares view of the periodogram, the periodogram is interpreted as an inverse measure of the goodness of  fit for a model.
    - If the sinusoidal model perfectly fits the data at some frequency $f_0$, then $\hat{\chi}^2(f_0) = 0$ and the periodogram is maximized at a value of $\hat{\chi}^2_0 / 2$.
    
- The minimum value of the periodogram can only be 0, and  therefore a possible normalization that keep the (unitless) values between 0 and 1 is:

```math
P_{norm}(f) = 1 - \frac{\hat{\chi}^2(f)}{\hat{\chi}^2_0}
```
"""

# ╔═╡ 78119e66-703b-4305-85cd-956f1053cc60
md"""
## Reference & Material

Material and papers related to the topics discussed in this lecture.

- [Vander Plas (2017) - "Understanding the Lomb-Scargle Periodogram”](https://ui.adsabs.harvard.edu/abs/2018ApJS..236...16V/abstract)
"""

# ╔═╡ 435fd658-c227-47ce-8dd6-cca20bff2a26
md"""
## Further Material

Papers for examining more closely some of the discussed topics.

- [Makarov et al. (2024) - "Robust 1-norm Periodograms for Analysis of Noisy Non-Gaussian Time Series with Irregular Cadences: Application to VLBI Astrometry of Quasars"](https://ui.adsabs.harvard.edu/abs/2024PASP..136e4503M/abstract)
- [Nilsson et al. (2018) - "Long-term optical monitoring of TeV emitting blazars. I. Data analysis"](https://ui.adsabs.harvard.edu/abs/2018A%26A...620A.185N/abstract)
- [Tarnopolski & Marchenko (2021) - "A Comprehensive Power Spectral Density Analysis of Astronomical Time Series. II. The Swift/BAT Long Gamma-Ray Bursts"](https://ui.adsabs.harvard.edu/abs/2021ApJ...911...20T/abstract)
- [Vio et al. (2010) - "Unevenly-sampled signals: a general formalism for the Lomb-Scargle periodogram"](https://ui.adsabs.harvard.edu/abs/2010A%26A...519A..85V/abstract)
- [Glynn et al. (2006) - "Detecting periodic patterns in unevenly spaced gene expression time series using Lomb–Scargle periodograms"](https://academic.oup.com/bioinformatics/article-pdf/22/3/310/48839224/bioinformatics_22_3_310.pdf)
- [Emmanoulopoulos et al. (2013) - "Generating artificial light curves: revisited and updated"](https://ui.adsabs.harvard.edu/abs/2013MNRAS.433..907E/abstract)
- [Perig et al. (2019) - "Periodicity in Volcanic Gas Plumes: A Review and Analysis"](https://ui.adsabs.harvard.edu/abs/2019Geosc...9..394P/abstract)
"""

# ╔═╡ 427b330c-e92a-4be3-8db2-4373f91f4d0e
md"""
### Credits
***

This notebook contains material obtained from [https://github.com/jakevdp/PracticalLombScargle/blob/master/figures/LombScargleVsClassical.ipynb](https://github.com/jakevdp/PracticalLombScargle/blob/master/figures/LombScargleVsClassical.ipynb), [https://github.com/jakevdp/PracticalLombScargle/blob/master/figures/Uncertainty.ipynb](https://github.com/jakevdp/PracticalLombScargle/blob/master/figures/Uncertainty.ipynb) and from [https://github.com/jakevdp/PracticalLombScargle/blob/master/figures/LINEAR_binary.ipynb](https://github.com/jakevdp/PracticalLombScargle/blob/master/figures/LINEAR_binary.ipynb).
"""

# ╔═╡ 6d25216d-ab50-4ce2-993b-84ce3070eb4e
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
    <td><a href="./open?path=Lectures/ScienceCase-X-RayBinaries/Lecture-X-RayBinaries.jl">Science case: X-ray binariess</a></td>
    <td><a href="./open?path=Lectures/ScienceCase-VariableStars/Lecture-VariableStars.jl">Science case about variable stars</a></td>
  </tr>
  <tr>
	<td>html</td>
    <td><a href="../../Lectures/ScienceCase-X-RayBinaries/Lecture-X-RayBinaries.html">Science case: X-ray binariess</a></td>
    <td><a href="../../Lectures/ScienceCase-VariableStars/Lecture-VariableStars.html">Science case about variable stars</a></td>
  </tr>
 </table>

"""

# ╔═╡ aa791077-751f-4542-8dc3-6cbb392c2302
md"""
**Copyright**

This notebook is provided as [Open Educational Resource](https://en.wikipedia.org/wiki/Open_educational_resources). Feel free to use the notebook for your own purposes. The text is licensed under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/), the code of the examples, unless obtained from other properly quoted sources, under the [MIT license](https://opensource.org/licenses/MIT). Please attribute the work as follows: *Stefano Covino, Time Domain Astrophysics - Lecture notes featuring computational examples, 2026*.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
CommonMark = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LombScargle = "fc60dff9-86e7-5f2f-a8a0-edeadbb75bd9"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.10.16"
CairoMakie = "~0.15.9"
CommonMark = "~1.0.1"
DataFrames = "~1.8.1"
FFTW = "~1.10.0"
LombScargle = "~1.0.3"
PlutoUI = "~0.7.79"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.5"
manifest_format = "2.0"
project_hash = "77adbe2868212e9cf160fdd50014e60a3af79d16"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "35ea197a51ce46fcd01c4a44befce0578a1aaeca"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.5.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "bca794632b8a9bbe159d56bf9e31c422671b35e0"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.3.2"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "8d8e0b0f350b8e1c91420b5e64e5de774c2f0f4d"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.16"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "fa072933899aae6dc61dde934febed8254e66c6a"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.15.9"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a21c5464519504e41e0cbc91f0188e8ca23d7440"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9cb23bbb1127eefb022b022481466c0f1127d430"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonMark]]
deps = ["PrecompileTools"]
git-tree-sha1 = "019ad9e55bb3549403f2d5a9b314fbb29a806ecb"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "1.0.1"

    [deps.CommonMark.extensions]
    CommonMarkMarkdownASTExt = "MarkdownAST"
    CommonMarkMarkdownExt = "Markdown"

    [deps.CommonMark.weakdeps]
    Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
    MarkdownAST = "d0879d2d-cac2-40c8-9cee-1863dc0c7391"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "3b4be73db165146d8a88e47924f464e55ab053cd"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.7"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d8928e9169ff76c6281f39a659f9bca3a573f24c"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.1"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "c55f5a9fd67bdbc8e089b5a3111fe4292986a8e8"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.6"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "fbcc7610f6d8348428f722ecbe0e6cfe22e672c6"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.123"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "66381d7059b5f3f6162f28831854008040a4e905"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+1"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6522cfb3b8fe97bec632252263057996cbd3de20"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.18.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.9.0"

    [deps.FilePaths.extensions]
    FilePathsGlobExt = "Glob"
    FilePathsURIParserExt = "URIParser"
    FilePathsURIsExt = "URIs"

    [deps.FilePaths.weakdeps]
    Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
    URIParser = "30578b45-9adc-5946-b283-645ec420af67"
    URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2f979084d1e13948a3352cf64a25df6bd3b4dca3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.16.0"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "1f5a80f4ed9f5a4aada88fc2db456e637676414b"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.10"

    [deps.GeometryBasics.extensions]
    GeometryBasicsGeoInterfaceExt = "GeoInterface"

    [deps.GeometryBasics.weakdeps]
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "93d5c27c8de51687a2c70ec0716e6e76f298416f"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dcc8d0cd653e55213df9b75ebc6fe4a8d3254c65"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.2.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
git-tree-sha1 = "8f3d257792a522b4601c24a577954b0a8cd7334d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.5"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "2cce1fed119ca7b6cc230c4a3b85202478af7924"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.3"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "d966f85b3b7a8e49d034d27a189e9a4874b4391a"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.13"

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

    [deps.IntervalSets.weakdeps]
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "4260cfc991b8885bf747801fb60dd4503250e478"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.11"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

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

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "97bbca976196f2a1eb9607131cb108c69ec3f8a6"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0205286d9eceadc518742860bf23f703779a3d6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.3+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LombScargle]]
deps = ["FFTW", "LinearAlgebra", "Measurements", "Random", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d64a0ce7539181136a85fd8fe4f42626387f0f26"
uuid = "fc60dff9-86e7-5f2f-a8a0-edeadbb75bd9"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "68af66ec16af8b152309310251ecb4fbfe39869f"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.9"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "7eb8cdaa6f0e8081616367c10b31b9d9b34bb02a"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.7"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf"]
git-tree-sha1 = "cb47f69a1cab9dcec7ff4a5d6e163410d6905866"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.14.1"

    [deps.Measurements.extensions]
    MeasurementsBaseTypeExt = "BaseType"
    MeasurementsJunoExt = "Juno"
    MeasurementsMakieExt = "Makie"
    MeasurementsRecipesBaseExt = "RecipesBase"
    MeasurementsSpecialFunctionsExt = "SpecialFunctions"
    MeasurementsUnitfulExt = "Unitful"

    [deps.Measurements.weakdeps]
    BaseType = "7fbed51b-1ef5-4d67-9085-a4a9b26f478c"
    Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f2b3b9e52a5eb6a3434c8cca67ad2dde011194f4"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.30+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "df9b7c88c2e7a2e77146223c526bf9e236d5f450"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.4.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e4cff168707d441cd6bf3ff7e4832bdf34278e4a"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.37"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3ac7038a98ef6977d44adeadc73cc6f596c08109"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.79"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "211530a7dc76ab59087f4d4d1fc3f086fbe87594"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.2.3"

    [deps.PrettyTables.extensions]
    PrettyTablesTypstryExt = "Typstry"

    [deps.PrettyTables.weakdeps]
    Typstry = "f0ed7684-a786-439e-b1e3-3b82803b501e"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "472daaa816895cb7aee81658d4e7aec901fa1106"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ebe7e59b37c400f694f52b58c93d26201387da70"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.9"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Statistics"]
git-tree-sha1 = "3949ad92e1c9d2ff0cd4a1317d5ecbba682f4b92"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5acc6a41b3082920f79ca3c759acbcecf18a8d78"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.7.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "d05693d339e37d6ab134c5ab53c29fce5ee5d7d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.4"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "a2c37d815bf00575332b7bd0389f771cb7987214"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.2"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "fa95b3b097bcef5845c142ea2e085f1b2591e92c"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.7.1"
weakdeps = ["Measurements", "StaticArraysCore", "Tables"]

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "08c10bc34f4e7743f530793d0985bf3c254e193d"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.8"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "57e1b2c9de4bd6f40ecb9de4ac1797b81970d008"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.28.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e015f211ebb898c8180887012b938f3851e719ac"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.55+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"
"""

# ╔═╡ Cell order:
# ╟─ff837227-c0bc-433c-b462-9e18ff3d947b
# ╟─53a01399-1f0b-4fbb-9966-97124e9eff42
# ╟─93f1b40e-8990-405d-bba1-9f98abd054ca
# ╟─8f9cd30e-5598-4e1b-b379-1473232ea867
# ╟─f7100fe5-2cfd-448d-abcc-66a623a8f11f
# ╟─31dd4ceb-8afd-4b55-aacb-2d888681ad1c
# ╟─e6386bbc-ea18-49a6-b64a-51b627d31565
# ╟─477704cc-f022-4cc5-a4a4-679387b819bd
# ╟─f95ec17d-9c56-4897-b5b2-65c3470fc294
# ╟─0f01200a-665e-4ae2-a6b9-b5c192955f3f
# ╟─a24629ad-7302-4aff-b7bf-fcc7635b9691
# ╟─533e41e2-08c3-450b-83eb-ae19664a6cfe
# ╟─ae7829b7-b652-4e54-9231-f7bc3c3b6093
# ╟─9a422815-783a-473c-8538-6ecf5c208bb0
# ╟─6fa31a81-cbbe-401b-acb0-0d44b0980363
# ╟─367c08ff-cedd-4e63-b53a-3178b0ec45cb
# ╟─ef919d43-ffb2-404f-94e4-a10b9a543e09
# ╟─ce5b3710-ac14-44f5-85d8-abc9756c1f11
# ╟─18958f04-3eb0-4e4e-9957-ca0c5540e5a5
# ╟─b985cda2-8242-468b-b5de-2299bfa50b33
# ╟─2e3e0747-48c2-43d1-abd4-b15f61a402bd
# ╟─ccb2f0e6-517b-4860-bac2-dac964e8d388
# ╟─0b098634-7890-4551-8343-34dbba3159ba
# ╟─2503892b-3bcf-41b8-a636-562d2613bd59
# ╟─f652aeb7-bcf0-481c-924c-412f406b3542
# ╟─637108d6-c977-4a2d-b312-e40bf375d5f1
# ╟─325a6b84-af05-4d78-94a2-3674a501c1c2
# ╟─50a1fb43-f72d-43bd-bea7-99bac76db19e
# ╟─7cd277fd-6b32-4ca8-9c8e-1eab55659c74
# ╟─46aea63c-b5a8-4977-8b87-f1cc479498e0
# ╟─1c023f88-4e2c-494b-9c03-002591955b60
# ╟─d5cd7643-6c82-4be1-a995-6e1f79547288
# ╟─ae2ac755-825b-4687-aeec-1de9a4ff24fb
# ╟─981eb94f-2bb9-4330-892b-9f0a5f929a00
# ╟─910a90a8-7b85-4d69-ba8d-3cabbfe89bc4
# ╟─8ca6fd17-566d-4aae-961c-cc4a4e6ab1d0
# ╟─2b0e312c-4285-4e63-9c68-b737646869d9
# ╟─171b2e3c-058e-4d82-9d0e-b5dd5f521404
# ╟─ec859add-e8e7-413b-b9d4-7d2b8694a46d
# ╟─eae2ed93-31fb-417a-9621-5586cf4f9732
# ╟─f38c7158-5206-467b-987d-38b4aa80e935
# ╟─ce90af7b-ef2c-4cdf-84cb-de76f5dd5c66
# ╟─0a342241-26b4-425e-b90f-fc555efb0b67
# ╟─17444051-6ee8-44af-8cfc-275bcb495d86
# ╟─467e6316-fa6d-4487-86cc-8ba408577fe7
# ╟─bd68480e-e2f6-4759-8333-066001c2ae19
# ╟─ec7fa3e8-2635-4af3-b2e7-6435e96ebcb2
# ╟─0dafbc5d-d383-44a8-8dc8-704f8a36f7fe
# ╟─fe1aa970-abc8-4720-9f09-5fe931100263
# ╟─9d413d79-c9fb-4c30-ac02-92380510b607
# ╟─26d6e75f-dd44-4876-ad72-e5606d27798f
# ╟─66da0ee4-4d50-47e6-ab43-3385b2a9f9fa
# ╟─e30cc0a2-ce7f-4e9f-9f3f-ebc14a76e936
# ╟─2ef299bf-e198-4bde-83cd-ecd86e566683
# ╟─6d2ccb03-4455-4249-9d3d-be4cb5d7a45f
# ╟─78119e66-703b-4305-85cd-956f1053cc60
# ╟─435fd658-c227-47ce-8dd6-cca20bff2a26
# ╟─427b330c-e92a-4be3-8db2-4373f91f4d0e
# ╟─6d25216d-ab50-4ce2-993b-84ce3070eb4e
# ╟─aa791077-751f-4542-8dc3-6cbb392c2302
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
