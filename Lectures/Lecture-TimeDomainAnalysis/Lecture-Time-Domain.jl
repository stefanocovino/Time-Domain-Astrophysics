### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 508fcfab-51c1-4dc5-992d-d5961517b9bf
begin
	using ARFIMA
	using CairoMakie
	using CommonMark
	using CSV
	using DataFrames
	using Distributions
	using HypothesisTests
	using LombScargle
	import Plots
	using PlutoUI
	using ShiftedArrays
	using StateSpaceModels
	using StatsBase
end

# ╔═╡ 14a0c11f-827d-46c8-bd5f-469a029f9afe
md"""
**What is this?**


*This jupyter notebook is part of a collection of notebooks on various topics discussed during the Time Domain Astrophysics course delivered by Stefano Covino at the [Università dell'Insubria](https://www.uninsubria.eu/) in Como (Italy). Please direct questions and suggestions to [stefano.covino@inaf.it](mailto:stefano.covino@inaf.it).*
"""

# ╔═╡ 2aeb4208-bc0c-4907-ab6d-da3f48210aa5
md"""
**This is a `Julia` notebook**
"""

# ╔═╡ fce942fc-7126-4e92-b758-30d36609117f
# ╠═╡ show_logs = false
md"""
$(LocalResource("Pics/TimeDomainBanner.jpg"))
"""

# ╔═╡ 97bebdb5-df0d-4eb2-b214-1a38dd9ffa85
# ╠═╡ show_logs = false
md"""
# Time-Domain analysis

***

- One can derive inferences about the spectral content of a time series by a spectral analysis.

- An analogous (i.e. complementary) analysis can be carried out in the time domain.

- Traditionally, analyses in the temporal domain are often aimed at deriving predictions, while in the spectral domain one often looks for periodicity, etc.

- There are indeed no fundamental technical reasons for that. 

$(LocalResource("Pics/timevsspectral.png"))

"""

# ╔═╡ 36377f43-7624-4113-9ea9-409fbd3ef079
# ╠═╡ show_logs = false
cm"""
## Time-Domain vs Spectral analysis: main tools
***


|  | **Time Domain**                  | **Frequency Domain**                            |
|:--:|:----------------------------------:|:-------------------------------------------------:|
|  | ``x_t`` linear combination of past | ``x_t`` linear combination of periodic components |
| Object of interest | population ACF | Spectral density |
| Data analysis tool | sample ACF | Periodogram |
| | | Identify dominant frequency(ies) |



"""

# ╔═╡ 91f7ef29-6570-44d1-a5e2-79658576d488
cm"""
### The Correlation Function
***

- One of the main statistical tools for the analysis of stochastic variability is the autocorrelation function. It represents a specialized case of the correlation function of two functions, ``f(t)`` and ``g(t)``, scaled by their standard deviations, and defined at time lag  as: 

```math
CF(\Delta t) = \frac{\lim_{T\to\infty} \int_{(T)} f(t) g(t+\Delta t)dt}{\sigma_f \sigma_g}
```

- where ``σ_f`` and ``σ_g`` are standard deviations of ``f(t)`` and ``g(t)``, respectively. 

- With this normalization, the correlation function is unity for ``Δt = 0``, without normalization by standard deviation, the above expression is equal to the covariance function.
"""

# ╔═╡ 2d403a10-ede5-4d72-a732-59d4715273bc
cm"""
### The Auto-Correlation Function
***

- It is assumed that both ``f`` and ``g`` are statistically weakly stationary functions (more on this topic later).

- With ``f(t)=g(t)=y(t)``, the autocorrelation of ``y(t)`` defined at time lag  is:

```math
ACF(\Delta t) = \frac{\lim_{T\to\infty} \int_{(T)} y(t) y(t+\Delta t)dt}{\sigma^2_y}
```

- The autocorrelation function and the PSD of function ``y(t)`` are Fourier pairs; this fact is known as the Wiener–Khinchin theorem and applies to stationary random processes.

- The sample auto-correlation function is defined as:

```math
ACF(k) = \frac{\sum_{t=1}^{n-k}(X_t - \bar{X})(X_{t+k} - \bar{X})}{\sum_{t=1}^k (X_t - \bar{X})^2}
```

- where the numerator is just the sample auto-covariance function and the denominator the sample variance. 

- ACF(k) is often called a correlogram.

- The ACF is fundamental measure of the serial correlation in a time series.

- It is defined for evenly spaced time-series, yet there are generalization to irregularly sampled data.

- Under the null hypothesis that the time series has no correlated structure and the population ACF is zero for all lags except for ACF(0) which is always unity.

- The distribution of the sample is asymptotically normal with mean ``−1/n`` and variance ``1/n``, i.e. the distribution of the null case ACF is ``ACF(k) = \mathcal{N}(−1/n,1/n)``. 

- However, this holds under the assumtpion one is testing whether a time-series is white noise, i.e. the residuals of a fit. Else, a different formula, based on the assumption that a time-series can be described as a moving average process (see later), is known as Bartlett's formula: `` \sigma_{\rm ACF(k)} ≈ (1 / \sqrt{N})  \sqrt{1 + 2\sum_{i=1}^{k-1} {\rm ACF(i)}^2}``.


    - Essentially, Bartlett's formula refines the confidence interval calculation for ACF by accounting for the potential dependence between ACF values in non-white noise processes, particularly those following a moving average structure.




- Simple quantities like the sample mean are valid estimates of the population mean for stationary processes, but its uncertainty is not the standard value when autocorrelation is present:

```math
\widehat{Var}(\bar{X}_n) = \frac{\sigma^2}{n} \left[1 + 2\sum_{k=1}^{n-1}(1-k/n){\rm ACF}(k) \right]
```

- Qualitatively, this can be understood as due to the decrease of independent measurements.

- This is indeed one very important result, although often (wrongly) ignored.

- Knowing the ACF one can compute the so-called integrated auto-correlation time ``\tau_{X,{\rm int}}``:

```math
\tau_{X,{\rm int}} = \int_0^\infty ACF (\tau) d\tau
```

- This allows one to define the *effective* number of samples in our dataset `` N_{\rm eff} = N / 2 \tau_{X,{\rm int}}`` so that the variance can be written as ``\widehat{Var}(\bar{X}_n) = \sigma^2 / N_{\rm eff} ``.

- The ACF is *non-negative definite*, i.e.: ``\sum_{i=1}^n \sum_{j=1}^n \alpha_j ACF(i-j) \alpha_j \ge 0`` for all positive integers ``n`` and vectors ``\alpha = (\alpha_1,\alpha_2,...,\alpha_n)' \in \mathbb{R}^n``. In fact, ``\sum_{i=1}^n \sum_{j=1}^n \alpha_j ACF(i-j) \alpha_j = Var \left( \sum_{i=1}^n \alpha_j x_j \right)``.
"""



# ╔═╡ 90b6cf29-44c0-425d-a814-910fb08009d2
cm"""
### The partial ACF (PACF)
***

- The PACF at lag k gives the autocorrelation at value k removing the effects of correlations at shorter lags.

- The value of the p-th coefficient PACF(p) is found by successively fitting autoregressive models (more on this topic later) with order ``1,2,...,p`` and setting the last coefficient of each model to the PACF parameter (this again will be clearer later).

- For instance, for a stationary AR(2) process:

```math
PACF(2) = \frac{ACF(2)-ACF(1)^2}{1-ACF(1)^2}
```

- The partial autocorrelation function of a stationary time series can be calculated by using the Durbin–Levinson Algorithm:

```math
\phi_{n,n} = \frac{\rho(n)-\sum_{k=1}^{n-1} \phi_{n-1,k}\rho(n-k)}{1-\sum_{k=1}^{n-1} \phi_{n-1,k}\rho(n-k)}
```

- where ``\phi_{n,k} = \phi_{n-1,k} - \phi_{n,n}\phi_{n-1,n-k}`` for ``1\le k\le n-1`` and ``\rho(n)`` is the ACF.
"""



# ╔═╡ 91c6247b-af12-47cf-811a-1440b6919576
cm"""
### Discrete Correlation Function
***

- The discrete correlation function is a procedure for computing the autocorrelation function that avoids interpolating an unevenly spaced dataset onto a regular grid.

- The method can treat both autocorrelation within one time series or cross-correlation between two unevenly spaced time series.

- Consider two datasets ``(x_i,t_{x_i})`` and ``(z_j,t_{z_j})`` with ``i = 1,2,...,n`` and ``j = 1,2,...,m`` points, respectively. For autocorrelation within a single dataset, let ``z = x``. 

- Construct two matrices, the **unbinned discrete correlation function** (UDCF) and its associated time lags:

```math
{\rm UDCF}_{i,j} = \frac{(x_i-\bar{x})(z_j-\bar{z})}{\sigma_x \sigma_z} \qquad \Delta_{ij}=t_j-t_i
```

- The UDCF is then grouped into a univariate function of lag-time ``τ`` by collecting the ``M(τ)`` data pairs with lags falling within the interval ``τ − Δτ /2 ≤ Δt_{ij} < τ + Δτ /2``. 

- The resulting discrete correlation function (DCF) and its variance are:

```math
{\rm DCF}(\tau) = \frac{1}{M(\tau)}\sum_{k=1}^{M(\tau)} {\rm UDCF}_{ij} \quad 
{\rm Var}(\tau) = \frac{1}{(M(\tau)-1)^2} \sum_{k=1}^{M(\tau)} [{\rm UDCF}_{ij} - {\rm DCF}(\tau) ]^2
```

"""

# ╔═╡ f1c10261-a2c7-473c-879b-005714127aee
cm"""
### Z-transformed DCF
***

- The drawback of the DFC is that its sample distribution of is known to be very skewed and far from normal. 

- If bins have an equal number of points (i.e. bins can have different length), The bin distribution becomes approximately binomial and, if a [Fisher z-transform](https://en.wikipedia.org/wiki/Fisher%27s_z-distribution) is applied:

```math
z(t) = \frac{1}{2} \ln \left( \frac{1+{\rm DCF}(\tau)}{1-{\rm DCF}(\tau)} \right)
```

- The ``z(τ)`` values are now normally distributed with known mean and variance (see [Alexander 1997](https://ui.adsabs.harvard.edu/abs/1997ASSL..218..163A/abstract)).
"""



# ╔═╡ b40d3fb9-f423-44d5-8a3e-32e15a6cc564
cm"""
### The Structure Function
***

- The structure function (sometimes called the Kolmogorov structure function) is again a measure of autocorrelation.

- The q-th order structure function is:

```math
D^q(\tau) = \langle|x(t)-x(t+\tau)|^q \rangle
```

- where the angular brackets ⟨ ⟩ indicate an average over the time series and q is the order of the structure function (not necessarily integer).

- The q = 2 structure function is also called the variogram. 

- When a dataset has a characteristic time-scale of variation ``τ_c``, the structure function is small at τ shorter than ``τ_c``, rises rapidly to a high level around ``τ_c``, and stays at this plateau for longer τ.

- If a structure function exhibits a power-law dependence on lag, ``D_q(\tau) \propto \tau^\alpha``, the time series has a multi-fractal behavior. In this case, the Wiener–Khinchin theorem links the power-law dependence of the q = 2 structure function to the power-law dependence of the Fourier power spectrum. 

- The structure function with q = 2 can also be seen as:

```math
SF(\Delta t) = SF_\infty[1-{\rm ACF}(\Delta t)]^{1/2}
```

- where ``SF_\infty`` is the standard deviation of the time series evaluated over an infinitely large time interval (or at least much longer than any characteristic timescale τ). 

- The structure function with q=2 is equal to the standard deviation of the distribution of the difference of ``y(t_2)-y(t_1)`` evaluated at many different ``t_1`` and ``t_2``. 

- As mentioned before, when the structure function ``SF \propto t^α``, then the ``{\rm PSD} \propto 1/f^{1+2α}``.

- The slope of the PSD directly affects the shape of the related time-series:

$(LocalResource("Pics/longmemory.png"))
"""

# ╔═╡ 40e140c2-a8fb-4e31-9636-a36104f60460
md"""
### Different stochastic processes can be categorized based on their ACF/PSD
***

- A stochastic process with $1/f^2$ spectrum is known as random walk (if discrete) or Brownian motion (or, more accurately, Wiener process) if continuous. These physically occur when the value being observed is subjected to a series of independent changes of similar size. It's also sometimes called as "red noise". Quasar variability (for instance) exhibits $1/f^2$ properties at high frequencies (that is, short time scales, below a year or so).

- A stochastic process with $1/f$ spectrum are sometimes called "long-term memory processes" (also sometimes know as "pink noise"). They have equal energy at all octaves (or over any other logarithmic frequency interval). This type of process has infinite variance and an undefined mean (similar to a Lorentzian distribution).

- A process with a constant PSD is frequently referred to as "white noise", i.e., it has equal intensity at all frequencies. This is a process with no memory, each measurement is independent of all others. White noise is identically and indepndently distributed (I.I.D.).
"""

# ╔═╡ 2f7623f2-359c-4c93-b6b2-ab48f606ec56
begin
	function GetACF(data, lags; sigma=1.96)
	    cc = StatsBase.autocor(data,0:lags)
	    return cc, -1/length(data)-sigma*sqrt(1/length(data)),-1/length(data)+sigma*sqrt(1/length(data))
	end
	
	function GetPACF(data, lags; sigma=1.96)
	    cc = StatsBase.pacf(data,0:lags)
	    return cc, -1/length(data)-sigma*sqrt(1/length(data)),-1/length(data)+sigma*sqrt(1/length(data))
	end

	function GetCrossCorr(x,y,lags)
		cc = StatsBase.crosscor(x, y, -lags:lags; demean=true)
		return cc
	end
end;

# ╔═╡ f3ed1ba4-849b-4713-8e48-4a537ee894d7
md"""
#### Exercize: cross-correlation of data from the Covid19 outbreak in Italy
***

- Italian data downloaded from: [https://github.com/pcm-dpc/COVID-19](https://github.com/pcm-dpc/COVID-19) and updated till April 07, 2022.
"""

# ╔═╡ f0e73dd6-2f4f-4dd6-9426-f0a739956602
begin
	nep = DataFrame(CSV.File("Data_Italy.csv"))
	first(nep,5)
end

# ╔═╡ dcb69463-c144-4726-8336-5f3851873d8a
md"""
- We need to slighly rearrange the data in order to compute the number of deaths per day.
"""

# ╔═╡ f8b03f2b-a241-4d77-8fc9-4c24a4c8ef90
begin
	dDead = diff(nep[!,"Dead"])
	nep[!,"dDead"] = pushfirst!(dDead,0)
end;

# ╔═╡ 577fa116-fa14-4077-8fc8-4d3ea544f410
md"""
- We need to slighly rearrange the data in order to compute the number of deaths per day.
"""

# ╔═╡ 20d1925b-ddb8-49ec-86a2-98f2efac3532
begin
	fg1 = Figure()
	
	ax1fg1 = Axis(fg1[1, 1],
	    xlabel="Time (days)",
	    ylabel="Log N",
	    yscale=log10,
	    title="Infected and deaths per day - Italy"
	    )
	
	flt = nep[!,"dDead"] .> 0
	
	lines!(nep[!,"Date"][flt],nep[!,"New Infected"][flt],label="Infected/day")
	lines!(nep[!,"Date"][flt],nep[!,"dDead"][flt],label="Deaths/day")
	
	
	axislegend(position=:lt)
	
	fg1
end

# ╔═╡ 12225d55-9c0a-482c-890b-8ca5f8139c5a
md"""
- There is a clear correlation between the number of infected people and deaths per day. Superposed to a seasonal, and shorter, variability together with "spikes" likely due to irregularities in recording the data.

- Lets compute a cross-correlation between these two datasets in order to try to quantify a delay.

> Disclaimer. This is just an exercize. Please, do not underestimate the complexity of a proper epidemiological study.
"""

# ╔═╡ db5d7b3e-66da-4498-a924-2f9b999fda56
begin
	fg2 = Figure()
	
	ax1fg2 = Axis(fg2[1, 1],
	    xlabel="Lags (days)",
	    ylabel="Cross-Correlation",
	    #yscale=log10,
	    #title="Infected and deaths per day - Italy"
	    )

	ccr = GetCrossCorr(nep[!,"dDead"],nep[!,"New Infected"],50)
	
	lines!(-50:50,ccr)
	
	
	#axislegend(position=:lt)
	
	fg2
end

# ╔═╡ 06a37524-c82e-4f68-8529-89759cdda2bd
md"""
- A $\sim$15-20 day delay appears, that is quite reasonable. 

- Of course, a much more meaningful analysis would have required to separate the infected peole basing upon age, sex, etc.
"""

# ╔═╡ 06f77be0-51c6-4ca2-82b1-eaa69d9481b7
md"""
##### Epidemic periodicity?
***

- The former plot might suggest some sort of quasi-periodicity of the outbreak. 

- Let's see whether a Lomb-Scargle analysis produce interestibg results:
"""

# ╔═╡ b1f35602-c627-41a7-8754-6634d5eff5a7
begin
	t = Float64.(collect(1:nrow(nep)))
	y = Float64.(nep[!,"New Infected"])
	lsres = lombscargle(t,y;samples_per_peak=10,minimum_frequency=1/700,maximum_frequency=1/2)
	
	fg3 = Figure()
	
	ax1fg3 = Axis(fg3[1, 1],
	    xlabel="Period (days)",
	    ylabel="Power",
	    #yscale=log10,
	    )
	
	
	lines!(1 ./lsres.freq,lsres.power,label="LS periodogram")
	vlines!(365.,linestyle=:dash)
	
	axislegend(position=:lt)
	
	fg3
end

# ╔═╡ c097e9b3-f9dc-473e-902a-cc5b03c15da3
md"""
- The LS periodogram shows a considerable power at long periods, with a large peak at about 400 days. This is likely the summer/winter cycle that did not reproduce exactly.

- The time-series shows features at shorter periods too.
"""

# ╔═╡ 08a6a975-b496-4d19-8e8b-5df7aae43ad9
begin
	t2 = Float64.(collect(1:nrow(nep)))
	y2 = Float64.(nep[!,"New Infected"])
	lsres2 = lombscargle(t2,y2;samples_per_peak=10,minimum_frequency=1/30,maximum_frequency=1/2)
	
	fg4 = Figure()
	
	ax1fg4 = Axis(fg4[1, 1],
	    xlabel="Period (days)",
	    ylabel="Power",
	    #yscale=log10,
	    )
	
	
	lines!(1 ./lsres2.freq,lsres2.power,label="LS periodogram")
	
	xlims!(1,10)
	
	axislegend(position=:lt)
	
	fg4
end

# ╔═╡ 2acb3d30-1b38-4255-b955-a0db9f77be72
md"""
- A distinct peat at a week is clearly visible, indicating the "periodicity" in the collection of data. 
"""

# ╔═╡ f34215b4-91fd-4175-96f1-524097ec7160
# ╠═╡ show_logs = false
md"""
## Stationarity
***

- This is a fundamental concept. A stationary time-series is a dataset where it is “meaningful” to draw predictions.

-  A stationary series is one in which the properties – mean, variance and covariance, do not vary with time.
    - Let us understand this using an intuitive example. Consider the three plots shown below:
    
$(LocalResource("Pics/threestationary.png"))

- In the first plot, we can clearly see that the mean varies (increases) with time which results in an upward trend. Thus, this is a non-stationary series. For a series to be classified as stationary, it should not exhibit a trend.
- Moving on to the second plot, there certainly isn't a trend in the series, but the variance of the series is a function of time. As mentioned previously, a stationary series must have a constant variance.
- The the third plot. The spread becomes closer as the time increases, which implies that the covariance is a function of time.

- The three examples shown above represent non-stationary time series.
    
$(LocalResource("Pics/truestationary.png"))

- In this case, the mean, variance and covariance are constant with time. This is what a stationary time series looks like.

- Most statistical models require the series to be stationary to make effective and precise predictions.
"""

# ╔═╡ 45f962f9-79be-492a-8add-737fc73925c2
cm"""
- Let's now try to define "stationarity" in a more formal way.

- A stochastic process ``X_t; t = 0, ±1,...`` is stationary if a joint distribution of ``X_t,..., X_{t+k}`` is same as a joint distribution of ``X_0,...,X_k`` for all t and all k. 

    - This is typically called strict stationarity. 

- A stochastic process ``X_t; t = 0, ±1,...`` is weakly (or second order) stationary if, for all t:

```math
\mathbb{E}(X_t) = \mu
```

- and ``{\rm Cov}(X_t, X_{t+h}) = {\rm Cov}(X_0, X_h) = C_X(h)`` is a function of h only, i.e. it does not depend on t. 
    - ``C_X`` is also known as autocovariance of X. ``C_X(h)/C_x(0)`` is known as autocorrelation.

- Let’s recall that: ``{\rm Cov}(X, Y) \equiv \mathbb{E}[(X − \mathbb{E}[X])(Y − \mathbb{E}[Y])]``

- It can be shown that:
    - If ``X_t`` has finite variance, strictly stationary implies ``X_t`` is weakly stationary.
    - If ``X_t`` is second order stationary and Gaussian it implies ``X_t`` is also strictly stationary.

- ``X_t`` is Gaussian if for each ``t_1,...,t_k`` the vector ``(X_{t_1},...,X_{t_k})^T`` has a Multivariate Normal Distribution.

- The process ``X_t`` is said to have a stationary covariance if: ``{\rm Cov}(X_t, X_s) = {\rm Cov}(X_{t+1}, X_{s+1}) = {\rm Cov}(X_{t+2}, X_{s+2})...``, etc.
"""

# ╔═╡ 5c1eff4e-8079-4568-a792-7528465db277
md"""
#### Exercize about time-series stationarity tests
***

- First read and pre-process a dataset reporting the [number of passengers on airplanes](https://www.analyticsvidhya.com/wp-content/uploads/2018/09/AirPassengers.csv) in the USA for several years. 

    - This is a simple yet interesting dataset for this kind of exercise. 
"""

# ╔═╡ 9057a011-76b9-475a-906c-81afd18e99ce
begin
	train = DataFrame(CSV.File("AirPassengers.csv"))
	first(train,5)
end

# ╔═╡ 3da138ae-cd0b-4f7b-ab5a-319d4f65dee2
md"""
- We want to determine whether a given series is stationary or not and deal with it accordingly. 

- Let's first try with a simple, yet always useful, visual test.
"""

# ╔═╡ 25a8201f-402b-4868-9cd6-262a03b980fb
begin
	fg5 = Figure()
	
	ax1fg5 = Axis(fg5[1, 1],
	    xlabel="Time",
	    ylabel="N",
	    )
	
	
	lines!(train[!,:Month],train[!,"#Passengers"])
	
	#xlims!(1,10)
	
	#axislegend(position=:lt)
	
	fg5
end

# ╔═╡ 698e9cc6-9de5-4ca6-b425-5bb038b08b2b
cm"""
- We see that there a clear trend, i.e. the mean is varying.

- It is anyway useulf to look for more formal tests.


##### Statistical tests
***

- We can use statistical tests like the unit root stationary tests. 
    - Unit root indicates that the statistical properties of a given series are not constant with time, which is the condition for stationary time series. 
    
- Suppose we have a time series:

```math
y_t = a*y_{t-1} + ε_t 
```

- where y``_t`` is the value at the time instant t and ε``_t`` is the error term. In order to calculate y``_t`` we need the value of y``_{t-1}``, which is:

```math
y_{t-1} = a*y_{t-2} + ε_{t-1} 
```

- If we do that for all observations, the value of ``y_t`` will come out to be:

```math
y_t = a^n*y_{t-n} + \sum ε_{t-i}*a^i 
```

- If the value of ``a`` is 1 (unit) in the above equation, then the predictions will be equal to the y``_{t-n}`` and sum of all errors from t-n to t, which means that the variance will increase with time. 

- This is knows as unit root in a time series. The unit root tests check the presence of unit root in the series by checking if value of a=1. 

- The ADF ([Augmented Dickey-Fuller](https://en.wikipedia.org/wiki/Augmented_Dickey%E2%80%93Fuller_test)) is one of the most commonly used unit root stationary tests:
  - It can be used to determine the presence of unit root in the series, and hence help us understand if the series is stationary or not. The null and alternate hypothesis of this test are:
    - Null Hypothesis: The series has a unit root (value of a=1). Alternate Hypothesis: The series has no unit root.
    - If we fail to reject the null hypothesis, we can say that the series is non-stationary. 
"""

# ╔═╡ f0996cba-03af-4247-99f7-d148c6e737a0
ADFTest(train[!,"#Passengers"],:none,1)

# ╔═╡ 0f220d9a-ffe5-4d84-ba4b-a65c5dbf4404
md"""
- The ADF tests gives the following results – test statistic, p value and the critical value at 1%, 5% , and 10% confidence intervals. 

- If the test statistic is less than the critical value, we can reject the null hypothesis and the series can be stationary. When the test statistic is greater than the critical value, we fail to reject the null hypothesis, which means the series is not stationary.

    - In our above example, the test statistic >> critical value, which implies that the series is not stationary, as expected.
"""

# ╔═╡ 59175381-65a4-45e7-ae85-a00467c4d181
md"""
### Types of Stationarity
***

- Strict Stationary: A strict stationary series satisfies the mathematical definition of a stationary process. For a strict stationary series, the mean, variance and covariance are not the function of time. The aim is to convert a non-stationary series into a strict stationary series for making predictions.
- Trend Stationary: A series that has no unit root but exhibits a trend is referred to as a trend stationary series. Once the trend is removed, the resulting series will be strict stationary. The KPSS test classifies a series as stationary on the absence of unit root. This means that the series can be strict stationary or trend stationary.
- Difference Stationary: A time series that can be made strict stationary by differencing falls under difference stationary. ADF test is also known as a difference stationarity test.
"""

# ╔═╡ 9a8ac940-8a2c-4d5a-809f-bbdb1dba12dc
md"""
#### Exercize about making a time-series stationary
***

- Is it possible to "stationarize" a time-series? 
    - In several cases it is. For instance if the non-stationarity is due to a trend, etc.
    
- Trends can be removed by differencation, i.e. we compute the difference of consecutive terms in the series. Differencing is typically performed to get rid of the varying mean: $ y_t‘ = y_t – y_{(t-1)}$, where y$_t$ is the value at a time t.

- Let's apply a differentiation on our series and plotting the results:
"""

# ╔═╡ 76917ce4-5249-4cbb-bdb7-da654ce2c75d
#transform!(train, "#Passengers" => ShiftedArrays.lag => "#Passengers_shift")
train[!,"#Passengers_diff"] = train[!,"#Passengers"] - ShiftedArrays.lag(train[!,"#Passengers"]);

# ╔═╡ a7c1c07b-23d7-4721-ac7d-70012483b056
begin
	fg6 = Figure()
	
	ax1fg6 = Axis(fg6[1, 1],
	    xlabel="Time",
	    ylabel="N",
	    )
	
	
	lines!(train[!,:Month],train[!,"#Passengers_diff"])
	
	fg6
end

# ╔═╡ cc765f72-cb9d-4b1b-9e11-f2e7b6a4c71e
md"""
- The average now seems to be fairly constant, but variance and probably covariance are not.

- It is possiboe to compute te differece with longer lags than 1. This is often called *seasonal differencing*.

- In seasonal differencing, instead of calculating the difference between consecutive values, we calculate the difference between an observation and a previous observation from the same "season": y$_t$‘ = y$_t$ – y$_{(t-n)}$.

"""

# ╔═╡ 7ecf612d-04ea-4e3b-b176-98476b443688
train[!,"#Passengers_diff"] = train[!,"#Passengers"] - ShiftedArrays.lag(train[!,"#Passengers"],6);

# ╔═╡ 7da9155f-27d8-48c0-93a6-13a3297fd302
begin
	fg7 = Figure()
	
	ax1fg7 = Axis(fg7[1, 1],
	    xlabel="Time",
	    ylabel="N",
	    )
	
	
	lines!(train[!,:Month],train[!,"#Passengers_diff"])
	
	fg7
end

# ╔═╡ 72700334-0d07-4486-ba3f-46fe0cfacb4f
md"""
- The curve is now cleaner, yet still far from being, even visually, stationary.

- Together in alternatve to differencing it is also possible to transform the data. 

- Transformations are used to stabilize the non-constant variance of a series. Common transformation methods include power transform, square root, and log transform. 
"""

# ╔═╡ 6021c695-2cd9-499c-bbbf-7d86b5a6b347
begin
	train[!,"#Passengers_log"] = log.(train[!,"#Passengers"])
	train[!,"#Passengers_log_diff"] = train[!,"#Passengers_log"] - ShiftedArrays.lag(train[!,"#Passengers_log"],6)
end;

# ╔═╡ 3d013635-b4ac-49ec-85f6-93c26cf7539c
begin
	fg8 = Figure()
	
	ax1fg8 = Axis(fg8[1, 1],
	    xlabel="Time",
	    ylabel="N",
	    )
	
	
	lines!(train[!,:Month],train[!,"#Passengers_log_diff"])
	
	fg8
end

# ╔═╡ c81ad1a5-b4ed-4056-a161-eab4f9e981d9
md"""
- The improvement seems to be significan over the previous plots. Let's apply the ADF stationarity test.

- In principle, a seasonality with one year (12 months) period could be preferred here. But we see that even a 6 month period does a good job in making the time-series stationary.
"""

# ╔═╡ f1bc344e-ac85-4815-928f-61f5f48c18a9
ADFTest(dropmissing(train)[!,"#Passengers_log_diff"],:none,1)

# ╔═╡ 381bc481-5c7e-4d9a-80c9-8f245718718d
md"""
- Now the criterion is reasonably satisfied. The difference of the log-transformed original time-series is stationary.
"""

# ╔═╡ ba0761ce-3fd9-4b28-a9b4-46f850dfbafe
md"""
### Covariance matrix
***

- Let's take the opportunity to formally define the covariance matrix $\Sigma$, whose components are $\sigma_{i,j} = {\rm Cov}(X_i,X_j)$:

```math
\Sigma=\left[
\begin{array}{ccc}
   \sigma_{11} & \cdots & \sigma_{1n} \\
   \vdots & \ddots & \vdots \\
   \sigma_{n1} & \cdots & \sigma_{nn}
\end{array}
\right]
```

- This kind of matrix is also a “Toeplitz” matrix (it has useful properties for massive computations).
"""

# ╔═╡ 9d1af19c-e1b3-48aa-9b9b-245347efd682
cm"""
- Let's now study some simple linear process of interest.

### White Noise
***

- ``\{w_t\}`` is a white noise process if ``w_t`` are uncorrelated and identically distributed random variables with:
    - ``\mathbb{E}[w_t] = 0`` and Var``[w_t] = σ^2``, for all t.
    
- If ``\{w_t\}`` are Normally (Gaussian) distributed, we call this Gaussian white noise.

- A white noise process is stationary.

- This is an example with ``σ^2 = 1``:
"""

# ╔═╡ e410b2af-6665-4ddb-bb89-4b38bec8bfc2
begin
	N9 = 10000
	
	d9 = Normal()
	x9 = rand(d9, N9)
end;

# ╔═╡ 2aca0524-9737-43e9-b456-9fa0b2f7939b
begin
	fg9 = Figure()
	
	ax1fg9 = Axis(fg9[1, 1],
	    )
	
	
	lines!(x9,label="sigma=1")
	
	axislegend()
	xlims!(0,N9)
	
	fg9
end

# ╔═╡ 0c8e81ce-7b48-427c-b72e-4fdad4015d51
md"""
- The ACF of a white noise shows a distinct and easy to recognize pattern: there should be no significant correlation for any lag but 0.
"""

# ╔═╡ 63a3b6f2-a46a-4a4c-ba9e-5639714cd93d
begin
	rs10 = GetACF(x9,40)
	fg10 = Figure()
	axfg10 = Axis(fg10[1, 1],)
	stem!(rs10[1])
	hlines!([rs10[2],rs10[3]],linestyle=:dash)
	fg10
end

# ╔═╡ d674149a-25a7-4c60-8770-e2b46cb2ed74
cm"""
### Random walk (with drift)
***

- A random walk, with drift, can be defined as:

```math
x_t = \delta + x_{t-1} + w_t 
```

- where ``δ`` is a constant, ``w_t`` is a white noise process, and we assume for simplicity ``x_0 = 0``.

- The equation can be rewritten as: `` x_t = t\delta + \sum_{j=1}^t w_j`` simply writing ``x_1 = \delta + w_1``, ``x_2 = \delta + x_1 + w_2 = \delta + \delta + w_1 + w_2``, etc.

- Is a random walk (with or without drift) stationary?

- The expectation value is simple to compute: `` \mathbb{E}[x_t] = \mathbb{E}[t\delta + \sum_{j=1}^t w_j] = t\delta ``, since a sum of random noise should give a zero result.  
    - Therefore, unless ``δ=0`` the expectation value depends on t.

- For the covariance we need some more steps: 
	- ``{\rm Cov}(x_t,x_{t+h}) = {\rm Cov}[t\delta + \sum_{j=1}^t w_j,(t+h)\delta + \sum_{j=1}^{t+h} w_j] ``. 
	- Recalling that ``{\rm Cov}(\sum_i x_, \sum_j x_j) = \sum_i \sum_j {\rm Cov}(x_i, x_j)`` and that ``{\rm Cov}(X,X) = {\rm Var}(X)`` we get that: 
	- ``{\rm Cov}(x_t,x_{t+h}) = {\rm Cov}[t\delta,(t+h)\delta] +  {\rm Cov}[\sum_{j=1}^t w_j,(t+h)\delta] + {\rm Cov}[t\delta,\sum_{j=1}^{t+h} w_j] + {\rm Cov}[\sum_{j=1}^t w_j,\sum_{j=1}^{t+h} w_j]``. 
	- All terms but the last are 0 since they include constant terms. We finally have: `` {\rm Cov}(x_t,x_{t+h}) = t{\rm Var}(w_t) = t\sigma^2``.
    - Clearly non-stationary (for any ``\delta``).
    
- And a random walk with drift (``\delta=0.1``) can be:
"""

# ╔═╡ ea14bb29-f6bf-4a0c-a6b1-d38e285207ee
begin
	N10 = 10000
	delta10 = 0.1
	
	d10 = Normal()
	sigma10 = rand(d10, N10)
	
	x10 = zeros(N10)
	for i in range(2,N10)
	    x10[i] = x10[i-1]+sigma10[i]+delta10
	end
end

# ╔═╡ 95047c44-9d05-4dc3-a50a-aa87cf9cd2c4
begin
	fg11 = Figure()
	
	ax1fg11 = Axis(fg11[1, 1],
	    )
	
	
	lines!(x10,label="delta="*string(delta10))
	
	axislegend(position=:lt)
	xlims!(0,N10)
	
	fg11
end

# ╔═╡ 61d7c290-6d4d-4da8-a4f6-f0c2c795a4f8
md"""
### Moving average (of the first order)
***

- A moving average process of the first order, or MA(1), is defined as:

```math
x_t = \beta_1 w_{t-1} + w_t
```

- where again $w_t$ is a white noise process.

- The expectation value of a MA(1) process is: $\mathbb{E}[x_t] = \mathbb{E}[\beta_1 w_{t-1} + w_t] = \beta_1\mathbb{E}[w_{t-1}] + \mathbb{E}[w_t] = 0$.

- The variance is: ${\rm Var}[x_t] = {\rm Var}[\beta_1 w_{t-1} + w_t] = \beta_1^2 {\rm Var}[w_{t-1}] + {\rm Var}[w_t] = \sigma^2 (1 + \beta_1^2) $. 

- The covariance is: ${\rm Cov}[x_t,x_{t+h}] = {\rm Cov}[\beta_1 w_{t-1} + w_t,\beta_1 w_{t+h-1} + w_{t+h}]$. 
    - Now, if $h=1$ we have: ${\rm Cov}[x_t,x_{t+1}] = {\rm Cov}[\beta_1 w_{t-1} + w_t,\beta_1 w_{t} + w_{t+1}] = {\rm Cov}[\beta_1 w_{t-1},\beta_1 w_{t}] + {\rm Cov}[\beta_1 w_{t-1},w_{t+1}] + {\rm Cov}[w_t,\beta_1 w_{t}] + {\rm Cov}[w_t,w_{t+1}]$. All terms with different index are zero (no covariance between independent white noise processes), and therefore: ${\rm Cov}[x_t,x_{t+1}] = \beta_1^2 \sigma^2$.  
    - If $h > 1$ all terms are zero and theerfore: ${\rm Cov}[x_t,x_{t+1}] = 0$.

- A MA(1) with $\beta_1 = 1$ can be:
"""

# ╔═╡ 962d9bb2-8fee-4d94-9054-a48f084127d7
begin
	N12 = 10000
	beta112 = 1.
	
	d12 = Normal()
	sigma12 = rand(d12, N12)
	
	x12 = zeros(N12)
	for i in range(2,N12)
	    x12[i] = sigma12[i]+beta112*sigma12[i-1]
	end
end

# ╔═╡ 33df98ca-16f8-4bc7-8848-494d920e60d9
begin
	fg12 = Figure()
	
	ax1fg12 = Axis(fg12[1, 1],
	    )
	
	
	lines!(x12,label="beta1="*string(beta112))
	
	axislegend()
	xlims!(0,N12)
	
	fg12
end

# ╔═╡ d4d93ec8-e268-4fc9-b5ee-cdd194636659
md"""
- Since the process is a sum of independent noise processes it is visually very close to a random noise.

- Let's analyse the ACF of the previous process and discover this is not fully the case:
"""

# ╔═╡ 9512aa2e-afaa-4ddd-aac4-826abae447ff
begin
	rs13 = GetACF(x12,40)
	fg13 = Figure()
	axfg13 = Axis(fg13[1, 1],)
	stem!(rs13[1])
	hlines!([rs13[2],rs13[3]],linestyle=:dash)
	fg13
end

# ╔═╡ d0dc8d40-75d7-4611-aab6-d8542f82947e
md"""
- Apart from lag=0, the only lag showing correlation significantly different from 0 is lag=1, the order of the process. This is not, as we will see, by chance.
"""

# ╔═╡ b0ae3114-1aca-4109-922a-02c25989e407
md"""
### Autoregressive (of the first order)
***

- An autoregressiv process of the first order, or AR(1), is defined as:

```math
x_t = \alpha_1 x_{t-1} + w_t

- and again $w_t$ is a white noise process.

- The process can be rewritten as follows: $ x_t = \alpha_1 x_{t-1} + w_t = \alpha_1 (\alpha_1 x_{t-2} + w_{t-1}) + w_t = \alpha_1^t x_0 + \sum_{i=o}^{t-1} \alpha_1^i w_{t-i}$

- The expectation value is: $\mathbb{E}[x_t] = \mathbb{E}[\alpha_1^t x_0 + \sum_{i=o}^{t-1} \alpha_1^i w_{t-i}] = \alpha_1^t \mathbb{E}[x_0]$. This is zero only if the starting position is zero.

- The variance is: ${\rm Var}[x_t] = {\rm Var}[\alpha_1 x_{t-1} + w_t] = \alpha_1^2 {\rm Var}[x_{t-1}] + {\rm Var}[w_t] \Rightarrow  {\rm Var}[x_t] = \frac{\sigma^2}{1-\alpha_1^2}$. This relation gives a valid variance only if $|\alpha_1| < 1$. 

- The covariance is: ${\rm Cov}[x_t,x_{t+h}] = {\rm Cov}[x_t,\alpha_1^{t+h} x_0 + \sum_{i=o}^{t+h-1} \alpha_1^i w_{t+h-i}] = {\rm Cov}[x_t,\alpha_1^h x_t + \sum_{i=o}^{h-1} \alpha_1^i w_{t+h-i}] = {\rm Cov}[x_t,\alpha_1^h x_t] = \alpha_1^{2h} {\rm Var}[x_t] = \frac{\alpha_1^{2h} \sigma^2}{1-\alpha_1^2}$, again valid if $|\alpha_1| < 1$.

- And, therefore, the correlation turns out to be: ${\rm Corr}[x_t,x_{t+h}] = \frac{{\rm Cov}[x_t,x_{t+h}]}{{\rm Std}[x_t] {\rm Std}[x_{t+h}]} = \alpha_1^{2h}$.

- In AR(1) process, the value of $\alpha_1$ determines whether the AR(1) process is stationary. 

- Let's see a AR(1) with $\alpha_1 = 0.9$.
"""

# ╔═╡ 924e216e-9ba7-4906-9d70-94f1813911e7
begin
	N14 = 10000
	alpha114 = 0.9
	
	
	d14 = Normal()
	sigma14 = rand(d14, N14)
	
	x14 = zeros(N14)
	for i in range(2,N14)
	    x14[i] = alpha114*x14[i-1]+sigma14[i]
	end
end

# ╔═╡ 30ece656-c5ba-42a6-8de2-e9827d72eabf
begin
	fg14 = Figure()
	
	ax1fg14 = Axis(fg14[1, 1],
	    )
	
	
	lines!(x14,label="alpha1="*string(alpha114))
	
	axislegend()
	xlims!(0,N14)
	
	fg14
end

# ╔═╡ b86d10a0-f3a0-4489-b595-64c3eafcfef2
md"""
- As one possibly infers looking at the plot the time-series shows fluctuations that do not apper totally random.

- Let's analyse the ACF of the previous process:
"""

# ╔═╡ 8a8d5e0e-7ccf-43a2-b18b-c389971252e9
begin
	rs15 = GetACF(x14,40)
	fg15 = Figure()
	axfg15 = Axis(fg15[1, 1],)
	stem!(rs15[1])
	hlines!([rs15[2],rs15[3]],linestyle=:dash)
	fg15
end

# ╔═╡ c6bab5e0-4a61-4a6c-aa54-34e0e4f2dcdc
md"""
- At variance with MA processes, a AR process shows a correlation decreases with increasing lag as a power-law.
"""

# ╔═╡ 2905155b-0de7-4f15-acb0-aaaf1b421c9b
# ╠═╡ show_logs = false
md"""
- The ACF is indeed a diagnostic tool to infer the nature and the order of a linear process:


$(LocalResource("Pics/acftab1.png"))

- Try to guess the kind of process looking at the plots below:

$(LocalResource("Pics/test1.jpg"))
$(LocalResource("Pics/test2.jpg"))


1. Upper left: The ACF shows a decreasing correlation with alternate signs. It is a MA(10) process. 
2. Upper right: The ACF shows a decreasing correlation as a power-law. Probably it is an AR process, altgough we cannot say the order. In principle it could also be a MA(4) process.
3. Bottom left: The ACF shows no correlation beyondd lag=0, it is a white noise process.
4. Bottom right: The ACF shows a 0 correlation after lag=1. It is thus a MA(1) process, although it might also be a AR process.

- As we have seenm the ACF is a powerful diagnostics but often we have ambiguous situations. The application of the PACF will help us to solve these ambiguities, as we are going to see later.
"""

# ╔═╡ fa1e8adf-969a-4476-81ae-b548cabcd8f6
md"""
#### The Autoregressive Moving Average scheme
***

- This approach were pioneered by [Yule](https://en.wikipedia.org/wiki/Udny_Yule) and [Slutsky](https://en.wikipedia.org/wiki/Eugen_Slutsky) in the 1920’s.
    - Yule’s researches led to the notion of the autoregressive scheme.
    - Slutsky’s researches led to the notion of a moving average scheme.

- We have an AutoRegressive scheme when a time series $x_t$ is assumed to be generated as a linear function of its past values, plus a random shock:

```math
x_t = φ_1 x_{t−1} + ... + φ_p x_{t−p} + u_t
```

-Conceptualy, an autoregressive scheme is one with a ‘memory’ in the sense thar each values is correlated with p preceding values. The constants $φ_1,...,φ_p$ are weights measuring the influence of preceding values $x_{t−1},...x_{t−p}$ on the value $x_t$.

- We have a Moving Average scheme when a  time series $x_t$ is assumed to be generated as a weighted linear sum of the last $q + 1$ random shocks:

```math
x_t = u_t + θ_1 u_{t−1} + ... + φ_p u_{t−q}
```

- If both schemes, the autoregressive and the moving average one, are used, we obtain the so-called Autoregressive Moving Average scheme:

```math
x_t = φ_1 x_{t−1} + ... + φ_p x_{t−p} + u_t + θ_1 u_{t−1} + ... + φ_p u_{t−q}
```
"""

# ╔═╡ 0edac0fb-f52e-4a2c-b3dc-54f8237ad325
md"""
## General linear processes
***

- The examples we have seen so far are all part of a general family

- We call “linear process” a linear combination of noise variates $w_t$:

```math
x_t = \sum_{i=0}^\infty \psi_i w_{t-i}
```

- For this kind of processes, if $\sum_{i=0}^\infty |\psi_i| < \infty$ then the process is stationary. 

- The autocovariance of this process can be written as: 

```math
\rho(h) = \sigma^2 \sum_{i=0}^\infty \psi_{i+h} \psi_i
```

- Now we provide a set of useful definitons.
"""

# ╔═╡ 2eb45539-ef03-4f6d-86be-52a1e9228088
md"""
### Backshift operator
***

- The backshift operator, **B**, is defined as:

```math
B x_t = x_{t-1}
```

- It can be extended to more powers as: $B^2 x_t = (B B) x_t = B (B x_t) = B x_{t-1} = x_{t-2}$, so that:

```math
B^k x_t = x_{t-k}
```


### Difference operator
***

- The difference operator, $\nabla$, is defined as:

```math
\nabla^d x_t = (1-B)^d x_t
```

- e.g., $\nabla^1 x_t = (1-B)^1 x_t = x_t - x_{t-1}$.

- The different operator is often used to convert non-stationary time-series to stationary.
"""

# ╔═╡ db519b68-35a3-4396-a302-e0bd804985c1
md"""
#### Exercise: let's convert MA(1) and AR(1) processes to the general linear form
***

- MA(1) can be expressed as: $x_t = \beta_1 w_{t-1} + w_t = \sum_{i=0}^\infty \psi_i w_{t-i}$, with $\psi_0 = 1, \psi_1 = \beta_1$ and $\psi_i = 0$ for $i \ge 2$.  

- AR(1) can be expressed as: $x_t = \alpha_1 x_{t-1} + w_t = w_t + \alpha_1 w_{t-1} + \alpha_1^2 w_{t-2}+...$ since $x_{t-1} = \alpha_1 x_{t-2} + w_{t-1}$ etc. 
    - We thus have $x_t = \sum_{i=0}^\infty \psi_i w_{t-i}$ with $\psi_i = \alpha_1^i$.
    
"""

# ╔═╡ a262674e-d56a-4252-8929-3d548150955c
md"""
### MA(q) and AR(p)
***

- MA(q): $x_t = w_t + \beta_1 w_{t-1} + \beta_2 w_{t-2} + ... + \beta_q w_{t-q}$
- AR(p): $x_t = \alpha_1 x_{t-1} + \alpha_2 x_{t-2} + ... + \alpha_p x_{t-p} + w_t$



- We can rewrite the above equations by the backshift operator:

- $\theta(B) = 1 + \beta_1 B + \beta_2 B^2 + ... + \beta_q B^q$
- $\phi(B) = 1 - \alpha_1 B - \alpha_2 B^2 - ... - \alpha_p B^p$



- So that:

```math
{\rm MA}(q): x_t = \theta(B) w_t \qquad {\rm AR}(p): \phi(B) x_t = w_t
```

- A very compact form making easier further manipulations.
"""

# ╔═╡ c0ca31e3-8026-46b6-aa90-ffe53b7c2305
md"""
### ARMA(p,q) processes
***

- A process, $x_t$, is said to be an ARMA(p,q) process if it has the form: 

```math
\phi(B) x_t = \theta(B) w_t
```

- It is possible to prove that:
    - ARMA(p,q) is stationary if and only if the roots of $\phi(B)w$ lie outside the unit circle (i.e. $\phi(B)w \ne 0$ for all $|w| \le 1$). 
    - ARMA(p,q) is invertible if and only if the roots of $\theta(B)w$ lie outside the unit circle (i.e. $\theta(B)w \ne 0$ for all $|w| \le 1$).

\

- It is also possible to convert an ARMA process, with $\phi(B)$ and $\theta(B)$ finite ordine polynomials, to infinite order MA or AR processes:

- MA: $x_t = \psi(B)w_t$, where $\psi(B) = \theta(B)/\phi(B)$.
- AR: $\pi(B)x_t = w_t$, where $\pi(B) = \phi(B)/\theta(B)$. 
"""

# ╔═╡ 44ab551f-b0cc-45cb-977e-3a4b77f15c9d
md"""
#### Exercize: the ACF of MA(p) and AR(q) processes
***

- MA(1), $\beta_1 = 0.5$.
"""

# ╔═╡ 7ad279de-35ed-4d38-bb78-9a097f5847d9
begin
	N16 = 10000
	beta116 = 0.5
	
	
	d16 = Normal()
	sigma16 = rand(d16, N16)
	
	x16 = zeros(N16)
	for i in range(2,N16)
	    x16[i] = sigma16[i]+beta116*sigma16[i-1]
	end
	
	rs16 = GetACF(x16,40)
	fg16 = Figure()
	axfg16 = Axis(fg16[1, 1],)
	stem!(rs16[1])
	hlines!([rs16[2],rs16[3]],linestyle=:dash)
	fg16
end

# ╔═╡ a7b6b298-49de-4a5f-a095-8306d96661a8
md"""
- MA(2), $\beta_1 = 1/6, \beta_2 = 1/2$.
"""

# ╔═╡ 973b9499-23be-4f59-a4a3-69bb3e1b2871
begin
	N17 = 10000
	beta117 = 1/6
	beta217 = 0.5
	
	
	d17 = Normal()
	sigma17 = rand(d17, N17)
	
	x17 = zeros(N17)
	for i in range(3,N17)
	    x17[i] = sigma17[i]+beta117*sigma17[i-1]+beta217*sigma17[i-2]
	end
	
	rs17 = GetACF(x17,40)
	fg17 = Figure()
	axfg17 = Axis(fg17[1, 1],)
	stem!(rs17[1])
	hlines!([rs17[2],rs17[3]],linestyle=:dash)
	fg17
end

# ╔═╡ d5caf4c8-59d3-488d-9fa2-f2e75eaa5b85
md"""
- MA(5), $\beta_{12} = -1/2, \beta_{345} = 1/4$.
"""

# ╔═╡ cc23ada7-044a-4620-b312-7e3ccf2eb25f
begin
	N18 = 10000
	beta1218 = -0.5
	beta34518 = 0.25
	
	
	d18 = Normal()
	sigma18 = rand(d18, N18)
	
	x18 = zeros(N18)
	for i in range(6,N18)
	    x18[i] = sigma18[i]+beta1218*sigma18[i-1]+beta1218*sigma18[i-2]+beta34518*sigma18[i-3]+beta34518*sigma18[i-4]+beta34518*sigma18[i-5]
	end
	
	rs18 = GetACF(x18,40)
	fg18 = Figure()
	axfg18 = Axis(fg18[1, 1],)
	stem!(rs18[1])
	hlines!([rs18[2],rs18[3]],linestyle=:dash)
	fg18
end

# ╔═╡ b679f761-1695-4710-afef-17b63b204e18
md"""
- MA(10), $\beta_j = 1/2$ for $j=1...10$.
"""

# ╔═╡ 86bbfeb7-ccb8-44e9-b9e5-0e7609a1d6db
begin
	N19 = 10000
	beta19 = 0.5
	
	
	d19 = Normal()
	sigma19 = rand(d19, N19)
	
	x19 = []
	
	x19 = zeros(N19)
	for i in range(11,N19)
	    xt = sigma19[i]
	    for j in range(1,10)
	        xt = xt + beta19*sigma19[i-j]
	    end
	    push!(x19,xt)
	end
	
	rs19 = GetACF(x19,40)
	fg19 = Figure()
	axfg19 = Axis(fg19[1, 1],)
	stem!(rs19[1])
	hlines!([rs19[2],rs19[3]],linestyle=:dash)
	fg19
end

# ╔═╡ 2f8bb6a9-e09d-4ff3-a158-dc00b26c0fa9
md"""
- AR(1): $\alpha_1 = 1/2$.
"""

# ╔═╡ 7bed93e6-0ee2-4374-8bc2-e4d988ad3c8c
begin
	N20 = 10000
	alpha20 = 0.5
	
	
	d20 = Normal()
	sigma20 = rand(d20, N20)
	
	x20 = [0.,]
	
	for i in range(2,N20)
	    push!(x20,alpha20*x20[i-1]+sigma20[i])
	end
	
	rs20 = GetACF(x20,40)
	fg20 = Figure()
	axfg20 = Axis(fg20[1, 1],)
	stem!(rs20[1])
	hlines!([rs20[2],rs20[3]],linestyle=:dash)
	fg20
end

# ╔═╡ 4457932a-49ab-449e-bd93-240e3458272b
md"""
- AR(2), $\alpha = 1/6, \alpha_2 = 1/2$.
"""

# ╔═╡ 99594bae-47cf-408b-aae7-ead6e783cbd3
begin
	N21 = 10000
	alpha121 = 1/6
	alpha221 = 0.5
	
	
	d21 = Normal()
	sigma21 = rand(d21, N21)
	
	x21 = [0.,0.]
	
	for i in range(3,N21)
	    push!(x21,alpha121*x21[i-1]+alpha221*x21[i-2]+sigma21[i])
	end
	
	rs21 = GetACF(x21,40)
	fg21 = Figure()
	axfg21 = Axis(fg21[1, 1],)
	stem!(rs21[1])
	hlines!([rs21[2],rs21[3]],linestyle=:dash)
	fg21
end

# ╔═╡ 9cd4d3a6-2c2c-467c-931f-cba17568046a
md"""
- AR(5), $\alpha_{12} = -1/2, \alpha_{345} = 1/4$.
"""

# ╔═╡ abb69c9b-2ccf-408e-8bed-75ded62f6c34
begin
	N22 = 10000
	alpha1222 = -0.5
	alpha34522 = 0.25
	
	
	d22 = Normal()
	sigma22 = rand(d22, N22)
	
	x22 = [0.,0.,0.,0.,0.]
	
	for i in range(6,N22)
	    push!(x22,alpha1222*x22[i-1]+alpha1222*x22[i-2]+alpha34522*x22[i-3]+alpha34522*x22[i-4]+alpha34522*x22[i-5]+sigma22[i])
	end
	
	rs22 = GetACF(x22,40)
	fg22 = Figure()
	axfg22 = Axis(fg22[1, 1],)
	stem!(rs22[1])
	hlines!([rs22[2],rs22[3]],linestyle=:dash)
	fg22
end

# ╔═╡ a7174765-5cb1-494c-9974-10010a900de3
md"""
- AR(8), $\alpha_j = 1/9$ for $j=1...8$.
"""

# ╔═╡ 07fe7c84-caa4-4d8e-a4db-f6e874ad7568
begin
	N23 = 10000
	alpha23 = 1/9
	
	
	d23 = Normal()
	sigma23 = rand(d23, N23)
	
	x23 = [0.,0.,0.,0.,0.,0.,0.,0.]
	
	for i in range(9,N23)
	    xt = sigma23[i]
	    for j in range(1,8)
	        xt = xt + alpha23*x23[i-j]
	    end
	    push!(x23,xt)
	end
	
	rs23 = GetACF(x23,40)
	fg23 = Figure()
	axfg23 = Axis(fg23[1, 1],)
	stem!(rs23[1])
	hlines!([rs23[2],rs23[3]],linestyle=:dash)
	fg23
end

# ╔═╡ d7d97154-4a7e-4e8f-a7e2-4eea6f9d1595
md"""
- ARMA(1,1): $\alpha_1 = 1/2, \beta_1 = 1/2$.
"""

# ╔═╡ ec56becc-5ba6-4a5c-9a30-cdfbaf823674
begin
	N24 = 10000
	alpha124 = 0.5
	beta124 = 0.5
	
	
	d24 = Normal()
	sigma24 = rand(d24, N24)
	
	x24 = [0.,]
	
	for i in range(2,N24)
	    xt = sigma24[i] + alpha124*x24[i-1] + beta124*sigma24[i-1]
	    push!(x24,xt)
	end
	
	rs24 = GetACF(x24,40)
	fg24 = Figure()
	axfg24 = Axis(fg24[1, 1],)
	stem!(rs24[1])
	hlines!([rs24[2],rs24[3]],linestyle=:dash)
	fg24
end

# ╔═╡ bf3e5cbf-6099-4bdd-966f-bcd59c3f09cf
md"""
- From now on we write our ARMA processes by means of a proper library.

- ARMA(2,1): $\alpha_1 = 1/6, \alpha_2=1/2, \beta_1 = 1/2$.
"""

# ╔═╡ d908b9c8-d4bd-47cd-b2f1-f662f734cd89
begin
	arma_rvs25 = arma(10000, 1., SVector(1/6,0.5),SVector(0.5))
	
	rs25 = GetACF(arma_rvs25,20)
	fg25 = Figure()
	axfg25 = Axis(fg25[1, 1],)
	stem!(rs25[1])
	hlines!([rs25[2],rs25[3]],linestyle=:dash)
	fg25
end

# ╔═╡ 5a27d099-c133-4845-b288-f88304287f73
md"""
- ARMA(2,2): $\alpha_1 = -1/2, \alpha_2=1/4, \beta_1 = -1/2, \beta_2 = 1/4$.
"""

# ╔═╡ cc15941d-5cf3-4c15-9f87-dcd40d0a2eac
begin
	arma_rvs26 = arma(10000, 1., SVector(-0.5,0.25),SVector(-0.5,0.25))
	
	rs26 = GetACF(arma_rvs26,20)
	fg26 = Figure()
	axfg26 = Axis(fg26[1, 1],)
	stem!(rs26[1])
	hlines!([rs26[2],rs26[3]],linestyle=:dash)
	fg26
end

# ╔═╡ bf84f03a-3c9f-4c6f-a886-8b40c24124b8
md"""
- ARMA(2,2): $\alpha_1 = 1/9, \alpha_2=1/9, \beta_1 = 1/2, \beta_2 = -1/4$.
"""

# ╔═╡ 796008cd-13ed-481c-ab43-7f9998e21ebd
begin
	arma_rvs27 = arma(10000, 1., SVector(1/9,1/9),SVector(0.5,-0.25))
	
	rs27 = GetACF(arma_rvs27,40)
	fg27 = Figure()
	axfg27 = Axis(fg27[1, 1],)
	stem!(rs27[1])
	hlines!([rs27[2],rs27[3]],linestyle=:dash)
	fg27
end

# ╔═╡ d94bbb24-b0fd-44d4-bedc-586941fa233c
# ╠═╡ show_logs = false
cm"""
- The ACF is a powerful method to derive the order of the process but as soon as the processes become complex cannot give unambiguos answers.

- A combined use of ACF and PACF can solve more cases.

$(LocalResource("Pics/acfpacf.png"))

- For instance, let's compare AR(1) with ``\alpha_1 = 0.6`` to ARMA(1,1) with ``\alpha_1 = 1/2, \beta_1 = 1/2``. 
"""

# ╔═╡ c1e80404-665c-4a68-978c-1df72a158059
begin
	
	arma_t = arma(10000, 1., SVector(0.5),SVector(0.5))
	ar_t = arma(10000, 1., SVector(0.6),nothing)
	
	armaacf = GetACF(arma_t,20)
	armapacf = GetPACF(arma_t,20)
	aracf = GetACF(ar_t,20)
	arpacf = GetPACF(ar_t,20)
	
	fg28 = Figure()
	ax1fg28 = Axis(fg28[1, 1],
	    title = "AR ACF")
	stem!(aracf[1])
	hlines!([aracf[2],aracf[3]],linestyle=:dash)
	
	ax2fg28 = Axis(fg28[1, 2],
	    title = "ARMA ACF")
	stem!(armaacf[1])
	hlines!([armaacf[2],armaacf[3]],linestyle=:dash)
	
	ax3fg28 = Axis(fg28[2, 1],
	    title = "AR PACF")
	stem!(arpacf[1])
	hlines!([arpacf[2],arpacf[3]],linestyle=:dash)
	
	ax4fg28 = Axis(fg28[2, 2],
	    title = "ARMA PACF")
	stem!(armapacf[1])
	hlines!([armapacf[2],armapacf[3]],linestyle=:dash)
	
	
	fg28
	
end

# ╔═╡ 73aecd0e-d6d0-4429-b63b-38f572115c50
md"""
## ARIMA(p,d,q)
***

- Autoregressive Integrated Moving Average

- A process $x_t$ is ARIMA(p,d,q) if $x_t$, differenced “d times" ($\nabla^d x_t$), is an ARMA(p,q) process.

```math
\phi(B) \nabla^d x_t = \theta(B) w_t \qquad \phi(B) (1-B)^d x_t = \theta(B) w_t
```

- Having chosen the right orders of the moving average autoregressive processes there are different strategies for determining the ARIMA parameters (Yule-Walker equations, maximum-likelihood, etc.). 

- Typically, one could follow a procedure known as [Box-Jenkins method](https://en.wikipedia.org/wiki/Box%E2%80%93Jenkins_method) that, with some simplification, essentially is:
    1. Plot the data
    2. Difference until series is stationary (find d)
    3. Examine ACF and PACF to guess p and q
    4. Fit ARIMA(p,d,q) to the original data
    5. Check model diagnostics 

- Or, even if time-consuming, a grid-search (or something more elaborated) can be a solution.
"""

# ╔═╡ 61cade13-3747-4481-aa8a-b59ff87ca05d
md"""
#### Exercise: fit a dataset by an ARIMA model
***

- Let's go back to the airline passenger dataset.

- Converting to log and differencing we got a stationary time-series.

- Let's study the ACF and PACF of the difference time-series.
"""

# ╔═╡ ca1ac1dd-b580-491b-94b6-8ece8c037298
begin
	tracf = GetACF(dropmissing(train)[!,"#Passengers_log_diff"],30)
	trpacf = GetPACF(dropmissing(train)[!,"#Passengers_log_diff"],30)
	
	fg29 = Figure()
	ax1fg29 = Axis(fg29[1, 1],
	    title = "ACF")
	stem!(tracf[1])
	hlines!([tracf[2],tracf[3]],linestyle=:dash)
	
	ax2fg29 = Axis(fg29[2, 1],
	    title = "PACF")
	stem!(trpacf[1])
	hlines!([trpacf[2],trpacf[3]],linestyle=:dash)
	
	fg29
end

# ╔═╡ f1e1cdde-a8b8-4ff5-9c49-6207de8ad125
md"""
- The ACF becomes 0 after 2 lags, and the seasonality is clearly visible, while the PACF is 0 after 8 or 9 lags.

- We might guess that the AR order could be 8 and the MA order 2.

- Diagnostic plots of the time series can be used along with heuristic rules to determine the hyperparameters of the ARIMA model.cThese are good in most, but perhaps not all, situations.

- As a matter if fact, we also try with a brute-force grid search:
"""

# ╔═╡ 6bf9db61-f56c-4e91-af6b-772a5066e709
begin
	bics = Dict()
	minbc = 1e6
	for p in 0:10
	    for q in 0:10
	        model_ARIMA = StateSpaceModels.SARIMA(collect(skipmissing(train[!,"#Passengers_log_diff"])); order = (p, 0, q), suppress_warns=true)
	        try
	            StateSpaceModels.fit!(model_ARIMA, save_hyperparameter_distribution=false, optimizer = Optimizer(StateSpaceModels.Optim.NelderMead()))
	            println("p: ", p, " q: ", q, " BIC: ", model_ARIMA.results.bic)
	            if model_ARIMA.results.bic < minbc
	                minbc = model_ARIMA.results.bic
	                bics["BIC"] = minbc
	                bics["p"] = p
	                bics["q"] = q
	            end
	        catch DomainError
	            print()
	        end
	    end
	end
	println(bics)
end

# ╔═╡ 435071db-e81b-4b3e-ab0e-a271030749e4
md"""
- All in all, the grid search did not provide results so different wrt to those suggested by the ACF/PACF plots.
"""

# ╔═╡ dc5a5f74-9289-4614-b1ae-553d43c33120
begin
	model_ARIMA = StateSpaceModels.SARIMA(collect(skipmissing(train[!,"#Passengers_log_diff"])); order = (bics["p"], 0, bics["q"]), suppress_warns=true)
	StateSpaceModels.fit!(model_ARIMA, save_hyperparameter_distribution=false, optimizer = Optimizer(StateSpaceModels.Optim.NelderMead()))
end

# ╔═╡ 61d51f9c-2368-478c-9d77-90c7d9a72610
begin
	kf = kalman_filter(model_ARIMA)
	
	plotdiagnostics(kf)
end

# ╔═╡ 15cf410a-13ee-4638-b91a-4dc286754c8c
md"""
- The fit is not perfect, in particular for the tail values, yet for this exercise is fine.
"""

# ╔═╡ e74e9fba-8062-45b0-aea9-b7ae1a646bf3
begin
	fg30 = Figure()
	
	ax1fg30 = Axis(fg30[1, 1],
	    title="ARMA fit",
	    )
	
	lines!(model_ARIMA.system.y,color=:blue,label="data")
	lines!(model_ARIMA.system.y .+ StateSpaceModels.get_innovations(kf)[:,1],color=:red,alpha=0.5)
	
	axislegend()
	
	fg30
end

# ╔═╡ 8d084ab0-3027-49d8-9b3b-24f06d0af830
md"""
- Indeed the modeling of the input dataset is rather satisfactory.

- Finally, we can plot our data with the original scale:
"""

# ╔═╡ e688750a-c6a9-4e08-ba7c-c1fd314405ec
begin
	origdata = collect(skipmissing(ShiftedArrays.lag(train[!,"#Passengers_log"],6))) .+ collect(skipmissing(train[!,"#Passengers_log_diff"]))
	arimamodel = collect(skipmissing(ShiftedArrays.lag(train[!,"#Passengers_log"],6))) .+ StateSpaceModels.get_innovations(kf)[:,1] .+ model_ARIMA.system.y
	#origdata = train['#Passengers_log'].shift(6) + train['#Passengers_log_diff']
	#arimamodel = train['#Passengers_log'].shift(6) + results_ARIMA.predict(0,typ='levels')
	
	fg31 = Figure()
	
	ax1fg31 = Axis(fg31[1, 1],
	    title="ARMA fit",
	    )
	
	lines!(exp.(origdata),color=:blue,label="data")
	lines!(exp.(arimamodel),color=:red,alpha=0.5,label="fit")
	
	axislegend(position = :lt)
	
	fg31
end

# ╔═╡ c9b352f2-016d-4d16-8053-48eedc39458d
md"""
- Nonstationarity is a complex issue, anyway. 

- Another form of nonstationarity occurs when the variance, rather than the local mean level, of  the time series changes during the observation. This is commonly called *volatility*, in econometric contexts.

- These kind of problems can be addressed by the autoregressive conditional heteroscedastic (ARCH) models. Here the variance is assumed to be a stochastic autoregressive process depending on previous values.

- For instance, in a ARCH(1) model: 

```math
{\rm Var}_{\rm ARCH(1)}(\epsilon_i) = w_t + \alpha_1 {\rm Var(\epsilon_{i-1})}
```

- Nevertheless, the zoo of general linear processes is rich...!
    - ARFIMA: autoregressive moving average fractionally integrated.
    - SARIMA: multiplicative seasonal autoregressive integrated moving average model.
    - CARFIMA: designed for irregularly sampled time series.
    - ARMAX: where X stands for “exogenous covariate”.

```math
x_t = \sum_{i=1}^p \alpha_i x_{t-i} + \sum_{j=1}^q \beta_j w_{t-j} + \sum_{k=1}^r c_k y_{t-k}
```

- where the exogenous covariates can be any function linear in the parameters: periodic, or some deterministic trend (e.g., polynomial), or even a tabulated variables without a simple mathematical form. 

"""

# ╔═╡ e3234616-d59f-43bb-91b3-bd2f727c68b1
md"""
### State space models
***

- SSM form a broader approach to parametric modeling of complicated time series, where autoregressive and other stochastic behaviors can be combined with linear or nonlinear deterministic trends or periodic components.

- Functional behaviors can be linear or nonlinear, and can involve derivatives of state variables. Noise can be Gaussian or non-Gaussian, white or autoregressive.

- State space models are defined hierarchically with one level describing relationships between “state variables” defining the temporal behavior of the system and other levels describing how the underlying system relates to  the observables. 

- A state-space model for a (possibly multivariate) time series {$Y_t, t = 1,2,...$} consists of two equations. The first, known as the observation equation, expresses the w-dimensional observation $Y_t$ as a linear function of a v-dimensional state variable $X_t$ plus noise. 

```math
Y_t = G_t X_t + W_t, \qquad t=1,2,...
```

- The second equation is called the state equation, and determines the state $X_{t+1}$ at time $t+1$ in terms of the previous state $X_t$ and a noise term:

```math
X_{t+1} = F_t X_t + V_t, \qquad t=1,2,...
```

- The noise terms are often (but not necessarily) supposed to be uncorrelated.

- It is possible to find a state-space representation for a large number of time-series.

- Neither ${X_t}$ nor ${Y_t}$ is necessarily stationary.

- The beauty of a state-space representation lies in the simple structure of the state equation, which permits relatively simple analysis of the process ${X_t}$.

- The behavior of ${Y_t}$ is then easy to determine from that of ${X_t}$ using the observation equation.

- The coefficients of the model are calculated by least squares or maximum likelihood estimation, often through a recursive algorithm such as, e.g., the *Kalman filter*. 
"""

# ╔═╡ 4a98a96f-a1b3-4144-b98a-9f65c0ca44aa
md"""
## Reference & Material

Material and papers related to the topics discussed in this lecture.

- [Ivezić et al. (2020) - "Statistics, Data Mining, and Machine Learning in Astronomy"](https://ui.adsabs.harvard.edu/abs/2020sdmm.book.....I/abstract)
- [Feigelson et al. (2016) - "Autoregressive Times Series Methods for Time Domain Astronomy”](https://ui.adsabs.harvard.edu/abs/2018FrP.....6...80F/abstract)
- [Feigelson & Babu (2013) - "Statistical Methods for Astronomy"](https://ui.adsabs.harvard.edu/abs/2013pss2.book..445F/abstract)
"""

# ╔═╡ e97eaccb-f03c-455e-8f15-b606f66704d9
md"""
## Further Material

Papers for examining more closely some of the discussed topics.

- [Alexander (1997) - "Is AGN Variability Correlated with Other AGN Properties? ZDCF Analysis of Small Samples of Sparse Light Curves"](https://ui.adsabs.harvard.edu/abs/1997ASSL..218..163A/abstract)
- [Kelly et al. (2014) - "Flexible and Scalable Methods for Quantifying Stochastic Variability in the Era of Massive Time-domain Astronomical Data Sets"](https://ui.adsabs.harvard.edu/abs/2014ApJ...788...33K/abstract)
- [Stone et al. (2022) - "Optical variability of quasars with 20-yr photometric light curves](https://ui.adsabs.harvard.edu/abs/2022MNRAS.514..164S/abstract)
- [Tarnopolski et al. (2020) - "A Comprehensive Power Spectral Density Analysis of Astronomical Time Series. I. The Fermi-LAT Gamma-Ray Light Curves of Selected Blazars"](https://ui.adsabs.harvard.edu/abs/2020ApJS..250....1T/abstract).
"""

# ╔═╡ ecfe0d5e-edc1-4557-8534-93c3970f366c
md"""
### Credits
***

This notebook contains material obtained from [https://www.analyticsvidhya.com/blog/2018/09/non-stationary-time-series-python/](https://www.analyticsvidhya.com/blog/2018/09/non-stationary-time-series-python/), [https://towardsdatascience.com/how-to-analyse-a-single-time-series-variable-11dcca7bf16c](https://towardsdatascience.com/how-to-analyse-a-single-time-series-variable-11dcca7bf16c), [https://machinelearningmastery.com/grid-search-arima-hyperparameters-with-python/](https://machinelearningmastery.com/grid-search-arima-hyperparameters-with-python/), [https://machinelearningmastery.com/arima-for-time-series-forecasting-with-python/](https://machinelearningmastery.com/arima-for-time-series-forecasting-with-python/), and from [https://www.analyticsvidhya.com/blog/2016/02/time-series-forecasting-codes-python/](https://www.analyticsvidhya.com/blog/2016/02/time-series-forecasting-codes-python/).
"""

# ╔═╡ eca5662e-1eee-46a0-82ce-ca3f5f70ae77
cm"""
## Course Flow

<table>
  <tr>
    <td>Previous lecture</td>
    <td>Next lecture</td>
  </tr>
  <tr>
    <td><a href="./open?path=Lectures/Science Case - Variable Stars/Lecture-VariableStars.jl">Science case about variable stars</a></td>
    <td><a href="./open?path=Lectures/ScienceCase-AGNandBlazars/Lecture-AGN-and-Blazars.jl">Science case about AGN and blazars</a></td>
  </tr>
 </table>


"""

# ╔═╡ bd9d24c9-3dc1-4759-b2ea-2663c6a49678
md"""
**Copyright**

This notebook is provided as [Open Educational Resource](https://en.wikipedia.org/wiki/Open_educational_resources). Feel free to use the notebook for your own purposes. The text is licensed under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/), the code of the examples, unless obtained from other properly quoted sources, under the [MIT license](https://opensource.org/licenses/MIT). Please attribute the work as follows: *Stefano Covino, Time Domain Astrophysics - Lecture notes featuring computational examples, 2025*.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ARFIMA = "9d0fb3db-ba49-4108-bc86-650b3813b6d5"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
CommonMark = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
HypothesisTests = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
LombScargle = "fc60dff9-86e7-5f2f-a8a0-edeadbb75bd9"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ShiftedArrays = "1277b4bf-5013-50f5-be3d-901d8477a67a"
StateSpaceModels = "99342f36-827c-5390-97c9-d7f9ee765c78"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
ARFIMA = "~0.4.0"
CSV = "~0.10.16"
CairoMakie = "~0.15.9"
CommonMark = "~1.0.1"
DataFrames = "~1.8.1"
Distributions = "~0.25.123"
HypothesisTests = "~0.11.6"
LombScargle = "~1.0.3"
Plots = "~1.41.6"
PlutoUI = "~0.7.79"
ShiftedArrays = "~1.0.0"
StateSpaceModels = "~0.6.7"
StatsBase = "~0.33.21"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.5"
manifest_format = "2.0"
project_hash = "cabde73850f21f3b11e59cd5c718c403e0029ef0"

[[deps.ADTypes]]
git-tree-sha1 = "f7304359109c768cf32dc5fa2d371565bb63b68a"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.21.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.ARFIMA]]
deps = ["LinearAlgebra", "Random", "StaticArrays", "Test"]
git-tree-sha1 = "e75e73b854b4f592466782c0d035f0e5b64ac83d"
uuid = "9d0fb3db-ba49-4108-bc86-650b3813b6d5"
version = "0.4.0"

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

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "856ecd7cebb68e5fc87abecd2326ad59f0f911f3"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.43"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

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

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "78b3a7a536b4b0a747a0f296ea77091ca0a9f9a3"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.23.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceAMDGPUExt = "AMDGPU"
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

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

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

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

[[deps.Combinatorics]]
git-tree-sha1 = "c761b00e7755700f9cdf5b02039939d1359330e1"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.1.0"

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

[[deps.CommonSolve]]
git-tree-sha1 = "78ea4ddbcf9c241827e7035c3a03e2e456711470"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.6"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

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

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "3b4be73db165146d8a88e47924f464e55ab053cd"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.7"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

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
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "c55f5a9fd67bdbc8e089b5a3111fe4292986a8e8"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.6"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "7ae99144ea44715402c6c882bfef2adbeadbc4ce"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.16"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

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

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

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
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

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

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "9340ca07ca27093ff68418b7558ca37b05f8aeb1"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.29.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

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

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "eef4c86803f47dcb61e9b8790ecaa96956fdd8ae"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.3.2"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

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

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "b7bfd56fa66616138dfe5237da4dc13bbd83c67f"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "44716a1a667cb867ee0e9ec8edc31c3e4aa5afdc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.24"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "be8a1b8065959e24fdc1b51402f39f3b6f0f6653"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.24+0"

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

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

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

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

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

[[deps.HypothesisTests]]
deps = ["Combinatorics", "Distributions", "LinearAlgebra", "Printf", "Random", "Roots", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "b67985bd11331ccef26109a6269dbaae01474a72"
uuid = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
version = "0.11.6"

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
weakdeps = ["ForwardDiff", "Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

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
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "ac0aaa807ed5eaf13f67afe188ebc07e828ff640"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.10.0"

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

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

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

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

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

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Printf"]
git-tree-sha1 = "9ea3422d03222c6de679934d1c08f0a99405aa03"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.5.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LinearMaps]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7f6be2e4cdaaf558623d93113d6ddade7b916209"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.11.4"
weakdeps = ["ChainRulesCore", "SparseArrays", "Statistics"]

    [deps.LinearMaps.extensions]
    LinearMapsChainRulesCoreExt = "ChainRulesCore"
    LinearMapsSparseArraysExt = "SparseArrays"
    LinearMapsStatisticsExt = "Statistics"

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

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

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

[[deps.MatrixEquations]]
deps = ["LinearAlgebra", "LinearMaps"]
git-tree-sha1 = "f765b4eda3ea9be8e644b9127809ca5151f3d9ea"
uuid = "99c1a7ee-ab34-5fd5-8076-27c950a045f4"
version = "2.4.2"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ff69a2b1330bcb730b9ac1ab7dd680176f5896b8"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.1010+0"

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

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "2c140d60d7cb82badf06d8783800d0bcd1a7daa2"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.8.1"

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

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "22df8573f8e7c593ac205455ca088989d0a2c7a0"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.7"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "25a6638571a902ecfb1ae2a18fc1575f86b1d4df"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.10.0"

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

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["Compat", "EnumX", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "48968edaf014f67e58fe4c8a4ce72d392aed3294"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.13.3"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

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
git-tree-sha1 = "f07c06228a1c670ae4c87d1276b92c7c597fdda0"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.35"

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

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "cb20a4eacda080e517e4deb9cfb6c7c518131265"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.6"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3ac7038a98ef6977d44adeadc73cc6f596c08109"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.79"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["Intervals", "LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "a1f7f4e41404bed760213ca01d7f384319f717a5"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "2.0.25"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

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

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "d7a4bff94f42208ce3cf6bc8e4e7d1d663e7ee8b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "d5b7dd0e226774cbd87e2790e34def09245c7eab"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

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

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

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

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "10a488dbecb88a9679c8f357d383d7d83dcc748d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.13"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

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

[[deps.SeasonalTrendLoess]]
deps = ["Statistics"]
git-tree-sha1 = "839dcd8152dc20663349781f7a7e8cf3d3009673"
uuid = "42fb36cb-998a-4034-bf40-4eee476c43a1"
version = "0.1.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ebe7e59b37c400f694f52b58c93d26201387da70"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.9"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.ShiftedArrays]]
git-tree-sha1 = "22395afdcf37d6709a5a0766cc4a5ca52cb85ea0"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "1.0.0"

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

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

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

[[deps.StateSpaceModels]]
deps = ["Distributions", "LinearAlgebra", "MatrixEquations", "Optim", "OrderedCollections", "Polynomials", "Printf", "RecipesBase", "SeasonalTrendLoess", "ShiftedArrays", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "1fca6ae24e606629e0701a347831a675a291e69b"
uuid = "99342f36-827c-5390-97c9-d7f9ee765c78"
version = "0.6.7"

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
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

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

[[deps.TZJData]]
deps = ["Artifacts"]
git-tree-sha1 = "7def47e953a91cdcebd08fbe76d69d2715499a7d"
uuid = "dc5dba14-91b3-4cab-a142-028a31da12f7"
version = "1.4.0+2025a"

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

[[deps.TimeZones]]
deps = ["Artifacts", "Dates", "Downloads", "InlineStrings", "Mocking", "Printf", "Scratch", "TZJData", "Unicode", "p7zip_jll"]
git-tree-sha1 = "38bb1023fb94bfbaf2a29e1e0de4bbba6fe0bf6d"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.21.2"
weakdeps = ["RecipesBase"]

    [deps.TimeZones.extensions]
    TimeZonesRecipesBaseExt = "RecipesBase"

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
weakdeps = ["ConstructionBase", "ForwardDiff", "InverseFunctions", "LaTeXStrings", "Latexify", "NaNMath", "Printf"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

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

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

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

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

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

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

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

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

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

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

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

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

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

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

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

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╟─14a0c11f-827d-46c8-bd5f-469a029f9afe
# ╟─2aeb4208-bc0c-4907-ab6d-da3f48210aa5
# ╠═508fcfab-51c1-4dc5-992d-d5961517b9bf
# ╟─fce942fc-7126-4e92-b758-30d36609117f
# ╟─97bebdb5-df0d-4eb2-b214-1a38dd9ffa85
# ╟─36377f43-7624-4113-9ea9-409fbd3ef079
# ╟─91f7ef29-6570-44d1-a5e2-79658576d488
# ╟─2d403a10-ede5-4d72-a732-59d4715273bc
# ╟─90b6cf29-44c0-425d-a814-910fb08009d2
# ╟─91c6247b-af12-47cf-811a-1440b6919576
# ╟─f1c10261-a2c7-473c-879b-005714127aee
# ╟─b40d3fb9-f423-44d5-8a3e-32e15a6cc564
# ╟─40e140c2-a8fb-4e31-9636-a36104f60460
# ╟─2f7623f2-359c-4c93-b6b2-ab48f606ec56
# ╟─f3ed1ba4-849b-4713-8e48-4a537ee894d7
# ╟─f0e73dd6-2f4f-4dd6-9426-f0a739956602
# ╟─dcb69463-c144-4726-8336-5f3851873d8a
# ╟─f8b03f2b-a241-4d77-8fc9-4c24a4c8ef90
# ╟─577fa116-fa14-4077-8fc8-4d3ea544f410
# ╟─20d1925b-ddb8-49ec-86a2-98f2efac3532
# ╟─12225d55-9c0a-482c-890b-8ca5f8139c5a
# ╟─db5d7b3e-66da-4498-a924-2f9b999fda56
# ╟─06a37524-c82e-4f68-8529-89759cdda2bd
# ╟─06f77be0-51c6-4ca2-82b1-eaa69d9481b7
# ╟─b1f35602-c627-41a7-8754-6634d5eff5a7
# ╟─c097e9b3-f9dc-473e-902a-cc5b03c15da3
# ╟─08a6a975-b496-4d19-8e8b-5df7aae43ad9
# ╟─2acb3d30-1b38-4255-b955-a0db9f77be72
# ╟─f34215b4-91fd-4175-96f1-524097ec7160
# ╟─45f962f9-79be-492a-8add-737fc73925c2
# ╟─5c1eff4e-8079-4568-a792-7528465db277
# ╠═9057a011-76b9-475a-906c-81afd18e99ce
# ╟─3da138ae-cd0b-4f7b-ab5a-319d4f65dee2
# ╠═25a8201f-402b-4868-9cd6-262a03b980fb
# ╟─698e9cc6-9de5-4ca6-b425-5bb038b08b2b
# ╠═f0996cba-03af-4247-99f7-d148c6e737a0
# ╟─0f220d9a-ffe5-4d84-ba4b-a65c5dbf4404
# ╟─59175381-65a4-45e7-ae85-a00467c4d181
# ╟─9a8ac940-8a2c-4d5a-809f-bbdb1dba12dc
# ╠═76917ce4-5249-4cbb-bdb7-da654ce2c75d
# ╠═a7c1c07b-23d7-4721-ac7d-70012483b056
# ╟─cc765f72-cb9d-4b1b-9e11-f2e7b6a4c71e
# ╠═7ecf612d-04ea-4e3b-b176-98476b443688
# ╠═7da9155f-27d8-48c0-93a6-13a3297fd302
# ╟─72700334-0d07-4486-ba3f-46fe0cfacb4f
# ╠═6021c695-2cd9-499c-bbbf-7d86b5a6b347
# ╠═3d013635-b4ac-49ec-85f6-93c26cf7539c
# ╟─c81ad1a5-b4ed-4056-a161-eab4f9e981d9
# ╠═f1bc344e-ac85-4815-928f-61f5f48c18a9
# ╟─381bc481-5c7e-4d9a-80c9-8f245718718d
# ╟─ba0761ce-3fd9-4b28-a9b4-46f850dfbafe
# ╟─9d1af19c-e1b3-48aa-9b9b-245347efd682
# ╠═e410b2af-6665-4ddb-bb89-4b38bec8bfc2
# ╠═2aca0524-9737-43e9-b456-9fa0b2f7939b
# ╟─0c8e81ce-7b48-427c-b72e-4fdad4015d51
# ╠═63a3b6f2-a46a-4a4c-ba9e-5639714cd93d
# ╟─d674149a-25a7-4c60-8770-e2b46cb2ed74
# ╠═ea14bb29-f6bf-4a0c-a6b1-d38e285207ee
# ╠═95047c44-9d05-4dc3-a50a-aa87cf9cd2c4
# ╟─61d7c290-6d4d-4da8-a4f6-f0c2c795a4f8
# ╠═962d9bb2-8fee-4d94-9054-a48f084127d7
# ╠═33df98ca-16f8-4bc7-8848-494d920e60d9
# ╟─d4d93ec8-e268-4fc9-b5ee-cdd194636659
# ╠═9512aa2e-afaa-4ddd-aac4-826abae447ff
# ╟─d0dc8d40-75d7-4611-aab6-d8542f82947e
# ╟─b0ae3114-1aca-4109-922a-02c25989e407
# ╠═924e216e-9ba7-4906-9d70-94f1813911e7
# ╠═30ece656-c5ba-42a6-8de2-e9827d72eabf
# ╟─b86d10a0-f3a0-4489-b595-64c3eafcfef2
# ╠═8a8d5e0e-7ccf-43a2-b18b-c389971252e9
# ╟─c6bab5e0-4a61-4a6c-aa54-34e0e4f2dcdc
# ╟─2905155b-0de7-4f15-acb0-aaaf1b421c9b
# ╟─fa1e8adf-969a-4476-81ae-b548cabcd8f6
# ╟─0edac0fb-f52e-4a2c-b3dc-54f8237ad325
# ╟─2eb45539-ef03-4f6d-86be-52a1e9228088
# ╟─db519b68-35a3-4396-a302-e0bd804985c1
# ╟─a262674e-d56a-4252-8929-3d548150955c
# ╟─c0ca31e3-8026-46b6-aa90-ffe53b7c2305
# ╟─44ab551f-b0cc-45cb-977e-3a4b77f15c9d
# ╠═7ad279de-35ed-4d38-bb78-9a097f5847d9
# ╟─a7b6b298-49de-4a5f-a095-8306d96661a8
# ╠═973b9499-23be-4f59-a4a3-69bb3e1b2871
# ╟─d5caf4c8-59d3-488d-9fa2-f2e75eaa5b85
# ╠═cc23ada7-044a-4620-b312-7e3ccf2eb25f
# ╟─b679f761-1695-4710-afef-17b63b204e18
# ╠═86bbfeb7-ccb8-44e9-b9e5-0e7609a1d6db
# ╟─2f8bb6a9-e09d-4ff3-a158-dc00b26c0fa9
# ╠═7bed93e6-0ee2-4374-8bc2-e4d988ad3c8c
# ╟─4457932a-49ab-449e-bd93-240e3458272b
# ╠═99594bae-47cf-408b-aae7-ead6e783cbd3
# ╟─9cd4d3a6-2c2c-467c-931f-cba17568046a
# ╠═abb69c9b-2ccf-408e-8bed-75ded62f6c34
# ╟─a7174765-5cb1-494c-9974-10010a900de3
# ╠═07fe7c84-caa4-4d8e-a4db-f6e874ad7568
# ╟─d7d97154-4a7e-4e8f-a7e2-4eea6f9d1595
# ╠═ec56becc-5ba6-4a5c-9a30-cdfbaf823674
# ╟─bf3e5cbf-6099-4bdd-966f-bcd59c3f09cf
# ╠═d908b9c8-d4bd-47cd-b2f1-f662f734cd89
# ╟─5a27d099-c133-4845-b288-f88304287f73
# ╠═cc15941d-5cf3-4c15-9f87-dcd40d0a2eac
# ╟─bf84f03a-3c9f-4c6f-a886-8b40c24124b8
# ╠═796008cd-13ed-481c-ab43-7f9998e21ebd
# ╟─d94bbb24-b0fd-44d4-bedc-586941fa233c
# ╠═c1e80404-665c-4a68-978c-1df72a158059
# ╟─73aecd0e-d6d0-4429-b63b-38f572115c50
# ╟─61cade13-3747-4481-aa8a-b59ff87ca05d
# ╠═ca1ac1dd-b580-491b-94b6-8ece8c037298
# ╟─f1e1cdde-a8b8-4ff5-9c49-6207de8ad125
# ╠═6bf9db61-f56c-4e91-af6b-772a5066e709
# ╟─435071db-e81b-4b3e-ab0e-a271030749e4
# ╠═dc5a5f74-9289-4614-b1ae-553d43c33120
# ╠═61d51f9c-2368-478c-9d77-90c7d9a72610
# ╟─15cf410a-13ee-4638-b91a-4dc286754c8c
# ╠═e74e9fba-8062-45b0-aea9-b7ae1a646bf3
# ╟─8d084ab0-3027-49d8-9b3b-24f06d0af830
# ╠═e688750a-c6a9-4e08-ba7c-c1fd314405ec
# ╟─c9b352f2-016d-4d16-8053-48eedc39458d
# ╟─e3234616-d59f-43bb-91b3-bd2f727c68b1
# ╟─4a98a96f-a1b3-4144-b98a-9f65c0ca44aa
# ╟─e97eaccb-f03c-455e-8f15-b606f66704d9
# ╟─ecfe0d5e-edc1-4557-8534-93c3970f366c
# ╟─eca5662e-1eee-46a0-82ce-ca3f5f70ae77
# ╟─bd9d24c9-3dc1-4759-b2ea-2663c6a49678
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
