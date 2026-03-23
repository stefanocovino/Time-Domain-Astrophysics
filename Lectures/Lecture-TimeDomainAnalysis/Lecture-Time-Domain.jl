### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 80a1ed5e-3cae-4cbc-9343-0326fd01ddc2
# ╠═╡ show_logs = false
import Pkg; Pkg.activate(".")

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

# ╔═╡ 3f22f35d-0b3b-4292-84ef-576490d13ca0
Pkg.instantiate()

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

$(LocalResource("Pics/tabtimevsspectral.png"))


### The Correlation Function
***

- One of the main statistical tools for the analysis of stochastic variability is the autocorrelation function. It represents a specialized case of the correlation function of two functions, ``f(t)`` and ``g(t)``, scaled by their standard deviations, and defined at time lag  as: 

```math
CF(\Delta t) = \frac{\lim_{T\to\infty} \int_{(T)} f(t) g(t+\Delta t)dt}{\sigma_f \sigma_g}
```

- where ``σ_f`` and ``σ_g`` are standard deviations of ``f(t)`` and ``g(t)``, respectively. 

- With this normalization, the correlation function is unity for ``Δt = 0``, without normalization by standard deviation, the above expression is equal to the covariance function.


### The Auto-Correlation Function
***

- It is assumed that both ``f`` and ``g`` are statistically weakly stationary functions (more on this topic later).

- With ``f(t)=g(t)=y(t)``, the autocorrelation of ``y(t)`` defined at time lag  is:

```math
ACF(\Delta t) = \frac{\lim_{T\to\infty} \int_{(T)} y(t) y(t+\Delta t)dt}{\sigma^2_y}
```

- he autocorrelation function and the PSD of function ``y(t)`` are Fourier pairs; this fact is known as the Wiener–Khinchin theorem and applies to stationary random processes.

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

- Simple quantities like the sample mean are valid estimates of the population mean for stationary processes, but its uncertainty is not the standard value when autocorrelation is present:

```math
\widehat{Var}(\bar{X}_n) = \frac{\sigma^2}{n} \left[1 + 2\sum_{k=1}^{n-1}(1-k/n){\rm ACF}(k) \right]
```

- Qualitatively, this can be understood as due to the decrease of independent measurements.

- This is indeed one very important result, although often practically ignored.


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


### Z-transformed DCF
***

- The drawback of the DFC is that its sample distribution of is known to be very skewed and far from normal. 

- If bins have an equal number of points (i.e. bins can have different length), The bin distribution becomes approximately binomial and, if a [Fisher z-transform](https://en.wikipedia.org/wiki/Fisher%27s_z-distribution) is applied:

```math
z(t) = \frac{1}{2} \ln \left( \frac{1+{\rm DCF}(\tau)}{1-{\rm DCF}(\tau)} \right)
```

- The ``z(τ)`` values are now normally distributed with known mean and variance (see [Alexander 1997](https://ui.adsabs.harvard.edu/abs/1997ASSL..218..163A/abstract)).


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
end

# ╔═╡ f3ed1ba4-849b-4713-8e48-4a537ee894d7
md"""
#### Exercize: cross-correlation of data from the Covid19 outbreak in Italy
***

- Italian data downloaded from: https://github.com/pcm-dpc/COVID-19 and updated till April 07, 2022.
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
- A $\sim$20 day delay appears, that is quite reasonable. 

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
- The LS periodogram shows a considerable power at long periods, with a large peak at about 400 days. This is likely the summer/winter cycle that did not reproduce exatcly.

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
    <td><a href="./open?path=Lectures/Science Case - AGN and Blazars/Lecture-AGN-and-Blazars.jl">Science case about AGN and blazars</a></td>
  </tr>
 </table>


"""

# ╔═╡ bd9d24c9-3dc1-4759-b2ea-2663c6a49678
md"""
**Copyright**

This notebook is provided as [Open Educational Resource](https://en.wikipedia.org/wiki/Open_educational_resources). Feel free to use the notebook for your own purposes. The text is licensed under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/), the code of the examples, unless obtained from other properly quoted sources, under the [MIT license](https://opensource.org/licenses/MIT). Please attribute the work as follows: *Stefano Covino, Time Domain Astrophysics - Lecture notes featuring computational examples, 2025*.
"""

# ╔═╡ Cell order:
# ╟─14a0c11f-827d-46c8-bd5f-469a029f9afe
# ╟─2aeb4208-bc0c-4907-ab6d-da3f48210aa5
# ╠═80a1ed5e-3cae-4cbc-9343-0326fd01ddc2
# ╠═3f22f35d-0b3b-4292-84ef-576490d13ca0
# ╠═508fcfab-51c1-4dc5-992d-d5961517b9bf
# ╟─fce942fc-7126-4e92-b758-30d36609117f
# ╟─97bebdb5-df0d-4eb2-b214-1a38dd9ffa85
# ╟─36377f43-7624-4113-9ea9-409fbd3ef079
# ╟─40e140c2-a8fb-4e31-9636-a36104f60460
# ╠═2f7623f2-359c-4c93-b6b2-ab48f606ec56
# ╟─f3ed1ba4-849b-4713-8e48-4a537ee894d7
# ╠═f0e73dd6-2f4f-4dd6-9426-f0a739956602
# ╟─dcb69463-c144-4726-8336-5f3851873d8a
# ╠═f8b03f2b-a241-4d77-8fc9-4c24a4c8ef90
# ╟─577fa116-fa14-4077-8fc8-4d3ea544f410
# ╠═20d1925b-ddb8-49ec-86a2-98f2efac3532
# ╟─12225d55-9c0a-482c-890b-8ca5f8139c5a
# ╠═db5d7b3e-66da-4498-a924-2f9b999fda56
# ╟─06a37524-c82e-4f68-8529-89759cdda2bd
# ╟─06f77be0-51c6-4ca2-82b1-eaa69d9481b7
# ╠═b1f35602-c627-41a7-8754-6634d5eff5a7
# ╟─c097e9b3-f9dc-473e-902a-cc5b03c15da3
# ╠═08a6a975-b496-4d19-8e8b-5df7aae43ad9
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
