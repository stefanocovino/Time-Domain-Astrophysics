### A Pluto.jl notebook ###
# v0.20.21

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

# ╔═╡ 6a59a8e6-39fc-44b3-9921-8976c677f4b1
begin
	using CairoMakie
	using CommonMark
	using CSV
	using DataFrames
	using Dates
	using Distributions
	using DSP
	using FFTW
	using Format
	using HTTP
	using LaTeXStrings
	using Latexify
	using PlutoUI
	using Statistics
end

# ╔═╡ d9329a8d-b07c-4207-93a6-1668da23e296
md"""
**What is this?**


*This jupyter notebook is part of a collection of notebooks on various topics discussed during the Time Domain Astrophysics course delivered by Stefano Covino at the [Università dell'Insubria](https://www.uninsubria.eu/) in Como (Italy). Please direct questions and suggestions to [stefano.covino@inaf.it](mailto:stefano.covino@inaf.it).*
"""

# ╔═╡ 3307644d-e875-4f15-a140-a41f7ca82a8f
md"""
**This is a `Pluto` notebook**
"""

# ╔═╡ 771d6af3-a5e2-4869-82c3-68540f71cb41
TableOfContents()

# ╔═╡ 95ee75d5-3112-42b6-83aa-29e639ac6eb0
# ╠═╡ show_logs = false
md"""
$(LocalResource("Pics/TimeDomainBanner.jpg"))
"""

# ╔═╡ b175b388-911f-4862-afdb-7449fec2bf9e
cm"""
# Time-Series
***

- A time-series is any sequene of observation such that the distribution of a given value depends on the previous values.
- Time is an exogeneous (outside the model) variable that is directional - measurements only depend on the past.
    - This is a statement of causality.

- However, the exogenous variable can be anything.

- Let's assume to have a set of data extracted from ``y(t) = A \sin(\omega t)`` with homoscedastic variance ``V = \sigma^2 + A^2/2``.
    - This is easy to prove if you compute the variance as ``\sum (y-\lt y \gt)^2 / N``. Since the average value is zero, this turns out to be ``V = \frac{A^2}{N} \sum \sin^2 (\omega t)`` giving the ``A^2/2`` term.

- We can compute the ``\chi^2`` for this toy model:
```math
\chi^2_{\rm dof} = \frac{1}{N} \sum (y/\sigma)^2 = \frac{1}{N \sigma^2} \sum (y)^2 = \frac{1}{\sigma^2} {VAR} = 1 + \frac{A^2}{2\sigma^2}
```

- With no variability (``A \sim 0``) the expectation value of the ``\chi^2_{\rm dof} \sim 1`` with standard deviation ``\sqrt{2/N}``, while it'll be larger in case of variability.
- Therefore, in order to have ``\chi^2_{\rm dof} > 1 + 3 \sqrt{2/N}`` we need ``A > \sigma \sqrt[4]{72/N}``, which shows that with ``N`` sufficiently large we can detect variability well below the uncetainty of the single points.
"""

# ╔═╡ 882fc480-b132-45a8-906a-db157059b92c
md"""
# Fourier Analysis
***

- The aim of Fourier analysis is to express any function as a sum of different sines and cosines, characterised by an angular frequency $ω$ or a corresponding time period $P = 2π/ω$.

- Such a decomposition is known as a Fourier series:

```math
f(t) = \frac{a(0)}{2}+\sum_{n=1}^\infty [ a(n)\sin n\omega_0 t + b(n)\cos n\omega_0 t]
```

- This expression can represent practically any periodic function with period $P_0$, with suitable adjustment of the coefficients $a(n)$ and $b(n)$, which are known as the Fourier coefficients.

- The conditions under which the Fourier decomposition is valid are that $f(t)$ has only a finite number of finite discontinuities and only a finite number of extreme values within a period.

    - These are known as Dirichlet conditions and the functions obeying them are called piece-wise regular.

- The Fourier coefficients $a(n)$ and $b(n)$ can be determined by performing the following integrations:

```math
a(n) = \frac{2}{P_0} \int_0^{P_0} f(t)\sin n\omega_0 t \, dt \;\; ; \;\;\; b(n) = \frac{2}{P_0} \int_0^{P_0} f(t)\cos n\omega_0 t \, dt, \;\; n=0,1,2...
```

- Alternatively, the Fourier decomposition may also be expressed in an equivalent exponential form:

```math
f(t) = \frac{1}{P_0}\sum_{n=-\infty}^{\infty} c(n) e^{i n\omega_0 t}; \;\;\mbox{\rm where}\;\; c(n)=\int_{0}^{P_0} f(t) e^{-in\omega_0 t} dt
```

- In the above expansions, the $n = 0$ term is often called the constant or the D.C. (Direct Current) component, the $n = 1$ term the fundamental and the terms with $n > 1$ the harmonics.
"""

# ╔═╡ a6681a0b-2cbb-4cbd-ad44-7f01953f875f
begin
	# Define the Square Wave and its Fourier Series approximation
	# Square wave with period 2π, oscillating between -1 and 1
	# Fourier Series: (4/π) * Σ [sin(n*x) / n] for odd n
	function fourier_square(x, n_terms)
	    val = 0.0
	    for i in 1:n_terms
	        n = 2i - 1  # Only odd harmonics
	        val += sin(n * x) / n
	    end
	    return (4 / π) * val
	end
	
	# Setup Data
	x_vals = range(-π, π, length=1000)
	# Define the number of components to visualize
	n_components = [1, 3, 10, 50]
	colors = [:blue, :orange, :green, :red]
	
	# Create the Visualization
	fig = Figure(size = (800, 600), font = "sans")
	ax = Axis(fig[1, 1], 
	    title = "Fourier Series Approximation of a Box-Car Function",
	    xlabel = "x", ylabel = "f(x)",
	    xticks = ([-π, 0, π], ["-π", "0", "π"]))
	
	# Plot the ideal box-car (square wave) for reference
	lines!(ax, x_vals, [sign(sin(x)) for x in x_vals], 
	    color = :black, linestyle = :dash, label = "Ideal Square Wave")
	
	# Loop through and plot different levels of approximation
	for (i, n) in enumerate(n_components)
	    y_vals = [fourier_square(x, n) for x in x_vals]
	    lines!(ax, x_vals, y_vals, color = (colors[i], 0.8), 
	        linewidth = 2, label = "n = $n harmonics")
	end
	
	axislegend(ax, position = :rt)
	
	# Display the plot
	fig
end

# ╔═╡ bc82a078-fea7-4a00-b38c-a849fb760594
md"""
- Example of a Fourier series decomposition. The original time domain function $f(t)$, shown in the dashed line, is a square wave , which equals 1 from $t=0$ to $\pi$ and -1 from $-\pi$ to $0$. Thin lines show the Fourier sum with different number of terms.
    - The gradual improvement in the approximation of the function with increasing number of terms in Fourier series is evident.
- At the point of discontinuity, all the reconstructions pass through the average of the left and right limits of the original function and the transition gets progressively sharper with larger number of terms.
"""

# ╔═╡ 3bdf3bea-7e3c-4f09-98f7-bdbc16ea2880
md"""
### Continuous Fourier Transform
***

- We define the Fourier transform (FT) of any function $f(t)$ as:

```math
F(\omega) = \int\limits_{-\infty}^{+\infty} f(t) e^{-i\omega t} dt
```

- This is a linear transformation and no information is lost. The representations of a function in time and frequency domains are equivalent.

- The original $f(t)$ can be recovered by applying the inverse Fourier transform:

```math
f(t) = {1\over 2\pi} \int\limits_{-\infty}^{+\infty} F(\omega) e^{i\omega t} d\omega
```

- The FT has a number of interesting properties. It is linear, not necessarily a real function, and its amplitude is invariant to time shift (but not its phase).

| $f(t)$ | F(ω) |
|:-------: |:--------:|
| Real    | H(-ω) = H$^*$(ω)   |
| Even    | Even   |
| Odd    | Odd   |
| Real and Even | Real and Even     |
| Real and Odd  | Imaginary and Odd |

- Unless the original function is even, an unlikely situation in the case of a time series, its Fourier transform is complex!
"""

# ╔═╡ 0dbb5f0e-292d-4d8a-bb83-7954e4a46f29
# ╠═╡ show_logs = false
md"""
$(LocalResource("Pics/FTexamples.png"))
"""

# ╔═╡ 875732fc-5644-4c06-8f8c-39dc30f1d3b4
md"""
- Scaling: $h(at) \longleftrightarrow H(f/a) / |a|$, "broad" $\longleftrightarrow$ "narrow"
- Shifting: $h(t-t_0) \longleftrightarrow H(f) * e^{2\pi i f t_0}$, "shift" $\longleftrightarrow$ "phase roll/gradient"
- Convolution: $h(t) \ast g(t) \longleftrightarrow H(f) G(f)$, "convolution" $\longleftrightarrow$ "multiplication"
"""

# ╔═╡ 56081488-4ade-436a-991e-d7138368edf1
md""" 
> How to compute the Fourier Transform in a few simple cases can be found [here](./open?path=Lectures/Lecture-SpectralAnalysis/Lecture-FT.jl).
"""

# ╔═╡ 6d88e888-8f08-47ab-8906-c7b4b86b2588
md"""
## Power density spectrum (PDS)
***

- The Power Density Spectrum (PSD) is defined as the Fourier Transform multiplied by its complex conjugate and therefore the square modulus of the Fourier Transform:

```math
P(\omega) = F(\omega)\cdot F^*(\omega) =  | F(\omega) |^2
```

- If the original function is real (usually the case for time series) the PSD is an even function and the values at negative frequencies are redundant.

- The FT is a linear function, the PDS is not. This means that while the FT of the sum of two signals is the sum of the FT of the signals, in the case of the PDS this is not true and there are cross-terms to be considered.

- E.g., if two signals are $f(t)$ and $g(t)$, the PSD of the sum of the two is:

```math
P[f(t)+g(t)] = |F[f(t)+g(t)]|^2 =  F[f(t)+g(t)]\cdot F^*[f(t)+g(t)] = 
```

```math
= P[f(t)] + P[g(t)] + 2 Re\{F[f(t)]\cdot F[g(t)] \}
```

- If the two signals are uncorrelated, the cross-term is zero and linearity applies.
"""

# ╔═╡ 99f8ea18-5807-482e-8516-e2bdaf1546e6
md"""
## Autocorrelation function (ACF)
***

- The autocorrelation (ACF) of a function $f(t)$ is defined as:

```math
A(t) = \int\limits_{-\infty}^{+\infty} f(\tau) f(t+\tau) d\tau \Longleftrightarrow F(f) F^*(f) \equiv |F(f)|^2
```

- The autocorrelation of a function is the Fourier transform of its PSD.

- It is also simple to derive, by the *Parseval’s theorem*, simply setting $t=0$, that:

```math
\int\limits_{-\infty}^{+\infty} |f(t)|^2 dt = {1\over 2\pi} \int\limits_{-\infty}^{+\infty} |f(\omega)|^2 d\omega
```

> For the interested readers, a proof that PDS and ACF are Fourier duals can be found [here](./open?path=Lectures/Lecture-SpectralAnalysis/Lecture-PDS-ACF.jl).

"""

# ╔═╡ e087bf1e-dd1a-41f3-8154-bfb39325d905
cm"""
## Discrete Fourier transform
***

- In the real world, we only have discrete measurements, the time series, commonly called, in astronomy, “light curves”. They consist of ``N`` measurements ``x_k`` taken (now) at equally-spaced times ``t_k`` from ``0`` to ``T``.

- In this case we can define the discrete Fourier transform (and its inverse) as:

```math
a_j = \sum\limits_{k=0}^{N-1} x_k e^{-2\pi ijk/N}   \quad   (j=-N/2,...,N/2-1)
```

```math
x_k = {1\over N} \sum\limits_{k=-N/2}^{N/2-1} a_j e^{2\pi ijk/N} \quad\quad (k=0,...,N-1)
```

- Since the data are equally spaced, the times are ``kT/N`` and the frequencies are ``j/T``.

- The time step is ``δt = T/N`` and the frequency step is ``δν = 1/T``.

- As the discrete time series has a time step ``δt`` and a duration ``T``, there are limitations to the frequencies that can be examined:
    - The lowest frequency is ``1/T``, corresponding to a sinusoid with a period equal to the signal duration.
    - The highest frequency that can be sampled, is called *Nyquist frequency*: ``\nu_\rm{Nyq} = \frac{1}{2\delta T} = \frac{1}{2}\frac{N}{T} ``.

- At the zero frequency, the FT value is just the sum of the signal values:

```math
a_0 = \sum\limits_{k=0}^{N-1} x_k e^{-2\pi i0k/N} = \sum\limits_{k=0}^{N-1} x_k
```

- *Parseval’s theorem* applies also to the discrete case and one can see that the variance of a signal is ``1/N`` times the sum of the ``a_j`` over all indices besides zero (also known as *Plancherel Theorem*):

```math
Var(x_k) = \sum_k (x_k - \bar{x})^2 = \sum_k x_k^2 + \sum_k \bar{x}^2 - \sum_k 2\bar{x}x_k = \sum_k x_k^2 + N \bar{x}^2 -  2N\bar{x}^2 = 
```

```math
= \sum_k x_k^2 - N\bar{x}^2 = \sum_k x_k^2 - \frac{1}{N}(\sum_k x_k)^2 = \frac{1}{N} \sum_j |a_j|^2 - \frac{1}{N}a_0^2 =
```

```math
\Longrightarrow Var(x_k) = \frac{1}{N}\sum_{j=-\frac{N}{2}}^{j=\frac{N}{2}-1} |a_j|^2,  j\ne0
```

"""

# ╔═╡ e2ca5449-f85f-44c8-adc3-eb257d2e3830
md"""

#### Exercise: let's compute the Fourier frequency grid for a simple case
***

- We assume to monitor a given phenomenum observing a quantity of interest for 100s and sampling it every 2s, i.e. we have 50 observations.
- Therefore, the lowest frequency (meaningful to analyse) turns out to be $1/T = 0.01$Hz, and the Nyquist frequency is  $\nu_\rm{Nyq} = \frac{1}{2\delta T} = \frac{1}{2 \times 2s} = 0.25$ Hz.
- The whole ( meaningful) Fourier frequency grid consists therefore of 50 entries, from $-\nu_\rm{Nyq}$ to $+\nu_\rm{Nyq}$, excluding the former or the latter value.
- Since for a real input function the power spectrum for negative or positive frequencies is identical, we restrict to the positive side.
- And the Fourier grid is finally, starting from 0 with 25 (i.e. $N/2$) $1/2s = 0.01$Hz steps, $0.0, 0.01, 0.02 \ldots 0.25$Hz.
"""

# ╔═╡ ec204c5e-7d37-4ffc-88de-3607f7f5fb07
begin
	function generate_noisy_signal(
	    duration::Float64 = 1000.0,
	    fs::Float64 = 100.0,
	    signal_freq::Float64 = 0.5,
	    amplitude::Float64 = 1.0,
	    noise_std::Float64 = 2.0
	)
	    t = 0:1/fs:duration-1/fs
	    signal = amplitude .* sin.(2π * signal_freq .* t)
	    noise = noise_std .* randn(length(t))
	    noisy_signal = signal .+ noise
	
	    return collect(t), signal, noisy_signal
	end
	
	
	# Genera la serie temporale
	t, clean_signal, noisy_signal = generate_noisy_signal()
end;

# ╔═╡ 50b7dca2-e81b-4982-93a4-31a2e13f11fa
function power_spectrum(signal::Vector{Float64}, fs::Float64)
    N = length(signal)

    # FFT e spettro di potenza (normalizzato)
    X = fft(signal)
    psd = (abs.(X).^2) ./ N

    # Solo frequenze positive (metà spettro)
    freqs = (0:N÷2) .* (fs / N)
    psd_one_sided = psd[1:N÷2+1]

    # Raddoppia le componenti (eccetto DC e Nyquist) per conservare la potenza
    psd_one_sided[2:end-1] .*= 2

    return freqs, psd_one_sided
end;

# ╔═╡ 9e780d0e-dcea-4308-bea2-aefbe212e60b
begin
	# --- Calcola spettri ---
	fs = 100.0
	freqs_noisy, psd_noisy = power_spectrum(noisy_signal, fs)
	freqs_clean, psd_clean = power_spectrum(clean_signal, fs)
end;

# ╔═╡ 4e3a40ee-3a09-4c6a-b2f4-43377d551483
begin
	fig2 = Figure(size = (1000, 750), fontsize = 13)
	
	# --- Pannello 1: Serie temporale ---
	ax1 = Axis(fig2[1, 1],
	    xlabel = "Time (s)",
	    ylabel = "Amplitude",
	    title  = "Time series",
	    xgridcolor = (:gray, 0.3),
	    ygridcolor = (:gray, 0.3),
	)
	
	lines!(ax1, t, noisy_signal,
	    color = (:steelblue, 0.6), linewidth = 0.8, label = "Signal + Noise")
	lines!(ax1, t, clean_signal,
	    color = :red, linewidth = 2.0, label = "Clean Signal (0.5 Hz)")
	axislegend(ax1, position = :rt)
	
	xlims!(10,50)
	
	# --- Pannello 2: Spettro di potenza ---
	ax2 = Axis(fig2[2, 1],
	    xlabel = "Frequency (Hz)",
	    ylabel = "PSD",
	    title  = "Power Spectrum (FFT)",
	    xgridcolor = (:gray, 0.3),
	    ygridcolor = (:gray, 0.3),
	    yscale = log10,          # scala logaritmica sull'asse Y
	    xticksvisible = true,
	)
	
	lines!(ax2, freqs_noisy, psd_noisy,
	    color = (:steelblue, 0.8), linewidth = 1.0, label = "Signal + Noise PSD")
	#lines!(ax2, freqs_clean, psd_clean,
	#    color = :red, linewidth = 1.5, label = "Clean Signal PSD")
	
	# Linea verticale sul picco a 0.5 Hz
	vlines!(ax2, [0.5], color = :orange, linewidth = 1.5, linestyle = :dash, label = "0.5 Hz")
	
	# Zoom sull'asse X fino a 1 Hz per vedere bene il picco
	xlims!(ax2, 0.45, 0.55)
	axislegend(ax2, position = :rt)
	
	fig2
end

# ╔═╡ 3f90207c-aff5-4c25-bb15-3900ec99ad78
md"""
- Top panel: the black line shows 40 seconds of a simulated time series 1000s long, consisting of a weak sinusoidal modulation (red line) "drowned" into a strong Gaussian noise. Bottom panel: corresponding PDS, zoomed on the relevant frequency range, where the modulation is clearly visible. A weak signal spread in time is collected into a single frequency bin at high significance.
"""

# ╔═╡ 58d43f9c-e7e9-4a7d-8471-1beca3568a1e
md"""
### Exercize about an analysis of weather data in France
***

- In the following exercize, we are going to analyse weather data spanning about 20 years in France obtained from the US National Climatic Data Center.

- Data are imported [http://www.ncdc.noaa.gov/cdo-web/datasets#GHCND](http://www.ncdc.noaa.gov/cdo-web/datasets#GHCND). The number "-9999" is used for N/A values. And we need to parse dates contained in the DATE column
"""

# ╔═╡ 3d75fecb-157b-408a-aa32-2d2096a766ba
begin
	url = "https://github.com/ipython-books/cookbook-2nd-data/blob/master/weather.csv?raw=true"
	#fname = "France_temperatures.csv"
	
	df = DataFrame(CSV.File(HTTP.get(url).body,missingstring="-9999"))
	#df = DataFrame(CSV.File(fname,missingstring="-9999"))
	
	df[!,:DateTime] = Date.(string.(df[!,:DATE]),DateFormat("yyyymmdd"))
end;

# ╔═╡ 2e2e09aa-6023-46a9-b816-bb9de2b321c1
md"""
- Let's now select only dates later than Jan 1994, and compute daily averages dropping the missing data and plot our dataset.

    - The temperature unit is in tenths of a degree, and we get the average value between the minimal and maximal temperature.
"""

# ╔═╡ 0148c42e-24fa-486e-8c61-eae25e772bd8
begin
	filter!(:DateTime => >=(Date("19940101",DateFormat("yyyymmdd"))),df)
	
	dropmissing!(df)
	
	gdf = combine(groupby(df, :DateTime), [:PRCP,:TMAX,:TMIN] .=> mean, renamecols=false)
	
	temp = (gdf[!,:TMAX] + gdf[!,:TMIN]) / 20.
	N = length(temp)
	
	
	fg1 = Figure()
	
	ax1fg1 = Axis(fg1[1, 1],
	    xlabel = "Date (yyyy/mm/dd)",
	    ylabel = L"Average daily temperature ($^\circ$C)"
	    )
	
	scatter!(gdf[!,:DateTime],temp,color=:blue)
	
	fg1
end

# ╔═╡ 7709a017-0c80-4399-9751-bf07b80f0b5e
md"""
- We now compute the Fourier transform and the spectral density of the signal using the **fft()** function.
- Once the FFT has been obtained, we take the square of its absolute value in order to get the **power spectral density (PSD)**.
- Then, we get the frequencies corresponding to the values of the PSD by the **fftfreq()** utility function.  
- Since the original unit is in days we change it to annual unit by the factor **1/365**.
"""

# ╔═╡ b49f1cda-582a-4d43-9e0a-a0e6adf85e78
begin
	# The FFTW library has been used for this exercize
	
	temp_fft = fft(temp)
	
	temp_psd = abs.(temp_fft).^2
	
	# 1 is the sampling frequency
	temp_freq = fftfreq(length(temp_psd), 1) * 365
end;

# ╔═╡ c601e37c-c08d-4ce4-8843-6fd52bd1e812
md"""
- We restrict to positive frequences only and let's check which is the maximum frequency for our dataset. Knowing the sampling rate we predict that it should correspond to a period of 2 days so that the Nyquist frequency turns out to be exactly 0.5 day$^{-1}$.
"""

# ╔═╡ aff12d14-4803-480e-aab0-6d33917304b3
begin
	posfr = temp_freq .> 0
	
	#println("Nyquist frequency: ", round(maximum(temp_freq) / 365, digits=1), L" day$^{-1}$")
	
end;

# ╔═╡ 6f29671a-048a-49e1-b7f9-e5d9cd733a0c
Markdown.parse("""
##### Nyquist frequency:  $(latexify(maximum(temp_freq) / 365,fmt="%.1f"))

""")

# ╔═╡ b8a2074a-490e-4e92-ba59-3735184d7c91
md"""
- Let's now plot the power spectral density of our signal, as a function of the frequency (in unit of **1/year**). We choose a logarithmic scale for the y axis (decibels).
"""

# ╔═╡ 40fe7f7d-8d30-40f1-b6fa-932374c93380
begin
	fg2 = Figure()
	
	ax1fg2 = Axis(fg2[1, 1],
	    xlabel = "Frequency (1/year)",
	    ylabel = "PSD (dB)",
	    )
	
	lines!(temp_freq[posfr],10 * log10.(temp_psd[posfr]),color=:blue)
	
	xlims!(0,maximum(temp_freq[posfr]))
	
	fg2
end

# ╔═╡ 1eb98eee-8373-4fce-b40f-3118d096f082
md"""
- Or with a better zoom in a region of our interest.
"""

# ╔═╡ 0f5adef3-5de0-4fa3-b4ca-8ee59968047e
md"""
Plot horizontal axis limit: $( @bind xmx PlutoUI.Slider(0.1:0.5:200, default=7) ) 
"""

# ╔═╡ 1b13be55-5153-4d53-ba09-92832f67a448
begin
	fg3 = Figure()
	
	ax1fg3 = Axis(fg3[1, 1],
	    xlabel = "Frequency (1/year)",
	    ylabel = "PSD (dB)",
	    )
	
	lines!(temp_freq[posfr],10 * log10.(temp_psd[posfr]),color=:blue)
	
	xlims!(0,xmx)
	ylims!(30,85)
	
	fg3
end

# ╔═╡ bf2f9df2-7df7-43d2-926d-7a561be79c0e
md"""
- Not surprisingly, the fundamental frequency of the signal is the yearly variation of the temperature at **f=1**.

- We can now "clean" our data cutting out frequencies higher than the fundamental frequency and by an **inverse FFT** we recover a signal that mainly contains the fundamental frequency.
"""

# ╔═╡ bea9ec97-57a7-45a4-9a29-ccca274fd5dc
begin
	temp_fft_bis = copy(temp_fft)
	temp_fft_bis[abs.(temp_freq) .> 1.1] .= 0
	
	temp_slow = real.(ifft(temp_fft_bis))
	
	
	l = gdf[!,:DateTime] .< Date(2000,1,1)
	
	
	fg4 = Figure()
	
	ax1fg4 = Axis(fg4[1, 1])
	
	scatter!(gdf[!,:DateTime][l],temp[l],color=:blue, label="Original data")
	lines!(gdf[!,:DateTime][l],temp_slow[l],color=:red, label="Filtered data")
	
	ylims!(-10,40)
	
	axislegend()
	
	fg4
end

# ╔═╡ b81793b6-dcec-4ab4-a806-4aa111d69f00
md"""
## What do we observe in the real world?
***
"""

# ╔═╡ c1363e2a-d8b4-4be9-8394-435231b62177
# ╠═╡ show_logs = false
cm"""
### Windowing and Sampling
***

- The CFT and the DFT can be connected easily taking into account that the FT of the product fo two functions is the convolution of the FT of the functions.

```math
F[x \cdot y] = F[x] \otimes F[y] = \int\limits_{-\infty}^{+\infty}F[x(\nu')]F[y(\nu-\nu')] d\nu'
```

- A discrete time series `x(t_k) ≡ x_k` can be seen as the product of a continuous function `f(t)` over `(−∞,∞)` and two additional functions: `w(t)` to limit it to the `(0,T)` interval and `s(t)` to sample it at times tk:

```math
x_k = h(t) \cdot w(t) \cdot s(t)
```

- ``w(t)`` is a boxcar window function, which is 1 in the ``(0,T) `` interval and zero outside. ``s(t)`` is a series of delta functions at ``t_k``, spaced by ``T/N``:

$(LocalResource("Pics/windowing.png"))
"""

# ╔═╡ f5fea622-f431-4da6-af3f-548fd10d906d
# ╠═╡ show_logs = false
md"""
### Windowing effects
***

- Let us consider a purely sinusoidal function $f(t) = sin(ωt)$, whose FT is a delta function at $ω$.

- The multiplication by the window function corresponds to the convolution of the delta function with the FT of the window.

- It is simple to calculate the FT of the window: we consider a window function that is unity in the $−T/2,T/2$ interval, as it is a real and even function, whose FT is also real and even:

```math
F(w(t)) = 2 {\sin(\pi\nu T)\over\pi\nu}
```

- This is the well known “sinc” function.

- An important general rule is that the FT peak is broader for shorter T. Something easily deducible from the general propertieds of FT.

    - The resolution of the signal FT is therefore higher the longer the observation is.
    
- In addition to the broadening, there is the formation of side lobes. They are much lower than the central peak, but cannot always be ignored.

$(LocalResource("Pics/window_ft.png"))


"""

# ╔═╡ 43caa1a2-3083-414d-9681-4153ade8a92b
md"""
### Sampling effects: aliasing
***

- The FT of a series of regularly spaced delta functions with spacing $T/N$ is itself a series of delta functions with spacing $N/T$:

```math
s(t)    = \sum\limits_{k=-\infty}^{+\infty} \delta(t - {kT\over N}) \Longleftrightarrow
    F(s(t)) = \sum\limits_{m=-\infty}^{+\infty} \delta(\nu - {mN\over T})
```
    
- Therefore, the effect of sampling on the FT of a sinusoidal signal with frequency $ν_0$ (a delta function at $ν_0$) is that of adding an infinite sequence of delta functions spaced by $N/T$, called *aliases*.

- Depending on the frequency of the original signal and the Nyquist frequency we can have different situations.
    - In fact, features at $ν = ν_{N/2} + ν_x$ also appear at $ν = ν_{N/2} − ν_x$.

- This happens because the transition from the CFT to the DFT involves two operations:
    - windowing, a convolution with the function $W(f)$, which is essentially a peak with a width $δf = 1/T$ plus sidelines,
    - and aliasing, a reflection of features above the Nyquist frequency back into the range ($0,ν_{N/2}$).
"""

# ╔═╡ dd78a7e2-f63e-4e84-af31-286d979074f5
# ╠═╡ show_logs = false
md"""
$(LocalResource("Pics/aliases.png"))
"""

# ╔═╡ 95095834-abd5-43e6-aad0-cdfa402b0366
md"""
- Top panel: Aliasing on a sinusoidal signal at 15 Hz sampled at 40 Hz. In black the true signal FT amplitude, in red the aliased one. The blue region marks the interval below $\nu_{Nyq}$. The signal is detected and the aliases are not.
- Bottom panel: same as the top panel, but with a signal at 35 Hz. Here the signal is not detected, but the 5 Hz alias is.
"""

# ╔═╡ 1cfc40d3-24a3-4a5c-8538-61c13542b403
# ╠═╡ show_logs = false
cm"""
- We have all experienced aliasing effects when looking at fast rotating objects like an air fan under fluorescent light.
    - The light provides a sampling at 50 Hz (or 60 Hz, depending on where you live), while the fan has a periodicity.
    - Depending on its angular speed, you can see it rotating apparently much slower, or even to stop and rotate in the opposite direction.

$(LocalResource("Pics/sampling.png"))

- Time domain example of aliasing. The blue signal has ``\nu_0 = 0.1`` Hz. If it is sampled with ``\nu_{Nyq} = 0.038`` Hz (red points) the red dashed alias is the best fit to the data, with ``\nu_a = 0.023`` Hz.
"""

# ╔═╡ e8786a24-4391-40ad-a9ba-a0cb3f7bb2b0
md"""
- In the real world we do not really sample signals, but integrate them over finite time bins, i.e. we convolve it with a binning function:

```math
b(t) = \begin{cases}
{N\over T} & t\in [-{T\over 2N},{T\over 2N}]\\
0                    & {\rm outside}
\end{cases}
```

- Therefore, the signal FT will be multiplied by that of the binning function, which is again a sinc function:

```math
B(\nu) = {\sin \pi\nu /2\nu_{Nyq}\over \pi\nu /2\nu_{Nyq}}
```

- $B(ν)$ is a broad function that reaches $0$ at $2\nu_{Nyq}$ and has the value of $2/π$ at $\nu_{Nyq}$.
"""

# ╔═╡ 39e1ee25-a2f4-4e16-9fec-4a4e9c2c653e
md"""
### Effects of binning on the harmonic analysis
***

- Binning is a common operation carried out for different reasons, e.g., for casting an irregularly sampled time series to a regular grid, for improving the S/N of each point, to reduce the computational burden, etc.

> Nevertheless, binning is never a price-free operation.

- Let's assume we have a dataset defined as $D = \{y_1, y_2,..., y_n \}$, with observations carried out at a regular time step. Its discrete Fourier transform is:

```math
Y(\omega) \equiv \sum_{t=1}^n y_t e^{i \omega t}
```

- This is defined for continuous values of $\omega$ but we know that no loss of data happens if we compute the Fourier transform at the Fourier grid $\omega_k \equiv 2 \pi k / n, \quad 0 \le k < n$.


#### Moving average
***

- Let's suppose our data are replaced by a moving average of past values (i.e. a rebinning):

```math
z_t \equiv \sum_{s=0}^{m-1} y_{t-s} \omega_s
```

- where $\omega_s$ is the weighting coefficient for lag $s$.

- Since this can be seen as a convolution of the input data we know that the Fourier transform of $z$ is the product of the Fourier transform of the original data and of the average function:

```math
Z(\omega) = W(\omega) Y(\omega)
```

- where $W(\omega) = \sum_{s=0}^{m-1} \omega_s e^{i \omega s}$.

- In particular, for uniform weighting, $w_s = 1/m, 0 \le s < m$, and we have the well known ["sinc"](https://en.wikipedia.org/wiki/Sinc_function) function:

```math
W(\omega) = \frac{1}{m} \sum_{s=0}^{m-1} e^{-i \omega s} = e^{-i\frac{\omega}{2}(m-1)} \left[\frac{\sin(m\omega/2)}{m \sin(\omega/s)}\right]
```

- It is clear, therefore, that together with changing the phase of the Fourier transform, any binning typically decreases the amplitude of the Fourier transform of the original data: i.e. it might be useful for plotting purposes and for saving computer power, yet binning should typically be avoided.
"""

# ╔═╡ 136b29c5-ac92-4e32-878e-3ad99216f78d
md"""
### Window carpentry
***

- Having a longer observation reduces the width of the main window peak, but does not change the possible spillover effects.

- In some cases it can be advantageous to multiply the data by another window, not boxcar-shaped.

- This results in a loss of signal, as some data are multiplied by a factor less than unity, but there are advantages, depending on the chosen window.

    - Many window functions with different characteristics have been designed and one can tailor them depending on what is needed.

- The main features that identify a window in its PDS are: the width of the main peak $∆ω$, the relative amplitude of the first side lobe $L$ (expressed in decibels) and the slope of the decay of side lobes n:
"""

# ╔═╡ d490adbd-376a-403b-adc3-3ab4d4e65bc5
# ╠═╡ show_logs = false
md"""
$(LocalResource("Pics/window_features.png"))
"""

# ╔═╡ 848e6f85-e3b9-4928-aed0-b316f97d18db
md"""
- The boxcar window is the one with the lowest $∆ω$, but with alternative windows it is possible to obtain a significant reduction of the amplitude of the side lobes.

| Window   | $\Delta\omega$ | L     | n  | Function                       |
|:--------:|:-------------: |:-----:|:--:|:------------------------------:|
| Boxcar   |     0.89       | -13db | 2  |         1                      |
| Hamming  |     1.36       | -43db | 2  | 0.54+0.46 cos(2πt)             |
| Hann     |     1.44       | -32db | 5  | 0.5(1-cos(2πt))                |
| Blackman |     1.68       | -58db | 5  | 0.42+0.5cos(2πt)+0.08cos(4πt)  |
| Gaussian |     1.55       | -56db | 2  | $\exp(-4.2 x^2$)               |
"""

# ╔═╡ 8d5c7719-647a-4e7e-905f-43051bb26716
# ╠═╡ show_logs = false
md"""
$(LocalResource("Pics/windows.png"))
"""

# ╔═╡ a9964f9d-1f3a-493b-9ab8-a8933ff51c72
md"""
- Left panel: the shape of five windows: boxcar, Hamming, Hann, Blackman and Gaussian. Right panel: the corresponding PDS, corresponding to (half) the shape of the PDS of a sinusoidal signal. The sidelobes of all but the boxcar window are too low to be seen.

"""

# ╔═╡ a0392d53-9e0d-4d70-a9a5-b99de7884f9d
# ╠═╡ show_logs = false
md"""
### Effect of gaps
***

$(LocalResource("Pics/wtot.png"))
"""

# ╔═╡ eae97444-2431-448c-abe6-c14aa64fa409
md"""
- Each of the four panels contains a signal (top) and its PDS (bottom). In all four, the signal consist of Gaussian noise plus a sinusoid at P=200 s, with 1-s binning.
    - Top left: continuous exposure of $10^5$ s.
    - Top right: continuous exposure of $10^4$ s.
    - Bottom left: exposure of $10^4$ s split into three intervals. The PDS is computed including the gaps as consisting of points at zero level.
    - Bottom right: exposure of $10^4$ s split into three intervals as in the previous case. The gaps are filled with Gaussian noise.
"""

# ╔═╡ 8d3eb863-5cb0-4221-8fc8-c48361615e95
md"""
### Fast Fourier transform
***

- Evaluation of the Discrete Fourier Transform (DFT) of $N$ samples involves $\sim N^2$ multiplication and addition operations -- for every Fourier component $a_j$, each of the $N$ samples $x_k$ needs to be multiplied by a phase factor $e^{-2\pi ijk/N}$ and then they have to be summed.  

- The Fast Fourier Transform algorithm has been devised to accomplish this computation in much fewer steps, typically with $\sim N\log_2 N$ multiplications and additions.  This provides enormous computational savings for large transforms, and has made Fourier analysis accessible to cases where it would have been otherwise prohibitive.
"""

# ╔═╡ 0b0b8578-4c67-463a-bd53-b3e6bd1f9712
md"""
### PSD normalization
***

- The FT is a linear transformation and the PSD is its squared modulus, the PSD then scales with the square of the intensity level of the signal.

- It is possible to normalize the PSD in different ways. One of the most common is the so-called *Leahy normalization*:

```math
P_j^{Leahy} = {2\over N_\gamma} |a_j|^2
```

- where $N_\gamma$ is the total number of photons in the signal. If instead of counts we have fluxes, and noise is $\mathcal{N}(0,σ)$, $N_\gamma$ is substituted by $N_\rm{data}σ^2$.

- The periodogram with this normalization is also known as *Classical* or *Schuster* periodogram.

- With the Leahy normalization the total variance can also be written as:

```math
Var(x_k) = \frac{1}{N}\sum_{j=-\frac{N}{2}}^{j=\frac{N}{2}-1} |a_j|^2 \ (j\ne0) =  \frac{N_\gamma}{N} \left( \sum_{j=1}^{j=\frac{N}{2}-1} P_j + \frac{1}{2} P_{N/2}\right)
```

- This normalization leads to a known statistical distribution of signal power:
    - if the signal is dominated by fluctuations due to Poisson statistics and if $N_\gamma$  is large, powers follow a $\chi^2$ distribution with 2 degrees of freedom, $<P>=2$ and $Var(P)=4$.
- The reason is that the periodogram is the sum of the squares of the real and imaginary parts of the FT and, for a stochastic process, the latter are normally distributed, so the sum of their squares is distributed as a $\chi^2$ with 2 degrees of freedom.

- For other noise distributions (Poisson, etc.), for large N, due to the central limit theorem the real and imaginary parts become still normal.


"""

# ╔═╡ 8e458457-f940-460a-8638-debe1f94defc
md"""
- Periodograms are intrinsically very noisy. If the signal is divided into $S$ segments and the resulting PDS are averaged and rebinned by a factor $M$, the powers will be distributed as a $\chi^2$ with $2SM$ degrees of freedom scaled by $1/SM$: therefore, the average power remains $<P>=2$, but the variance is now $Var(P)=4/SM$.
    - The technique of dividing the time series into equal-duration intervals and averaging the corresponding PSD is called *Bartlett’s method*.

- The reduction in time duration $T$ increases the minimum frequency in the PSD $ν_\rm{min} = 1/T$.

- This method, in principle, also allows one to skip over (short) data gaps, which have dramatic effects on the PSD.


"""

# ╔═╡ b2e2bace-7e50-41f0-87fb-780737ff6616
md"""
### Exercise about PSD manipulation
***

- Let's generate two arrays of relative timestamps, one 8 seconds long and one 1600 seconds long, with dt = 0.03125 s, and make two signals in units of counts. The signal is a sine wave with amplitude = 300 cts/s, frequency = 2 Hz, phase offset = 0 radians, and mean = 1000 cts/s. We then add Poisson noise to the light curve and plot the shortest of the two.
"""

# ╔═╡ 1112b278-802d-4611-8954-2795db178e88
begin
	pf = @bind pl_fr NumberField(1:10, default=2)
	af = @bind am_fr NumberField(100:100:1000, default=300)
end;

# ╔═╡ 19a6c251-adc1-422d-a702-20ea9b971449
cm"Frequency for the test signal:"

# ╔═╡ a8b99222-18a7-4b47-8f98-8961a708ec57
pf

# ╔═╡ 77b0596d-8c4a-4d3d-9bd8-54d3fde618c8
cm"Amplitude for the test signal:"

# ╔═╡ d542c639-1c93-4edd-9132-facd18fc89fb
af

# ╔═╡ b96525bf-90e5-47b3-b305-1db0ad675c2c
begin
	dt = 0.03125  # seconds
	exposure = 8.  # seconds
	long_exposure = 1600. # seconds
	times = range(start=0, stop=exposure-dt, step=dt)  # seconds
	long_times = range(start=0, stop=long_exposure-dt, step=dt)  # seconds
	
	signal = am_fr .* sin.(pl_fr .* pi .* times ./ 0.5) .+ 1000  # counts/s
	long_signal = am_fr .* sin.(pl_fr .* pi .* long_times ./ 0.5) .+ 1000  # counts/s
	
	noisy = [rand(Poisson(theta)) for theta in signal .* dt]  # counts
	long_noisy = [rand(Poisson(theta)) for theta in long_signal .* dt]  # counts
	
	fg5 = Figure()
	
	ax1fg5 = Axis(fg5[1, 1],
	    xlabel = "Time (s)",
	    ylabel = "Counts (cts)",
	    )
	
	lines!(times,noisy,color=:blue,label="Noisy signal")
	lines!(times,signal .* dt,color=:green,label="Signal")
	
	axislegend()
	
	fg5
end

# ╔═╡ ae8a4ebf-da8d-4778-ad9a-87e3aab72399
md"""
- Now let's compute the periodogram
"""

# ╔═╡ b7213b38-21ac-4fb5-ab53-694d3ac2fda9
pf

# ╔═╡ 78dc57ce-1cbd-4f9c-9741-54bf430b9186
af

# ╔═╡ 25fcf199-7da4-4632-a2fa-d0131f58d533
begin
	# We use the DSP package
	
	psd = periodogram(noisy; fs=1/dt)
	
	fg6 = Figure()
	
	ax1fg6 = Axis(fg6[1, 1],
	    yscale = log10,
	    xlabel = "Frequency (Hz)",
	    ylabel = "Power",
	    )
	
	lines!(psd.freq,psd.power,color=:blue,label="PSD")
	
	xlims!(1,15)
	
	axislegend()
	
	
	fg6
end

# ╔═╡ 3e27b7d8-9dec-4e54-9050-5ca9275d8499
md"""
- The computed periodogram (with this function) is normalized so that the area under the periodogram is equal to the uncentered variance (or average power) of the original signal.

- Since the negative Fourier frequencies (and their associated powers) are discarded, the number of time bins per segment `n` is twice the length of `freq` and `power`.

- The zero frequency is the sum of the signal value, and it not typically used.
"""

# ╔═╡ 817993c0-a1ce-46d2-a213-ed6448d5159f
Markdown.parse("""
Number of data points:   $(latexify(length(noisy)))

Number of data points:   $(latexify(length(psd.freq[2:end])))
""")

# ╔═╡ 7fedf685-66f5-4d46-b52e-99c7305d3f95
md"""
- The power spectrum is a bit noisy. Let's try averaging together power spectra from multiple segments of data using the long time series.

- We want to average periodograms computed each 8 seconds, averaging therefore 1600/8 = 200 periodograms of 256 elements each.
"""

# ╔═╡ 9985012a-ad66-434b-84c2-13eead6a1af7
pf

# ╔═╡ 39538b6b-414d-436a-a1c4-19c6f6754dc7
begin
	apsd = welch_pgram(long_noisy, 256, 0, fs=1/dt)
	
	fg7 = Figure()
	
	ax1fg7 = Axis(fg7[1, 1],
	    yscale = log10,
	    xlabel = "Frequency (Hz)",
	    ylabel = "Power",
	    )
	
	lines!(apsd.freq,apsd.power,color=:blue,label="Averaged PSD")
	
	xlims!(1,15)
	
	axislegend()
	
	fg7
end

# ╔═╡ 549cc9b7-01ff-4639-b8cb-66df220279a9
md"""
- With a clear increase of the S/N.

- Let's try now to compute periodograms using different windows functions.
"""

# ╔═╡ cef2b9a7-701e-4975-8e71-f8a0afff1641
pf

# ╔═╡ b8c5baed-5a58-4380-b941-a7d4035c8d01
begin
	psd_rect = periodogram(noisy; fs=1/dt, window=nothing)
	psd_ham = periodogram(noisy; fs=1/dt, window=hamming)
	psd_tri = periodogram(noisy; fs=1/dt, window=triang)
	psd_cos = periodogram(noisy; fs=1/dt, window=cosine)
	
	
	
	fg8 = Figure()
	
	ax1fg8 = Axis(fg8[1, 1],
	    yscale = log10,
	    xlabel = "Frequency (Hz)",
	    ylabel = "Power",
	    )
	
	lines!(psd_rect.freq,psd_rect.power,label="Rectangular")
	lines!(psd_ham.freq,psd_ham.power,label="Hamming")
	lines!(psd_tri.freq,psd_tri.power,label="Triangular")
	lines!(psd_cos.freq,psd_cos.power,label="Cosine")
	
	
	xlims!(1,5)
	
	axislegend()
	
	fg8
end

# ╔═╡ ee646903-dc54-46d0-b2f5-b30e41a68853
md"""
### Auto and Cross-Correlation
***

- The Cross-correlation of two functions $f(t)$ and $g(t)$ is defined as:

```math
C(\tau)=f\star g=\int\limits_{-\infty}^{\infty}f^*(t)g(t+\tau)dt
```

- The result is a function of the *lag* $τ$ introduced between the two functions, and is often used to estimate the similarity between two different time series, as a function of lag.
- If a common underlying process causes the time variation of intensity at two different electromagnetic bands with differential delays while propagating to the observer, then the cross correlation function of the two time series will exhibit a peak at the corresponding lag, namely the relative delay between the two bands.

- The autocorrelation function is a special case where a function is correlated with itself, which would always show a peak at zero lag.

> Convolution is an operation akin to the Cross-correlation, but the function $g(t)$ in the integrand is inverted to $g(−t)$ before adding the shift.

- A few important properties of cross-correlation include:

```math
[f\star g](\tau) = [g^*\star f^*](-\tau)
```

```math
[f\star g]\star[f\star g] = [f\star f]\star[g\star g]
```

```math
g\star(f\otimes h)  =  [g\star f]\otimes h \\
```

```math
F[f\star g]  =  F(f)\cdot F^*(g)
```

- where $F$ represents the Fourier transform.
"""

# ╔═╡ 454fc15d-ea34-4938-b86e-76f2cb0d5125
md"""
> We now propose a few definitions for completeness, but without extensive discussions.

### Cross-spectra, phase lag spectra, coherence
***

- Given two signals $f(t)$ and $g(t)$ and their respective FTs, $F(ω)$ and $G(ω)$, we define the cross spectrum as:

```math
CS(\omega) = F(\omega)\cdot G^*(\omega)
```

- Analogous to the PSD and the autocorrelation, the cross spectrum between two signals is the Fourier transform of their cross-correlation (and vice-versa).

- In its essence, the cross spectrum of two signals at each frequency is a complex number whose argument represents the phase delay between the signals at that frequency.

- The autocorrelation function may be thought of as the second order correlation:

```math
c_2(\tau)=\langle f(t)f(t+\tau) \rangle
```

- where the angular brackets denote an ensemble average.

- Bispectrum is an extension of the above concept to triple correlations. The third order correlation function:

```math
c_3(\tau_1,\tau_2)=\langle f(t)f(t+\tau_1)f(t+\tau_2)\rangle
```

- that allows one to define the Bispectrum:

```math
B(\omega_1,\omega_2)=\int\limits_{-\infty}^{\infty}\int\limits_{-\infty}^{\infty} c_3(\tau_1,\tau_2) e^{-i(\omega_1\tau_1+\omega_2\tau_2)}\,d\tau_1\,d\tau_2
```
"""

# ╔═╡ 6e178c7b-1846-4390-b842-b4d303688bd8
md"""
## Power Spectrum statistics
***

- The first thing we need to know is the probability distribution of the noise power.

- We now assume that it is additive and independent of the frequency (i.e. noise is “white”):

```math
P_j = P_{j,noise} + P_{j,signal}
```

- The “null hypothesis” is that the periodogram is consistent with pure noise.

- Let’s remind us that If we have $x_k ≡ y_k + z_k$ and $b_j$ and $c_j$ are the FT of $y_k$ and $z_k$, respectively, we have $a_j = b_j + c_j$.

- This does not hold for power spectra unless the signal are uncorrelated random noise:

```math
|a_j|^2 = |b_j + c_j|^2 = |b_j|^2 + |c_j|^2 + {\rm cross\ terms}
```

- For a wide range of types of noise, $P_{j,noise}$ follows a $\chi^2$ distribution with 2 degrees of freedom (but at the Nyquist frequency, where dof is 1).

- For other noise distributions (Poisson, etc.), for large N, due to the central limit theorem the Fourier coefficients $A_j$ and $B_j$ become still normal.

- This suggests a simple consistency test: compute the standard deviation in each frequency bin and divide by the mean power. Results should be (if the hypotheses hold) close to 1.

- In practice one finds that noise powers are nearly always $\chi^2$ distributed, not only for Poisson noise, but also for many other types of noise.

- With the Leahy normalization, the probability for $P_{j,noise}$ to exceed a given threshold is given by:

```math
Prob(P_{j,noise} > P_{threshold}) = Q(P_{threshold} | 2) \quad (j = 1, N/2-1)
```

```math
Q(\chi^2|\nu) \equiv \left[2^{\nu/2} \Gamma(\frac{\nu}{2})\right]^{-1} \int_{\chi^2}^\infty t^{\frac{\nu}{2}-1} e^{-\frac{t}{2}} dt
```

- where $\nu$ is the number of dof and $\Gamma$ is the gamma function, the generalization of factorial to real and complex numbers [$\Gamma(n) = (n-1)!$]

- For $\nu = 2$:

```math
Q(\chi^2|2) = \frac{1}{2} \int_{\chi^2}^\infty e^{-\frac{t}{2}} dt = e^{-\frac{\chi^2}{2}}
```

"""

# ╔═╡ 83aece6e-cc50-4ce4-b9b5-f518342620e5
md"""
- Power spectra are unavoidably very noisy.

    - Standard deviation of noise power is equal to the mean value ($σ_{P_j} = < P_j > = 2$).

- More interesting, this cannot be improved increasing the number of data points (i.e. the length of the time series). This merely increases the number of powers.

- One can decrease the large variance rebinning the power spectrum and/or dividing the data in multiple segments of equal length. This of course degrade the frequency resolution.  

    - If the number of segments is large the power statistics tends to become Normal.

- Let’s define  confidence detection level as the power with only $ε$ probability to be exceeded by noise.

- This holds for a single frequency. If your spectrum consists of $N_{\rm trial}$ (independent) frequencies, the confidence detection level decreases to take into account the multiple trials:

$$(1 - \epsilon)^{N_{\rm trial}} \sim 1 - \epsilon N_{\rm trial} \quad {\rm for\ } \epsilon << 1$$

- In general if noise is not Poissonian or Gaussian noise power spectrum will not be flat anymore.

    - However, often noise powers still follow a $\chi^2$ distribution with 2 dof, but with a different normalisation (in general depending on $j$).
"""

# ╔═╡ 06df9913-445f-4170-b110-13ec4208d940
cm"""
### The Likelihood for Periodograms
***

- The ``\chi^2`` distribution defines a sampling distribution or likelihood, i.e. the probability distribution of observing a given data set given some underlying (true, unknown) power spectrum.

- If we define a model power ``S_j(θ)`` at frequency ``ν_j``, specified by a set of parameters ``θ``, we can then compute the probability of having observed periodogram power ``P_j`` at that same frequency:

```math
p(P_j|S_j(θ)) = \frac{1}{S_j(θ)} e^{-\frac{P_j}{S_j(θ)}}
```

- The likelihood (also known as *Whittle likelihood*) for a periodogram over ``N/2`` observed powers ``P_j`` is then defined as the product of individual probabilities for each frequency ``ν_j``.

- One generally (and equivalently), defines the logarithm of the likelihood as the sum of logarithm of all probabilities, such that:

```math
\log(\mathcal{L}(\theta)) = \sum_{j=1}^{N/2} \log(p(P_j|S_j(θ))) =  -\sum_{j=1}^{N/2} ( \log(S_j(θ)) + \frac{P_j}{S_j(θ)})
```

- A likelihood for averaged periodograms can be derived from the ``\chi^2_{2LM}/2ML`` sampling distribution for periodograms averaged over ``L`` independent segments and ``M`` independent neighbouring frequencies by:

```math
\log(\mathcal{L}_{avg}(\theta)) = -2ML\sum_{j=1}^{N/2} \left[ \frac{P_j}{S_j(θ)}) + \log(S_j(θ)) + (\frac{1}{ML} - 1) \log(P_j) + c(2ML) \right]
```

- where ``c(2ML)`` is a factor independent of ``P_j`` or ``S_j``, and thus unimportant to the parameter estimation problem considered here.

"""

# ╔═╡ 9a17f40a-2aad-4b52-96a4-7c275fda9cee
# ╠═╡ show_logs = false
md"""
### Noise color
***

- Fourier analysis (and related techniques) has proven to be very effective in identifying periodic behaviours.

- In astrophysics, however, we often have to deal with phenomena too long for having a reliable coverage (decades or more) and with only approximately cyclical behaviours (quasi-periodicities).

- More important, in order to compute the statistical significance of any possible periodicity, one needs to properly model the noise affecting a time series.

$(LocalResource("Pics/colnoise.jpg"))

- "Coloured" noise directly affect a time-series shape:

$(LocalResource("Pics/redblunoise.png"))
"""

# ╔═╡ a76024f5-7355-4d57-9167-ed522a76ff51
md"""
## Final considerations
***

- Fourier analysis reveals nothing about the evolution in time, but rather reveals the variance of the signal at different frequencies.

- The classical periodogram is an estimator of the spectral density, i.e. the Fourier transform of the autocovariance function.

- Fourier analysis has restrictive assumptions: an infinitely long data of equally-spaced observations; homoscedastic Gaussian noise with purely periodic signal of sinusoidal shape.

- The classical periodogram is not a good estimator, it is “inconsistent” because the number of parameters grows with the number of datapoints.

- The DFT and its probabilities depends on several strong assumptions that are rarely achieved in real astronomical data: evenly spaced data of infinite duration with a high sampling rate (Nyquist frequency), Gaussian noise, single frequency periodicity with sinusoidal shape and stationary behavior.

- Each of these constraints is often violated in various astronomical problems. Data spacing may be affected by daily/monthly/orbital cycles. Periods may be comparable to the sampling time. Several periods may be present (e.g. helioseismology). Shape may be non-sinusoidal (e.g. elliptical orbits, eclipses, recurrent flares). Periods may not be constant (e.d. QPOs in accretion disks).
"""

# ╔═╡ ddfd682b-56a7-49c9-883a-0d348c9cb8c8
md"""
## Reference & Material

Material and papers related to the topics discussed in this lecture.

- [Belloni & Bhattacharya (2022) - "Basics of Fourier Analysis for High-Energy Astronomy”](https://ui.adsabs.harvard.edu/abs/2022hxga.book....7B/abstract)
- [van der Klis (1988) - "Fourier techniques in X-ray timing"](https://ui.adsabs.harvard.edu/abs/1989ASIC..262...27V/abstract)
"""

# ╔═╡ 95ec0443-95e3-4986-b170-7589328246b2
md"""
## Further Material

Papers for examining more closely some of the discussed topics.

- [Vaughan (2010) - "A Bayesian test for periodic signals in red noise"](https://ui.adsabs.harvard.edu/abs/2010MNRAS.402..307V/abstract)
- [Barret & Vaughan (2012) - "Maximum Likelihood Fitting of X-Ray Power Density Spectra: Application to High-frequency Quasi-periodic Oscillations from the Neutron Star X-Ray Binary 4U1608-522](https://ui.adsabs.harvard.edu/abs/2012ApJ...746..131B/abstract)
- [Covino eta al. (2019) - "Gamma-ray quasi-periodicities of blazars. A cautious approach"](https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.1270C/abstract)
"""

# ╔═╡ 22c6e22a-4b5d-47cd-b4b1-97543a7c9d74
md"""
### Credits
***

This notebook contains material obtained from [https://www.tutorialspoint.com/power-spectral-density-psd-and-autocorrelation-function#](https://www.tutorialspoint.com/power-spectral-density-psd-and-autocorrelation-function#) and from [https://ipython-books.github.io/101-analyzing-the-frequency-components-of-a-signal-with-a-fast-fourier-transform/](https://ipython-books.github.io/101-analyzing-the-frequency-components-of-a-signal-with-a-fast-fourier-transform/).
"""

# ╔═╡ 27d77a9f-afa4-4b79-87f8-f8a92b87381e
cm"""
## Course Flow

<table>
  <tr>
    <td>Previous lecture</td>
    <td>Next lecture</td>
  </tr>
  <tr>
    <td><a href="./open?path=Lectures/Lecture-StatisticsReminder/Lecture-BayesianReminder.jl">Lecture about Bayesian statistics</a></td>
    <td><a href="./open?path=Lectures/ScienceCase-SunspotNumber/Lecture-SunspotNumber.jl">Science case about Sunspot number</a></td>
  </tr>
 </table>


"""

# ╔═╡ 321154fd-5092-4d93-ba82-1bcde10efcb5
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
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Format = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.10.16"
CairoMakie = "~0.15.9"
CommonMark = "~1.0.1"
DSP = "~0.8.4"
DataFrames = "~1.8.1"
Distributions = "~0.25.123"
FFTW = "~1.10.0"
Format = "~1.3.7"
HTTP = "~1.11.0"
LaTeXStrings = "~1.4.0"
Latexify = "~0.16.10"
PlutoUI = "~0.7.79"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.5"
manifest_format = "2.0"
project_hash = "a6abc2fde11c66c8cd35d8fa95d1ecc3cd7242c2"

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

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

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

[[deps.DSP]]
deps = ["Bessels", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "5989debfc3b38f736e69724818210c67ffee4352"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.8.4"
weakdeps = ["OffsetArrays"]

    [deps.DSP.extensions]
    OffsetArraysExt = "OffsetArrays"

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

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

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

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "Setfield", "SparseArrays"]
git-tree-sha1 = "2d99b4c8a7845ab1342921733fa29366dae28b24"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.1"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"
    PolynomialsRecipesBaseExt = "RecipesBase"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

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

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

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
# ╟─d9329a8d-b07c-4207-93a6-1668da23e296
# ╟─3307644d-e875-4f15-a140-a41f7ca82a8f
# ╟─6a59a8e6-39fc-44b3-9921-8976c677f4b1
# ╟─771d6af3-a5e2-4869-82c3-68540f71cb41
# ╟─95ee75d5-3112-42b6-83aa-29e639ac6eb0
# ╟─b175b388-911f-4862-afdb-7449fec2bf9e
# ╟─882fc480-b132-45a8-906a-db157059b92c
# ╟─a6681a0b-2cbb-4cbd-ad44-7f01953f875f
# ╟─bc82a078-fea7-4a00-b38c-a849fb760594
# ╟─3bdf3bea-7e3c-4f09-98f7-bdbc16ea2880
# ╟─0dbb5f0e-292d-4d8a-bb83-7954e4a46f29
# ╟─875732fc-5644-4c06-8f8c-39dc30f1d3b4
# ╟─56081488-4ade-436a-991e-d7138368edf1
# ╟─6d88e888-8f08-47ab-8906-c7b4b86b2588
# ╟─99f8ea18-5807-482e-8516-e2bdaf1546e6
# ╟─e087bf1e-dd1a-41f3-8154-bfb39325d905
# ╟─e2ca5449-f85f-44c8-adc3-eb257d2e3830
# ╠═ec204c5e-7d37-4ffc-88de-3607f7f5fb07
# ╠═50b7dca2-e81b-4982-93a4-31a2e13f11fa
# ╠═9e780d0e-dcea-4308-bea2-aefbe212e60b
# ╠═4e3a40ee-3a09-4c6a-b2f4-43377d551483
# ╟─3f90207c-aff5-4c25-bb15-3900ec99ad78
# ╟─58d43f9c-e7e9-4a7d-8471-1beca3568a1e
# ╟─3d75fecb-157b-408a-aa32-2d2096a766ba
# ╟─2e2e09aa-6023-46a9-b816-bb9de2b321c1
# ╟─0148c42e-24fa-486e-8c61-eae25e772bd8
# ╟─7709a017-0c80-4399-9751-bf07b80f0b5e
# ╟─b49f1cda-582a-4d43-9e0a-a0e6adf85e78
# ╟─c601e37c-c08d-4ce4-8843-6fd52bd1e812
# ╟─aff12d14-4803-480e-aab0-6d33917304b3
# ╟─6f29671a-048a-49e1-b7f9-e5d9cd733a0c
# ╟─b8a2074a-490e-4e92-ba59-3735184d7c91
# ╟─40fe7f7d-8d30-40f1-b6fa-932374c93380
# ╟─1eb98eee-8373-4fce-b40f-3118d096f082
# ╟─0f5adef3-5de0-4fa3-b4ca-8ee59968047e
# ╟─1b13be55-5153-4d53-ba09-92832f67a448
# ╟─bf2f9df2-7df7-43d2-926d-7a561be79c0e
# ╟─bea9ec97-57a7-45a4-9a29-ccca274fd5dc
# ╟─b81793b6-dcec-4ab4-a806-4aa111d69f00
# ╟─c1363e2a-d8b4-4be9-8394-435231b62177
# ╟─f5fea622-f431-4da6-af3f-548fd10d906d
# ╟─43caa1a2-3083-414d-9681-4153ade8a92b
# ╟─dd78a7e2-f63e-4e84-af31-286d979074f5
# ╟─95095834-abd5-43e6-aad0-cdfa402b0366
# ╟─1cfc40d3-24a3-4a5c-8538-61c13542b403
# ╟─e8786a24-4391-40ad-a9ba-a0cb3f7bb2b0
# ╟─39e1ee25-a2f4-4e16-9fec-4a4e9c2c653e
# ╟─136b29c5-ac92-4e32-878e-3ad99216f78d
# ╟─d490adbd-376a-403b-adc3-3ab4d4e65bc5
# ╟─848e6f85-e3b9-4928-aed0-b316f97d18db
# ╟─8d5c7719-647a-4e7e-905f-43051bb26716
# ╟─a9964f9d-1f3a-493b-9ab8-a8933ff51c72
# ╟─a0392d53-9e0d-4d70-a9a5-b99de7884f9d
# ╟─eae97444-2431-448c-abe6-c14aa64fa409
# ╟─8d3eb863-5cb0-4221-8fc8-c48361615e95
# ╟─0b0b8578-4c67-463a-bd53-b3e6bd1f9712
# ╟─8e458457-f940-460a-8638-debe1f94defc
# ╟─b2e2bace-7e50-41f0-87fb-780737ff6616
# ╟─1112b278-802d-4611-8954-2795db178e88
# ╟─19a6c251-adc1-422d-a702-20ea9b971449
# ╟─a8b99222-18a7-4b47-8f98-8961a708ec57
# ╟─77b0596d-8c4a-4d3d-9bd8-54d3fde618c8
# ╟─d542c639-1c93-4edd-9132-facd18fc89fb
# ╟─b96525bf-90e5-47b3-b305-1db0ad675c2c
# ╟─ae8a4ebf-da8d-4778-ad9a-87e3aab72399
# ╟─b7213b38-21ac-4fb5-ab53-694d3ac2fda9
# ╟─78dc57ce-1cbd-4f9c-9741-54bf430b9186
# ╟─25fcf199-7da4-4632-a2fa-d0131f58d533
# ╟─3e27b7d8-9dec-4e54-9050-5ca9275d8499
# ╟─817993c0-a1ce-46d2-a213-ed6448d5159f
# ╟─7fedf685-66f5-4d46-b52e-99c7305d3f95
# ╟─9985012a-ad66-434b-84c2-13eead6a1af7
# ╟─39538b6b-414d-436a-a1c4-19c6f6754dc7
# ╟─549cc9b7-01ff-4639-b8cb-66df220279a9
# ╟─cef2b9a7-701e-4975-8e71-f8a0afff1641
# ╟─b8c5baed-5a58-4380-b941-a7d4035c8d01
# ╟─ee646903-dc54-46d0-b2f5-b30e41a68853
# ╟─454fc15d-ea34-4938-b86e-76f2cb0d5125
# ╟─6e178c7b-1846-4390-b842-b4d303688bd8
# ╟─83aece6e-cc50-4ce4-b9b5-f518342620e5
# ╟─06df9913-445f-4170-b110-13ec4208d940
# ╟─9a17f40a-2aad-4b52-96a4-7c275fda9cee
# ╟─a76024f5-7355-4d57-9167-ed522a76ff51
# ╟─ddfd682b-56a7-49c9-883a-0d348c9cb8c8
# ╟─95ec0443-95e3-4986-b170-7589328246b2
# ╟─22c6e22a-4b5d-47cd-b4b1-97543a7c9d74
# ╟─27d77a9f-afa4-4b79-87f8-f8a92b87381e
# ╟─321154fd-5092-4d93-ba82-1bcde10efcb5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
