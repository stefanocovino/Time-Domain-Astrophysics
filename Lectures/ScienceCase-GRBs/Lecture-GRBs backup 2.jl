### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 6a1315d1-9a6d-4ce0-b1c0-3fe22beb1ec2
begin
	using CairoMakie
	using CommonMark
	using Distributions
	using FITSIO
	using LsqFit
	using LombScargle
	using PlutoUI
	using ProgressLogging
	using Random
	using StatsBase
end

# ╔═╡ 7f04786b-a1c6-42ac-ac41-2cbe68878149
include("src/P4J.jl")

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
# Science Case: gamma-ray bursts
***


- GRBs are bright, rapid and short-duration gamma-ray signals appearing randomly in space and time. This is the light-curve of the first GRB ever detected, by one if the [Vela satellites](https://en.wikipedia.org/wiki/Vela_(satellite)):

$(LocalResource("Pics/vela.png"))

- After a few years to collect a sufficient number of events results were published and the GRB era began with this paper by [Klebesadel et al. (1973)](https://ui.adsabs.harvard.edu/abs/1973ApJ...182L..85K/abstract).

$(LocalResource("Pics/klebesadel.png"))

- And since the very beginning the main question was to determine where do GRBs come from?

. A substantial, although indirect, breakthrough came after the launch of the [Compton Gamma-Ray Observatory](https://en.wikipedia.org/wiki/Compton_Gamma_Ray_Observatory).

$(LocalResource("Pics/compton.png"))

- Compton-GRO was Launched in 1991. One of the main goals of the mission was to study GRBs. In a few years of operation thousands of them were detected.
"""

# ╔═╡ 72cf6fcf-2e2f-4dee-994b-ce2bfef51802
md"""

## Extreme variability
***

$(LocalResource("Pics/variability.png"))

- Rapid variability allows one to derive an estimate of the source size. It turns out to be very compact: of the order of tens of km.
"""

# ╔═╡ e1f2b5ca-49e4-4f8f-b69c-c3f621ffc068
cm"""
## The problem of GRB localization
***

- High-energy satellites can frequently detect GRBs thanks to wide-field of view capabilitied. 
    - However, at the times, the localization capabilites yielded an uncertainty of the order of a degree radius. 
    - In such an area typically you have billons of optical sources, Galactic and/or extra-Galactic.

- The problem can be addressed statistically.
    - First of all, the distribution of Galactic sources is far from isotropic in wky, as shown by the following all-sky view obtained by the [Fermi](https://en.wikipedia.org/wiki/Fermi_Gamma-ray_Space_Telescope) satellite.

$(LocalResource("Pics/fermi.png"))

- On the contrary, the distribution of galaxies farther than ``\sim 100`` Mpc is isotropic.

- Results based upon the observations of thousands of GRBs is impressive:

$(LocalResource("Pics/batse.png"))

- The distribution is highly isotropic! GRBs form a population of cosmological sources.
"""

# ╔═╡ 6bbb867f-4166-4ee9-9049-039dd46773d4
cm"""

## The [*Beppo*SAX](https://en.wikipedia.org/wiki/BeppoSAX) contribution
***

$(LocalResource("Pics/sax.jpg"))

- *Beppo*SAX was a Dutch-Italian satellite launched in 1996. 
    - It had a revolutionary capability at the time: the ability to repoint after an alert in a few hours.
    - After that a GRB was located by the high-energy instrument with an error below a few tens of arcmin, the satellte repointed and the target location could be observed by small field of vire but better resolution soft X-ray telescopes.

- This quickly brought to the identification of the first low-energy counterpart of a GRB: the afterglow era was beginning.

$(LocalResource("Pics/afterglow.png"))

<br>

- Only a few monnths later, for [GRB090726](https://www.mpe.mpg.de/~jcg/grb090726.html), the identified of the the optical counterpart was rapid enough to allow researchers to obtain a spectrum of the source that revealed to be a distant galaxy at a redshift ``z \sim 2.7``.

$(LocalResource("Pics/grb960726.png"))

- The energy output of these events is huge, comparable to a SN but in a timescale of seconds!

- The cosmological origin of GRBs and their “extreme” features have been one of the hottest problems for astrophysics in the 2000s.


"""

# ╔═╡ c1637e39-fcbe-434c-9197-6426d6743d1a
cm"""

## The need to be *Swift*
***

- The next step starts the golden age of GRB research, with the launch of the [Neil Gehrels *Swift*](https://en.wikipedia.org/wiki/Neil_Gehrels_Swift_Observatory) observatory.

$(LocalResource("Pics/swift.png"))

- Launched on 2004, Nov 20, and still operational.

- Swift has been designed to point very rapidly after an alert, typically with a time-scale of a minute.

"""

# ╔═╡ df340319-d931-4906-893f-bc0b0c043be7
cm"""

## So, what is a GRB?
***

- Let's summarize the results of decades of hard work in a plot:

$(LocalResource("Pics/grb.png"))

- There are at least two families of GRBs: the long and short duration, with possible further subdivisions.

$(LocalResource("Pics/progenitors.png"))

- We have observational evidence for both cases, although ther could be more complex scenarios.

"""

# ╔═╡ b83b11be-cbc2-45d1-aa27-92b90dec6506
cm"""

### Exercise: the proposed QPO for GRB211211A
***

- [GRB211211A](https://gcn.gsfc.nasa.gov/other/211211A.gcn3) was a bright GRB detected by several high-energy satellite that attracted a considerable attention since ``\sim 20``Hz QPOs have been reported in various phases of its prompt emission, e.g., [Xiao et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024ApJ...970....6X/abstract) and [Chirenti et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024ApJ...967...26C/abstract).

- This is a peculiarly difficult case: the GRB light-curve was short and the superposed oscillation could have lasted only a few cycles.
    - A non-stationary time-series with a transient behavior. 
   


- Within this exercize we analyse [Fermi-GBM](https://fermi.gsfc.nasa.gov/science/instruments/gbm.html) data for this event.
- Let's see the data:

"""

# ╔═╡ c813bd45-a811-4f84-bb5d-55175a895f96
# n0,n1,n2,n5,n9 and na, 8-200 keV
tte_file="glg_tte_na_bn211211549_v00.fit";

# ╔═╡ 46aa7e89-4b15-4149-8065-afb603ec1b06
begin
	
	# Open the FITS file
	f = FITS(tte_file)
	
	# Find corresponding CH1, CH2 from 'EBOUNDS' extension
	# In Julia, we access the HDU by name or index directly
	ebound_hdu = f["EBOUNDS"]
	emin = read(ebound_hdu, "E_MIN")
	emax = read(ebound_hdu, "E_MAX")
	
	# Get data from 'EVENTS' extension
	events_hdu = f["EVENTS"]
	time = read(events_hdu, "TIME")
	ch = read(events_hdu, "PHA")  # Or "DETCHANS", depending on your specific TTE file header
	
	# Close the file when done
	#close(f)
end

# ╔═╡ 2d9914d0-b7dc-4364-9374-3b3f668a8005
cm"- Let's set the energy range (keV):"

# ╔═╡ 37973ee8-e45f-48bc-834a-28a2deebab0c
begin
	# Define energy bounds
	E1 = 8.0
	E2 = 200.0
	
	# Create a boolean mask for the energy range
	# The '.' ensures the comparison happens for every element in the array
	index_for_EminEmax = (emin .>= E1) .& (emax .< E2)
	
	# Find the indices where the condition is true
	# findall returns a vector of indices, equivalent to np.where()[0]
	ch_arr = findall(index_for_EminEmax)
	
	# Get the first and last channel index
	ch1 = ch_arr[1]
	ch2 = ch_arr[end]
	
	# Filter the events based on the channel range
	ch_index = (ch .>= ch1) .& (ch .<= ch2)
	
	# Apply the filter to time and ch arrays
	time_sel = time[ch_index]
	ch_sel = ch[ch_index]
end;

# ╔═╡ 55e51bb5-211d-444a-be27-3c441212a310
cm"- And also the time interval (seconds before and after the GRB time):"

# ╔═╡ a4bb200d-9f3a-413c-a0b3-78a351f6cc42
begin
	# Define time window
	t1 = -0.5
	t2 = 1.0
	# Use the read_header() function instead of the .header field
	header_data = read_header(f[1])
	trigtime = header_data["TRIGTIME"]
	
	# Or in a single line:
	trigtime = read_header(f[1])["TRIGTIME"]
	
	# The rest of your logic remains the same:
	time_sel_sel = time_sel .- trigtime
	time_index = (time_sel_sel .>= t1) .& (time_sel_sel .<= t2)
	time_sel_T = time_sel_sel[time_index]
	ch_sel_T = ch_sel[time_index]
end

# ╔═╡ 8e06fdfb-6a39-4131-9a82-e178fc060a6c
cm"- Finally, let's choose the binsize (0.01s):"

# ╔═╡ 6be178b8-f38b-4f64-983b-f4ade0d40e55
begin
	binsize = 0.01
	
	# Create the bin ranges
	# Julia's range syntax is start:step:stop
	tbins = t1:binsize:t2
	
	# Calculate the histogram
	# fit(Histogram, data, bins) is the standard approach
	h = fit(Histogram, time_sel_T, tbins)
	
	# Access the counts and the bin edges
	histvalue = h.weights
	histbin = h.edges[1]
end

# ╔═╡ 0e7f2c00-f012-45dc-98ec-a36f8e653c2e
begin
	# Calculate bin centers for the x-axis
	# In Julia, histbin[1:end-1] is the same as histbin[:-1]
	plottime = histbin[1:end-1] .+ binsize/2.0
	plotrate = histvalue
	
	# Create the figure and axis
	fig = Figure()
	ax = Axis(fig[1, 1], 
	    xlabel = "Time (s, 0=GRB time)", 
	    ylabel = "Rate"
	)
	
	# Plot using the 'stairs' function (equivalent to plt.step)
	# Use step = :center to match your plottime calculation
	stairs!(ax, plottime, plotrate, step = :center, label = "Fermi-GBM")
	
	# Add the legend
	axislegend(ax)
	
	# Display the figure
	fig
end

# ╔═╡ 8d1d5f66-22bc-4dc7-83d0-f663b0ca2024
cm"""
- This episode of GRB activity is very short, and the pulse follows a typical profile with a verey fast rise and a slower, maybe exponential, decay.

- Let's now compute a LS periodogram for these data. In principle a simple DFT could be applied since the input light-curve is evenly sampled.
"""

# ╔═╡ 5e16d0ab-1f08-461a-92aa-8b135e125c0d
begin
	# Define frequency range
	# np.linspace(start, stop, num) -> range(start, stop, length=num)
	freq = range(1, 0.9/0.02, length=2000)

	# Force to 1D Floats
	final_t = convert(Vector{Float64}, vec(plottime))
	final_s = convert(Vector{Float64}, vec(plotrate))

	ls = lombscargle(final_t, final_s, frequencies=freq, normalization=:psd)
	
	# Access the power values
	pwrGBM = power(ls)
end

# ╔═╡ d4011a10-eac2-4f98-b57a-5ede937b4182
begin
	# 1. Create the Figure
	figls = Figure(size = (800, 600))
	
	# 2. Create the Axis with log-log scaling
	# xscale and yscale handle the logarithmic transformation
	axls = Axis(figls[1, 1], 
	    xlabel = "Frequency (Hz)", 
	    ylabel = "Power",
	    xscale = log10, 
	    yscale = log10,
	    title = "Lomb-Scargle Periodogram"
	)

	maskls = pwrGBM .> 0.
	
	# 3. Plot the data
	# lines! is the equivalent of plt.plot
	lines!(axls, freq[maskls], pwrGBM[maskls], label = "Fermi-GBM", color = :blue)
	
	# 4. Add the legend
	axislegend(axls)
	
	# 5. Show the result
	figls
end

# ╔═╡ b87886f6-865b-47a9-8318-7a08913db0c6
cm"""
- We see a peak at ``\sim`` 20Hz, but we also see that the periodogram is characterized by a red-noise behavior with power quickly growing toward the lowest frequencies.

- This kind of behavior is rather typical, and it is possible to model the periodogram in order to test the significance of the detected period as we did, e.g., for the science case devoted to AGN and blazars ([notebook](./open?path=Lectures/ScienceCase-AGNandBlazars/Lecture-AGN-and-Blazars.jl), [html](Lectures/ScienceCase-AGNandBlazars/Lecture-AGN-and-Blazars.html)).

- However, if we recall that a LS (or DFT) periodogram can be interpreted as a measure of the goodness of fit with sinusoids with different frequencies, it is likely that the red-noise behavior is (mainly) due to the modulation of the light-curve that is temptativey fit with low frequency sinusoids.

- We can try a different approach, i.e. fitting the pulse and removing it from the data and then analyze a "whitened" version of the light-curve. 
    - Please, pay attention that in temporal analysis there are no free lunches. Removing the results of a fit implies that the residuals now suffer also from the accumulated uncertainties of the fit too. Resulting data are noisier. 
    - In addition, unavoidably, after the fit removal, data are no more independent of each other since the operation introduces some correlation.
"""

# ╔═╡ 432e78d5-5d73-4833-bd14-f0f4c62e93fb
cm"""

> Quite beyond the purpose of this exercize, let's stress that any analysis of non-stationary time-series involves several additional complexities that must be properly evaluated. See, for instance, [Hübner et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJS..259...32H/abstract).

"""

# ╔═╡ 02b09798-fb21-4b2b-80a3-15f0006a9b7b
cm"""
- In spite of these difficulties, we now model the pulse and remove it from the data and analyse the results.

- In order to model the pulse, we adopt a formula defined in [Norris et al. (1996)](https://ui.adsabs.harvard.edu/abs/1996ApJ...459..393N/abstract):
"""

# ╔═╡ b58dd13b-3e1f-4d02-8cac-f24cbb3e45d3
function Norris_1996_fast(x, A, tmax, sigma_r, sigma_d, mu, B)
    @. ifelse(x < tmax, 
        B + A * exp(-(abs(x - tmax) / sigma_r)^mu), 
        B + A * exp(-(abs(x - tmax) / sigma_d)^mu)
    )
end

# ╔═╡ da981eda-aa53-4751-b498-a5f0e27d9167
cm"- And carry out a fit with this function:"

# ╔═╡ dececc3e-5676-416e-94dd-76e0a3b8c1f9
function norris_model(x, p)
    A, tmax, sigma_r, sigma_d, mu, B = p
    
    # Ensure sigma and mu are treated as positive to avoid DomainErrors
    # We use abs() here so the solver doesn't crash if it 'probes' negative values
    sr = abs(sigma_r)
    sd = abs(sigma_d)
    m  = abs(mu)

    return @. ifelse(x < tmax, 
        B + A * exp(-(abs(x - tmax) / sr)^m), 
        B + A * exp(-(abs(x - tmax) / sd)^m)
    )
end

# ╔═╡ e1d78c1a-61ac-4364-8e5b-d8a4d3a84e76
begin
	# Your initial guesses (p0)
	p0 = [30.0, -0.01, 0.06, 0.3, 0.087, 1.0]
	
	# Lower bounds for [A, tmax, sigma_r, sigma_d, mu, B]
	# We set sigmas and mu to a tiny positive number (1e-6)
	lb = [-Inf, -Inf, 1e-6, 1e-6, 1e-6, -Inf]
	ub = [Inf, Inf, Inf, Inf, Inf, Inf]

	fit_result = curve_fit(norris_model, plottime, plotrate, p0, lower=lb, upper=ub)
	
	# Extract optimized parameters (popt)
	popt0GBM = fit_result.param
	
	# Extract the covariance matrix (pcov)
	pcov0GBM = estimate_covar(fit_result)
	
	# Bonus: Get standard errors (the square root of the diagonal of pcov)
	errors = stderror(fit_result)
end

# ╔═╡ e55d34c8-8e97-4f04-b339-406068f862cf
begin
	# 1. Calculate the model and residuals
	# Using the dot . to broadcast the model over plottime
	fit_values = norris_model(plottime, popt0GBM)
	sub_GBM = plotrate .- fit_values
	
	# 2. Create the Figure
	figsub = Figure(size = (800, 800))
	
	# 3. Top Panel: Data vs Fit
	ax1sub = Axis(figsub[1, 1], 
	    ylabel = "Rate", 
	    title = "Norris 1996 Pulse Fit",
	    xticksvisible = false, xticklabelsvisible = false # Hide x-labels for the top plot
	)
	
	# scatter! or errorbars! can be used for the data
	scatter!(ax1sub, plottime, plotrate, color = (:black, 0.5), markersize = 8, label = "8-200 keV (GBM)")
	lines!(ax1sub, plottime, fit_values, color = :red, linewidth = 2, label = "pulse (GBM)")
	axislegend(ax1sub, position = :rt)
	
	# 4. Bottom Panel: Residuals
	ax2sub = Axis(figsub[2, 1], 
	    xlabel = "Time (s, 0=GRB time)", 
	    ylabel = "Residuals"
	)
	
	# Plot the pulse-subtracted data
	lines!(ax2sub, plottime, sub_GBM, color = (:blue, 0.5), label = "pulse-subtracted (GBM)")
	hlines!(ax2sub, [0], color = :black, linestyle = :dash) # Add a zero-line for reference
	axislegend(ax2sub, position = :rt)
	
	# 5. Fine-tune layout (removes gap between subplots)
	rowgap!(figsub.layout, 10)
	
	figsub
end

# ╔═╡ d73fb514-e125-439c-afcd-4b2b0c019541
cm"""
- The fit looks like satisfactory, and therefore we now analyze the subtracted curve.

- Let's than compute a LS periodogram and also conpute FAP levels by means of a boostrapping techniques.

- Please, be aware that for *independent and identically distributed* (IID) bootstrap the data points are randomly sampled with replacement. The IID bootstrap destroys not only the periodicity but also any time correlation or structure in the time series. This results in underestimation of the confidence bars. 
    - For data with serial correlations it is better to use moving block (MB) bootstrap. In MB bootstrap blocks of data of a given length are patched together to create a new time series. The block length is a parameter. Because light curves are irregularly sampled we set a block length in days rather than number of points. The ideal is to set the length so that it destroys the periodicity and preserves most of the serial correlation

- Bootstrap applied to time series is discussed in [Bühlmann (2002) - "Bootstraps for time series."](https://www.jstor.org/stable/3182810).

- Here, however, for didactic purposes, we follow a simple approach with a standard bootstrap technique.
"""

# ╔═╡ 87e5b3d7-bcde-46e5-88c0-6c649adf99a8
begin	
	# 2. Compute the periodogram on the residuals (sub_GBM)
	# We ensure sub_GBM is a 1D Vector of Floats to avoid type errors
	ls_sub = lombscargle(plottime, Float64.(vec(sub_GBM)), frequencies=freq, normalization=:psd)
	
	# 3. Extract the power
	pwrGBM_sub = power(ls_sub)
end

# ╔═╡ e6f7338c-5ae8-471e-bc89-8864e2b8f107
function bootstrap_fap(times, signal, freq, probs; n_iterations=1000)
    max_powers = zeros(n_iterations)
    
    for i in 1:n_iterations
        # Shuffle the signal to create "null hypothesis" data
        shuffled_signal = shuffle(signal)
        
        # Compute periodogram on shuffled data
        ls_boot = lombscargle(times, shuffled_signal, frequencies=freq, normalization=:psd)
        max_powers[i] = maximum(power(ls_boot))
    end
    
    # Find the power levels corresponding to the requested quantiles
    # If FAP = 0.1, we want the 90th percentile (0.9)
    return [quantile(max_powers, 1 - p) for p in probs]
end

# ╔═╡ d12af425-6a65-476d-915f-4548628f4d2d
begin
	# Define the probabilities (1 - confidence level)
	# 1-0.9997 is the 3-sigma level
	probabilities = [0.1, 0.05]
	
	# Run the bootstrap
	resprGBM_boot = bootstrap_fap(plottime, Float64.(vec(sub_GBM)), freq, probabilities, n_iterations=1000)
	
	println("FAP Levels (Bootstrap): ", resprGBM_boot)
end

# ╔═╡ 459347f8-da91-47d1-ab14-ec4af1c3c650
begin
	# 1. Create the Figure
	figls2 = Figure(size = (800, 600))
	
	# 2. Create the Axis with log-log scaling
	# xscale and yscale handle the logarithmic transformation
	axls2 = Axis(figls2[1, 1], 
	    xlabel = "Frequency (Hz)", 
	    ylabel = "Power",
	    #xscale = log10, 
	    #yscale = log10,
	    title = "Lomb-Scargle Periodogram"
	)

	maskls_sub = pwrGBM_sub .> 0.
	
	# 3. Plot the data
	# lines! is the equivalent of plt.plot
	lines!(axls2, freq[maskls_sub], pwrGBM_sub[maskls_sub], label = "Fermi-GBM", color = :blue)

	colors = [:green, :orange, :red]
	labels = ["90%", "95%", "3-sigma"]

	for (i, level) in enumerate(resprGBM_boot)
    	hlines!(axls2, [level], color = colors[i], linestyle = :dash, label = "$(labels[i]) FAP")
	end

	
	# 4. Add the legend
	axislegend(axls2)
	
	# 5. Show the result
	figls2
end

# ╔═╡ 4ed94405-688f-4f72-80d6-e09d632f907d
cm"""
- Given the assumptions, and quite interestingly, the periodicity at ``\sim 20``Hz does not show a strong significance. 

- However, again, please consider that this is only an exemplificatory analysis, and the considerations discussed in [Hübner et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJS..259...32H/abstract) definitely hold.
"""

# ╔═╡ 907b2aea-ff2b-426f-a336-215467adc25e
cm"""
- We have mentioned that red noise bahavior visible in the LS periodogram is likely due to the attempt of the algorithm to model the envelope of the light-curve evolution.

- We might wonder what the periodogram could become if we adopt a non-parametric method, where the power in the periodogram is a measure if the similarity without invoving any modeling of the data.

- Coverting a light-curve from the time-domain to the phase-domain, given a period, typically destroyes the correlation affecting the data.
"""

# ╔═╡ 3d81e826-977c-401c-aec1-0e63c1c8db8b
begin
	pg_qmi = Periodogram("qmieu")
	ffmin = 5
	ffmax = 1/0.02
	ffres = 1e-1
	fit!(pg_qmi, plottime, plotrate; dy = plotrate/10, fmin = ffmin, fmax = ffmax, resolution = ffres, n_local_optima=0);
end

# ╔═╡ edf2856d-5c77-456c-bde5-db3bc4de8ed5
begin
	figqmi = Figure(size = (800, 600))
	
	axqmi = Axis(figqmi[1, 1],
	    xlabel = "Frequency (Hz)",
	    ylabel = "Power",
	    title = "QMI periodogram",
	    xgridvisible = true,
	    ygridvisible = true
	)
	
	# Linea principale del periodogramma
	lines!(axqmi, pg_qmi.frequencies, pg_qmi.scores, color = :dodgerblue, linewidth = 2)	
	
	vlines!(axqmi, pg_qmi.best_frequency, ymin = 0, ymax = 1, color = :orange, linewidth = 8, alpha = 0.25)
	
	figqmi
end

# ╔═╡ f72f2ffc-e826-4102-b849-a95a3f305627
cm"""
- We now see that the periodogram looks flat and the main peak is at ``\sim 20Hz`` and higher harmonics.

- It is not, anyway, truly dominant the periodogram. This sugegsts that it should not be highly significant.
"""

# ╔═╡ 97cf17ef-694a-4902-9e55-4c2a9eaee802
cm"""


- Standard bootstrapping assumes every data point is independent. In GRB lightcurves or power-law noise scenarios, points close in time are often correlated. 
	- By picking a "block" of length block_length and moving it to the new time series, you preserve those correlations, making your significance tests (like FAP) much more robust against "fake" peaks caused by colored noise.

- We use a moving block bootstrap (or "overlapping block bootstrap") algorithm. 
	- It's particularly useful for time-series data where the noise might be correlated over short timescales, as it preserves the local structure of the data better than a simple shuffle.

"""

# ╔═╡ 7e025fee-891b-4b0e-bca5-1e4dd513b4d2
function block_bootstrap(mjd, mag, err; block_length=10.0, rseed=nothing)
    # Set the random seed if provided
    !isnothing(rseed) && Random.seed!(rseed)
    
    N = length(mjd)
    mjd_boot = zeros(Float64, N)
    mag_boot = zeros(Float64, N)
    err_boot = zeros(Float64, N)
    
    k = 1  # Julia is 1-indexed
    last_time = 0.0
    
    # Equivalent to finding max_idx: 
    # How many points (roughly) fit into a block_length from the end
    max_idx = 2
    for i in 2:N
        max_idx = i
        if mjd[end] - mjd[end - i + 1] > block_length
            break
        end
    end

    while k <= N
        # Pick a random starting index
        # rand(min:max) is inclusive on both ends
        idx_start = rand(1:(N - max_idx))
        
        idx_end = idx_start + 1
        for j in (idx_start + 1):N
            idx_end = j
            # Check if block exceeds length or if we're filling the last spots in N
            if mjd[idx_end] - mjd[idx_start] > block_length || k + (idx_end - idx_start) >= N
                break
            end
        end

        # Calculate the length of the slice to copy
        len = idx_end - idx_start
        range_target = k:(k + len - 1)
        range_source = idx_start:(idx_end - 1)

        # Assign values
        mjd_boot[range_target] .= mjd[range_source] .- mjd[idx_start] .+ last_time
        mag_boot[range_target] .= mag[range_source]
        err_boot[range_target] .= err[range_source]

        # Update tracking variables
        last_time = (mjd[idx_end] - mjd[idx_start]) + last_time
        k += len
    end

    return mjd_boot, mag_boot, err_boot
end

# ╔═╡ 8fd8d5b1-9201-4b72-bd8c-dd6fbc09fb54
begin
	niter = 200
	nmaxima = 20
	pbest_bootstrap = zeros(Float64, niter, nmaxima)
end;

# ╔═╡ f6a030e2-83df-433b-b2a8-8f02afae4d33
@progress for i in 1:niter
    tt_b, cc_b, ee_b = block_bootstrap(plottime, plotrate, plotrate/10, block_length=0.3)
    fit!(pg_qmi, plottime, plotrate; dy = plotrate/10, fmin = ffmin, fmax = ffmax, resolution = ffres, n_local_optima=0)
    pbest_bootstrap[i, :] = reverse(sort(pg_qmi.scores))[1:nmaxima]
end

# ╔═╡ d0bba050-d6a8-4b7f-a187-d488ff6ba7f6
# 1. Flatten the 2D array into a 1D vector (equivalent to .ravel())
# vec() creates a view, whereas pbest_bootstrap[:] creates a copy
data_flat = vec(pbest_bootstrap)

# 2. Fit the Generalized Extreme Value distribution
# fit(Type, data) automatically estimates parameters
d_gev = fit(GeneralizedExtremeValue, data_flat)

# 3. Access parameters (location μ, scale σ, shape ξ)
# Note: Julia's ξ (xi) corresponds to SciPy's -c (usually)
μ, σ, ξ = params(d_gev)

# 4. 'rv' in SciPy is the distribution object. 
# In Julia, d_gev is already that object.
# You can now call functions on it:
p_val = pdf(d_gev, 0.5)      # Probability Density Function
c_val = cdf(d_gev, 0.5)      # Cumulative Distribution Function
quantile_val95 = quantile(d_gev, 0.95) # 95th percentile
quantile_val90 = quantile(d_gev, 0.90) # 90th percentile

# ╔═╡ 2d596c28-74bc-4ff8-a030-fbac18dbceb0
cm"""
## Reference & Material

Material and papers related to the topics discussed in this lecture.

- [Piran (2004) - "The physics of gamma-ray bursts"](https://ui.adsabs.harvard.edu/abs/2004RvMP...76.1143P/abstract)
- [Hübner et al. (2022) - "Pitfalls of Periodograms: The Nonstationarity Bias in the Analysis of Quasiperiodic Oscillations"](https://ui.adsabs.harvard.edu/abs/2022ApJS..259...32H/abstract)
- [Huppenkothen et al. (2025) - "Searching for quasi-periodicities in short transients: the curious case of GRB 230307A"](https://ui.adsabs.harvard.edu/abs/2025arXiv250410153H/abstract)
"""

# ╔═╡ 3adba216-1df6-43f0-8000-05084c4d58c2
cm"""
## Further Material

Papers for examining more closely some of the discussed topics.

- [Klebesadel et al. (1973) - "Observations of Gamma-Ray Bursts of Cosmic Origin"](https://ui.adsabs.harvard.edu/abs/1973ApJ...182L..85K/abstract)
- [Xiao et al. (2024) - "The Peculiar Precursor of a Gamma-Ray Burst from a Binary Merger Involving a Magnetar"](https://ui.adsabs.harvard.edu/abs/2024ApJ...970....6X/abstract)
- [Chirenti et al. (2024) - "Evidence of a Strong 19.5 Hz Flux Oscillation in Swift BAT and Fermi GBM Gamma-Ray Data from GRB 211211A"](https://ui.adsabs.harvard.edu/abs/2024ApJ...967...26C/abstract)
- [Norris et al. (1996) - "Attributes of Pulses in Long Bright Gamma-Ray Bursts"](https://ui.adsabs.harvard.edu/abs/1996ApJ...459..393N/abstract)
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
    <td><a href="./open?path=Lectures/Lecture-NonParametricAnalysis/Lecture-NonParametricPeriodograms.jl">Lecture about non-parametric periodograms</a></td>
    <td><a href="./open?path=Lectures/Lecture-SingularSpectrumAnalysis/Lecture-SSA.jl">Lecture about singlular spectrum analysis"</a></td>
  </tr>
  <tr>
	<td>html</td>
    <td><a href="Lectures/Lecture-NonParametricAnalysis/Lecture-NonParametricPeriodograms.html">Lecture about non-parametric periodograms</a></td>
<td><a href="Lectures/Lecture-SingularSpectrumAnalysis/Lecture-SSA.html">Lecture about singlular spectrum analysis</a></td>
  </tr>
</table>


"""

# ╔═╡ 206474b8-0811-4785-8a71-acdcfd20b76c
md"""
**Copyright**

This notebook is provided as [Open Educational Resource](https://en.wikipedia.org/wiki/Open_educational_resources). Feel free to use the notebook for your own purposes. The text is licensed under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/), the code of the examples, unless obtained from other properly quoted sources, under the [MIT license](https://opensource.org/licenses/MIT). Please attribute the work as follows: *Stefano Covino, Time Domain Astrophysics - Lecture notes featuring computational examples, 2026*.
"""

# ╔═╡ Cell order:
# ╟─4d477519-c44f-434c-b7e0-8daaa5009358
# ╟─ec67de24-d88a-46fa-ae46-c0cd7b797adc
# ╠═6a1315d1-9a6d-4ce0-b1c0-3fe22beb1ec2
# ╠═7f04786b-a1c6-42ac-ac41-2cbe68878149
# ╟─3ddd0f61-79d0-473c-8fec-a0e0c3fc72bf
# ╟─5029a214-0841-40fb-b397-4a2e1047bfb7
# ╟─404060d3-23ec-400b-84cf-779e63b90293
# ╟─72cf6fcf-2e2f-4dee-994b-ce2bfef51802
# ╟─e1f2b5ca-49e4-4f8f-b69c-c3f621ffc068
# ╟─6bbb867f-4166-4ee9-9049-039dd46773d4
# ╟─c1637e39-fcbe-434c-9197-6426d6743d1a
# ╟─df340319-d931-4906-893f-bc0b0c043be7
# ╟─b83b11be-cbc2-45d1-aa27-92b90dec6506
# ╠═c813bd45-a811-4f84-bb5d-55175a895f96
# ╠═46aa7e89-4b15-4149-8065-afb603ec1b06
# ╟─2d9914d0-b7dc-4364-9374-3b3f668a8005
# ╠═37973ee8-e45f-48bc-834a-28a2deebab0c
# ╟─55e51bb5-211d-444a-be27-3c441212a310
# ╠═a4bb200d-9f3a-413c-a0b3-78a351f6cc42
# ╟─8e06fdfb-6a39-4131-9a82-e178fc060a6c
# ╠═6be178b8-f38b-4f64-983b-f4ade0d40e55
# ╠═0e7f2c00-f012-45dc-98ec-a36f8e653c2e
# ╟─8d1d5f66-22bc-4dc7-83d0-f663b0ca2024
# ╠═5e16d0ab-1f08-461a-92aa-8b135e125c0d
# ╠═d4011a10-eac2-4f98-b57a-5ede937b4182
# ╟─b87886f6-865b-47a9-8318-7a08913db0c6
# ╟─432e78d5-5d73-4833-bd14-f0f4c62e93fb
# ╟─02b09798-fb21-4b2b-80a3-15f0006a9b7b
# ╠═b58dd13b-3e1f-4d02-8cac-f24cbb3e45d3
# ╟─da981eda-aa53-4751-b498-a5f0e27d9167
# ╠═dececc3e-5676-416e-94dd-76e0a3b8c1f9
# ╠═e1d78c1a-61ac-4364-8e5b-d8a4d3a84e76
# ╠═e55d34c8-8e97-4f04-b339-406068f862cf
# ╟─d73fb514-e125-439c-afcd-4b2b0c019541
# ╠═87e5b3d7-bcde-46e5-88c0-6c649adf99a8
# ╠═e6f7338c-5ae8-471e-bc89-8864e2b8f107
# ╠═d12af425-6a65-476d-915f-4548628f4d2d
# ╠═459347f8-da91-47d1-ab14-ec4af1c3c650
# ╟─4ed94405-688f-4f72-80d6-e09d632f907d
# ╟─907b2aea-ff2b-426f-a336-215467adc25e
# ╠═3d81e826-977c-401c-aec1-0e63c1c8db8b
# ╠═edf2856d-5c77-456c-bde5-db3bc4de8ed5
# ╟─f72f2ffc-e826-4102-b849-a95a3f305627
# ╟─97cf17ef-694a-4902-9e55-4c2a9eaee802
# ╠═7e025fee-891b-4b0e-bca5-1e4dd513b4d2
# ╠═8fd8d5b1-9201-4b72-bd8c-dd6fbc09fb54
# ╠═f6a030e2-83df-433b-b2a8-8f02afae4d33
# ╠═d0bba050-d6a8-4b7f-a187-d488ff6ba7f6
# ╟─2d596c28-74bc-4ff8-a030-fbac18dbceb0
# ╟─3adba216-1df6-43f0-8000-05084c4d58c2
# ╟─b36fd613-95c8-44bf-876d-4eb345c26f08
# ╟─206474b8-0811-4785-8a71-acdcfd20b76c
