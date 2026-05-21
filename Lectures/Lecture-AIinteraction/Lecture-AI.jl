### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 2f82cbd2-4c34-11f1-901e-01c419d48099
begin
	using CairoMakie
	using CSV
	using DataFrames
	using DSP
	using FFTW
	using HDF5
	using LombScargle
	using PlutoUI
	using PlutoTeachingTools
	using Statistics
end

# ╔═╡ d5b55ac6-0f38-4603-aa72-bd964d0604a8
md"""
**What is this?**


*This notebook is part of a collection of `pluto` notebooks on various topics discussed during the Time Domain Astrophysics course delivered by Stefano Covino at the [Università dell'Insubria](https://www.uninsubria.eu/) in Como (Italy). Please direct questions and suggestions to [stefano.covino@inaf.it](mailto:stefano.covino@inaf.it).*
"""

# ╔═╡ efc84c48-a322-4c98-869d-2fe1e32151d5
TableOfContents()

# ╔═╡ 72168e90-1440-42ac-a859-834ec6dcc703
md"""
$(LocalResource("Pics/TimeDomainBanner.jpg"))
"""

# ╔═╡ 9b66dbc6-7dd9-4762-a791-7fe1e433961b
md"""
# Interacting with an AI
***
"""

# ╔═╡ fbe81ddf-836e-4e2e-a6e1-475257301c41
md"""
- Or, more correctly, with a [LLM](https://en.wikipedia.org/wiki/Large_language_model)!

$br

- This exercize has the goal to show a possible (very simple) interaction with a common LLM to solve time-domain astrophysics (exploratory) propblems.

- These algorithms can be very useful if correctly used, and it should become quickly clear that they can make things quicker, and often better, if we have specific knowlewdge about the problem we are addressing.

- However, worlking blindly and if an excess of trust interacting with these software can be a real waste of time, or occasionally worse.
"""

# ╔═╡ e1e5bb07-e4a8-4048-8878-9f7c49d57afd
md"""

> We are working now with *Claude Sonnet 4.6*, part of *Anthropic's Claude 4.6* family (which currently includes *Claude Opus 4.6* and *Claude Sonnet 4.6*). 

"""

# ╔═╡ 71bcd5bc-afb3-48cc-9185-41e44d433ff9
md"""

## [GW150914](https://en.wikipedia.org/wiki/First_observation_of_gravitational_waves): the first gravitational wave signal ever detected
***

- We have downloaded the strain recorded by the [LIGO](https://en.wikipedia.org/wiki/LIGO) observatory, in particular by the interferometer at the [Hanford site](https://en.wikipedia.org/wiki/Hanford_Site).
"""

# ╔═╡ fb8ce658-f20a-410e-aadc-b06e4eba6c51
question_box("Hi! I have the file with the GW strain for the first even ever GW detected.  How can I show a plot to see the GW signal? I work with Julia and CairoMakie.")

# ╔═╡ 7066d72c-6f84-4311-9a12-a05f897731e7
answer_box("Could you share the file so I can see its format (HDF5, CSV, txt…)? That'll let me write the exact Julia + CairoMakie code to load and plot it correctly.")

# ╔═╡ 2e9b26ef-78eb-4df9-a419-a099cf285719
md" - Once having uploaded the file, this is the system reply:"

# ╔═╡ 22d009b8-995f-4150-9efc-4a55d02b82ab
answer_box(md"""
		   Got everything I need. The structure is clear:

- `strain/Strain` — 131072 samples of h(t)
- `Xstart` = 1126259446 (GPS time)
- `Xspacing` = 0.000244140625 s → sample rate = **4096 Hz**
- Duration = 32 s

Here's the Julia + CairoMakie code:

```julia
using HDF5
using CairoMakie

# --- Load data ---
fid = h5open("GW150914_H1_32s.hdf5", "r")

strain  = read(fid["strain/Strain"])          # h(t) array
t_start = read(fid["strain/Strain"]).data     # use attrs instead
attrs   = HDF5.attrs(fid["strain/Strain"])
t0      = read(attrs["Xstart"])               # GPS start
dt      = read(attrs["Xspacing"])             # 1/4096 s
n       = length(strain)

close(fid)

# Time axis relative to GPS start (seconds)
t = range(0.0, step=dt, length=n)

# The merger is at GPS 1126259462.4 → ~16.4 s into the 32-s segment
t_merger = 1126259462.4 - t0   # ≈ 16.4 s

# --- Plot ---
fig = Figure(size = (900, 400), figure_padding = 20)

ax = Axis(fig[1, 1],
    xlabel = "Time from GPS $(t0) (s)",
    ylabel = "Strain  h(t)  [×10⁻²¹]",
    title  = "GW150914 — LIGO Hanford (H1)  |  Raw strain")

lines!(ax, collect(t), strain .* 1e21,
    color     = :royalblue,
    linewidth = 0.8)

# Mark the merger time
vlines!(ax, [t_merger],
    color     = (:red, 0.7),
    linewidth = 1.5,
    linestyle = :dash)

text!(ax, t_merger + 0.3, maximum(strain .* 1e21) * 0.85,
    text     = "merger",
    color    = :red,
    fontsize = 13)

save("GW150914_H1_strain.png", fig, px_per_unit = 2)
fig
```

A few notes:

- **`HDF5.attrs(dataset)`** gives you access to the metadata attributes — that's how you get `Xstart` and `Xspacing` without hardcoding them.
- The strain values are of order 10⁻²¹, so multiplying by `1e21` makes the y-axis readable.
- The merger event (the famous chirp) falls at GPS **1126259462.4**, which is about **16.4 s** into this 32-s segment — that's where you'll see the amplitude suddenly grow.
- The raw strain is quite noisy; if you want to see the chirp more clearly, the next step would be bandpass filtering (e.g. 35–350 Hz) using `DSP.jl`. Let me know if you'd like that too!
""")

# ╔═╡ 47ec8223-3159-433e-b91f-10a540edb5bb
md"- However, it turns out that the code is not correct, and I needed a few interactions to get to a working code"

# ╔═╡ 120af7fd-32e9-4936-950d-e55955380ad0
question_box("I have this error: FieldError: type Array has no field `data`, available fields: `ref`, `size`")

# ╔═╡ d6330ed7-457e-40e0-9f92-94bd3a75c46d
answer_box(md"""
That's a simple fix — `read()` already returns the array directly, so you don't need `.data`. The attribute reading is the right approach though. Here's the corrected loading block:

```julia
using HDF5
using CairoMakie

# --- Load data ---
fid = h5open("GW150914_H1_32s.hdf5", "r")

strain_ds = fid["strain/Strain"]
strain    = read(strain_ds)
attrs     = HDF5.attrs(strain_ds)
t0        = read(attrs["Xstart"])    # GPS start (Int)
dt        = read(attrs["Xspacing"])  # 0.000244140625 s

close(fid)

# Time axis relative to GPS start (seconds)
n = length(strain)
t = range(0.0, step=dt, length=n)

t_merger = 1126259462.4 - t0   # ≈ 16.4 s into the segment
```

The key fix: remove the intermediate `read(fid["strain/Strain"]).data` line entirely — just read the dataset reference **before** closing the file, then call `read()` on it. The rest of the plotting code stays exactly the same.
""")

# ╔═╡ 64ed9b01-8c08-4a4c-8d12-ad8b4ab3895f
question_box(md"""Again an error: MethodError: no method matching read(::Int64)
The function `read` exists, but no method is defined for this combination of argument types.

Closest candidates are:
  read(::Any, ::Type{Union{}}, ::Any...; kwargs...)
   @ Base __io.jl:238__
  read(::FilePathsBase.FileBuffer)
   @ FilePathsBase ~/.julia/packages/FilePathsBase/NV2We/src/__buffer.jl:55__
  read(::FilePathsBase.FileBuffer, ::Type{String})
   @ FilePathsBase ~/.julia/packages/FilePathsBase/NV2We/src/__buffer.jl:62__
  ...""")

# ╔═╡ b65614f7-60a5-4ff7-a8ad-a90e2c71bc99
answer_box(md"""
The issue is that `HDF5.attrs` returns a dict-like object where values are already plain Julia values, not datasets — so you shouldn't call `read()` on them. Just index directly:

```julia
fid = h5open("GW150914_H1_32s.hdf5", "r")

strain_ds = fid["strain/Strain"]
strain    = read(strain_ds)
attrs     = HDF5.attrs(strain_ds)
t0        = attrs["Xstart"]    # already an Int64, no read() needed
dt        = attrs["Xspacing"]  # already a Float64, no read() needed

close(fid)

n = length(strain)
t = range(0.0, step=dt, length=n)

t_merger = 1126259462.4 - t0
```

So the rule is:
- `read(dataset)` → to get the **array data** from a dataset
- `attrs["key"]` → to get **attribute values** (already decoded, no `read()` needed)
""")

# ╔═╡ c9f20dcb-3bbd-4fd0-a672-b86ce2f14ec5
md"- Now, everything seems fine, and I copy the code adding the right filename:"

# ╔═╡ a4e6c2f2-7c1d-4bb9-b1dc-92803bc8d908
begin
	fid = h5open("GW150914_H1_32s.hdf5", "r")
	
	strain_ds = fid["strain/Strain"]
	strain    = read(strain_ds)
	attrs     = HDF5.attrs(strain_ds)
	t0        = attrs["Xstart"]    # already an Int64, no read() needed
	dt        = attrs["Xspacing"]  # already a Float64, no read() needed
	
	close(fid)
end

# ╔═╡ eee3d157-7e61-4a9b-bd3c-c5ba96480962
begin
	n = length(strain)
	t = range(0.0, step=dt, length=n)
	
	t_merger = 1126259462.4 - t0
end;

# ╔═╡ d9ac1d30-f90a-4df7-a863-aa62f2154dce
begin
	figst = Figure(size = (900, 400), figure_padding = 20)
	
	axst = Axis(figst[1, 1],
	    xlabel = "Time from GPS $(t0) (s)",
	    ylabel = "Strain  h(t)  [×10⁻²¹]",
	    title  = "GW150914 — LIGO Hanford (H1)  |  Raw strain")
	
	lines!(axst, collect(t), strain .* 1e21,
	    color     = :royalblue,
	    linewidth = 0.8)
	
	# Mark the merger time
	vlines!(axst, [t_merger],
	    color     = (:red, 0.7),
	    linewidth = 1.5,
	    linestyle = :dash)
	
	text!(axst, t_merger + 0.3, maximum(strain .* 1e21) * 0.85,
	    text     = "merger",
	    color    = :red,
	    fontsize = 13)
	
	figst
end

# ╔═╡ 23650722-09c9-4df1-b277-fa388d3fe160
md"- The strain signal needs to be cleaned to show the chirp of the GW event. One pass is to reduce its frequency range:"

# ╔═╡ 0ff6019f-c988-489d-834a-ec3b58a76b9d
question_box("I would like, now, to bandpass the dataset.")

# ╔═╡ c6006389-9924-47c1-8134-49bd68cd840f
answer_box(md"""
Here's the code to bandpass filter the strain using `DSP.jl`:

```julia


# --- Bandpass filter (35–350 Hz) ---
fs      = 1.0 / dt          # sample rate = 4096 Hz
bpf     = Bandpass(35.0, 350.0; fs=fs)
design  = Butterworth(8)    # 8th-order Butterworth
filt_c  = digitalfilter(bpf, design)

# Zero-phase filtering (forward + backward pass, like scipy's filtfilt)
strain_filt = filtfilt(filt_c, strain)
```

Then just swap `strain` for `strain_filt` in the plot:

```julia
lines!(ax, collect(t), strain_filt .* 1e21,
    color     = :royalblue,
    linewidth = 1.0)
```

A few notes:

- **35–350 Hz** is the standard band used by the LIGO team for GW150914 — it cuts low-frequency seismic noise and high-frequency instrumental noise while preserving the chirp signal (which sweeps roughly 35→150 Hz at merger).
- **`filtfilt`** applies the filter forward and backward, giving zero phase distortion — this is important so the chirp isn't time-shifted.
- **8th-order Butterworth** gives a fairly sharp rolloff without ringing. You can try 4 if you see edge artifacts.
- The edges of the time series (~0.5 s on each side) will have filter transients — you can crop them with `t[1025:end-1024]` and `strain_filt[1025:end-1024]` if needed.
"""
)

# ╔═╡ df96e17c-d8cd-4d93-b74a-7b8bef848e9f
md"- Unfortunately, the code is again not running correctly."

# ╔═╡ 18d6ffd2-5d52-44c8-86dd-b07739ca7f69
question_box(md"""
I have an error: MethodError: no method matching DSP.Filters.Bandpass(::Float64, ::Float64; fs::Float64)
This method does not support all of the given keyword arguments (and may not support any).

Closest candidates are:
  DSP.Filters.Bandpass(::T, ::V) where {T<:Real, V<:Real} got unsupported keyword argument "fs"
   @ DSP ~/.julia/packages/DSP/b3mah/src/Filters/__design.jl:281__
""")

# ╔═╡ 61101ce8-99ec-4736-a328-3234ac93d716
answer_box(md"""
Easy fix — in DSP.jl, `fs` is passed to `digitalfilter`, not to `Bandpass`:

```julia
fs     = 1.0 / dt          # 4096 Hz

bpf    = Bandpass(35.0, 350.0)                        # no fs here
design = Butterworth(8)
filt_c = digitalfilter(bpf, design; fs=fs)            # fs goes here

strain_filt = filtfilt(filt_c, strain)
```

That's the correct DSP.jl API — the filter type (`Bandpass`, `Lowpass`, etc.) just holds the cutoff frequencies, while `digitalfilter` handles the frequency normalisation using `fs`.
""")

# ╔═╡ 081d09ac-692a-46c1-a87e-3350c3c2baad
begin
	fs     = 1.0 / dt          # 4096 Hz
	
	bpf    = Bandpass(35.0, 350.0)                        # no fs here
	design = Butterworth(8)
	filt_c = digitalfilter(bpf, design; fs=fs)            # fs goes here
	
	strain_filt = filtfilt(filt_c, strain)
end;

# ╔═╡ 57ba7408-bd56-4227-a58d-effa843d64e4
md"- In the plot below I added limits on the `y`-axis to better center the bandpassed strain data."

# ╔═╡ d267408c-246b-4f92-bbbb-410157356a1d
begin
	figst2 = Figure(size = (900, 400), figure_padding = 20)
	
	axst2 = Axis(figst2[1, 1],
	    xlabel = "Time from GPS $(t0) (s)",
	    ylabel = "Strain  h(t)  [×10⁻²¹]",
	    title  = "GW150914 — LIGO Hanford (H1)  |  Raw strain")
	
	lines!(axst2, collect(t), strain_filt .* 1e21,
	    color     = :royalblue,
	    linewidth = 0.8)
	
	# Mark the merger time
	vlines!(axst2, [t_merger],
	    color     = (:red, 0.7),
	    linewidth = 1.5,
	    linestyle = :dash)
	
	text!(axst2, t_merger + 0.3, maximum(strain_filt .* 1e21) * 0.85,
	    text     = "merger",
	    color    = :red,
	    fontsize = 13)

	ylims!(-7.5,7.5)
	
	figst2
end

# ╔═╡ 5a0cc1ac-1ecb-43fa-a687-048ae9f6cc56
md"- Now, I propose to the AI to better clean the data by a *whitening* procedure, i.e. removing (oart of) the noise to make the signal more evident."

# ╔═╡ cd3f0947-78b7-4593-b7b7-ef82d7a15273
question_box("It works, yet the signal is still not visible. Might we tray to whiten it?")

# ╔═╡ 5fa720ff-19da-4e22-8972-eb87e3c3d382
answer_box(md"""
Yes, whitening is exactly what the LIGO team does to make the chirp visible — it flattens the noise power spectral density so no frequency band dominates. Here's the code:

```julia
using FFTW

# --- Whiten the bandpassed strain ---
function whiten(data, dt)
    n     = length(data)
    freqs = rfftfreq(n, 1.0/dt)       # frequency array
    
    # FFT of the data
    data_fft = rfft(data)
    
    # Estimate the ASD (amplitude spectral density) by smoothing |FFT|
    asd = abs.(data_fft)
    
    # Smooth the ASD with a running median to avoid dividing by noise spikes
    win = 64   # smoothing window (bins)
    asd_smooth = similar(asd)
    for i in 1:length(asd)
        lo = max(1, i - win÷2)
        hi = min(length(asd), i + win÷2)
        asd_smooth[i] = median(asd[lo:hi])
    end
    
    # Divide FFT by smoothed ASD → flat spectrum
    white_fft = data_fft ./ (asd_smooth .+ 1e-50)   # tiny epsilon avoids /0
    
    # Back to time domain
    white = irfft(white_fft, n)
    
    # Normalise
    white ./= std(white)
    return white
end

strain_white = whiten(strain_filt, dt)
```

Then plot it:

```julia
lines!(ax, collect(t), strain_white,
    color     = :royalblue,
    linewidth = 0.8)
```

A few notes:

- **Whiten first on `strain_filt`** (the bandpassed signal) so the very low and very high frequency noise doesn't dominate the ASD estimate.
- The y-axis after whitening is in units of **σ** (standard deviations of the noise) — the merger should show up as an excursion of several σ around t ≈ 16.4 s.
- The smoothing window `win = 64` is a reasonable starting point; increase it (e.g. 128 or 256) for a smoother ASD estimate and stronger whitening.
- Edge transients will be more visible after whitening — cropping ~0.5–1 s from each end is a good idea:

```julia
crop = round(Int, 1.0/dt)   # 1 second worth of samples
t_crop      = collect(t)[crop:end-crop]
white_crop  = strain_white[crop:end-crop]

lines!(ax, t_crop, white_crop, color = :royalblue, linewidth = 0.8)
```
""")

# ╔═╡ 05a82acf-b130-4642-a6e3-bd910bad9b01
# --- Whiten the bandpassed strain ---
function whiten(data, dt)
    n     = length(data)
    freqs = rfftfreq(n, 1.0/dt)       # frequency array
    
    # FFT of the data
    data_fft = rfft(data)
    
    # Estimate the ASD (amplitude spectral density) by smoothing |FFT|
    asd = abs.(data_fft)
    
    # Smooth the ASD with a running median to avoid dividing by noise spikes
    win = 64   # smoothing window (bins)
    asd_smooth = similar(asd)
    for i in 1:length(asd)
        lo = max(1, i - win÷2)
        hi = min(length(asd), i + win÷2)
        asd_smooth[i] = median(asd[lo:hi])
    end
    
    # Divide FFT by smoothed ASD → flat spectrum
    white_fft = data_fft ./ (asd_smooth .+ 1e-50)   # tiny epsilon avoids /0
    
    # Back to time domain
    white = irfft(white_fft, n)
    
    # Normalise
    white ./= std(white)
    return white
end

# ╔═╡ 208b8cb2-b34f-44ba-97c1-83db0da6fb33
strain_white = whiten(strain_filt, dt);

# ╔═╡ 3f4a1fe0-80cb-46e9-acc9-d800dea51285
md"- Here, again, I fine tuned the plot axis to magnify the strain data."

# ╔═╡ 589f5ed7-084d-48ad-9162-03df4c1a9609
begin
	figst3 = Figure(size = (900, 400), figure_padding = 20)
	
	axst3 = Axis(figst3[1, 1],
	    xlabel = "Time from GPS $(t0) (s)",
	    ylabel = "Strain  h(t)  [×10⁻²¹]",
	    title  = "GW150914 — LIGO Hanford (H1)  |  Raw strain")
	
	crop = round(Int, 1.0/dt)   # 1 second worth of samples
	t_crop      = collect(t)[crop:end-crop]
	white_crop  = strain_white[crop:end-crop]

	lines!(axst3, t_crop, white_crop, color = :royalblue, linewidth = 0.8)
	
	# Mark the merger time
	vlines!(axst3, [t_merger],
	    color     = (:red, 0.7),
	    linewidth = 1.5,
	    linestyle = :dash)
	
	text!(axst3, t_merger + 0.3, maximum(strain_white) * 0.85,
	    text     = "merger",
	    color    = :red,
	    fontsize = 13)
	xlims!(16.0,16.8)
	ylims!(-5,5)
	figst3
end

# ╔═╡ a985cfbd-dc02-4f8f-af40-f55ef584a055
md"""
- This plot above, while truly if interest, does not fully show how the event evolved. 

- Probably, a wavelet, i.e. time-frequency-power, plot could be more expressive.
"""

# ╔═╡ 164f7527-6434-4df1-8cab-d25f6363ecb1
question_box("This is good. Do you think you could write a code to have a plot of the chirp by, e.g., a wavelet?")

# ╔═╡ 8ac36e25-05ac-4477-958e-a35c9fb4646a
answer_box(md"""
Yes! A wavelet spectrogram (scalogram) is exactly how LIGO visualises the chirp — you can see the frequency sweeping upward in time right at the merger. The standard choice is the **Morlet wavelet**, which gives good time-frequency resolution.

DSP.jl doesn't have a built-in CWT, but it's straightforward to implement in the frequency domain:

```julia
using FFTW
using CairoMakie
using Statistics

# --- Morlet CWT in the frequency domain ---
function morlet_cwt(data, dt, freqs; omega0=6.0)
    n       = length(data)
    data_fft = rfft(data)
    ang_freqs = rfftfreq(n, 1.0/dt) .* 2π   # angular frequencies

    scalogram = zeros(Float64, length(freqs), n)

    for (i, f) in enumerate(freqs)
        # Scale corresponding to this frequency
        scale = omega0 / (2π * f)

        # Morlet wavelet in frequency domain
        psi_hat = @. (π^(-0.25)) * exp(-0.5 * (scale * ang_freqs - omega0)^2)
        psi_hat .*= sqrt(2π * scale / dt)

        # Convolve via multiplication in frequency domain, then iFFT
        col = irfft(data_fft .* psi_hat, n)
        scalogram[i, :] = abs.(col)
    end
    return scalogram
end

# --- Frequency axis (log-spaced, 30–400 Hz) ---
freqs = exp10.(range(log10(30.0), log10(400.0), length=200))

# --- Compute scalogram on whitened+bandpassed strain ---
# Use a cropped version to avoid edge effects
crop  = round(Int, 1.0/dt)
t_crop     = collect(t)[crop:end-crop]
white_crop = strain_white[crop:end-crop]

println("Computing CWT... (this may take a few seconds)")
scalogram = morlet_cwt(white_crop, dt, freqs)

# --- Plot ---
fig = Figure(size=(1000, 500), figure_padding=20)

ax = Axis(fig[1, 1],
    xlabel      = "Time from GPS $(t0) (s)",
    ylabel      = "Frequency (Hz)",
    title       = "GW150914 — H1 Wavelet Spectrogram (Morlet)",
    yscale      = log10,
    yticks      = ([32, 64, 128, 256], ["32", "64", "128", "256"]),
    yminorticksvisible = true)

# Clip dynamic range for better contrast (like LIGO plots)
slog = log10.(scalogram .+ 1e-10)
vmin, vmax = quantile(vec(slog), 0.05), quantile(vec(slog), 0.9995)

hm = heatmap!(ax, t_crop, freqs, slog',
    colorrange  = (vmin, vmax),
    colormap    = :inferno)

# Merger line
vlines!(ax, [t_merger],
    color     = (:white, 0.6),
    linewidth = 1.5,
    linestyle = :dash)

text!(ax, t_merger + 0.2, 300.0,
    text     = "merger",
    color    = :white,
    fontsize = 13)

Colorbar(fig[1, 2], hm, label="log₁₀ |CWT|")

save("GW150914_H1_wavelet.png", fig, px_per_unit=2)
fig
```

A few notes:

- **Morlet with ω₀=6** is the standard choice for GW data — it balances time and frequency resolution well. Increasing ω₀ gives better frequency resolution but smears in time.
- The frequency axis is **log-spaced** (30–400 Hz) because the chirp sweeps logarithmically and this matches LIGO's standard plots.
- The **`quantile` clipping** of the dynamic range (`0.05`–`0.9995`) is key for contrast — without it the loud merger drowns out everything else.
- The chirp should appear as a bright curved track sweeping from ~35 Hz up to ~150 Hz, peaking at `t_merger`. If the track isn't visible, try tightening `vmax` to `quantile(..., 0.999)` or even `0.99`.
- The CWT loop over 200 frequencies on 131072 samples takes a handful of seconds — perfectly fine for a 32-s segment.
""")

# ╔═╡ 9d73d318-6cdf-49e5-b997-b6ac204ec19d
# --- Morlet CWT in the frequency domain ---
function morlet_cwt(data, dt, freqs; omega0=6.0)
    n       = length(data)
    data_fft = rfft(data)
    ang_freqs = rfftfreq(n, 1.0/dt) .* 2π   # angular frequencies

    scalogram = zeros(Float64, length(freqs), n)

    for (i, f) in enumerate(freqs)
        # Scale corresponding to this frequency
        scale = omega0 / (2π * f)

        # Morlet wavelet in frequency domain
        psi_hat = @. (π^(-0.25)) * exp(-0.5 * (scale * ang_freqs - omega0)^2)
        psi_hat .*= sqrt(2π * scale / dt)

        # Convolve via multiplication in frequency domain, then iFFT
        col = irfft(data_fft .* psi_hat, n)
        scalogram[i, :] = abs.(col)
    end
    return scalogram
end

# ╔═╡ 3cf92273-1a38-413a-9ed9-87572007aa54
begin
	# --- Frequency axis (log-spaced, 30–400 Hz) ---
	freqs = exp10.(range(log10(30.0), log10(400.0), length=200))
	
	# --- Compute scalogram on whitened+bandpassed strain ---
	# Use a cropped version to avoid edge effects
	cropsc  = round(Int, 1.0/dt)
	t_cropsc     = collect(t)[cropsc:end-cropsc]
	white_cropsc = strain_white[cropsc:end-cropsc]
	
	println("Computing CWT... (this may take a few seconds)")
	scalogram = morlet_cwt(white_cropsc, dt, freqs)
	
end;

# ╔═╡ 56f6c7d9-086e-41cb-9618-69ac9d103fc5
md"- Again, I needed to improve the plot changing the axis limits, etc."

# ╔═╡ ac23e27e-6ad4-405e-b360-c13e674c7528
begin
	# --- Plot ---
	figsc = Figure(size=(1000, 500), figure_padding=20)
	
	axsc = Axis(figsc[1, 1],
	    xlabel      = "Time from GPS $(t0) (s)",
	    ylabel      = "Frequency (Hz)",
	    title       = "GW150914 — H1 Wavelet Spectrogram (Morlet)",
	    yscale      = log10,
	    yticks      = ([32, 64, 128, 256], ["32", "64", "128", "256"]),
	    yminorticksvisible = true)
	
	# Clip dynamic range for better contrast (like LIGO plots)
	slogsc = log10.(scalogram .+ 1e-10)
	vminsc, vmaxsc = quantile(vec(slogsc), 0.05), quantile(vec(slogsc), 0.9995)
	
	hmsc = heatmap!(axsc, t_cropsc, freqs, slogsc',
	    colorrange  = (vminsc, vmaxsc),
	    colormap    = :inferno)
	
	# Merger line
	vlines!(axsc, [t_merger],
	    color     = (:white, 0.6),
	    linewidth = 1.5,
	    linestyle = :dash)
	
	text!(axsc, t_merger + 0.05, 300.0,
	    text     = "merger",
	    color    = :white,
	    fontsize = 13)

	xlims!(16.2,16.6)
	
	Colorbar(figsc[1, 2], hmsc, label="log₁₀ |CWT|")
	
	figsc
end

# ╔═╡ da665df8-9724-412a-acd9-9bd7e7994e7f
md"""
- Now, the *chirp* is clear!
- The plot could be definitely improved, by smoothing it, etc. in order to make it publication ready.
- Even a better subtraction of the background could help.
- However, for this exercize we can be satisfied.
"""

# ╔═╡ 8e4d77a6-6f77-4621-a7c7-f5b444fd2411
md"""
- As commented above, this is only an exploratory exercize. Yet, it could show how a blind interaction, for avdanced topics, with an AI can be improductive.
- On the other hand, with a sufficient knowledge of the topic, one can drive the AI with the right questions and getting truly interesting results.
"""

# ╔═╡ 845b8722-0d66-42c9-b1b8-eddcfebce6a6
md"""
- A possible step forward could be to really analyse the strain data in order to identify the event, rather than locating and plotting it.
- This is possible, but definitely more difficult, althoug the many tutorials about this event available on the web could allow the AI to be properly trained.
"""

# ╔═╡ 791f7bc6-95e5-4beb-a908-563947d1609d
md"""
## [Delta Cephei](https://en.wikipedia.org/wiki/Delta_Cephei): one of the most important variable stars
***
"""

# ╔═╡ 1555efa5-cf2c-449d-819a-9d7258c7c24b
md"""
- We have downloaded the photoemtric data from the American Association of Variable Star Observers, [AAVSO](https://www.aavso.org/), site.
- It is a standard [CSV](https://it.wikipedia.org/wiki/Comma-separated_values), comma-separeted, text file.
"""

# ╔═╡ 1998e1eb-6c63-48de-b33b-aa248151d166
question_box("Hi, I have a CSV file with photometry of the delta cephei star. How can I read it? I work with julia and CairoMakie.")

# ╔═╡ 1a530ba0-c19c-4056-8408-b75956810571
md"- The AI complains I have not uploaded the file, so I do it and get a more specific response."

# ╔═╡ dc85b4d0-5dcd-4fa6-ba93-be27991145a0
answer_box(md"""

**The file structure:** A comma-separated AAVSO file with 24 columns and ~69,400 rows. The useful ones are `JD`, `Magnitude`, `Uncertainty`, and `Band`. Many `Uncertainty` values are missing (especially older visual observations), which is expected.

**Reading it in Julia:**
```julia
using CSV, DataFrames

df = CSV.read("DeltaCep.csv", DataFrame,
    select = ["JD", "Magnitude", "Uncertainty", "Band"],
    types  = Dict("JD" => Float64, "Magnitude" => Float64, "Uncertainty" => Float64),
    silencewarnings = true
)
dropmissing!(df, :Magnitude)
```

**Plotting the V-band light curve with CairoMakie:**
```julia
using CairoMakie

df_V = filter(row -> row.Band == "V" && !ismissing(row.Uncertainty), df)

fig = Figure(size = (1000, 400))
ax  = Axis(fig[1,1],
    xlabel = "Julian Date", ylabel = "Magnitude",
    title = "δ Cephei — V band", yreversed = true)

scatter!(ax, df_V.JD, df_V.Magnitude,
    color = :steelblue, markersize = 3, alpha = 0.5)

save("delta_cep_lc.png", fig)
```

Shall I move on to the phase-folded light curve?
""")

# ╔═╡ b7f2b3ea-0643-4385-bdc0-8fa9263ba7a8
begin
	df = CSV.read("DeltaCep.csv", DataFrame,
	select = ["JD", "Magnitude", "Uncertainty", "Band"],
	types = Dict("JD" => Float64, "Magnitude" => Float64, "Uncertainty" => Float64), silencewarnings = true)
	dropmissing!(df, :Magnitude)
end

# ╔═╡ 805bf93e-4723-4556-9fca-1931e8611f0e
begin
	df_V = filter(row -> row.Band == "V" && !ismissing(row.Uncertainty), df)

	figlc = Figure(size = (1000, 400))

	axlc = Axis(figlc[1,1],
		xlabel = "Julian Date", ylabel = "Magnitude",
		title = "δ Cephei — V band", yreversed = true)

	scatter!(axlc, df_V.JD, df_V.Magnitude, color = :steelblue, markersize = 3, alpha = 0.5)

	figlc
end

# ╔═╡ 7db81808-236e-4b93-af75-d2db036b7873
question_box("Yes, I'd like to see the phase-folded light-curve")

# ╔═╡ 887a46fe-c37f-452a-8504-18ce0164b29c
answer_box(md"""
Sure! Here's the code for the phase-folded light curve. The standard ephemeris for δ Cephei uses a period of **5.366341 days** and a well-known reference epoch.

```julia
using CSV, DataFrames, CairoMakie

# --- Read data ---
df = CSV.read("DeltaCep.csv", DataFrame,
    select = ["JD", "Magnitude", "Uncertainty", "Band"],
    types  = Dict("JD" => Float64, "Magnitude" => Float64, "Uncertainty" => Float64),
    silencewarnings = true
)
dropmissing!(df, :Magnitude)

# --- Filter V-band with uncertainty ---
df_V = filter(row -> row.Band == "V" && !ismissing(row.Uncertainty), df)

# --- Phase folding ---
period  = 5.366341          # days
t0      = 2436084.0         # reference epoch of maximum light (HJD)

df_V.phase = mod.((df_V.JD .- t0) ./ period, 1.0)

# Sort by phase for a clean look
sort!(df_V, :phase)

# --- Plot ---
fig = Figure(size = (900, 450), figure_padding = 30)
ax  = Axis(fig[1,1],
    xlabel     = "Phase",
    ylabel     = "Magnitude",
    title      = "δ Cephei — Phase-folded light curve (V band)",
    yreversed  = true,
    xticks     = 0:0.1:1.0,
)

scatter!(ax, df_V.phase, df_V.Magnitude,
    color      = :steelblue,
    markersize = 3,
    alpha      = 0.4,
    label      = "V-band observations"
)

# Duplicate points shifted by +1 to show continuity across the cycle boundary
scatter!(ax, df_V.phase .+ 1.0, df_V.Magnitude,
    color      = :steelblue,
    markersize = 3,
    alpha      = 0.15,
)

xlims!(ax, 0.0, 2.0)
axislegend(ax, position = :rb)

save("delta_cep_phased.png", fig)
fig
```

A few notes on the choices made here:

- **`t0 = 2436084.0`** is a standard epoch of maximum light for δ Cep. Phase 0 corresponds to brightness maximum (the fast rise).
- **Duplicating points shifted by +1** in a faded style is a common trick to show the cycle wraps around continuously — it makes the steep rise near phase 0 much easier to see.
- The **asymmetric shape** should be very clear: a rapid brightening over ~0.1–0.2 in phase, followed by a slow decline — the classic Cepheid sawtooth.

Would you like to also add a **binned mean curve** on top of the scatter, to highlight the average pulse shape more cleanly?
""")

# ╔═╡ 98a5702a-8ab1-41a7-b0f3-1e7ef035b786
begin
	# --- Phase folding ---
	period_dc = 5.366341 # days
	t0_dc = 2436084.0 # reference epoch of maximum light (HJD)
	
	df_V.phase = mod.((df_V.JD .- t0) ./ period_dc, 1.0)
	
	# Sort by phase for a clean look
	sort!(df_V, :phase)
end

# ╔═╡ db1e7d54-a93b-4185-955e-9a2cd29a4f73
begin
	figpp = Figure(size = (900, 450), figure_padding = 30)

	axpp = Axis(figpp[1,1],
		xlabel = "Phase",
		ylabel = "Magnitude",
		title = "δ Cephei — Phase-folded light curve (V band)",
		yreversed = true,
		xticks = 0:0.1:1.0,
		)
	
	scatter!(axpp, df_V.phase, df_V.Magnitude, color = :steelblue, markersize = 3, alpha = 0.4, label = "V-band observations")
	
	# Duplicate points shifted by +1 to show continuity across the cycle boundary
	scatter!(axpp, df_V.phase .+ 1.0, df_V.Magnitude, color = :steelblue, markersize = 3, alpha = 0.15,)
	
	xlims!(axpp, 0.0, 2.0)
	axislegend(axpp, position = :rb)

	figpp
end

# ╔═╡ f817a31a-26f5-4fbb-ae4a-bc2adf3c9b1a
question_box("No, but I would like to find the best period by myself")

# ╔═╡ 40a2c488-c0b6-491a-8061-bed6a536ba93
answer_box(md"""
Great idea — this is a classic period-finding problem, perfect for a Cepheid. The natural tool here is the **Lomb-Scargle periodogram**, since the data are unevenly sampled and span 145 years.

Here's the code:

```julia
using CSV, DataFrames, LombScargle, CairoMakie

# --- Read and filter ---
df = CSV.read("DeltaCep.csv", DataFrame,
    select = ["JD", "Magnitude", "Uncertainty", "Band"],
    types  = Dict("JD" => Float64, "Magnitude" => Float64, "Uncertainty" => Float64),
    silencewarnings = true
)
dropmissing!(df, :Magnitude)

df_V = filter(row -> row.Band == "V" && !ismissing(row.Uncertainty), df)
sort!(df_V, :JD)

# --- Lomb-Scargle ---
t   = df_V.JD
mag = df_V.Magnitude
err = df_V.Uncertainty

plan = LombScargle.plan(t, mag, err,
    minimum_frequency = 1.0 / 6.0,   # periods shorter than 6 days
    maximum_frequency = 1.0 / 4.0,   # periods longer than 4 days
    samples_per_peak  = 10
)

pgram = lombscargle(plan)

# --- Find best period ---
best_freq  = freqat(pgram, argmaxpower(pgram))
best_period = 1.0 / best_freq
println("Best period: ", round(best_period, digits=6), " days")

# --- Plot periodogram ---
fig = Figure(size = (900, 400), figure_padding = 30)
ax  = Axis(fig[1,1],
    xlabel = "Period (days)",
    ylabel = "Lomb-Scargle Power",
    title  = "δ Cephei — Lomb-Scargle Periodogram (V band)"
)

periods = 1.0 ./ pgram.freq
lines!(ax, periods, pgram.power, color = :steelblue, linewidth = 0.8)

# Mark the best period
vlines!(ax, best_period, color = :crimson, linewidth = 1.5, linestyle = :dash,
    label = "Best period: $(round(best_period, digits=4)) d")

axislegend(ax, position = :rt)
save("delta_cep_periodogram.png", fig)
fig
```

A couple of things worth noting:

- I've **restricted the frequency range** to periods between 4 and 6 days, since we already have a good prior on δ Cep's period (~5.37 d). You can widen it if you want a blind search.
- **`samples_per_peak = 10`** gives good frequency resolution without being too slow. Increase it to 20–50 for finer sampling near the peak.
- The package is `LombScargle.jl` — if you don't have it yet, just run `] add LombScargle` in the REPL.

Once you have the best period, you can plug it straight into the phase-folding code from before. How close does it come to the canonical 5.366341 days?
""")

# ╔═╡ 954f2bdc-973d-4bc0-8402-0107fb3382b7
md"- But, as it happens frequently, the code is not correct. I need to interact a bit to have a correct (i.e. running) answer."

# ╔═╡ 26ddf772-d989-46e2-952a-ccaa9661921f
question_box("I have this error: MethodError: no method matching plan(::Vector{Float64}, ::Vector{Float64}, ::Vector{Union{Missing, Float64}}; minimum_frequency::Float64, maximum_frequency::Float64, samples_per_peak::Int64)
The function `plan` exists, but no method is defined for this combination of argument types.

Closest candidates are:
  plan(::AbstractVector{<:Real}, ::AbstractVector{<:Real}; kwargs...)
   @ LombScargle ~/.julia/packages/LombScargle/kygQJ/src/__planning.jl:130__
  plan(::AbstractVector{<:Real}, ::AbstractVector{<:Real}, ::AbstractVector{<:Real}; kwargs...)
   @ LombScargle ~/.julia/packages/LombScargle/kygQJ/src/__planning.jl:134__
  plan(::AbstractVector{<:Real}, ::AbstractVector{<:Measurements.Measurement}; kwargs...)
   @ LombScargle ~/.julia/packages/LombScargle/kygQJ/src/__planning.jl:144__")

# ╔═╡ e88c6a1c-7d31-43d8-a35e-bfa0104a93f4
answer_box(md"""
The issue is that `err` still has type `Vector{Union{Missing, Float64}}` even after `dropmissing!` on `:Magnitude` — because we never dropped rows with missing *uncertainty*. The fix is to either drop those rows too, or simply not pass the errors to `plan` (which is fine for a period search).

**Option 1 — drop missing uncertainties (recommended, keeps error weighting):**

```julia
dropmissing!(df_V, :Uncertainty)

t   = df_V.JD
mag = df_V.Magnitude
err = Float64.(df_V.Uncertainty)   # convert to plain Vector{Float64}
```

**Option 2 — ignore uncertainties entirely (simpler, slightly less accurate):**

```julia
t   = df_V.JD
mag = df_V.Magnitude

plan = LombScargle.plan(t, mag,   # no err argument
    minimum_frequency = 1.0 / 6.0,
    maximum_frequency = 1.0 / 4.0,
    samples_per_peak  = 10
)
```

I'd go with **Option 1** — δ Cep has many observations of varying quality, so weighting by uncertainty gives a cleaner periodogram. Just add the `dropmissing!` line right after your V-band filter:

```julia
df_V = filter(row -> row.Band == "V" && !ismissing(row.Uncertainty), df)
dropmissing!(df_V, :Uncertainty)   # ← add this
sort!(df_V, :JD)

t   = df_V.JD
mag = df_V.Magnitude
err = Float64.(df_V.Uncertainty)
```

The `Float64.()` broadcast conversion is the key step — it strips the `Union{Missing, Float64}` type even when no values are actually missing.
""")

# ╔═╡ de7d5ad8-d6f1-412a-a2c9-069733e47d8b
md"- I agree **option 1** is the best one. But there is another problem:"

# ╔═╡ 6a299349-4b1b-40e4-b875-b795367b6cd1
question_box("Another problem: UndefVarError: `freqat` not defined in this notebook.
Suggestion: check for spelling errors or missing imports.")

# ╔═╡ fe31c076-8e78-4c4b-bfa7-ccdb50d4f647
answer_box(md"""
The function is called `findmaxfreq` in `LombScargle.jl`. Replace this part:

```julia
best_freq   = freqat(pgram, argmaxpower(pgram))
```

with:

```julia
best_freq   = findmaxfreq(pgram)[1]
best_period = 1.0 / best_freq
println("Best period: ", round(best_period, digits=6), " days")
```

`findmaxfreq` returns a vector (in case of ties), so we just take the first element with `[1]`. Everything else in the code stays the same.
""")

# ╔═╡ 7d16e530-80ad-4169-8c23-664153184305
md"""
- However, the code it is still wrong, and even after a few interactions the AI cannot find the rigtht solution.
- Therefore, I work directly on the code, fixing a few problems. Actually, really minor problems due to the specific syntax of the findmax command and to the filtering of the data frame.
- In addition, several uncertainties were exactly ``0``, not allowing the Lomb-Scargle algorithm to work correctly.
"""

# ╔═╡ b3d42d14-9ff7-43d1-a0ee-3d69c04acbfb
begin
	# --- Lomb-Scargle ---
	
	df_V_dt = filter(row -> row.Band == "V", df)
	dropmissing!(df_V_dt)
	
	sort!(df_V_dt, :JD)

	t_dc = df_V_dt.JD[df_V_dt.Uncertainty .> 0]
	mag_dc = df_V_dt.Magnitude[df_V_dt.Uncertainty .> 0]
	err_dc = df_V_dt.Uncertainty[df_V_dt.Uncertainty .> 0]

	plan = LombScargle.plan(t_dc, mag_dc, err_dc,
		minimum_frequency = 1.0 / 6.0, # periods shorter than 6 days
		maximum_frequency = 1.0 / 4.0, # periods longer than 4 days
		samples_per_peak = 10)
	
	pgram = lombscargle(plan)

	# --- Find best period ---
	best_freq_dc = pgram.freq[findmax(pgram.power)[2]]
	best_period_dc = 1.0 ./ best_freq_dc
	println("Best period: ", round(best_period_dc, digits=6), " days")
end

# ╔═╡ 377d78bf-d250-42b6-8d4a-96e634b6e44f
md"- The right period!"

# ╔═╡ fc7331b5-5cd8-44ad-b2c2-49b240e55087
md"""
> While `julia` offers better performances for scientific computing, most of the tutorials are still written in `python`. This would probably make working with `python`, aided by an AI, a simpler task.
"""

# ╔═╡ ddb1955a-d1e7-4452-8bd6-a80913e39487
md"""
## Course Flow
"""

# ╔═╡ 24385054-0edb-4b78-8cce-c89426964dfb
html"""
<table>
  <tr>
	<td></td>
    <td>Previous lecture</td>
    <td>Next lecture</td>
	<td>Course Summary</td>	
  </tr>
  <tr>
    <td>notebook</td>
	<td><a href="./open?path=Lectures/Lecture-GaussianProcesses/Lecture-CO2.jl">Science case about CO₂ content in atmosphere</a></td>    
    <td><a href="./open?path=Lectures/Lecture-AstrostatisticsFuture/Lecture-AstrostatisticsFuture.jl">Lecture about astrostatistics future</a></td>
	<td><a href="./open?path=Course.jl">Course Summary</a></td>    
  </tr>
  <tr>
    <td>html</td>
	<td><a href="../../Lectures/Lecture-GaussianProcesses/Lecture-CO2.html">Science case about CO₂ content in atmosphere</a></td>    
    <td><a href="../../Lectures/Lecture-AstrostatisticsFuture/Lecture-AstrostatisticsFuture.html">Lecture about astrostatistics future</a></td>
	<td><a href="../../Course.html">Course Summary</a></td>    
  </tr>

 </table>
"""

# ╔═╡ e44397bf-7be0-4c6f-bf38-057df876fb92
md"""
**Copyright**

This notebook is provided as [Open Educational Resource](https://en.wikipedia.org/wiki/Open_educational_resources). Feel free to use the notebook for your own purposes. The text is licensed under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/), the code of the examples, unless obtained from other properly quoted sources, under the [MIT license](https://opensource.org/licenses/MIT). Please attribute the work as follows: *Stefano Covino, Time Domain Astrophysics - Lecture notes featuring computational examples, 2026*.
"""

# ╔═╡ 8ed11f8c-4810-433c-9c4f-bbd8cc5f23e9
md"Notebook v1.0.0 - 21 May 2026"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
LombScargle = "fc60dff9-86e7-5f2f-a8a0-edeadbb75bd9"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.10.16"
CairoMakie = "~0.15.9"
DSP = "~0.8.4"
DataFrames = "~1.8.2"
FFTW = "~1.10.0"
HDF5 = "~0.17.2"
LombScargle = "~1.0.3"
PlutoTeachingTools = "~0.4.7"
PlutoUI = "~0.7.80"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.6"
manifest_format = "2.0"
project_hash = "7191a3db40f8556b1d2e0050747da92d1522cfbe"

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
git-tree-sha1 = "0761717147821d696c9470a7a86364b2fbd22fd8"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.5.2"
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
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d0efe2c6fdcdaa1c161d206aa8b933788397ec71"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.6+0"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9cb23bbb1127eefb022b022481466c0f1127d430"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "12177ad6b3cad7fd50c8b3825ce24a99ad61c18f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.CodecZstd]]
deps = ["TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "da54a6cd93c54950c15adf1d336cfd7d71f51a56"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.8.7"

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

[[deps.CoreMath]]
deps = ["CoreMath_jll"]
git-tree-sha1 = "8c0480f92b1b1796239156a1b9b1bfb1b39499b4"
uuid = "b7a15901-be09-4a0e-87d2-2e66b0e09b5a"
version = "0.1.0"

[[deps.CoreMath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a692a4c1dc59a4b8bc0b6403876eb3250fde2bc3"
uuid = "a38c48d9-6df1-5ac9-9223-b6ada3b5572b"
version = "0.1.0+0"

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
git-tree-sha1 = "5fab31e2e01e70ad66e3e24c968c264d1cf166d6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.2"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e86f4a2805f7f19bec5129bc9150c38208e5dc23"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.4"

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
git-tree-sha1 = "12184a8cf11c7cbd90a4db8b2cb2f7b6f057cc46"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.124"

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
git-tree-sha1 = "8f05e9a2e7c2e3eb524102bb2926c5743c07fbe1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.8.0+0"

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
git-tree-sha1 = "6866aec60ef98e3164cd8d6855225684207e9dff"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.12+0"

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
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

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

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "MPIPreferences", "Mmap", "Preferences", "Printf", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "e856eef26cf5bf2b0f95f8f4fc37553c72c8641c"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.17.2"

    [deps.HDF5.extensions]
    MPIExt = "MPI"

    [deps.HDF5.weakdeps]
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"

[[deps.HDF5_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "LibCURL_jll", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "OpenSSL_jll", "TOML", "Zlib_jll", "libaec_jll"]
git-tree-sha1 = "e94f84da9af7ce9c6be049e9067e511e17ff89ec"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.14.6+0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XML2_jll", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "baaaebd42ed9ee1bd9173cfd56910e55a8622ee1"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.13.0+1"

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
deps = ["CRlibm", "CoreMath", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "f1c42fcaca2d8034fe392f3e86c2e0809f75b2a1"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.6"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticIrrationalConstantsExt = "IrrationalConstants"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    IrrationalConstants = "92d709cd-6900-40b7-9082-c6be49f344b6"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "79d6bd28c8d9bccc2229784f1bd637689b256377"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.14"

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
git-tree-sha1 = "3e846e18560a65dcef26febd2ede0160c6831c1c"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.5.1"

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
git-tree-sha1 = "c0c9b76f3520863909825cbecdef58cd63de705a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.5+0"

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
git-tree-sha1 = "17b94ecafcfa45e8360a4fc9ca6b583b049e4e37"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.1.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

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
git-tree-sha1 = "cc3ad4faf30015a3e8094c9b5b7f19e85bdf2386"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d620582b1f0cbe2c72dd1d5bd195a9ce73370ab1"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.0+0"

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

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "9341048b9f723f2ae2a72a5269ac2f15f80534dc"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.3.2+0"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "8e98d5d80b87403c311fd51e8455d4546ba7a5f8"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.12"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "675df097f8eeb28998b2cfe3b25655af73d5f7df"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.5.6+0"

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

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bc95bf4149bf535c09602e3acdf950d9b4376227"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+3"

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
git-tree-sha1 = "9ac7c730c53b3b5d9a73fb900ac4b4fc263774db"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.4.9+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML", "Zlib_jll"]
git-tree-sha1 = "2f3d05e419b6125ffe06e55784102e99325bdbe2"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "5.0.10+0"

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
git-tree-sha1 = "58e5ed5e386e156bd93e86b305ebd21ac63d2d04"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "5d5e0a78e971354b1c7bff0655d11fdc1b0e12c8"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.4"

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

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "90b41ced6bacd8c01bd05da8aed35c5458891749"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.7"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "fbc875044d82c113a9dee6fc14e16cf01fd48872"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.80"

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
git-tree-sha1 = "624de6279ab7d94fc9f672f0068107eb6619732c"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.3.2"

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
git-tree-sha1 = "5e8e8b0ab68215d7a2b14b9921a946fee794749e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.3"

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
git-tree-sha1 = "2700b235561b0335d5bef7097a111dc513b8655e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.7.2"
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
git-tree-sha1 = "ad8002667372439f2e3611cfd14097e03fa4bccd"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.3"

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
git-tree-sha1 = "86f5831495301b2a1387476cb30f86af7ab99194"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.0"
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
deps = ["CodecZstd", "ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "9ca5f1f2d42f80df4b8c9f6ab5a64f438bbd9976"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.9"

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
git-tree-sha1 = "0716e01c3b40413de5dedbc9c5c69f27cddfddfc"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.3"

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

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "80d3930c6347cfce7ccf96bd3bafdf079d9c0390"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.9+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

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

[[deps.libaec_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1411bc34c180946d3cef591de1384012afa6edee"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.1.6+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "850b06095ee71f0135d644ffd8a52850699581ed"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.3+0"

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
git-tree-sha1 = "e51150d5ab85cee6fc36726850f0e627ad2e4aba"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.58+0"

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
git-tree-sha1 = "da8c1f6eee04831f14edcfa5dae611d309807e57"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.3.0+0"

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
# ╟─d5b55ac6-0f38-4603-aa72-bd964d0604a8
# ╟─2f82cbd2-4c34-11f1-901e-01c419d48099
# ╟─efc84c48-a322-4c98-869d-2fe1e32151d5
# ╟─72168e90-1440-42ac-a859-834ec6dcc703
# ╟─9b66dbc6-7dd9-4762-a791-7fe1e433961b
# ╟─fbe81ddf-836e-4e2e-a6e1-475257301c41
# ╟─e1e5bb07-e4a8-4048-8878-9f7c49d57afd
# ╟─71bcd5bc-afb3-48cc-9185-41e44d433ff9
# ╟─fb8ce658-f20a-410e-aadc-b06e4eba6c51
# ╟─7066d72c-6f84-4311-9a12-a05f897731e7
# ╟─2e9b26ef-78eb-4df9-a419-a099cf285719
# ╟─22d009b8-995f-4150-9efc-4a55d02b82ab
# ╟─47ec8223-3159-433e-b91f-10a540edb5bb
# ╟─120af7fd-32e9-4936-950d-e55955380ad0
# ╟─d6330ed7-457e-40e0-9f92-94bd3a75c46d
# ╟─64ed9b01-8c08-4a4c-8d12-ad8b4ab3895f
# ╟─b65614f7-60a5-4ff7-a8ad-a90e2c71bc99
# ╟─c9f20dcb-3bbd-4fd0-a672-b86ce2f14ec5
# ╠═a4e6c2f2-7c1d-4bb9-b1dc-92803bc8d908
# ╠═eee3d157-7e61-4a9b-bd3c-c5ba96480962
# ╟─d9ac1d30-f90a-4df7-a863-aa62f2154dce
# ╟─23650722-09c9-4df1-b277-fa388d3fe160
# ╟─0ff6019f-c988-489d-834a-ec3b58a76b9d
# ╟─c6006389-9924-47c1-8134-49bd68cd840f
# ╟─df96e17c-d8cd-4d93-b74a-7b8bef848e9f
# ╟─18d6ffd2-5d52-44c8-86dd-b07739ca7f69
# ╟─61101ce8-99ec-4736-a328-3234ac93d716
# ╠═081d09ac-692a-46c1-a87e-3350c3c2baad
# ╟─57ba7408-bd56-4227-a58d-effa843d64e4
# ╟─d267408c-246b-4f92-bbbb-410157356a1d
# ╟─5a0cc1ac-1ecb-43fa-a687-048ae9f6cc56
# ╟─cd3f0947-78b7-4593-b7b7-ef82d7a15273
# ╟─5fa720ff-19da-4e22-8972-eb87e3c3d382
# ╠═05a82acf-b130-4642-a6e3-bd910bad9b01
# ╠═208b8cb2-b34f-44ba-97c1-83db0da6fb33
# ╟─3f4a1fe0-80cb-46e9-acc9-d800dea51285
# ╟─589f5ed7-084d-48ad-9162-03df4c1a9609
# ╟─a985cfbd-dc02-4f8f-af40-f55ef584a055
# ╟─164f7527-6434-4df1-8cab-d25f6363ecb1
# ╟─8ac36e25-05ac-4477-958e-a35c9fb4646a
# ╠═9d73d318-6cdf-49e5-b997-b6ac204ec19d
# ╠═3cf92273-1a38-413a-9ed9-87572007aa54
# ╟─56f6c7d9-086e-41cb-9618-69ac9d103fc5
# ╟─ac23e27e-6ad4-405e-b360-c13e674c7528
# ╟─da665df8-9724-412a-acd9-9bd7e7994e7f
# ╟─8e4d77a6-6f77-4621-a7c7-f5b444fd2411
# ╟─845b8722-0d66-42c9-b1b8-eddcfebce6a6
# ╟─791f7bc6-95e5-4beb-a908-563947d1609d
# ╟─1555efa5-cf2c-449d-819a-9d7258c7c24b
# ╟─1998e1eb-6c63-48de-b33b-aa248151d166
# ╟─1a530ba0-c19c-4056-8408-b75956810571
# ╟─dc85b4d0-5dcd-4fa6-ba93-be27991145a0
# ╠═b7f2b3ea-0643-4385-bdc0-8fa9263ba7a8
# ╟─805bf93e-4723-4556-9fca-1931e8611f0e
# ╟─7db81808-236e-4b93-af75-d2db036b7873
# ╟─887a46fe-c37f-452a-8504-18ce0164b29c
# ╠═98a5702a-8ab1-41a7-b0f3-1e7ef035b786
# ╟─db1e7d54-a93b-4185-955e-9a2cd29a4f73
# ╟─f817a31a-26f5-4fbb-ae4a-bc2adf3c9b1a
# ╟─40a2c488-c0b6-491a-8061-bed6a536ba93
# ╟─954f2bdc-973d-4bc0-8402-0107fb3382b7
# ╟─26ddf772-d989-46e2-952a-ccaa9661921f
# ╟─e88c6a1c-7d31-43d8-a35e-bfa0104a93f4
# ╟─de7d5ad8-d6f1-412a-a2c9-069733e47d8b
# ╟─6a299349-4b1b-40e4-b875-b795367b6cd1
# ╟─fe31c076-8e78-4c4b-bfa7-ccdb50d4f647
# ╟─7d16e530-80ad-4169-8c23-664153184305
# ╠═b3d42d14-9ff7-43d1-a0ee-3d69c04acbfb
# ╟─377d78bf-d250-42b6-8d4a-96e634b6e44f
# ╟─fc7331b5-5cd8-44ad-b2c2-49b240e55087
# ╟─ddb1955a-d1e7-4452-8bd6-a80913e39487
# ╟─24385054-0edb-4b78-8cce-c89426964dfb
# ╟─e44397bf-7be0-4c6f-bf38-057df876fb92
# ╟─8ed11f8c-4810-433c-9c4f-bbd8cc5f23e9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
