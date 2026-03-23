### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 1a24d083-b62a-4c89-856b-06cf921c3e06
import Pkg; Pkg.activate(".")

# ╔═╡ 0bf5e933-da14-4023-b11a-dccbc5ee80c1
begin
	using AbstractFFTs
	using CairoMakie
	using CommonMark
	using DataFrames
	using Distributions
	using DSP
	using FITSIO
	using Format
	using LaTeXStrings
	using PlutoUI
	using Statistics
end

# ╔═╡ eab8931d-6c92-4125-9fa4-d77e4f88e245
md"""
**What is this?**


*This jupyter notebook is part of a collection of notebooks on various topics discussed during the Time Domain Astrophysics course delivered by Stefano Covino at the [Università dell'Insubria](https://www.uninsubria.eu/) in Como (Italy). Please direct questions and suggestions to [stefano.covino@inaf.it](mailto:stefano.covino@inaf.it).*
"""

# ╔═╡ 8090d68e-690f-44c5-91ff-424b5a7f584d
md"""
**This is a `Julia` notebook**
"""

# ╔═╡ fb40ee32-4d46-46db-9bd6-0d44665c7f50
Pkg.instantiate()

# ╔═╡ 93e2f620-df83-4fc2-9199-8d320ce07551
# ╠═╡ show_logs = false
md"""
$(LocalResource("Pics/TimeDomainBanner.jpg"))
"""

# ╔═╡ bd58f0cf-4008-4ff1-a831-767beab264b4
md"""
# Science Case: X-Ray Binaries
***
"""

# ╔═╡ 77027e5f-1b87-4d6e-bd4f-c12c49648ac6
# ╠═╡ show_logs = false
cm"""
## X-ray Binaries (XRBs)
***

- X-ray binaries are a class of binary stars that are luminous in X-rays. 

- The X-rays are produced by matter falling from one component, called the “donor” (usually a relatively normal star), to the other component, called the “accretor”.

- The latter is very compact: a neutron star or a black hole. 

- The infalling matter releases gravitational potential energy, up to several tenths of its rest mass, as X-rays.

    - ``L_X \sim 10^{35-39}\ {\rm erg\ s}^{-1}`` in the X-ray band.
    
- The donor star nature affects the time evolution of the binary.

$(LocalResource("Pics/xrbs.png"))
"""

# ╔═╡ 2a8bc9e7-e7bb-484c-bfe9-0a03e1ef9591
# ╠═╡ show_logs = false
cm"""
- Emission can be observed under different “states”.

$(LocalResource("Pics/xrbstates.png"))

$(LocalResource("Pics/hardsoftstates.png"))

- The lifetime and the mass-transfer rate in an X-ray binary depends on the evolutionary status of the donor star, the mass ratio between the stellar components, and their orbital separation.

- X-ray binaries are further subdivided into several (sometimes overlapping) subclasses. Note that the classification by mass (high, intermediate, low) refers to the optically visible donor, not to the compact X-ray emitting accretor.

$(LocalResource("Pics/xrbzoo.jpg"))
"""

# ╔═╡ 5e73c451-6b1e-4913-b782-ed5631ccf460
# ╠═╡ show_logs = false
cm"""
### XRB evolution
***

- This is a very simple sketch of binary evolution

$(LocalResource("Pics/xrbevolution.png"))

- Evolutionary time scales:
    - XMXB: ``t \sim 10-50`` Myrs, comparabke to duration of star forming event so that these systems are star formation tracers.
    - LMXB: ``t \sim 1-10`` Gyrs, comparable to live time of the host galaxies so that these systems are stellar mass tracers.
"""

# ╔═╡ dd834c1d-962f-492e-a563-a01f2d7d6304
md"""
## XRB Zoo
***

- Low-mass X-ray binaries (LMXBs)
    - Soft X-ray transients (SXTs)
    - Symbiotic X-ray binaries
    - Super soft X-ray sources or Super soft sources (SSXs), (SSXB)

- Intermediate-mass X-ray binaries (IMXBs)
    - Ultracompact X-ray binaries (UCXBs)

- High-mass X-ray binaries (HMXBs)
    - Be/X-ray binaries (BeXRBs)
    - Supergiant X-ray binaries (SGXBs)
    - Supergiant Fast X-ray Transients (SFXTs)

- Others
    - X-ray bursters
    - X-ray pulsars
    - Microquasars (radio-jet X-ray binaries that can house either a neutron star or a black hole)

"""

# ╔═╡ cc0cb8bf-292e-4bb4-8dbb-320a3ae6e188
# ╠═╡ show_logs = false
md"""
## GX339-4 / V821 Ara
***

- GX-339-4 is a Galactic Low Mass X-ray Binary (LMXB), and candidate black-hole.

- It is a variable source showing occasionally a flaring activity.

- During the outbursts GX 339-4 shows evolution of quasi-periodic oscillations (QPOs).

- A strong, variable relativistic jet, emitting from radio to infrared wavelengths was observed by several studies.

- An artistic view of the GX339-4 system:

$(LocalResource("Pics/gx339-4.png"))
"""

# ╔═╡ 7b5a229a-1da1-495e-bcd2-a5da91ae73d6
md"""
## Exercize: the QPO in GX339
***

- We are going to analyse an observation of [this source](https://simbad.u-strasbg.fr/simbad/sim-id?Ident=GX+339-4) carried out by the [Rossi X-ray Timing Explorer (RXTE)](https://en.wikipedia.org/wiki/Rossi_X-ray_Timing_Explorer) satellite. 
"""

# ╔═╡ f18e7672-11df-4647-b2bc-2d5c61f1ab1c
begin
	hdu = FITS("gx339_qpo.lc")
	
	tbl = DataFrame(hdu[2])
end

# ╔═╡ 600d1ce7-df15-4c36-9f1f-de0926ece2cf
md"""
- **FRACEXP** is the fraction of a bin of length $\Delta t$ efffectively covered by the observations. We don't need it, but this information can be useful for a careful estimate of the statistics of the light curve.

- Let's look at the spacing of the observations.
"""

# ╔═╡ c43762df-6c71-46e3-97e4-068ce5337ae8
begin
	dt = mean(diff(tbl[!,:TIME]))
	dtvar = var(diff(tbl[!,:TIME]))
	        
	printfmtln("Sampling, mean: {:.5f}s variance: {:.5f}", dt, dtvar)
end

# ╔═╡ 33633521-9c34-4bbf-916c-4e3f73dda3ec
md"""
- Variance "0" means that all time bins are of the same length and the sampling time turns out to be 3.9ms. Therefore the Nyquist frequency is about 128Hz.

- Let's rename, rescale and plot our input variables:
"""

# ╔═╡ e48eb192-244a-49db-aac1-62d3bd6a0bf2
begin
	t = tbl[!,:TIME] .- minimum(tbl[!,:TIME])
	y = tbl[!,:RATE] 
	ey = tbl[!,:ERROR]
	
	
	fg1 = Figure()
	
	ax1fg1 = Axis(fg1[1, 1],
	    xlabel="Time (s)",
	    ylabel=L"Rate (counts s$^{-1}$)",
	    title="RXTE observation of GX339-4"
	    )
	
	scatter!(t,y,label="Data",color=:blue)
	errorbars!(t,y,ey,label="",color=:blue)
	
	axislegend()
	
	fg1
end

# ╔═╡ 88c27d21-788f-4fac-adaa-358a2b171a9d
md"""
- It is difficult to see the details. With such a rapid sampling we have a very large number of points (almost 180000!). Let's try a zoom showing the first 0.5s of data.
"""

# ╔═╡ 48425a2a-812e-4b12-925e-7ee8c1e1db07
begin
	fg2 = Figure()
	
	ax1fg2 = Axis(fg2[1, 1],
	    xlabel="Time (s)",
	    ylabel=L"Rate (counts s$^{-1}$)",
	    title="RXTE observation of GX339-4"
	    )
	
	scatter!(t,y,label="Data",color=:blue)
	errorbars!(t,y,ey,label="",color=:blue)
	
	xlims!(0,0.5)
	
	axislegend()
	
	fg2
end

# ╔═╡ 15fbb4b4-cb89-4547-8209-79d1b34f4abe
md"""
- Let's now compute and show the periodogram
"""

# ╔═╡ 893ef88e-2e28-4216-b8b6-65b29d91f5fa
begin
	# We are using the AbstractFFTs package
	
	
	function FourierPeriodogram(signal,fs)
	    N = length(signal)
	    freqs = fftfreq(N,fs)
	    positive = freqs .> 0  
	    ft = fft(signal)
	    powers = abs.(ft).^2
	    return freqs[positive], powers[positive]
	end
	
	rxtefr,rxtepw = FourierPeriodogram(y,1/dt)
	
	rxtepwleahy = 2*rxtepw/(var(y)*length(y))    # Convert to Leahy normalization
	
	fg3 = Figure()
	
	ax1fg3 = Axis(fg3[1, 1],
	    xlabel="Frequency (Hz)",
	    ylabel="Power",
	    )
	
	lines!(rxtefr,rxtepwleahy,label="DFT")
	
	xlims!(1,120)
	ylims!(0,2e2)
	
	axislegend()
	
	fg3
	
end

# ╔═╡ 67b86998-3a60-40c6-a7b7-8163b1db8bb3
md"""
- A poweful peak is clearly visible. Again let's zoom in.
"""

# ╔═╡ bef6b948-944e-4f94-be0d-c636c76c5b73
begin
	fg4 = Figure()
	
	ax1fg4 = Axis(fg4[1, 1],
	    xlabel="Frequency (Hz)",
	    ylabel="Power",
	    )
	
	lines!(rxtefr,rxtepwleahy,label="DFT")
	
	xlims!(1,10)
	ylims!(0,2e2)
	
	axislegend()
	
	fg4
end

# ╔═╡ 087b056d-31bd-4949-b91a-afd75d33552d
md"""
- And we can see a subtructure in the peak. Indeed, it is firmed by a large number of peaks with similar frequency (i.e. a QPO!). The peak centroid is close to $\sim$ 5.6 Hz. 

- We might want to test whether the identified peak can be due to noise fluctuation. 
    - A real analysis would ask to consider all the substructure visible in the periodogram. However, to simplify the discussion, let's pretend that the highest peak is the only one visible in the spectrum.
"""

# ╔═╡ 9d182ef8-8a54-4023-86d9-b7ce729b76fb
begin
	flt = (rxtefr .> 1) .& (rxtefr .< 10)
	pmax, idmax = findmax(rxtepwleahy[flt])
	fmax =  rxtefr[flt][idmax]
	printfmtln("Maximum power: {:.2f} at frequency {:.2f}Hz in a {:d} frequency bin long periodogram.", pmax, fmax, length(rxtefr))
end

# ╔═╡ 06346f7c-115d-4f3a-bc9e-7a8ba5c28f1f
md"""
- Let's not forget that any significance has to be computed considering also the number of indipendent frequncies considered for the analysis. This is often the leading a factor.

- Let's assume to have a peak with probability 0.9997, i.e. $3\sigma$. This is generally a rather promoment feature.
- And let's also assume our input time-serien consists of 200 points, geberating therefore a 100 frequency bin periodogram.
- Ledt's compute the final significance:
"""

# ╔═╡ 7362ae7a-64e5-404b-ba31-381f7d1905c9
begin
	eps = 1-0.9997 # 3sigma
	printfmtln("FAP: {:.2f}%", eps*100)
	Ppre = (1-eps)
	Ppost = (1-eps)^100
	printfmtln("Single-trial probability: {:.2f}%",Ppre*100)
	printfmtln("Post-trial probability: {:.2f}%", Ppost*100)
end

# ╔═╡ f4106f3c-9d48-4016-8249-c5bef66f198d
md"""
- Quite interestingly, the final probability can justufy a claim of a hint of periodocity, but nothing more: high sigificance periodicities need to be of very high single-trial significancd.

- Coming back to our exericize, the measure power is yields a very low probability to be due to noise. Yet, the number of frequencies in the periodogram is also very high. 

- Nevertheless, the identified features, within the noise assumptions (white noise), are highly significant.
"""

# ╔═╡ d78e47ad-4102-4514-b130-777378036e3e
begin
	cs = Chisq(2)
	
	probmax = ccdf(cs,pmax)
	printfmtln("Single trial probability: {:.2g}%", 100*(1-probmax))
	postprobmax = (1-probmax)^length(rxtefr)
	printfmtln("Post trial probability: {:.2g}%", 100*postprobmax)
end

# ╔═╡ f2254265-cedc-4efa-b303-b9798743b17a
md"""
- The large number of peaks wit similar but not identical frequency might be due to some sort of periodicity with period drifting in time.

- we can explore this possibility computing the DFT in several time-slot, in order to obtain a time-resolved DFT (more on this and other related techniques in Lecture 6.
"""

# ╔═╡ 5a175c4f-0bb7-4d5a-8a47-25e0478019f7
begin
	tresper = spectrogram(y,div(length(y),70),0,fs=1/dt)
	
	
	fg5 = Figure()
	
	ax1fg5 = Axis(fg5[1, 1],
	    xlabel="Time (s)",
	    ylabel="Frequency (Hz)",
	    )
	
	p = heatmap!(tresper.time,tresper.freq,tresper.power',colorscale=log10,colorrange=(1e4,5e5))
	Colorbar(fg5[:, end+1], p)
	
	
	#xlims!(1,10)
	ylims!(1,10)
	
	#axislegend()
	
	fg5
	
end

# ╔═╡ 00588dbb-867a-44f8-9cb8-4d6563456620
md"""
- We see that computing the peripdogram in 70 10s bins the drifting of trhe QPO frequency is clearly visible.
"""

# ╔═╡ 1d16be3d-0ed8-4ee7-be96-29ac7a87d46e
md"""
## Reference & Material

Material and papers related to the topics discussed in this lecture.

- [Nespoli et al. (2003) - "A transient variable 6 Hz QPO from GX 339-4”](https://ui.adsabs.harvard.edu/abs/2003A%26A...412..235N/abstract)
"""

# ╔═╡ baa5e1a6-fe23-4a48-9ffd-b66be4253418
md"""
## Further Material

Papers or sites for examining more closely some of the discussed topics.

- [Bahramian & Degenaar (2023) - "Low-Mass X-ray Binaries"](https://ui.adsabs.harvard.edu/abs/2023hxga.book..120B/abstract)
- [Tan C. (2021) - "High-Mass X-ray binary: Classification, Formation, and Evolution"](https://iopscience.iop.org/article/10.1088/1742-6596/2012/1/012119)
"""

# ╔═╡ f74f9d82-6079-45af-9798-e4e0ad9600ee
md"""
### Credits
***

This notebook contains RXTE observations kindly provided by dr. Tomaso Belloni (INAF / Brera Astronomical Observatory).
"""

# ╔═╡ 545551de-0e43-4fee-9870-ac595f5cbf19
cm"""
## Course Flow

<table>
  <tr>
    <td>Previous lecture</td>
    <td>Next lecture</td>
  </tr>
  <tr>
    <td><a href="../Science%20Case%20-%20Sunspot%20Number/Lecture-SunspotNumber.ipynb">Science case about Sunspot number</a></td>
    <td><a href="../Lecture%20-%20Lomb-Scargle/Lecture-Lomb-Scargle.ipynb">Lecture about irregular sampling</a></td>
  </tr>
 </table>


"""

# ╔═╡ 8382d0c7-1e8b-41b0-a3e5-9a2514ad355e
md"""
**Copyright**

This notebook is provided as [Open Educational Resource](https://en.wikipedia.org/wiki/Open_educational_resources). Feel free to use the notebook for your own purposes. The text is licensed under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/), the code of the examples, unless obtained from other properly quoted sources, under the [MIT license](https://opensource.org/licenses/MIT). Please attribute the work as follows: *Stefano Covino, Time Domain Astrophysics - Lecture notes featuring computational examples, 2024*.
"""

# ╔═╡ Cell order:
# ╟─eab8931d-6c92-4125-9fa4-d77e4f88e245
# ╟─8090d68e-690f-44c5-91ff-424b5a7f584d
# ╠═1a24d083-b62a-4c89-856b-06cf921c3e06
# ╠═fb40ee32-4d46-46db-9bd6-0d44665c7f50
# ╠═0bf5e933-da14-4023-b11a-dccbc5ee80c1
# ╟─93e2f620-df83-4fc2-9199-8d320ce07551
# ╟─bd58f0cf-4008-4ff1-a831-767beab264b4
# ╟─77027e5f-1b87-4d6e-bd4f-c12c49648ac6
# ╟─2a8bc9e7-e7bb-484c-bfe9-0a03e1ef9591
# ╟─5e73c451-6b1e-4913-b782-ed5631ccf460
# ╟─dd834c1d-962f-492e-a563-a01f2d7d6304
# ╟─cc0cb8bf-292e-4bb4-8dbb-320a3ae6e188
# ╟─7b5a229a-1da1-495e-bcd2-a5da91ae73d6
# ╠═f18e7672-11df-4647-b2bc-2d5c61f1ab1c
# ╟─600d1ce7-df15-4c36-9f1f-de0926ece2cf
# ╠═c43762df-6c71-46e3-97e4-068ce5337ae8
# ╟─33633521-9c34-4bbf-916c-4e3f73dda3ec
# ╠═e48eb192-244a-49db-aac1-62d3bd6a0bf2
# ╟─88c27d21-788f-4fac-adaa-358a2b171a9d
# ╠═48425a2a-812e-4b12-925e-7ee8c1e1db07
# ╟─15fbb4b4-cb89-4547-8209-79d1b34f4abe
# ╠═893ef88e-2e28-4216-b8b6-65b29d91f5fa
# ╟─67b86998-3a60-40c6-a7b7-8163b1db8bb3
# ╠═bef6b948-944e-4f94-be0d-c636c76c5b73
# ╟─087b056d-31bd-4949-b91a-afd75d33552d
# ╠═9d182ef8-8a54-4023-86d9-b7ce729b76fb
# ╟─06346f7c-115d-4f3a-bc9e-7a8ba5c28f1f
# ╠═7362ae7a-64e5-404b-ba31-381f7d1905c9
# ╟─f4106f3c-9d48-4016-8249-c5bef66f198d
# ╠═d78e47ad-4102-4514-b130-777378036e3e
# ╟─f2254265-cedc-4efa-b303-b9798743b17a
# ╠═5a175c4f-0bb7-4d5a-8a47-25e0478019f7
# ╟─00588dbb-867a-44f8-9cb8-4d6563456620
# ╟─1d16be3d-0ed8-4ee7-be96-29ac7a87d46e
# ╟─baa5e1a6-fe23-4a48-9ffd-b66be4253418
# ╟─f74f9d82-6079-45af-9798-e4e0ad9600ee
# ╟─545551de-0e43-4fee-9870-ac595f5cbf19
# ╟─8382d0c7-1e8b-41b0-a3e5-9a2514ad355e
