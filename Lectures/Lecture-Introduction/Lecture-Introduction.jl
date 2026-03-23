### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ a9389fd0-72fd-4d3c-94a7-37d856ccf58b
using PlutoUI

# ╔═╡ b8bbafd8-25bd-4681-96bd-607baabfa138
md"""
**What is this?**


*This jupyter notebook is part of a collection of notebooks on various topics discussed during the Time Domain Astrophysics course delivered by Stefano Covino at the [Università dell'Insubria](https://www.uninsubria.eu/) in Como (Italy). Please direct questions and suggestions to [stefano.covino@inaf.it](mailto:stefano.covino@inaf.it).*
"""

# ╔═╡ 0ad7e253-35e7-4683-ba47-e2bef62bc491
md"""
**This is a `Pluto` notebook**
"""

# ╔═╡ 5f37e9b0-0f58-4499-9063-b7765e422f3e
TableOfContents()

# ╔═╡ 0d0c5ef4-0e29-44cb-8d9b-91a0b62f0a58
md"""
$(LocalResource("Pics/TimeDomainBanner.jpg"))
"""

# ╔═╡ 28b9ac44-044b-41bd-af41-a21ad5fe0b5e
md"""
# Introduction
***
"""

# ╔═╡ f7af2d56-3d0a-46d7-8859-a200189fd159
md"""
## Contacts
***

$(LocalResource("Pics/Stefano.png"))

- Stefano Covino
- INAF / Brera Astronomical Observatory
- +39 02 72320475
- +39 3316748534 (if urgent…)
- Emails: [stefano.covino@inaf.it](mailto:stefano.covino@inaf.it)  - [stefano.covino@uninsubria.it](mailto:stefano.covino@uninsubria.it)
- Web: [https://sites.google.com/a/inaf.it/stefano-s-site/](https://sites.google.com/a/inaf.it/stefano-s-site/)

$(LocalResource("Pics/Banner.png"))
"""

# ╔═╡ 26bd7a49-bdda-415f-8220-211016784c75
md"""
## Main Goal of the course: Have fun!
***

 $(LocalResource("Pics/data.jpg", :width => 300)) $(LocalResource("Pics/regression.jpg", :width => 300))

"""

# ╔═╡ 7d83b80c-4721-4454-997a-d981ccb3eca1
md"""
#### - This is serious. Specialized academic courses, in a manner of speaking, should be enjoyed...!

"""

# ╔═╡ 1ef75ad8-5f7f-4146-af0f-66982f2e7d62
md"""
## Time-Series are ubiquitous
***

- Anytime we have a measurement repetated multiple times we have a time-series.

$(LocalResource("Pics/CO2T.png",:width=>600))

$(LocalResource("Pics/CO2.png",:width=>600))

$(LocalResource("Pics/Neptune.png",:width=>600))

- As a matter of fact, a time-series does not need to have "time" as index!

$(LocalResource("Pics/PAMELA.png",:width=>600))

$(LocalResource("Pics/satellite.png",:width=>600))
"""

# ╔═╡ d3443a4d-62b5-4ca9-84dc-66e19c284ac2
md"""
## Temptative program (it may change…)
***

1. Introduction
2. Statistics reminder - part I
3. Statistics reminder - part II
4. Spectral analysis - part I
5. Spectral analysis - part II
6. Science cases: Sunspots Number - X-ray Binaries
7. Irregularly sampled time series - part I
8. Irregularly sampled time series - part II
9. Science Cases - Variable Stars  
10. Time domain analysis - part I
11. Time domain analysis - part II
12. Science Cases - AGN and blazars
13. Wavelet analysis and climatology science case
14. Time of arrival analysis and paleo-climatology science case
15. Science case: FRBs
16. Non-parametric methods - part I
17. Non-parametric methods - part II
18. Singular spectrum analysis - part I
19. Singular spectrum analysis - part II
20. Gaussian processes - part I
21. Gaussian processes - part II
22. Science case: GRBs
23. Astrostatistics: final considerations
"""

# ╔═╡ 482cfb39-4f17-4138-b68e-317bad0b325a
md"""
## How is the course managed?
***

### Frontal lectures

- These are the traditional university lectures.

- Although this increases the organizational complexity substantially, I am availbale to stream and record my lectures, if needed.

- There are contraindications. As a matter of fact, this is one of few cases where a remote access is not even close as effective as being in presence.

$(LocalResource("Pics/FrontalLectures.jpg"))

### Real research life examples…

- Scientists working in the field will deliver "didactic lectures", allowing one to see most of ideas developed during the course applied in a real research environment.

$(LocalResource("Pics/Paperino.jpg"))

### (Optional) papers to deepen our knowledge…

- Most of the topics discussd during the course can be investigated thoroughly and papers from astrophysical (mainly) literature are presented for particularly concerned readers.

$(LocalResource("Pics/Papersetal.jpg"))

### Question time

- The course is divided in several main sections. At the end of each of them, some time will be devoted to open discussions and questions.

$(LocalResource("Pics/Questions.gif"))

### Lectures from specialists in the field

- Together with regular lectures, a few specialists in the field, i.e. scientist carrying out researches by time-domain tools and techniques, are invited to describe their works.

$(LocalResource("Pics/Nilus.jpg"))

### Language

- According to university guidelines, lectures will be delivered in English. Of course, a fair evaluation of the context might ask some flexibility.

$(LocalResource("Pics/language.jpg"))

### Statistical framework

- During this course we are going to work in a Bayesian framework.

- Bayesian statistics is an approach to inferential statistics based on Bayes' theorem, where available knowledge about parameters in a statistical model is updated with the information in observed data. The background knowledge is expressed as a prior distribution and combined with observational data in the form of a likelihood function to determine the posterior distribution. The posterior can also be used for making predictions about future events.

- Nevertheless, we are not dogmatic and mentions or applications based on familiar "frequentist" approaches are presented and discussed, when we deem it opportune.

$(LocalResource("Pics/Bayesians.png"))

### Programming languages

- Most of the examples we are going to analyze during the course are based on some sort of computer analysis.

- `Python` is *de-facto* the standard language in data science.
    - Yet, while this language is definitely truly amazing, well designed and worth mastering, for the specific needs of scientific computing there are alternatives of growing popularity.

- We threfore provide examples mainly with `Julia`, and encourage the students to get some confidence with this programming language too.

- Indeed, we provide examples mainly with `Julia`, and encourage the students to get some confidence with this programming language too.

- Notebooks are written by the [markdown language](https://www.markdownguide.org/basic-syntax/), a simple language integrating features of the HTML and latex languages.


 $(LocalResource("Pics/python.png"))  $(LocalResource("Pics/julia.png"))

"""

# ╔═╡ 3c960134-e9d9-4b66-891f-6b3522c691c4
md"""
- A remarkable introducti0n to the `julia` language for a scientist is available online, e.g. [Julia data science](https://github.com/tirthajyoti/Julia-data-science). 
"""

# ╔═╡ 6827c9f6-b98a-4a88-a192-5c61f375f1d5
md"""
## Warning! The course is not only for astrophysicists!

- It is indeed part of the set of courses for future astrophysicists. Nevertheles, almost nothing we are going to discuss is truly only for astrophysics. In reality, several applications and ideas are taken from other fields, i.e. economics, social sciences, climatology, etc.

$(LocalResource("Pics/astrophysics.jpg"))
"""

# ╔═╡ d5e83a80-2c3a-48a5-937f-95b50e745bc4
md"""
## Final assessment

- The final examination is an oral one.

- *Students* must interact with the teacher in advance of the examination and a science case obtained by the modern literature will be selected.

- The *student* will be asked to properly describe the main formal aspects of the study and discuss critically the reliability and limits of the presented results.


> As a general rule, in order to take the exam, attending $\sim$ 50% of the lectures is required.

"""

# ╔═╡ c5257f84-999b-4655-9ade-2bcbc7eb324c
md"""
## Gitlab repository

- Slides, notebooks, papers, etc. are available on [gitlab](https://www.ict.inaf.it/gitlab/stefano.covino/TimeDomainAstrophysics.git)
- Check the repository frequently since is (rather often) updated  during the course.

 $(LocalResource("Pics/gitlab.jpg", :width => 200))  $(LocalResource("Pics/gitlabcourse.png", :width => 200)) 


### How to use the repository:

- Just surf the site with your preferred web browser.
    - It is probably enough although you cannot have real interactions (you just read...).
    - There are also *git* graphical clients, for essentially any OS, that can make the surfing esasier.

<br>

- There is also a "nerdier" solution!
- Open a terminal and move to a directory where material will be stored.
```
cd mydir
git clone https://www.ict.inaf.it/gitlab/stefano.covino/TimeDomainAstrophysics.git
cd TimeDomainAstrophysics
git pull
```

- This clones the whole tree (i.e. the course material, about 3GB).

> Repeating frequently the last command (`git pull`) you will always have the tree fully updated and you notebooks, data, papers, etc. ready to be used on your computer.

"""

# ╔═╡ 4520ebfa-8d97-4f1c-83d4-27dfc62e54f5
md"""
## Course calendar and remote connection

The course calendar, notes, advised, topics discussed during a given lecture, and remote connection details can be found at this [url](https://calendar.google.com/calendar/u/0?cid=Y19iMmFmN2RiNjQ0OWNjMTdjY2ZmMzJlMzE3ZjVhZWQ4N2FkYzliN2FkZWFmMzY1YjllYmMwOWFkODA4MjhlNzZjQGdyb3VwLmNhbGVuZGFyLmdvb2dsZS5jb20).
"""

# ╔═╡ 6bce4cd0-67bc-408b-94c0-6c8de26755c0
md"""
## Relaxing time(-series...)

$(LocalResource("Pics/relaxing.png"))
"""

# ╔═╡ 87f81f95-2f33-496a-a9df-cbc711e51e3c
md"""
## Reference & Material

- The course is based on published scientific papers distributed by the teacher before any main topic is addressed.

- Science cases are based on actual scientific papers as well.

- Slides prepared by the teacher will also be distributed.

    - A general introductory text to time series analysis as: [“Introduction to Time Series and Forecasting”, by P.J. Brockwell and R.A Davis](https://link.springer.com/book/10.1007/978-3-319-29854-2) might be useful. However, any other analogous text easily obtainable by the student will be fine as well.

- Two textbooks more strictly related to the topics discussed during the course mainly, but not only, for astrophysical applications are:
    - [“Modern Statistical Methods for Astronomy”, by E.D. Feigelson and G.J. Babu](https://www.cambridge.org/core/books/modern-statistical-methods-for-astronomy/941AE392A553D68DD7B02491BB66DDEC)
    - [“Statistics, data Mining and Machine Learning in Astronomy”, by Ivezić et al.](https://press.princeton.edu/books/hardcover/9780691198309/statistics-data-mining-and-machine-learning-in-astronomy)
"""

# ╔═╡ 2b826fbe-9dda-483f-9af6-8b89d7a6837d
md"""
## Further Material

Papers for examining more closely some of the discussed topics.

- [Voughan et al. (2013) - "Random Time Series in Astronomy"](https://royalsocietypublishing.org/doi/10.1098/rsta.2011.0549)
- [Storopoli et al. (2021) - "Julia Data Science"](https://juliadatascience.io/)
"""

# ╔═╡ c99cc57a-d012-4d11-a3e6-1e2ba35ce92e
md"""
## Course Flow
"""

# ╔═╡ 59bfe64e-000e-467f-859e-07249d0c9273
html"""
<table>
  <tr>
    <td>Previous lecture</td>
    <td>Next lecture</td>
  </tr>
  <tr>
    <td><a href="./open?path=Course.jl">Course Summary</a></td>    
    <td><a href="./open?path=Lectures/Lecture-StatisticsReminder/Lecture-StatisticsReminder.jl">Statistics Reminder</a></td>
  </tr>
 </table>
"""

# ╔═╡ 0fc5c982-e5e1-4a43-a0d6-86fe706e2601
md"""
**Copyright**

This notebook is provided as [Open Educational Resource](https://en.wikipedia.org/wiki/Open_educational_resources). Feel free to use the notebook for your own purposes. The text is licensed under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/), the code of the examples, unless obtained from other properly quoted sources, under the [MIT license](https://opensource.org/licenses/MIT). Please attribute the work as follows: *Stefano Covino, Time Domain Astrophysics - Lecture notes featuring computational examples, 2026*.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.61"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.5"
manifest_format = "2.0"
project_hash = "95f3f934b7e2c5249ef6e1068e98d6029094e806"

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

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

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
# ╟─b8bbafd8-25bd-4681-96bd-607baabfa138
# ╟─0ad7e253-35e7-4683-ba47-e2bef62bc491
# ╟─a9389fd0-72fd-4d3c-94a7-37d856ccf58b
# ╟─5f37e9b0-0f58-4499-9063-b7765e422f3e
# ╟─0d0c5ef4-0e29-44cb-8d9b-91a0b62f0a58
# ╟─28b9ac44-044b-41bd-af41-a21ad5fe0b5e
# ╟─f7af2d56-3d0a-46d7-8859-a200189fd159
# ╟─26bd7a49-bdda-415f-8220-211016784c75
# ╟─7d83b80c-4721-4454-997a-d981ccb3eca1
# ╟─1ef75ad8-5f7f-4146-af0f-66982f2e7d62
# ╟─d3443a4d-62b5-4ca9-84dc-66e19c284ac2
# ╟─482cfb39-4f17-4138-b68e-317bad0b325a
# ╟─3c960134-e9d9-4b66-891f-6b3522c691c4
# ╟─6827c9f6-b98a-4a88-a192-5c61f375f1d5
# ╟─d5e83a80-2c3a-48a5-937f-95b50e745bc4
# ╟─c5257f84-999b-4655-9ade-2bcbc7eb324c
# ╟─4520ebfa-8d97-4f1c-83d4-27dfc62e54f5
# ╟─6bce4cd0-67bc-408b-94c0-6c8de26755c0
# ╟─87f81f95-2f33-496a-a9df-cbc711e51e3c
# ╟─2b826fbe-9dda-483f-9af6-8b89d7a6837d
# ╟─c99cc57a-d012-4d11-a3e6-1e2ba35ce92e
# ╟─59bfe64e-000e-467f-859e-07249d0c9273
# ╟─0fc5c982-e5e1-4a43-a0d6-86fe706e2601
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
