### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 2f82cbd2-4c34-11f1-901e-01c419d48099
begin
	using PlutoUI
	using PlutoTeachingTools
end

# ╔═╡ d5b55ac6-0f38-4603-aa72-bd964d0604a8
md"""
**What is this?**


*This notebook is part of a collection of `pluto` notebooks on various topics discussed during the Time Domain Astrophysics course delivered by Stefano Covino at the [Università dell'Insubria](https://www.uninsubria.eu/) in Como (Italy). Please direct questions and suggestions to [stefano.covino@inaf.it](mailto:stefano.covino@inaf.it).*
"""

# ╔═╡ 355d961d-ad28-479c-a0a1-1865bc809c99
WidthOverDocs()

# ╔═╡ efc84c48-a322-4c98-869d-2fe1e32151d5
TableOfContents()

# ╔═╡ 72168e90-1440-42ac-a859-834ec6dcc703
md"""
$(LocalResource("Pics/TDA-banner.jpeg"))
"""

# ╔═╡ 9b66dbc6-7dd9-4762-a791-7fe1e433961b
md"""
# Machine Learning Approaches for Time Series
***
"""

# ╔═╡ fbe81ddf-836e-4e2e-a6e1-475257301c41
md"""
- Machine learning methods for time series aim to overcome some limitations of classical models by either automating the feature extraction or by using flexible nonlinear models that can capture complex patterns in data. 

- We can loosely divide ML approaches into two groups: 
   1. techniques that transform the time series into features for use with generic ML algorithms (feature-based approaches);
   2. specialized ML models that directly handle sequential data, including deep learning models.
"""

# ╔═╡ 71bcd5bc-afb3-48cc-9185-41e44d433ff9
md"""

## Feature-Based and Ensemble Methods
***

- One straightforward approach to apply machine learning is to convert the time series prediction problem into a supervised learning problem by extracting features from past observations. 
- For example, to forecast a value ``y_{t+h}`` (h steps ahead), one can take the recent history ``(y_t,y_{t−1},...,y_{t−k})`` and various summary statistics as input features to a regression model. These features might include lags, moving averages, recent trend estimates, seasonal indicators (to encode the time of year or day), etc. 
- Once the time series is represented in a feature vector form, any regression algorithm can be used: [linear regression](https://en.wikipedia.org/wiki/Linear_regression), [support vector machines](https://en.wikipedia.org/wiki/Support_vector_machine), [decision trees](https://en.wikipedia.org/wiki/Decision_tree), [random forests](https://en.wikipedia.org/wiki/Random_forest), [gradient boosting](https://en.wikipedia.org/wiki/Gradient_boosting) machines, etc. 
  - This approach does not explicitly use a sequential model, but rather relies on the feature engineering to capture temporal patterns.
"""

# ╔═╡ fb8ce658-f20a-410e-aadc-b06e4eba6c51
md"""
- Ensemble tree-based methods, such as [Random Forests](https://en.wikipedia.org/wiki/Random_forest) and [Gradient Boosting](https://en.wikipedia.org/wiki/Gradient_boosting) Machines (e.g. [XGBoost](https://en.wikipedia.org/wiki/XGBoost), [LightGBM](https://en.wikipedia.org/wiki/LightGBM), [CatBoost](https://en.wikipedia.org/wiki/CatBoost)), have proven extremely effective in many time series forecasting competitions and applications. 

- These tree-based models can naturally handle nonlinear relationships and interactions between features, and they are fairly robust to outliers and missing data. 

- One should note, however, that these models do not inherently account for the temporal order beyond what is coded in the features, so one must be careful to include enough lag features and relevant context.
"""

# ╔═╡ 7066d72c-6f84-4311-9a12-a05f897731e7
md"""
- A critical practical consideration when applying any ML model to time series is the evaluation methodology. 

- Standard k-fold [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_(statistics)), which randomly shuffles data into training and test folds, is inappropriate for time series because it breaks the temporal ordering and introduces data leakage: the model may train on future data to predict the past, producing overly optimistic performance estimates. 

- The correct approach is time series cross-validation (also called walk-forward validation or expanding/sliding window validation). In this scheme, the training set always consists of observations before the test set in time. For example, one trains on months 1–12, tests on month 13; then trains on months 1–13, tests on month 14; and so on. This respects the causal structure of the data and produces realistic estimates of out-of-sample performance. 

> Failure to use proper temporal validation is one of the most common pitfalls when applying machine learning to time series in practice.
"""

# ╔═╡ 2e9b26ef-78eb-4df9-a419-a099cf285719
md"""
- A complementary strategy is automated feature extraction. [Feature engineering](https://en.wikipedia.org/wiki/Feature_engineering) can systematically compute hundreds of statistical features from a time series (e.g. entropy, number of peaks, autoregressive coefficients, wavelet energies) and then use relevance filtering to select the most informative ones. 

- This approach makes it possible to apply standard classifiers or regressors to time series data with minimal manual feature engineering, and has been adopted in fields ranging from predictive maintenance to medical diagnostics.


- Another traditional approach related to features is using distance-based or instance-based learning for time series. 

- For classification tasks, a classic method is the [nearest-neighbor classifier](https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm) under a time-series similarity measure like [Dynamic Time Warping (DTW)](https://en.wikipedia.org/wiki/Dynamic_time_warping). DTW measures similarity between time series that may be stretched or misaligned in time by optimally warping the time axis. 

- A [1-NN classifier](https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm) with DTW distance was a strong baseline for time series classification for many years. Some approaches learn shape-based features (so-called shapelets) that are small subsequences particularly discriminative of classes.

- More recently, the ROCKET (RandOm [Convolutional KErnel Transform](https://en.wikipedia.org/wiki/Convolutional_neural_network)) family of methods has emerged as a state-of-the-art approach for time series classification. ROCKET applies a large number of random convolutional kernels to the input series, extracts simple summary statistics (max value and proportion of positive values) from each convolution output, and feeds the resulting feature vector into a linear classifier.

- In astronomy, these methods are promising for the rapid classification of transients in survey pipelines where computational speed is critical.

"""

# ╔═╡ 47ec8223-3159-433e-b91f-10a540edb5bb
md"""

## Hybrid and Specialized Models
***

- Many real-world applications benefit from combining the strengths of classical and ML approaches. For example, one can use an [ARIMA model](https://en.wikipedia.org/wiki/Autoregressive_integrated_moving_average) to handle seasonal patterns and trend, and feed its residuals (which ideally contain only more complex variations) into a machine learning model to capture any remaining nonlinear structure. This is sometimes called ARIMA residual learning or hybrid modeling. 

- In a different vein, another hybrid approach is to use machine learning to choose among or combine forecasts from multiple classical models (sometimes via stacking or [meta-learning](https://en.wikipedia.org/wiki/Meta-learning_(computer_science))).

- Another area of development has been in tailored ML models for time series. 
  - For example, shapelet transformation methods explicitly search for small subsequences that best differentiate classes, and then use those distances as features in a classifier. 
  - Methods like time series forest create random interval features from the series (e.g., mean, standard deviation over random time intervals) and build an ensemble of trees, which has proven to be a competitive approach for classification. 

- These methods often blend the line between manual feature engineering and automated feature learning.

- It is also worth mentioning the increasing role of unsupervised learning on time series. 

  - Techniques such as clustering or anomaly detection often proceed by defining a distance measure between series (like [DTW](https://en.wikipedia.org/wiki/Dynamic_time_warping) or [correlation-based distances](https://en.wikipedia.org/wiki/Distance_correlation)) and then using algorithms like [k-means](https://en.wikipedia.org/wiki/K-means_clustering) or [DBSCAN](https://en.wikipedia.org/wiki/DBSCAN) to cluster similar time series or to identify outliers. 

  - Anomaly detection can also be tackled by one-class classification methods or by building forecasting models and flagging large prediction residuals as anomalies. 

- In astronomy, anomaly detection is crucial for finding novel transient events or unusual variable stars in surveys; ML approaches like [autoencoders](https://en.wikipedia.org/wiki/Autoencoder) (see later) and clustering in feature space have been used to let the data itself indicate what’s ”normal” versus ”odd” without having to manually label anomalies.

"""

# ╔═╡ 5c01fb12-0d8d-4978-b080-defa6e881060
md"""
- Overall, classical ML approaches (feature-based, distance-based, ensembles) have enriched the time series toolkit by enabling flexible nonlinear modeling and data-driven feature extraction. 

- However, they often still rely on human expertise to choose the right features or transformations of the time axis. 

- The next major leap in time series analysis came with the rise of deep learning, which aims to automatically learn the relevant features and complex patterns directly from raw sequential data.
"""

# ╔═╡ 120af7fd-32e9-4936-950d-e55955380ad0
md"""

## Deep Learning for Time Series
***

- In recent years, deep learning has revolutionized many domains by enabling end-to-end learning of complex patterns from raw data. Time series are no exception.

- Deep neural networks, with their ability to approximate highly nonlinear functions, have been applied to time series tasks with considerable success. The most common deep learning architectures for sequential data are [Recurrent Neural Networks](https://en.wikipedia.org/wiki/Recurrent_neural_network) (RNNs), [Convolutional Neural Networks](https://en.wikipedia.org/wiki/Convolutional_neural_network) (CNNs), and more recently [Transformer](https://en.wikipedia.org/wiki/Transformer_(deep_learning)) networks with [attention](https://en.wikipedia.org/wiki/Attention_(machine_learning)) mechanisms. Each offers different advantages for modeling time series, and often they are used in combination.
"""

# ╔═╡ d6330ed7-457e-40e0-9f92-94bd3a75c46d
md"""

### Recurrent Neural Networks and LSTMs
***

- Recurrent neural networks are specifically designed to handle sequential inputs by maintaining an internal state (memory) that is updated at each time step. 

  - Unlike a standard feed-forward network that assumes inputs are independent, an RNN processes data one step at a time, using the output (hidden state) from the previous step as an additional input to the next. 

  - This creates a chain-like dependency that allows the network to retain information from the past. Conceptually, the RNN maintains a “memory” through a feedback loop: its output at each step is fed back as part of the input at the next step. This recurrent structure makes RNNs a natural architecture for time series data and sequential data more broadly.

- Formally, a basic RNN can be described as follows: at each time ``t``, the network takes the current input ``x_t`` (which could be the value of the time series at ``t`` or a vector of observations at time ``t``) and the previous hidden state ``h_{t−1}``, and computes a new hidden state ``h_t = f(Wx_t +Uh_{t−1} +b)`` for some nonlinear activation function f (like ``\tanh``) and weight matrices ``W,U``. 

  - The hidden state ``h_t`` can be thought of as an encoding of all relevant information seen up to time ``t``. The network may also produce an output ``y_t``, i.e., typically a prediction, based on ``h_t``.

- Through training on sequential data, the RNN learns to update its hidden state in a way that captures patterns in the sequence and to use that to make predictions.
"""

# ╔═╡ 64ed9b01-8c08-4a4c-8d12-ad8b4ab3895f
md"""
- RNNs were actually introduced decades ago, but only became truly effective with the advent of better training techniques and more data. 

  - One key development was the [Long Short-Term Memory](https://en.wikipedia.org/wiki/Long_short-term_memory) (LSTM) network, introduced to address the ”vanishing gradient” problem that made training basic RNNs difficult for long sequences. 

- As a matter of fact, standard RNNs tend to have trouble learning long-term dependencies. LSTMs (and variants as the [Gated Recurrent Unit](https://en.wikipedia.org/wiki/Gated_recurrent_unit), mitigate this by using a special gated architecture with separate mechanisms (gates) to control what information is kept, forgotten, or output at each time step. 

"""

# ╔═╡ b65614f7-60a5-4ff7-a8ad-a90e2c71bc99
md"""
- In practice, RNNs (especially LSTMs/GRUs) became the workhorse for many time series and sequence modeling tasks in the 2010s. 

  - For example, in forecasting, an LSTM can be trained to ingest a sequence of past values and directly output a sequence of future values (sequence-to-sequence forecasting). 

  - In anomaly detection, an LSTM-based autoencoder can be trained to reconstruct normal time series sequences and flag those with large reconstruction error as anomalies.

   - In astronomy, RNNs have been used for classifying light curves of variable stars and transients. 

   - RNNs do not assume linearity or stationarity and can, in principle, approximate very complex functions, they can model phenomena that classical models cannot. 

       - For instance, an LSTM could learn the irregular pattern of a star that has multiple outburst states (something a single ARIMA model would struggle with). 

   - RNNs have also been applied in cosmology, for example to emulate expensive physics simulations by learning from time series of simulation data, or to predict the evolution of cosmological parameters.
"""

# ╔═╡ c9f20dcb-3bbd-4fd0-a672-b86ce2f14ec5
md"""

#### Single-step and multi-step forecasting.

- An important practical consideration for any forecasting model, but especially
relevant for deep learning, is the distinction between *single-step* and *multi-step* forecasting. 

- Given observations up to time ``t``, the forecast horizon ``H`` is the number of future time steps we wish to predict. 

  - In single-step forecasting (``H = 1``) the model predicts only the next time step ``y_{t+1}``. 

- For predictions further into the future (``H > 1``) there are three main strategies: 

    1. the recursive (or iterated) strategy, where the model’s own prediction is fed back as input to generate the next, cascading errors as we move further along the horizon; 
    2. the direct strategy, where a separate model is trained for each lead time ``h = 1,...,H``, avoiding error accumulation but ignoring dependencies between different lead times; 
    3. the multi-output (or sequence-to-sequence) strategy, where a single model outputs the entire forecast vector (``y_{t+1},..., y_{t+H}``) at once. The multi-output approach, naturally suited to encoder-decoder architectures, has become the dominant paradigm in deep learning for time series, as it balances computational efficiency with the ability to capture inter-horizon dependencies.

"""

# ╔═╡ f9ae178e-9219-4efa-9764-dcb670f299c5
md"""

### Convolutional Networks for Time Series
***

- Convolutional neural networks are well-known for image and signal processing tasks, but they can also be applied to one-dimensional time series. A 1D CNN slides convolutional filters over the sequence to detect local patterns. 

- It was showed that a [Temporal Convolutional Network](https://en.wikipedia.org/wiki/Temporal_network) (TCN), which is a 1D CNN architecture with causal convolutions (no leaking from future to past) and dilations (skipping inputs to exponentially increase receptive field), outperformed canonical RNNs like LSTMs on a suite of sequence modeling tasks. 

  - The TCN was able to achieve longer effective memory and better accuracy, suggesting that convolutional architectures can serve as a viable alternative to RNNs for many problems.

- The key aspects of the TCN are: 

   1. Causal convolutions ensure the model is feed-forward in time (outputs at time ``t`` only depend on inputs up to ``t``), 
   2. Dilated convolutions allow the receptive field to grow exponentially with depth, meaning a relatively shallow network can model very long sequences. Residual connections are also used to ease training of deep networks. 

- CNN-based models have beed used for instance in healthcare time series such as patient vital-sign monitoring to detect patterns of disease onset. This technique has also been applied to traffic flow forecasting or electricity load forecasting. 

- CNNs, sometimes in combination with RNNs, have achieved strong results by capturing local trends and repeating patterns.

- In astronomy, convolutional networks have proven particularly effective. 1D CNNs applied directly to LIGO strain time series can detect gravitational-wave signals from compact binary coalescences with accuracy comparable to matched filtering, but orders of magnitude faster, a result with significant implications for real-time detection pipelines. 

- CNNs have also been employed for classifying periodic variable star light curves by learning shape features automatically, and for real-time transient detection in survey data streams.

- The success of convolutional approaches shows that explicitly sequential processing (like RNN) is not the only way to model sequences – learnable filters can capture temporal structure efficiently, and their parallelizability makes them well-suited to the high data rates of modern astronomical surveys.
"""

# ╔═╡ 23650722-09c9-4df1-b277-fa388d3fe160
md"""

### Attention Mechanisms and Transformers
***

- The latest development in sequence modeling has been the rise of [Transformers](https://en.wikipedia.org/wiki/Transformer_(deep_learning)), which forgo both recurrence and convolution in favor of a mechanism called self-[attention](https://en.wikipedia.org/wiki/Attention_(machine_learning)). 

- The Transformer architecture, introduced in the context of natural language processing, has since become state-of-the-art in Natural Language Processing (NLP) and is making inroads in time series analysis as well.

- The core idea of self-attention is to allow each element of a sequence to attend to (i.e., compute weighted interactions with) every other element, capturing long-range dependencies directly. 

- Given an input sequence, each element is projected into three vectors – a query ``\mathbf{q}``, a key ``\mathbf{k}``, and a value ``\mathbf{v}`` – and the attention output is computed as:

```math
{\rm Attention}(\mathbf{Q},\mathbf{K},\mathbf{V}) = {\rm softmax}\left( \frac{\mathbf{Q} \mathbf{K}^T}{\sqrt{d_k}} \right) \mathbf{V})
```

- where ``d_k`` is the dimension of the key vectors and the softmax ensures that the attention weights sum to one. The scaling factor ``\sqrt{d_k}`` prevents the dot products from becoming too large, which would push the softmax into regions of very small gradients.

- In a Transformer model, the sequence is processed as a whole (not one step at a time as in an RNN), and multiple attention heads operate in parallel to capture different types of relationships.

- The position in the sequence is encoded via positional encodings since the model itself has no inherent notion of order. 
  - For irregularly sampled series, time-aware positional encodings and continuous-time attention variants extend these models to non-uniform grids.

- For time series forecasting, Transformers offer the appealing ability to look at very long input histories and identify which past time points are most relevant to predicting the future. This can be very useful for data with long-term seasonal effects or irregular but long-range dependencies.

- Transformers have been successfully applied to a variety of time series tasks, including univariate and multivariate forecasting, anomaly detection (where a Transformer can learn normal patterns and identify deviant behavior), and classification (mapping an entire sequence to a class, useful in medicine for classifying EEG signals, or in astronomy for classifying light curves).

"""

# ╔═╡ 57ba7408-bd56-4227-a58d-effa843d64e4
md"- In the plot below I added limits on the `y`-axis to better center the bandpassed strain data."

# ╔═╡ 5a0cc1ac-1ecb-43fa-a687-048ae9f6cc56
md"- Now, I propose to the AI to better clean the data by a *whitening* procedure, i.e. removing (oart of) the noise to make the signal more evident."

# ╔═╡ cd3f0947-78b7-4593-b7b7-ef82d7a15273
md"""

- Deep learning models, while powerful, also bring challenges: they are often less interpretable than statistical models (though attention weights can sometimes be interpreted, and there is work on explaining neural forecasts), and they require careful tuning and significantly more data to avoid overfitting.

  - In scientific domains like astrophysics, where understanding the model can be (at least) as important as predictive accuracy, this is an important consideration.
"""

# ╔═╡ 5fa720ff-19da-4e22-8972-eb87e3c3d382
md"""

### Foundation Models for Time Series
***

- A rapidly emerging trend is the development of foundation models for time series, i.e., large, pre-trained models that can generalize across diverse tasks and domains with minimal or no task-specific fine-tuning. 

- Inspired by the success of [Large Language Models](https://en.wikipedia.org/wiki/Large_language_model) (LLMs) in natural language processing, several groups have recently proposed pre-trained time series models.
  - These models can perform zero-shot forecasting on unseen time series, often achieving competitive accuracy with task-specific models that were trained on the target data. 

- The appeal of foundation models lies in their potential to democratize time series analysis: a single pre-trained model could be applied out-of-the-box to forecasting tasks across astronomy, finance, healthcare, and engineering without requiring domain-specific training data or feature engineering. 

- However, important open questions remain regarding their robustness to distribution shifts, their ability to handle highly irregular or multivariate scientific data, and whether they can match the performance of carefully tuned domain-specific models. 

  - In astrophysics, where data have unique noise properties and cadences, the applicability of these general-purpose models is an active area of investigation.
"""

# ╔═╡ 1f2ad00a-fc62-4314-bede-7b94c0ad1f1a
md"""
## Reference & Material
***

Material and papers related to the topics discussed in this lecture.

- [Pagliaro & Anzalone (2026) - "Time Series Analysis in Machine Learning”](https://ui.adsabs.harvard.edu/abs/2026arXiv260611746P/abstract)
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
    <td><a href="./open?path=Lectures/Lecture-AIinteraction/Lecture-AI.jl">Lecture about AI interaction</a></td>
	<td><a href="./open?path=Course.jl">Course Summary</a></td>    
  </tr>
  <tr>
    <td>html</td>
	<td><a href="../../Lectures/Lecture-GaussianProcesses/Lecture-CO2.html">Science case about CO₂ content in atmosphere</a></td>    
    <td><a href="../../Lectures/Lecture-AIinteraction/Lecture-AI.html">Lecture about  AI interaction</a></td>
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
md"Notebook v1.0.0 - 4 September 2026"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoTeachingTools = "~0.4.7"
PlutoUI = "~0.7.80"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.0"
manifest_format = "2.1"
project_hash = "b4c18bbc242930204ae5cca25c12954846038fac"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
registries = "General"
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
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

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
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
registries = "General"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
registries = "General"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
registries = "General"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
registries = "General"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
registries = "General"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c0c9b76f3520863909825cbecdef58cd63de705a"
registries = "General"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.5+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
registries = "General"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
registries = "General"
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

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.18.0+1"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.1+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.103+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
registries = "General"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
registries = "General"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.8.13"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.46.0+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "90b41ced6bacd8c01bd05da8aed35c5458891749"
registries = "General"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.7"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "fbc875044d82c113a9dee6fc14e16cf01fd48872"
registries = "General"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.80"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
registries = "General"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

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
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
registries = "General"
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
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
registries = "General"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

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

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.67.1+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.2+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╟─d5b55ac6-0f38-4603-aa72-bd964d0604a8
# ╟─355d961d-ad28-479c-a0a1-1865bc809c99
# ╟─2f82cbd2-4c34-11f1-901e-01c419d48099
# ╟─efc84c48-a322-4c98-869d-2fe1e32151d5
# ╟─72168e90-1440-42ac-a859-834ec6dcc703
# ╟─9b66dbc6-7dd9-4762-a791-7fe1e433961b
# ╟─fbe81ddf-836e-4e2e-a6e1-475257301c41
# ╟─71bcd5bc-afb3-48cc-9185-41e44d433ff9
# ╟─fb8ce658-f20a-410e-aadc-b06e4eba6c51
# ╟─7066d72c-6f84-4311-9a12-a05f897731e7
# ╟─2e9b26ef-78eb-4df9-a419-a099cf285719
# ╟─47ec8223-3159-433e-b91f-10a540edb5bb
# ╟─5c01fb12-0d8d-4978-b080-defa6e881060
# ╟─120af7fd-32e9-4936-950d-e55955380ad0
# ╟─d6330ed7-457e-40e0-9f92-94bd3a75c46d
# ╟─64ed9b01-8c08-4a4c-8d12-ad8b4ab3895f
# ╟─b65614f7-60a5-4ff7-a8ad-a90e2c71bc99
# ╟─c9f20dcb-3bbd-4fd0-a672-b86ce2f14ec5
# ╟─f9ae178e-9219-4efa-9764-dcb670f299c5
# ╟─23650722-09c9-4df1-b277-fa388d3fe160
# ╟─57ba7408-bd56-4227-a58d-effa843d64e4
# ╟─5a0cc1ac-1ecb-43fa-a687-048ae9f6cc56
# ╟─cd3f0947-78b7-4593-b7b7-ef82d7a15273
# ╟─5fa720ff-19da-4e22-8972-eb87e3c3d382
# ╟─1f2ad00a-fc62-4314-bede-7b94c0ad1f1a
# ╟─ddb1955a-d1e7-4452-8bd6-a80913e39487
# ╟─24385054-0edb-4b78-8cce-c89426964dfb
# ╟─e44397bf-7be0-4c6f-bf38-057df876fb92
# ╟─8ed11f8c-4810-433c-9c4f-bbd8cc5f23e9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
