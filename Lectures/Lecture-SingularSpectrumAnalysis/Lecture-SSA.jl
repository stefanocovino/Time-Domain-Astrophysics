### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 6a1315d1-9a6d-4ce0-b1c0-3fe22beb1ec2
begin
	using CairoMakie
	using CommonMark
	using LaTeXStrings
	using LinearAlgebra
	using PlutoUI
	using Random
	using Statistics
end

# ╔═╡ cd347311-bb1e-4404-a93c-3ac7a37e04e2
module SSAJL

using LinearAlgebra
using Statistics
using CairoMakie
# We will use Vector{<:Real} as the general type for the time series
# for simplicity, similar to a NumPy array or list.

export SSA, components_to_matrix, reconstruct, calc_wcorr!, plot_wcorr

# --- Struct for SSA Data ---
# Corresponds to the SSA class's instance variables
struct SSA{T<:Real}
    # User inputs
    L::Int            # Window length
    orig_TS::Vector{T} # Original time series

    # Derived dimensions
    N::Int            # Length of time series
    K::Int            # Number of column vectors in trajectory matrix (N - L + 1)
    d::Int            # Rank of the trajectory matrix

    # Decomposition results
    X::Matrix{T}      # Trajectory Matrix
    U::Matrix{T}      # Left singular vectors
    Sigma::Vector{T}  # Singular values
    V::Matrix{T}      # Right singular vectors (transposed V in Python's SVD is VT)

    # Time Series Components (Diagonally Averaged)
    TS_comps::Matrix{T} # Array where each column is an elementary component F_i

    # W-Correlation Matrix
    Wcorr::Matrix{Float64}

    # Elementary Matrices (only if save_mem=false)
    X_elem::Union{Array{T, 3}, String}
end

"""
    SSA(tseries::Vector{T}, L::Int; save_mem::Bool=true) where T<:Real

Decomposes the given time series using Singular Spectrum Analysis.

# Arguments
- `tseries`: The original time series as a `Vector{<:Real}`.
- `L`: The window length. Must be an integer `2 <= L <= N/2`, where `N` is the length of the time series.
- `save_mem`: Conserve memory by not retaining the elementary matrices. Defaults to `true`.

# Returns
- An `SSA` struct containing the decomposition results.
"""
function SSA(tseries::Vector{T}, L::Int; save_mem::Bool=true) where T<:Real
    N = length(tseries)

    # Input checks
    if !(2 <= L <= N/2)
        throw(DomainError(L, "The window length L must be in the interval [2, N/2]."))
    end

    K = N - L + 1
    orig_TS = tseries

    # 1. Embedding / Trajectory Matrix Construction (X)
    # The Python code does: X = np.array([orig_TS.values[i:L+i] for i in range(0, K)]).T
    # This creates an L x K matrix where column j is the vector (tseries[j], tseries[j+1], ..., tseries[j+L-1])
    X = Matrix{T}(undef, L, K)
    for j in 1:K # Column index (Julia is 1-based)
        # The segment starts at index j and ends at j + L - 1
        X[:, j] = orig_TS[j:j+L-1]
    end

    # 2. Singular Value Decomposition (SVD)
    # Julia's svd(X) returns SVD(U, S, V) such that X = U * Diagonal(S) * V'
    # Python's svd returns U, Sigma, VT (V transpose). So Julia's V is Python's V.T
    F = svd(X)
    U, Sigma, V_py_VT = F.U, F.S, F.Vt # V_py_VT is V' (V-transpose) in Julia.

    d = rank(X) # The effective rank of X

    TS_comps = zeros(T, N, d)
    X_elem = "Re-run with save_mem=false to retain the elementary matrices."

    # 3. Diagonal Averaging (Reconstruction)
    # Diagonally average the elementary matrices to get the time series components

    if !save_mem
        # Construct and save all the elementary matrices
        X_elem_array = Array{T, 3}(undef, L, K, d)

        for i in 1:d
            # Elementary matrix X_i = sigma_i * u_i * v_i'
            # Note: Julia's F.Vt[i, :] is the i-th row, which is v_i' (Python's VT[i,:])
            X_i = Sigma[i] * (U[:, i] * V_py_VT[i, :]')
            X_elem_array[:, :, i] = X_i

            # Diagonal Averaging
            # Python's X_rev = X_elem[::-1] reverses rows.
            # Julia's equivalent: reverse(X_i; dims=1)
            X_rev = reverse(X_i; dims=1)

            # The indices for the diagonals range from -(L-1) to (K-1) in Julia's diag
            # In Python, the range is -X_rev.shape[0]+1 to X_rev.shape[1]-1, which is -(L-1) to K-1
            for j in -(L - 1):(K - 1)
                TS_comps[j + L, i] = mean(diag(X_rev, j)) # j+L is used to map the diagonal index j to 1-based index 1 to N
            end
        end
        X_elem = X_elem_array
    else # save_mem = true
        for i in 1:d
            # Reconstruct the elementary matrix without storing it
            X_i = Sigma[i] * (U[:, i] * V_py_VT[i, :]')

            # Diagonal Averaging
            X_rev = reverse(X_i; dims=1)

            for j in -(L - 1):(K - 1)
                # The mean of the diagonal is assigned to the corresponding element of the reconstructed series.
                # The index mapping in the Python code is somewhat complex, let's use the standard formula.
                # The diagonal index j maps to the time index t = j + L (1-based) for -L+1 <= j <= K-1.
                TS_comps[j + L, i] = mean(diag(X_rev, j))
            end
        end
    end

    # 4. Calculate W-Correlation Matrix (Wcorr)
    # The calc_wcorr logic is implemented as a separate function, then called.
    initial_wcorr = calc_wcorr(L, K, d, TS_comps)

    return SSA{T}(L, orig_TS, N, K, d, X, U, Sigma, V_py_VT', TS_comps, initial_wcorr, X_elem)
end

# --- W-Correlation Logic ---

function calc_wcorr(L::Int, K::Int, d::Int, TS_comps::Matrix{T}) where T<:Real
    # Calculate the weights w
    # Python: list(np.arange(self.L)+1) + [self.L]*(self.K-self.L-1) + list(np.arange(self.L)+1)[::-1]
    #println(L," ",K," ",d," ",K - L + 1 - L)
    w1 = collect(1:L)
    #w2 = fill(L, K - L + 1 - L) # (K - L - 1) is the Python length, which means K - L is the Julia length.
    if K - L >= 1
        w2 = fill(L, K - L) # Center part
    else
        w2 = Int[]
    end
    w3 = reverse(collect(1:L-1)) # The center part takes care of L

    # Julia's indexing for time series of length N
    w_julia = vcat(collect(1:L), fill(L, max(0, K-L)), reverse(collect(1:L-1)))

    # The length of the weights vector must be N (L + K - 1)
    N = L + K - 1

    # The weights vector w has length N
    w_final = zeros(Int, N)
    for i in 1:N
        if i <= L
            w_final[i] = i
        elseif i <= K
            w_final[i] = L
        else
            w_final[i] = N - i + 1
        end
    end

    # Weighted Inner Product function
    function w_inner(F_i::Vector{T}, F_j::Vector{T}) where T<:Real
        return dot(w_final, F_i .* F_j)
    end

    # Calculated weighted norms, ||F_i||_w, then invert.
    F_wnorms = [w_inner(TS_comps[:, i], TS_comps[:, i]) for i in 1:d]
    F_wnorms_inv_sqrt = F_wnorms .^ -0.5

    # Calculate Wcorr.
    Wcorr = Matrix{Float64}(I, d, d)
    for i in 1:d
        for j in i+1:d
            val = abs(w_inner(TS_comps[:, i], TS_comps[:, j]) * F_wnorms_inv_sqrt[i] * F_wnorms_inv_sqrt[j])
            Wcorr[i, j] = val
            Wcorr[j, i] = val
        end
    end
    return Wcorr
end

# This function should be used *after* the SSA struct is created, to update Wcorr if needed.
function calc_wcorr!(ssa::SSA)
    ssa.Wcorr = calc_wcorr(ssa.L, ssa.K, ssa.d, ssa.TS_comps)
    return ssa.Wcorr
end

# --- Methods (External Functions) ---

"""
    components_to_matrix(ssa::SSA, n::Int=0)

Returns a matrix of the first `n` time series components.
If `n` is 0, all components are returned.
"""
function components_to_matrix(ssa::SSA, n::Int=0)
    n = (n > 0) ? min(n, ssa.d) : ssa.d
    return ssa.TS_comps[:, 1:n]
end


"""
    reconstruct(ssa::SSA, indices::Union{Int, AbstractVector{Int}, AbstractRange{Int}})

Reconstructs the time series from its elementary components, using the given indices.
"""
function reconstruct(ssa::SSA, indices::Union{Int, AbstractVector{Int}, AbstractRange{Int}})
    if isa(indices, Int)
        indices = [indices]
    end

    # Indices are 1-based in Julia
    ts_vals = sum(ssa.TS_comps[:, indices], dims=2)[:]
    return ts_vals
end

# --- Plotting with CairoMakie ---

"""
    plot_wcorr(ssa::SSA; min_idx::Int=1, max_idx::Int=ssa.d)

Plots the w-correlation matrix for the decomposed time series using CairoMakie.
"""
function plot_wcorr(ssa::SSA; min_idx::Int=1, max_idx::Int=ssa.d, ptitle="")

    # Ensure indices are within bounds
    min_idx = max(1, min_idx)
    max_idx = min(ssa.d, max_idx)

    # Extract the relevant part of the Wcorr matrix
    Wcorr_subset = ssa.Wcorr[min_idx:max_idx, min_idx:max_idx]

    f = Figure()
    ax = Axis(f[1, 1],
              xlabel = L"$\tilde{F}_i$",
              ylabel = L"$\tilde{F}_j$",
              xticks = (1:size(Wcorr_subset, 1), ["$i" for i in min_idx:max_idx]),
              yticks = (1:size(Wcorr_subset, 2), ["$j" for j in min_idx:max_idx]),
              title = ptitle,
              aspect = 1
              )

    # Plot the heatmap
    CairoMakie.heatmap!(ax, Wcorr_subset, colorrange=(0, 1))

    # Add a color bar
    Colorbar(f[1, 2],
             #colormap = Makie.default_colormap[], # Use default colormap
             limits = (0, 1),
             label = L"$W_{i,j}$"
             )

    ax.yreversed = true

    return f
end

end # module SSAJL

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
# Singular-Spectrum Analysis
***

- We introduce now technique of **singular-spectrum analysis (SSA)**. In simple words, SSA decomposes a time series into a set of summable components that are grouped together and interpreted as *trend*, *periodicity* and *noise*. 

- SSA emphasises **separability** of the underlying components, and can readily separate periodicities that occur on different time scales, even in very noisy time series data. The original time series is recovered by summing together all of its components.

- SSA can be used to analyse and reconstruct a time series with or without different components as desired. For example, you could apply SSA to:
    1. construct a smoothed version of a time series using a small subset of its components; 
    2. investigate a time series' periodic components to understand the underlying processes that generated the time series;
    3. reconstruct the original time series without its periodic components;
    4. remove all trend and periodic components from the series, leaving just the 'noise', which may be meaningful in and of itself...

- Unlike the commonly used autoregressive integrated moving average (ARIMA) method, SSA makes no assumptions about the nature of the time series, and has just a single adjustable (and easily interpretable) parameter.

- Singular-Spectrum Analysis is connected to the popular exploratory and dimensionality reduction technique [Principal Component Analysis](https://en.wikipedia.org/wiki/Principal_component_analysis). An introduction to the PCA can be found here ([notebook](./open?path=Lectures/Lecture-SingularSpectrumAnalysis/Lecture-PCA.jl), [html](Lectures/Lecture-SingularSpectrumAnalysis/Lecture-PCA.html)).

"""

# ╔═╡ 72cf6fcf-2e2f-4dee-994b-ce2bfef51802
md"""
## A Toy Time Series
***

- We'll first define an arbitrary, toy time series, $F = \{f_0, f_1, \ldots, f_{N-1}\}$, containing trend, periodic and noise components:

$$f_t = 0.001 \times (t - 100)^2 + 2\sin(\frac{2\pi t}{p_1}) + 0.75\sin(\frac{2\pi t}{p_2}) + \text{Rand}\{-1,1\}$$

- where $t = \{0, 1,\ldots, N-1\}$ is a discrete time moment, $p_1$ and $p_2$ are set time periods and $\text{Rand}\{-1,1\}$ is a random number uniformly distributed between –1 and 1.

- The $0.001 \times (t - 100)^2$ term defines the (parabolic) trend of the series,  $2\sin(\frac{2\pi t}{p_1})$ and $0.75\sin(\frac{2\pi t}{p_2})$ are two periodic components with differing periodicities and amplitudes, while $\text{Rand}\{-1,1\}$ introduces noise.


"""

# ╔═╡ 0646411b-bbbe-4b6d-bf94-33c4e9522ed5
begin
	N = 200 # The number of time 'moments' in our toy series
	t = 0:N-1
	trend = 0.001 * (t .- 100).^2
	p1, p2 = 20, 30
	periodic1 = 2 * sin.(2π*t/p1)
	periodic2 = 0.75 * sin.(2π*t/p2)
	
	Random.seed!(123)
	noise = 2 .* (rand(N) .- 0.5)
	F = trend .+ periodic1 .+ periodic2 .+ noise
	
	# Plot everything
	fg1 = Figure()
	
	ax1 = Axis(fg1[1, 1],
	    xlabel = "t",
	    ylabel = "F(t)",
	    title = "The Toy Time Series and its Components",
	)
	lines!(t, F, linewidth=2.5, label="Toy Series (F)")
	lines!(t, trend, alpha=0.75, label="Trend")
	lines!(t, periodic1, alpha=0.75, label="Periodic #1")
	lines!(t, periodic2, alpha=0.75, label="Periodic #2")
	lines!(t, noise, alpha=0.5, label="Noise")
	    
	axislegend(ax1,framevisible = false,position = :ct)
	
	
	fg1
end

# ╔═╡ 2f2430a4-ea1f-4373-b5a7-10db4756b29b
md"- Based on visual inspection, the presence of trend (yellow), periodic (green and magenta) and noise (cyab) components in the time series is clear. However, the second periodic component (magenta) in the toy series is not apparent. Could we use SSA to recover the trend, periodic and noise components from the toy time series?"

# ╔═╡ e210c43a-1205-40b8-8989-7975ea85e648
cm"""
## Introducing the SSA Method
***

### From a Time Series to a Trajectory Matrix
***

The first step of SSA is to map the time series ``F`` to a sequence of multi-dimensional lagged vectors. Let an integer ``L`` be the **window length**, ``2 \le L \le N/2``. We form a 'window', given by the subseries ``\{f_i, \ f_{i+1}, \ldots , \ f_{i+L-1}\}``, for  ``i=0,\ldots,N-L``. We slide this window along the time series, forming a column vector, ``X_i``, for each subseries. That is, we have

```math
\begin{align*}
X_0 & = (f_0, \ f_1, \ f_2,  \ldots, \ f_{L-1} )^{\text{T}} \\
X_1 & = (f_1, \ f_2, \ f_3,  \ldots, \ f_L )^{\text{T}} \\
X_2 & = (f_2, \ f_3, \ f_4,  \ldots, \ f_{L+1} )^{\text{T}} \\
X_3 & = (f_3, \ f_4, \ f_5,  \ldots, \ f_{L+2} )^{\text{T}} \\
& \quad \quad \quad  \vdots \\
X_{N-L} & = (f_{N-L}, \ f_{N-L+1}, \ f_{N-L+2}, \ \ldots, \ f_{N-1} )^{\text{T}}.
\end{align*}
```

These column vectors form the ``L``-*trajectory matrix*, ``\mathbf{X}``, of the time series (hereafter just *trajectory matrix*):

```math
\mathbf{X} = \begin{bmatrix}
f_0 & f_1 & f_2 & f_3 &\ldots & f_{N-L} \\ 
f_1 & f_2 & f_3 & f_4 &\ldots & f_{N-L+1} \\
f_2 & f_3 & f_4 & f_5 &\ldots & f_{N-L+2} \\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\
f_{L-1} & f_{L} & f_{L+1} & f_{L+2} & \ldots & f_{N-1} \\ 
\end{bmatrix}
```

From writing out the matrix above, it is clear that the elements of the *anti-diagonals* (that is, the diagonals running from bottom-left to top-right) are equal. This type of matrix is known as a ***Hankel*** matrix.

For our toy time series, we'll set the window length to 70, and defer discussion on how to select an appropriate window length. Let ``K = N - L + 1`` represent the number of columns in the trajectory matrix. **We'll refer to the columns of ``\mathbf{X}`` as the ``L``-lagged vectors, and the rows as ``K``-lagged vectors.**
"""

# ╔═╡ a759da24-db69-402a-bde3-25edbae59834
begin
	L = 70 # The window length.
	K = N - L + 1 # The number of columns in the trajectory matrix.
	# Create the trajectory matrix by pulling the relevant subseries of F, and stacking them as columns.
	as = [F[i:i+L-1] for i in range(1,K)]
	X = reduce(hcat, as)
end;

# ╔═╡ 37551c89-38d4-46e6-9109-52b3a065ed73
begin
	fg2 = Figure()
	
	ax2 = Axis(fg2[1, 1],
	    xlabel = "L-Lagged Vectors",
	    ylabel = "K-Lagged Vectors",
	    title = "The Trajectory Matrix for the Toy Time Series",
	)
	hm = heatmap!(X)
	Colorbar(fg2[:, end+1], hm)
	
	ax2.yreversed = true
	
	fg2
end

# ╔═╡ 99fb1b5c-2857-43a3-b69b-fd55a7e35e71
md"""
### Decomposing the Trajectory Matrix
***

- The second step is decomposing the trajectory matrix with a [singular-value decomposition (SVD)](https://en.wikipedia.org/wiki/Singular-value_decomposition).

```math
\mathbf{X} = \mathbf{U\Sigma V}^{\text{T}}
``` 

- where:
    * ``\mathbf{U}`` is an ``L \times L`` unitary matrix containing the orthonormal set of ***left singular vectors*** of ``\mathbf{X}`` as columns;
    * ``\mathbf{\Sigma}`` is an ``L \times K`` rectangular diagonal matrix containing ``L`` ***singular values*** of ``\mathbf{X}``, ordered from largest to smallest; and
    * ``\mathbf{V}`` is a ``K \times K`` unitary matrix containing the orthonormal set of ***right singular vectors*** of ``\mathbf{X}`` as columns.





- The SVD of the trajectory matrix can be rewritten as 
```math
\begin{align*}
    \mathbf{X} & = \sum_{i=0}^{d-1}\sigma_i U_i V_i^{\text{T}} \\
               & \equiv \sum_{i=0}^{d-1}\mathbf{X}_i
\end{align*}
```
- where ``\sigma_i`` is the ``i``th singular value, ``U_i`` and ``V_i`` are vectors representing the ``i``th columns of ``\mathbf{U}`` and ``\mathbf{V}``, respectively, ``d \le L`` is the *rank* of the trajectory matrix, and ``\mathbf{X}_i = \sigma_i U_i V_i^{\text{T}}`` is the ``i``th **elementary matrix** of ``\mathbf{X}``. The collection ``\{U_i, \sigma_i, V_i\}`` will be denoted the ``i``th **eigentriple** of the SVD.

- To build a picture of what all of this means, let's inspect the ``\mathbf{U}``, ``\mathbf{V}`` and ``\mathbf{\Sigma}`` matrices in turn. 
"""

# ╔═╡ b393afda-d52b-4b48-aee0-cf1a42054afe
md"""
#### The ``\mathbf{U}`` Matrix
***

- ``\mathbf{U}`` is an ``L \times L`` matrix whose columns are orthonormal, that is:
```math
    U_i \cdot U_j = \left\{
  \begin{array}{lr}
    1 \ & i = j \\
    0 \ & i \ne j
  \end{array}
\right.
```

- This means that ``\mathbf{UU}^{\text{T}} = \mathbf{U}^{\text{T}}\mathbf{U} = \mathbf{1}``, making ``\mathbf{U}`` a unitary matrix. 

- To elucidate the role that ``\mathbf{U}`` plays in the expansion for ``\mathbf{X}`` above, let ``Z_i = \sigma_i V_i`` be a column vector, so that;
```math
\mathbf{X} = \sum_{i=0}^{d-1} U_i Z_i^{\text{T}}
```

- and each ``L``-lagged column vector, ``X_j``, is then given by:

```math
	X_j = \sum_{i=0}^{d-1}z_{j,i}U_i
```

- where ``z_{j,i}`` is the ``j``th component of the vector ``Z_i``. The expression for ``X_j`` suggests that ``\mathcal{U} = \{U_0, \ldots, U_{d-1} \}`` is a basis set spanning the *column space* of the trajectory matrix, and ``z_{j,i}`` is the ``i``th coefficient of the lagged vector ``X_j`` represented in the basis ``\mathcal{U}``.

- In other words, the columns of the ``\mathbf{U}`` matrix form an orthonormal basis set that describes the time subseries ``\left\{ f_i, \ldots, f_{i+L-1}\right\}_{i=0}^{N-L}`` in the columns of the trajectory matrix. 
"""

# ╔═╡ 61850790-88da-4cf9-8ed3-f96282f04e2d
md"""
#### The ``\mathbf{V}`` Matrix
***

- The matrix ``\mathbf{V}`` is a ``K \times K`` matrix with orthonormal columns, which, like the ``\mathbf{U}`` matrix, makes it unitary.

- To interpret the columns of ``\mathbf{V}`` in the SVD of the trajectory matrix, we first note that for any appropriately shaped matrices ``\mathbf{A}`` and ``\mathbf{B}``, ``\left(\mathbf{AB}\right)^{\text{T}} = \mathbf{B}^{\text{T}}\mathbf{A}^{\text{T}}``. Taking the transpose of ``\mathbf{X}``, we therefore have
```math
\begin{align*}
\mathbf{X}^{\text{T}} & = \mathbf{V \Sigma}^{\text{T}}\mathbf{U}^{\text{T}} \\
                      & = \sum_{i=0}^{d-1}V_i Y_i^{\text{T}}
\end{align*}
```
- where we have set ``Y_i = \sigma_i U_i``. Then,

```math
X^{(\text{T})}_j = \sum_{i=0}^{d-1}y_{j,i}V_i
``` 

- where ``X^{(\text{T})}_j`` is the ``j``th column of ``\mathbf{X}^{\text{T}}``, and ``y_{j,i}`` is the ``j``th component of the vector ``Y_i``. This expression suggests that the ``\mathcal{V} = \{V_0, \ldots, V_{d-1}\}`` is a basis set spanning the column space of ``\mathbf{X}^{\text{T}}``, and ``y_{j,i}`` is the ``i``th coefficient of the lagged vector ``X^{(\text{T})}_j`` represented in the basis ``\mathcal{V}``.

- Equivalently, ``\mathcal{V}`` is a basis set spanning the *row space* of ``\mathbf{X}``. That is, the columns of the ``\mathbf{V}`` matrix form an orthonormal basis set that describe the time subseries ``\{ f_i, \ldots, f_{i+N-L}\}_{i=0}^{L-1}`` in the rows of the trajectory matrix.
"""

# ╔═╡ 35390351-c3d9-4218-aa8a-02cabda3c192
md"""
#### The ``\mathbf{\Sigma}`` Matrix
***

- The ``\mathbf{\Sigma}`` matrix is an ``L \times K`` rectangular diagonal matrix containing the *singular values* of ``\mathbf{X}``. The singular values are ordered from largest to smallest, i.e. ``\sigma_0 \ge \sigma_1 \ge \ldots \ge \sigma_{L-1} \ge 0``.

> We can interpret ``\sigma_i`` as a scaling factor that determines the relative importance of the eigentriple ``(U_i, \sigma_i, V_i)`` in the expansion ``\mathbf{X} = \sum_{i=0}^{d-1}\sigma_i U_i V_i^{\text{T}}``.

- The *Frobenius norm* of ``\mathbf{X}``, ``\lvert\lvert \mathbf{X} \rvert\rvert_{\text{F}}``, is given by
```math
\lvert\lvert \mathbf{X} \rvert\rvert_{\text{F}} = \sqrt{\sum_{j=0}^{L-1}\sum_{k=0}^{K-1} \lvert x_{j,k}\rvert^2}
```
- where ``x_{j,k}`` denotes the element in the ``j``th row and ``k``th column of ``\mathbf{X}``.

- Let's turn our attention to the elementary matrices ``\mathbf{X}_i = \sigma_i U_i V_i^{\text{T}}``. Now, for an outer product such as ``U_i V_i^{\text{T}}``, we have  ``\lvert \lvert U_i V_i^{\text{T}} \rvert \rvert_{\text{F}} = \lvert \lvert U_i \rvert \rvert_{\text{F}} \lvert \lvert V_i \rvert \rvert_{\text{F}}``, which is simply equal to 1 due to ``U_i`` and ``V_i`` being normalised.

    - From this result, it is then clear that ``\lvert\lvert \mathbf{X}_i \rvert\rvert_{\text{F}} = \sigma_i``.

- It also turns out that 
```math
\lvert\lvert \mathbf{X} \rvert\rvert_{\text{F}}^2 = \sum_{i=0}^{d-1} \sigma_i^2
```
- i.e. the squared Frobenius norm of the trajectory matrix is equal to the sum of the squared singular values. This suggests that we can take the ratio ``\sigma_i^2 / \lvert\lvert \mathbf{X} \rvert\rvert_{\text{F}}^2`` as a measure of the contribution that the elementary matrix ``\mathbf{X}_i`` makes in the expansion of the trajectory matrix.

- Further, if we right-multiply the original SVD of ``\mathbf{X}`` by ``\mathbf{X}^{\text{T}}``:
```math
\begin{align*}
    \mathbf{XX}^{\text{T}} & = \mathbf{U\Sigma V}^{\text{T}}\mathbf{X}^{\text{T}} \\
               & = \mathbf{U\Sigma V}^{\text{T}} \mathbf{V \Sigma}^{\text{T}}\mathbf{U}^{\text{T}} \\
               & = \mathbf{U\Sigma} \mathbf{\Sigma}^{\text{T}}\mathbf{U}^{\text{T}}
\end{align*}
```

- Letting the square diagonal matrix ``\mathbf{\Sigma}^2 = \mathbf{\Sigma \Sigma}^{\text{T}}``, and multiplying on the right by ``\mathbf{U}``, gives:

```math
(\mathbf{XX}^{\text{T}})\mathbf{U} = \mathbf{U}\mathbf{\Sigma}^2
```

- which, given that ``\mathbf{\Sigma}^2`` is a diagonal matrix with elements ``\sigma_i^2``, demonstrates that the columns of ``\mathbf{U}`` are eigenvectors of the matrix ``\mathbf{XX}^{\text{T}}``, with eigenvalues ``\{\sigma_0^2, \ldots , \sigma_{L-1}^2\}``.

- Following a similar argument, multiplying ``\mathbf{X}`` on the left by ``\mathbf{X}^{\text{T}}`` shows that the columns of ``\mathbf{V}`` are eigenvectors of the matrix ``\mathbf{X}^{\text{T}}\mathbf{X}``, also with eigenvalues ``\{\sigma_0^2, \ldots , \sigma_{L-1}^2\}``.

"""

# ╔═╡ 86d2e08b-bfab-46a3-9db6-6ef8a22147b2
md"""
#### Putting it all Together
***

- Let's quickly recap everything so far: we have mapped a time series ``F = \{f_0, \ldots, f_{N-1}\}`` to a collection of multi-dimensional lagged vectors, ``X_i = (f_i, f_{i+1}, \ldots, f_{i+L-1})^{\text{T}}, i = 0, \ldots, N-L``, which together comprise the columns of the trajectory matrix ``\mathbf{X}``.

- We then decomposed this matrix with an SVD; in doing so, we found two orthonormal basis sets, ``\mathcal{U}`` and ``\mathcal{V}``, which span the column- and row-space, respectively, of the trajectory matrix. The SVD of ``\mathbf{X}`` can be written as
```math
\begin{align*}
    \mathbf{X} & = \sum_{i=0}^{d-1}\sigma_i U_i V_i^{\text{T}} \\
               & \equiv \sum_{i=0}^{d-1}\mathbf{X}_i
\end{align*}
```
- where ``\mathbf{X}_i`` is the ``i``th elementary matrix of ``\mathbf{X}``, determined by the eigentriple  ``\{U_i, \sigma_i, V_i\}``. The ``i``th singular value, ``\sigma_i``, determines the relative contribution of ``\mathbf{X}_i`` in the expansion of ``\mathbf{X}`` above.

- The integer ``d \le L`` is the intrinsic dimensionality of the time series' trajectory space, and we may choose to obtain a lower-dimensional approximation of ``\mathbf{X}`` by summing only the first ``r < d`` elementary matrices.

> Much of what has been covered so far is general to the SVD of *any* matrix, not just the trajectory matrix of a time series.

- From now on, we are going to focus on reconstructing the components of a time series from its elementary matrices.
"""

# ╔═╡ 490e3a2a-cf23-42e3-9530-3a6fad82bb30
begin
	d = rank(X)   # The intrinsic dimensionality of the trajectory space.
	
	U, Σ, V = svd(X)
	
	# Calculate the elementary matrices of X, storing them in a multidimensional array.
	# This requires calculating sigma_i * U_i * (V_i)^T for each i, or sigma_i * outer_product(U_i, V_i). 
	# Note that Sigma is a 1D array of singular values, instead of the full L x K diagonal matrix.
	X_elem = [Σ[i] * U[:, i] * (V[:, i])' for i in 1:d]
	
	# Quick sanity check: the sum of all elementary matrices in X_elm should be equal to X, to within a 
	# *very small* tolerance:
	if !isapprox(X, sum(X_elem), atol=1e-10)
	    println("WARNING: The sum of X's elementary matrices is not equal to X!")
	end
end

# ╔═╡ aaebbbcd-6bb6-46e0-8564-8029a4f2d1cf
md"- Let's take a peak at the first 12 elementary matrices:"

# ╔═╡ 2b3f952d-db20-4790-8482-9c6e981f60e0
begin
	function plot_2d(fig_position::Makie.GridPosition, m::AbstractMatrix, title::LaTeXString)
	    # 1. Create the Axis at the specified grid position
	    ax = Axis(fig_position, 
	              title = title,
	              xticksvisible = false, 
	              yticksvisible = false,
	              xticklabelsvisible = false,
	              yticklabelsvisible = false,
	              aspect = 1
	             )
	
	    # 2. Plot the matrix as a heatmap
	    Makie.heatmap!(ax, m)
	
	    ax.yreversed = true
	    
	    return ax
	end
	
	
	function plot_series_in_grid(X_elem, n::Int)
	    # X_elem is assumed to be an array/vector of 2D matrices/arrays
	    # n is the number of elements to plot (e.g., n = min(12, d))
	    
	    f = Figure(size = (800, 800)) # Set a reasonable size for a 4x4 grid
	
	    # Grid dimensions for subplot(4, 4, i)
	    rows = 4
	    cols = 4
	
	    for i in 1:n
	        r = cld(i, cols) # Ceiling division for row
	        c = i - (r - 1) * cols
	        
	        #title_text = L"\mathbf{X}_{%$(i-1)}$" 
	        #title_text = "X($i)"
	        title_text = L"\mathbf{X}_{%$i}"
	        
	        plot_2d(f[r, c], X_elem[i], title_text)
	    end
	    
	    rowgap!(f.layout, 5) # Small vertical gap
	    colgap!(f.layout, 5) # Small horizontal gap
	
	    return f
	end
end

# ╔═╡ ef9e9f9a-7598-4e22-9add-f1dd2a74c516
begin
	n = min(12, d) # In case d is less than 12 for the toy series. Say, if we were to exclude the noise component.
	plot_series_in_grid(X_elem, n)
end

# ╔═╡ 3d354f7c-ea01-4465-9b0e-b911801d742b
md"""
- From visual inspection of the ``\mathbf{X}_i`` above, it is obvious that the elementary matrices lack the anti-diagonal structure of the trajectory matrix.

- Even without inspecting the ``U_i`` vector associated with each ``\mathbf{X}_i``, or reconstructing a time series of each component, the appearance of the ``\mathbf{X}_i`` hints at the nature of each component, be it trend, periodicity or noise.

    - For example, the ``L``- and ``K``-lagged vectors in ``\mathbf{X}_0`` and ``\mathbf{X}_1`` vary relatively slowly across the matrix, suggesting that ``\mathbf{X}_0`` and ``\mathbf{X}_1`` may be associated with the overall trend in the time series.
    - The matrices ``\mathbf{X}_2`` to ``\mathbf{X}_5`` show large checkerboard patterns, suggesting periodicity.
    - ``\mathbf{X}_6`` may lie somewhere between periodicity and trend.
    - The matrices of ``\mathbf{X}_7`` onwards (and all the way up to ``\mathbf{X}_{69}``) appear to alternate quickly between a few values; these elementary matrices are likely to be associated with the noise in the original time series.

- Let's plot now the relative contributions, ``\dfrac{\sigma_i^2}{\sum_{k=0}^{d-1} \sigma_k^2}``, and the cumulative contributions, ``\dfrac{\sum_{j=0}^i \sigma_j^2}{\sum_{k=0}^{d-1} \sigma_k^2}``, of the first 12 elementary matrices to the trajectory matrix of the toy time series:
"""

# ╔═╡ 54273740-9fb0-46b6-a928-c74b9077eaa3
begin
	sigma_sumsq = sum(Σ.^2)
	fg3 = Figure(size=(1000,600))
	
	ax13 = Axis(fg3[1, 1],
	    xlabel = L"$i$",
	    ylabel = "Contribution (%)",
	    title = L"Relative Contribution of $\mathbf{X}_i$ to Trajectory Matrix",
	    titlesize=20
	)
	lines!(Σ.^2 / sigma_sumsq * 100, linewidth=2.5)
	xlims!(0,11) 
	
	
	ax23 = Axis(fg3[1, 2],
	    xlabel = L"$i$",
	    ylabel = "Contribution (%)",
	    title = L"Cumulative Contribution of $\mathbf{X}_i$ to Trajectory Matrix",
	    titlesize=20
	)
	xlims!(0,11) 
	lines!(cumsum(Σ.^2) / sigma_sumsq * 100, linewidth=2.5)
	
	fg3
end

# ╔═╡ 018667ce-b84d-46f2-9fce-dd890fc9a557
md"""
- The plots above depict the relative and cumulative contributions of the first 12 ``\mathbf{X}_i`` in the expansion ``\mathbf{X} = \sum_{i=0}^{d-1}\mathbf{X}_i``.

- The elementary matrices ``\mathbf{X}_0`` and ``\mathbf{X}_1`` contribute 52% and 22%, respectively, to the expansion of  ``\mathbf{X}``. Together, the first seven elementary matrices contribute 97%.

- Elementary matrices that make equal contributions to the expansion (that is, ``\sigma_i \approx \sigma_{i+1}``) are likely to be grouped together when reconstructing the time series, and appear as "breaks" in the plot of relative contributions. For example, the "breaks" in the plot above suggest that ``\mathbf{X}_2`` and ``\mathbf{X}_3``, and ``\mathbf{X}_4`` and ``\mathbf{X}_5``, should be grouped together. 

> It is important to note that the elementary matrices represent an optimal (although possibly non-unique) separation of components in trajectory space: by definition, the rows and columns of one elementary matrix are orthogonal to the rows and columns of the other elementary matrices. However, this separation may not coincide with what we would consider a useful, interpretable 'component' of the time series. In fact, there are restrictions on the types of time series components that are exactly separable under this formalism.
"""

# ╔═╡ f05b56b9-801c-42e1-9062-4aa3234bc482
md"""
### Reconstructing the Time Series
***

- So far, we have mapped a time series ``F`` to a series of ``L``-lagged vectors, forming the *trajectory matrix* of ``F``.

- We then decomposed this matrix with a singular-value decomposition, and constructed a set of *elementary matrices* which comprise the trajectory matrix.

- We then gave a bit of a hand-waving explanation to classify these elementary matrices as *trend*, *periodicity* and *noise*.

- In a perfect world, all the components of a time series ``F = \sum_j F^{(j)}`` would be separable, and we would have grouped the the resulting elementary matrices ``\mathbf{X}_i`` appropriately, such that:

```math
\begin{align*}
\mathbf{X} &  = \sum_{k \in \mathcal{S}}\mathbf{X}_k + \sum_{l \in \mathcal{T}}\mathbf{X}_l + \ldots \\
             &  = \sum_j \mathbf{X}^{(j)}
\end{align*}
```

- where ``\mathcal{S}`` and ``\mathcal{T}`` are disjoint (i.e. non-overlapping) sets of indices, and ``\mathbf{X}^{(j)}`` is the trajectory matrix of the time series component ``F^{(j)}``.

- In this case, each ``\mathbf{X}^{(j)}`` would have a [*Hankel structure*](https://en.wikipedia.org/wiki/Hankel_matrix) like the original trajectory matrix, and construction of each ``F^{(j)}`` would be simple.

- However, in the real world, no component trajectory matrices will have equal values on their anti-diagonals. Therefore, we seek a process to transform an elementary matrix to a Hankel matrix, and then into a time series.


"""

# ╔═╡ 911757a4-a582-469d-8a2a-9273065e8d83
md"""

- To extract a time series from the elementary matrices, we'll employ ***diagonal averaging***, which defines the values of the reconstructed time series ``\tilde{F}^{(j)}`` as averages of the corresponding anti-diagonals of the matrices ``\mathbf{X}^{(j)}``.

- Formally, this is represented by introducing the *Hankelisation* operator, ``\hat{\mathcal{H}}``, that acts on the ``L \times K`` matrix ``\mathbf{X}^{(j)}`` to give a Hankel matrix ``\mathbf{\tilde{X}}^{(j)}``; that is:
 
```math
\mathbf{\tilde{X}}^{(j)} = \hat{\mathcal{H}}\mathbf{X}^{(j)}
```

- The element ``\tilde{x}_{m,n}`` in ``\mathbf{\tilde{X}}^{(j)}``, for ``s = m+n``, is given by:

```math
\tilde{x}_{m,n} = \left\{
  \begin{array}{lr}
    \frac{1}{s+1}\sum_{l=0}^{s} x_{l, s-l} & \ 0 \le s \le L-1 \\
    \frac{1}{L-1}\sum_{l=0}^{L-1} x_{l, s-l} & \ L \le s \le K-1 \\
    \frac{1}{K+L-s-1}\sum_{l=s-K+1}^{L} x_{l, s-l} & \ K \le s \le K+L-2 \\
  \end{array}
\right.
```

- At first glance, the above looks like an impenetrable soup of matrix indices. However, all it is doing is calculating the given ``\tilde{x}_{m,n}`` by averaging the rest of the elements of the anti-diagonal wherein ``\tilde{x}_{m,n}`` belongs.

- The number of anti-diagonal elements to sum depends on the location of ``m`` and ``n`` in the matrix, and hence the index soup.

> In practice, we don't need the full Hankel matrix ``\mathbf{\tilde{X}}^{(j)}``, and can cut straight to the construction of the time series ``\tilde{F}^{(j)}``. However, we have included the definition of ``\hat{\mathcal{H}}\mathbf{X}^{(j)}`` above to complete the mathematical exposition of SSA.

- It is important to note that ``\hat{\mathcal{H}}`` is a linear operator, i.e.  ``\hat{\mathcal{H}}(\mathbf{A} + \mathbf{B}) = \hat{\mathcal{H}}\mathbf{A} + \hat{\mathcal{H}}\mathbf{B}``. Then, for a trajectory matrix ``\mathbf{X}``:
   
```math
\hat{\mathcal{H}}\mathbf{X} = \hat{\mathcal{H}} \left( \sum_{i=0}^{d-1} \mathbf{X}_i \right) = \sum_{i=0}^{d-1} \hat{\mathcal{H}} \mathbf{X}_i \equiv \sum_{i=0}^{d-1} \tilde{\mathbf{X}_i}
```

- As ``\mathbf{X}`` is already a Hankel matrix, then by definition ``\hat{\mathcal{H}}\mathbf{X} = \mathbf{X}``. Therefore, the trajectory matrix can be expressed in terms of its Hankelised elementary matrices:

```math
\mathbf{X} = \sum_{i=0}^{d-1} \tilde{\mathbf{X}_i}
```

- As a time series is uniquely determined from a Hankel matrix, the expression above also defines the time series ``F`` as a sum of its components ``\tilde{F}_i``. It is up to us to group these components together, and classify them as trend, periodicity or noise, and then we're free to decide how we use them. 
"""

# ╔═╡ 1d406a03-8158-4b3f-a734-d1b80c259847
"""
    Hankelise(X::AbstractMatrix)

Hankelises the matrix X, returning H(X).
"""
function Hankelise(X::AbstractMatrix)
    # Get dimensions. Julia uses size(X) and multiple assignment (similar to Python)
    L, K = size(X)
    transpose_flag = false

    if L > K
        # Transpose the matrix. In Julia, the non-mutating transpose is X'
        # The Hankelisation below only works for matrices where L < K.
        # To Hankelise a L > K matrix, first swap L and K and tranpose X.
        # Set flag for HX to be transposed before returning. 
        X = X'
        L, K = K, L # Swap dimensions
        transpose_flag = true
    end

    # Initialize HX. Julia uses `zeros(type, rows, cols)`. Float64 is the common default.
    HX = zeros(Float64, L, K)

    # Note: Julia uses 1-based indexing for arrays/matrices.
    # The Python range(N) goes from 0 to N-1.
    # The Julia equivalent for iterating indices 1 to N is 1:N or 0:(N-1) + 1.

    # L is rows, K is columns
    for m in 1:L
        for n in 1:K
            # Indices m and n are 1-based here.
            # The sum s from the formula uses 0-based indexing (m-1) + (n-1) in Python.
            # In Julia, m and n are the indices, so the sum of the 0-based indices is:
            s = (m - 1) + (n - 1)
            
            # Case 1: 0 <= s <= L-1 (Python indices)
            if 0 <= s <= L - 1
                # Loop variable l (0-based) ranges from 0 to s
                # Julia 1-based loop: 1 to s+1
                for l_idx in 1:(s + 1) 
                    # l_0 = l_idx - 1 (the 0-based index)
                    l_0 = l_idx - 1
                    
                    # Original indices: X[l, s-l] (Python 0-based)
                    # New indices (Julia 1-based): X[l_0 + 1, (s - l_0) + 1]
                    # l_0 + 1 = l_idx
                    # s - l_0 + 1 = s - (l_idx - 1) + 1 = s - l_idx + 2
                    
                    # Division: s+1 is the number of terms
                    # Note: / in Julia always performs floating-point division
                    HX[m, n] += (1 / (s + 1)) * X[l_idx, s - l_0 + 1]
                end

            # Case 2: L <= s <= K-1 (Python indices)
            elseif L <= s <= K - 1
                # Loop variable l (0-based) ranges from 0 to L-1
                # Julia 1-based loop: 1 to L
                for l_idx in 1:L 
                    l_0 = l_idx - 1 # 0-based index
                    
                    # Division: L-1 is the number of terms
                    HX[m, n] += (1 / (L - 1)) * X[l_idx, s - l_0 + 1]
                end

            # Case 3: K <= s <= K+L-2 (Python indices)
            elseif K <= s <= K + L - 2
                # Loop variable l (0-based) ranges from s-K+1 to L-1
                # Julia 1-based loop: (s-K+1)+1 to (L-1)+1
                start_l = s - K + 2 # (s - K + 1) + 1
                end_l = L
                
                for l_idx in start_l:end_l
                    l_0 = l_idx - 1 # 0-based index
                    
                    # Division: K+L-s-1 is the number of terms
                    HX[m, n] += (1 / (K + L - s - 1)) * X[l_idx, s - l_0 + 1]
                end
            end
        end
    end

    if transpose_flag
        # Julia's non-mutating transpose is the apostrophe operator: '
        return HX'
    else
        return HX
    end
end

# ╔═╡ 38fc4698-f890-4104-b37b-0d938a02f012
begin
	nh = min(12, d) # In case d is less than 12 for the toy series. Say, if we were to exclude the noise component.
	els = [Hankelise(X_elem[i]) for i in 1:nh]
	plot_series_in_grid(els, nh)
end

# ╔═╡ d45f4fbd-8ef0-47db-8278-72cd1230ab5e
cm"""
- Inspection of the Hankelised elementary matrices of the toy time series confirms our suspicians about the elementary matrices:

	- ``\tilde{\mathbf{X}}_0`` and ``\tilde{\mathbf{X}}_1`` vary slowly over the whole time series, and can be grouped together as the trend component.
	- ``\tilde{\mathbf{X}}_2`` and ``\tilde{\mathbf{X}}_3`` are both periodic, with the same frequency, and can be grouped as the first periodic component.
	- ``\tilde{\mathbf{X}}_4`` and ``\tilde{\mathbf{X}}_5`` are also periodic, with a different frequency to ``\tilde{\mathbf{X}}_2`` and ``\tilde{\mathbf{X}}_3``, and will be grouped as the second periodic component.
	- ``\tilde{\mathbf{X}}_6``, lacking obvious periodicity, will be grouped with the trend components.

- We'll lump together all components from ``\tilde{\mathbf{X}}_7`` and beyond as noise.
    
- To summarise:
```math
\begin{align*}
\tilde{\mathbf{X}}^{\text{(trend)}} & = \tilde{\mathbf{X}}_0 + \tilde{\mathbf{X}}_1 + \tilde{\mathbf{X}}_6 
    & \implies &  \tilde{F}^{\text{(trend)}} = \tilde{F}_0 + \tilde{F}_1 + \tilde{F}_6 \\
\tilde{\mathbf{X}}^{\text{(periodic 1)}} & = \tilde{\mathbf{X}}_2 + \tilde{\mathbf{X}}_3 
    & \implies & \tilde{F}^{\text{(periodic 1)}} = \tilde{F}_2 + \tilde{F}_3  \\
\tilde{\mathbf{X}}^{\text{(periodic 2)}} & = \tilde{\mathbf{X}}_4 + \tilde{\mathbf{X}}_5 
    & \implies & \tilde{F}^{\text{(periodic 2)}} = \tilde{F}_4 + \tilde{F}_5\\
\tilde{\mathbf{X}}^{\text{(noise)}} & = \tilde{\mathbf{X}}_7 + \ldots + \tilde{\mathbf{X}}_{69}
    & \implies & \tilde{F}^{\text{(noise)}} = \tilde{F}_7 + \ldots + \tilde{F}_{69}
\end{align*}
```

- While we have defined the time series component grouping in terms of Hankelised elementary matrices, we will no longer calculate the full Hankel matrix ``\tilde{\mathbf{X}}_i``, and instead calculate ``\tilde{F}_i`` directly from ``\mathbf{X}_i``. 
"""

# ╔═╡ 5385afbd-33a9-4dba-af14-dfe98c841611
"""
    X_to_TS(X_i)

Averages the anti-diagonals of the given elementary matrix, X_i, and returns a time series.

This function is a core step in the reconstruction phase of Singular Spectrum Analysis (SSA).
"""
function X_to_TS(X_i::AbstractMatrix)
    # 1. Reverse the row ordering of X_i.
    X_rev = X_i[end:-1:1, :]

    # 2. Iterate through all possible diagonals (k) of the reversed matrix.
    # The diagonals of the reversed matrix correspond to the anti-diagonals of X_i.
    
    # M, N = size(X_i)
    # The range for k is -(M-1) to N-1.
    M = size(X_i, 1)
    N = size(X_i, 2)
    k_start = -(M - 1)
    k_end = N - 1

    # 3. Use `diag(X_rev, k)` to extract the diagonal and `mean()` to average it.
    # The result is collected into a vector.
    return [mean(diag(X_rev, k)) for k in k_start:k_end]
end

# ╔═╡ 18745681-9e0c-46f1-b9fe-d00b61fbafb8
md"- Let's go ahead and construct the first 12 elementary components, $\tilde{F}_i$, for the toy time series."

# ╔═╡ 4c55c72d-3c4e-4fb9-8fa0-4dd6d6259637
begin
	n4 = min(12,d) # In case of noiseless time series with d < 12.
	
	fg4 = Figure()
	
	ax14 = Axis(fg4[1, 1],
	    xlabel = "t",
	    ylabel = L"\tilde{F}_i(t)",
	    title = "The First 12 Components of the Toy Time Series",
	)
	
	pl = []
	pll = []
	for i in 1:n4
	    F_i = X_to_TS(X_elem[i])
	    pli = lines!(t, F_i, linewidth=2,label=L"\tilde{F}_i(t)")
	    push!(pll,L"\tilde{F}_{%$i}(t)")
	    push!(pl,pli)
	end
	
	Legend(fg4[1, 2],pl,pll)
	
	fg4
end

# ╔═╡ 9993023e-588a-4ab8-999a-377702cc93a3
md"""
- As mentioned earlier, the elementary components separated in the time series' trajectory space may not coincide with a single, interpretable component in the time series.
    - For example, ``\tilde{F}_0$ and $\tilde{F}_1`` both look vaguely like the trend—are they *really* separate components? Similarly, ``\tilde{F}_2`` and ``\tilde{F}_3`` are almost identical, except near the boundaries of the time series.

- We'll introduce a way to quantify which ``\tilde{F_i}`` should be grouped together, but for a moment, let's follow our instincts and apply our earlier grouping for ``\tilde{F}^{\text{(trend)}}``,  ``\tilde{F}^{\text{(periodic 1)}}``,  ``\tilde{F}^{\text{(periodic 2)}}`` and  ``\tilde{F}^{\text{(noise)}}``, and see how the SSA-separated components compare with the original components that compose the toy time series:
"""

# ╔═╡ 7c912664-b822-46e6-9ff8-7743c6569447
begin
	# Assemble the grouped components of the time series.
	F_trend = X_to_TS(sum(X_elem[[1, 2, 7], :], dims=1)[1])
	F_periodic1 = X_to_TS(sum(X_elem[[3, 4], :], dims=1)[1])
	F_periodic2 = X_to_TS(sum(X_elem[[5, 6], :], dims=1)[1])
	F_noise = X_to_TS(sum(X_elem[8:end, :], dims=1)[1]);
end;

# ╔═╡ af222a9a-cf1f-4926-9ace-65723b03fc69
begin
	fg5 = Figure()
	
	ax15 = Axis(fg5[1, 1],
	    xlabel = "t",
	    ylabel = L"\tilde{F}_i(t)",
	    title = "Grouped Time Series Components",
	    )
	
	lines!(t,F, linewidth=1,label="F")
	lines!(t, F_trend,label=L"\tilde{F}(trend)")
	lines!(t, F_periodic1,label=L"\tilde{F}(periodic1)")
	lines!(t, F_periodic2,label=L"\tilde{F}(periodic2)")
	lines!(t, F_noise, alpha=0.5,label=L"\tilde{F}(noise)")
	
	axislegend(position=:ct)
	
	
	fg5
end

# ╔═╡ 8060ad99-5a09-4c77-8686-eab07dd08204
begin
	components = [("Trend", trend, F_trend),
	              ("Periodic 1", periodic1, F_periodic1),
	              ("Periodic 2", periodic2, F_periodic2),
	              ("Noise", noise, F_noise)]
	
	fg = Figure()
	
	cols = 2
	
	
	for e in enumerate(components)
	    r = cld(e[1], cols)
	    c = e[1] - (r - 1) * cols
	    ax = Axis(fg[r, c],title=e[2][1])
	    lines!(ax, t, e[2][2], linewidth=2.5, alpha=0.7, linestyle=:dash)
	    lines!(ax, t, e[2][3])
	end
	
	fg
end

# ╔═╡ ef8aa983-c67b-432a-84a6-1c1e45f92e22
md"""
- After grouping the elementary components together, it looks like SSA has done a great job separating the original components of the toy time series—especially the two periodic components with differing frequencies and amplitudes.
    - However, the separation isn't perfect: all the components deteriorate near the boundaries, especially the trend and second periodic component.
    - This is common in SSA, and arises from the fact that, under the SSA formalism, most types of series (i.e. polynomial, sine, exponential, etc.) are not exactly separable. Therefore, our attempt to recover the *exact* parabolic trend and periodic components from the toy series was always doomed to fail.

- However, that does not mean components of a time series cannot be *approximately separable*, as we witnessed above. (There is also the concept of *asymptotic separability* when the length of the time series approaches infinity.
"""

# ╔═╡ 099bc6ae-ba2d-45f8-be0d-40fb83326b3f
md"""
### Time Series Component Separation and Grouping
***

- So far, we have grouped the eigentriples/components of the toy time series together by visual inspection; that is, we decided which components belonged together by their appearance.

- This is fine for a short and simple time series, however, for longer and more complicated time series, we seek a method that quantifies whether a reconstructed component ``\tilde{F}_i`` can be considered separate from another component ``\tilde{F}_j``, so we don't need to make grouping decisions by visually inspecting each ``\tilde{F}_i``.

- For two reconstructed time series, ``\tilde{F}_i`` and ``\tilde{F}_j``, of length ``N``, and a window length ``L``, we define the *weighted inner product*, ``(\tilde{F}_i, \tilde{F}_j)_w`` as:

```math
(\tilde{F}_i, \tilde{F}_j)_w = \sum_{k=0}^{N-1} w_k \tilde{f}_{i,k} \tilde{f}_{j,k}
```
- where ``\tilde{f}_{i,k}`` and ``\tilde{f}_{j,k}`` are the ``k``th values of ``\tilde{F}_i`` and ``\tilde{F}_j``, respectively, and ``w_k`` is given by
```math
w_{k} = \left\{
  \begin{array}{lr}
    k+1 & \ 0 \le k \le L-1 \\
    L & \ L \le k \le K-1 \\
    N - k & \ K \le k \le N-1 \\
  \end{array}
\right.
```

- remembering that ``K = N - L + 1``.

- The weight ``w_k`` simply reflects the number of times ``\tilde{f}_{i,k}`` and ``\tilde{f}_{j,k}`` appear in the Hankelised matrices ``\mathbf{\tilde{X}}_i`` and ``\mathbf{\tilde{X}}_j``, from which the time series ``\tilde{F}_i`` and ``\tilde{F}_j`` have been obtained.

- Put simply, if ``(\tilde{F}_i, \tilde{F}_j)_w = 0``, ``\tilde{F}_i`` and ``\tilde{F}_j`` are *w-orthogonal* and the time series components are separable.
    - Of course, total w-orthogonality does not occur in real life, so instead we define a ``d \times d`` ***weighted correlation*** matrix, ``\mathbf{W}_{\text{corr}}``, which measures the deviation of the components ``\tilde{F}_i`` and ``\tilde{F}_j`` from w-orthogonality.
    
- The elements of ``\mathbf{W}_{\text{corr}}`` are given by:

```math
W_{i,j} = \frac{(\tilde{F}_i, \tilde{F}_j)_w}{\lVert \tilde{F}_i \rVert_w \lVert \tilde{F}_j \rVert_w}
```

- where ``\lVert \tilde{F}_k \rVert_w = \sqrt{(\tilde{F}_k, \tilde{F}_k)_w}`` for ``k = i,j``.

- The interpretation of ``W_{i,j}`` is straightforward: if ``\tilde{F}_i`` and ``\tilde{F}_j`` are arbitrarily close together (but not identical), then ``(\tilde{F}_i, \tilde{F}_j)_w \rightarrow \lVert \tilde{F}_i \rVert_w \lVert \tilde{F}_j \rVert_w`` and therefore ``W_{i,j} \rightarrow 1``.
    - Of course, if ``\tilde{F}_i`` and ``\tilde{F}_j`` are w-orthogonal, then ``W_{i,j} = 0``.
    
- Moderate values of ``W_{i,j}`` between 0 and 1, say ``W_{i,j} \ge 0.3``, indicate components that may need to be grouped together.




- Let's now construct the w-correlation matrix for the toy time series:
"""

# ╔═╡ 4c70cf9d-54f7-42ae-86bf-1df2230d986b
begin
	# Get the weights w first, as they'll be reused a lot.
	# Julia's 1:L creates the range 1 to L.
	# repeat([L], K-L-1) repeats L K-L-1 times.
	# The `vcat` function concatenates the arrays vertically (into a single vector).
	# The use of `Float64` ensures consistent floating-point arithmetic.
	w = vcat(
	    collect(1:L),
	    repeat([L], K - L - 1),
	    collect(L:-1:1)
	) |> x -> convert(Vector{Float64}, x)
	
	
	# Get all the components of the toy series, store them as columns in F_elem array.
	# In Julia, a comprehension is used to create an array of results from X_to_TS.
	# The resulting F_elem will be a Vector of Vectors (or a Vector of Arrays),
	# which is the idiomatic translation of the NumPy array where columns are time series.
	F_elem = [X_to_TS(X_elem[i]) for i in 1:d]
	
	# Calculate the individual weighted norms, ||F_i||_w, first, then take inverse square-root.
	# The dot product `w' * (F_elem[i] .^ 2)` calculates the weighted norm squared.
	# `.^` is the element-wise power operator.
	# The result is stored as a vector.
	F_wnorms_sq = [w' * (F_elem[i] .^ 2) for i in 1:d]
	F_wnorms = F_wnorms_sq .^ -0.5
	
	# Calculate the w-corr matrix.
	# Initialize with an identity matrix of size d x d.
	Wcorr = Matrix{Float64}(I, d, d)
	
	for i in 1:d
	    # Start j from i+1 because W[i,j] = W[j,i] and the diagonal W[i,i] is 1.
	    for j in (i+1):d
	        # w-dot product: w' * (F_elem[i] .* F_elem[j])
	        # `.*` is the element-wise multiplication operator.
	        # Note: Julia uses 1-based indexing, so F_wnorms[i] and F_wnorms[j] are used directly.
	
	        dot_product = w' * (F_elem[i] .* F_elem[j])
	
	        Wcorr[i, j] = abs(dot_product * F_wnorms[i] * F_wnorms[j])
	        Wcorr[j, i] = Wcorr[i, j]
	    end
	end
end

# ╔═╡ 5d5ff0d5-f24c-4260-8c86-b9f5162940b9
md"- Plot the w-correlation matrix."

# ╔═╡ 72a030bd-518c-47e2-a020-d460f4062f6d
begin
	fig = Figure()
	ax = Axis(fig[1, 1],
	    title = "The W-Correlation Matrix for the Toy Time Series",
	    xlabel = L"\tilde{F}_i",
	    ylabel = L"\tilde{F}_j"
	)
	
	hm6 = CairoMakie.heatmap!(ax, Wcorr,
	    colorrange = (0, 1)
	)
	
	Colorbar(fig[1, 2],
	    hm6,
	    label = L"W_{ij}",
	)
	
	ax.yreversed = true
	
	fig
end

# ╔═╡ fd32e9e5-546a-4be3-90f9-9d36738152b1
md"""
- The structure of $\mathbf{W}_{\text{corr}}$ shows a lot of correlation between the time series components, particularly in the range $7 \le i,j \le 69$. As these were the components we classified as belonging to the noise in the time series, it is no surprise that there are non-negligible correlations between all of them; this is a natural result of the noise having no underlying structural component that can be further separated.

- It is important to note that $\mathbf{W}_{\text{corr}}$ is roughly split into two 'blocks': $0 \le i,j \le 6$, and $7 \le i,j \le 69$. This corresponds to two main groupings: a smoothed time series (i.e. the trend plus the two periodic components), and the residual noise. Zooming into the first seven components in $\mathbf{W}_{\text{corr}}$:
"""

# ╔═╡ 7fcb3776-fa5c-46f7-abee-688d4ea28466
begin
	fig2 = Figure()
	axf2 = Axis(fig2[1, 1],
	    title = "The W-Correlation Matrix for the Toy Time Series",
	    xlabel = L"\tilde{F}_i",
	    ylabel = L"\tilde{F}_j",
	)
	
	xlims!(axf2,0.5,6.5),
	ylims!(axf2,6.5,0.5)
	
	hm2 = CairoMakie.heatmap!(axf2, Wcorr,
	    colorrange = (0, 1)
	)
	
	Colorbar(fig2[1, 2],
	    hm2,
	    label = L"W_{ij}",
	    width = 30, # width in pixels
	)
	
	axf2.yreversed = true
	
	fig2
end

# ╔═╡ 919a747c-c657-45b1-b578-cf9f532e54bf
md"""
- The initial appearance-based groupings we made for the first six components are supported by the corresponding w-correlation values. ``\tilde{F}_0`` and ``\tilde{F}_1`` have ``W_{0,1} = 0.40``, suggesting they should be paired. ``\tilde{F}_1`` and ``\tilde{F}_6`` also have ``W_{1,6} = 0.39``, suggesting ``\tilde{F}_6`` should also be grouped with ``\tilde{F}_0`` and ``\tilde{F}_1`` as a trend component. However, ``\tilde{F}_6`` also has a slight w-correlation with ``\tilde{F}_4`` and ``\tilde{F}_5``, but since ``\tilde{F}_4`` and ``\tilde{F}_5`` have no w-correlation with ``\tilde{F}_0`` and ``\tilde{F}_1``, we choose to keep ``\tilde{F}_6`` with ``\tilde{F}_0`` and ``\tilde{F}_1``.

- Our prior groupings of  ``\tilde{F}^{\text{(periodic 1)}} = \tilde{F}_2 + \tilde{F}_3`` and ``\tilde{F}^{\text{(periodic 2)}} = \tilde{F}_4 + \tilde{F}_5`` are clearly justified by the w-correlation matrix, with ``W_{2,3} = 0.98`` and ``W_{4,5} = 0.98``.
"""

# ╔═╡ 097f54d6-2a9c-4eeb-8008-07f087d592fc
md"""
## A julia module for SSA
***

- It is time to collect the SSA code into a handy module, imaginatively named *`SSAJL`*, which will form the basis for the rest of this notebook. Each instance of the class will contain the decomposition of a time series for some window length ``L``, and provide useful methods to analyse, plot and reconstruct the time series.

- To summarise the SSA algorithm:
    1. For a time series ``F = (f_0, \ f_1, \ldots, \ f_{N-1})``, and a window length ``L``, form the trajectory matrix ``\mathbf{X}``, with columns given by the vectors ``(f_i, \ldots, f_{L+i-1})^{\text{T}}``, ``0 \le i \le N-L``.
    2. Decompose ``\mathbf{X}`` with the singular value decomposition, ``\mathbf{X} = \sum_{i=0}^{d-1}\sigma_i U_i V^{\text{T}}_i``.
    3. Construct the ``d`` elementary matrices ``\mathbf{X}_i = \sigma_i U_i V^{\text{T}}_i``.
    4. Diagonally average the ``\mathbf{X}_i`` to form the elementary time series components ``\tilde{F}_i``, such that ``F = \sum_{i=0}^{d-1} \tilde{F}_i``.
    5. Calculate and store the weighted correlation matrix, ``\mathbf{W}_{\text{corr}}``, for the ``\tilde{F}_i``.

- The task of grouping and classifying the elementary components ``\tilde{F}_i`` is left to the user.
"""

# ╔═╡ 692bd29f-8949-42a6-885a-cada66630eac
md"""
## The Window Length
***

- We have now established the machinery to easily investigate the effect of the window length parameter, $L$, on the decomposition of our toy time series.

#### ``L = 2``

- A window length of 2 may seem like a useless choice, but it's a good place to start and watch the time series get decomposed into more and more components.
"""

# ╔═╡ bd0ac540-1ca1-42f7-b0df-7c41b315b1f5
begin

	F_ssa_L2 = SSAJL.SSA(F, 2)
	
	cpmm = SSAJL.components_to_matrix(F_ssa_L2)
	orts = F_ssa_L2.orig_TS
	
	fg7 = Figure()
	
	ax7 = Axis(fg7[1, 1],
	    xlabel = "t",
	    ylabel = L"\tilde{F}_i(t)",
	    title = "L=2 for the Toy Time Series",
	)
	
	lines!(t, orts, alpha=0.4)
	lines!(t, cpmm[:,1], label="F1")
	lines!(t, cpmm[:,2], label="F2")
	
	
	axislegend(ax7,framevisible = false,position = :ct)
	
	fg7
end

# ╔═╡ 128845f9-9def-4f79-8562-e8ec14120594
md"""
- For ``L=2`` we can only expect two elementary components to be returned. Even for such a small window length, the SSA algorithm has started to separate the high-frequency noise from the series, giving us a somewhat-denoised version of the original series in the component ``\tilde{F}_1``.

#### ``L = 5``
- Let's go up to a window length of 5, and see what happens to the elementary components:
"""

# ╔═╡ 88b4bbbe-90df-434e-8ab9-3c2783733db4
begin
	F_ssa_L5 = SSAJL.SSA(F, 5)
	
	cpmm8 = SSAJL.components_to_matrix(F_ssa_L5)
	orts8 = F_ssa_L5.orig_TS
	
	fg8 = Figure()
	
	ax81 = Axis(fg8[1, 1],
	    xlabel = "t",
	    ylabel = L"\tilde{F}_i(t)",
	    title = "L=5 for the Toy Time Series",
	)
	
	lines!(t, orts8, alpha=0.4)
	lines!(t, cpmm8[:,1], label="F1")
	lines!(t, cpmm8[:,2], label="F2")
	lines!(t, cpmm8[:,3], label="F3")
	lines!(t, cpmm8[:,4], label="F4")
	lines!(t, cpmm8[:,5], label="F5")
	
	
	
	axislegend(ax81,framevisible = false,position = :ct)
	
	fg8
end

# ╔═╡ d0ecf62d-4f12-4921-a43b-392fa41f039f
md"""
- We see that ``\tilde{F}_1`` is now a well-and-truly denoised version of the original series. ``\tilde{F}_2`` is a poorly resolved periodic component, while ``\tilde{F}_3`` to ``\tilde{F}_5`` are just noise.

#### ``L = 20``

- Let's quadruple the window length, and instead of inspecting elementary components, we'll look at the resulting w-correlation matrix and make some grouping decisions first.
"""

# ╔═╡ ff15d461-5923-4bf0-9ce7-739f9040375c
begin
	F_ssa_L20 = SSAJL.SSA(F, 20)
	
	SSAJL.plot_wcorr(F_ssa_L20,ptitle=L"W-Correlation for Toy Time Series, $L=20$")
	
end

# ╔═╡ 4177b9f8-957f-4b2e-9ae3-60b370ff5fc4
md"""
- The w-correlation matrix for ``L=20`` is split (roughly) into two blocks: ``\tilde{F}_1`` to ``\tilde{F}_4``, and ``\tilde{F}_5`` to ``\tilde{F}_{20}``.

- Within those blocks, the size of the ``W_{i,j}`` values suggest that we need to group ``\tilde{F}_2``, ``\tilde{F}_3`` and ``\tilde{F}_4``, and group all ``\tilde{F}_5, \cdots, \tilde{F}_{20}``.

- This grouping is certainly not ideal, as ``\tilde{F}_4`` has non-negligible w-correlation with components in the second block. We'll plot our chosen component groupings, along with ``\tilde{F}_4`` on its own, and see if we're justified in our choice of grouping:
"""

# ╔═╡ 0b935d48-7d19-49ad-a907-da757ecdec27
begin
	fg9 = Figure()
	
	ax91 = Axis(fg9[1, 1],
	    xlabel = "t",
	    ylabel = L"\tilde{F}_i(t)",
	    title = "Component Groupings for Toy Time Series, L=20",
	)
	
	lines!(t, SSAJL.reconstruct(F_ssa_L20,1), label=L"\tilde{F}_1")
	lines!(t, SSAJL.reconstruct(F_ssa_L20,[2,3,4]), label=L"\tilde{F}_2+\tilde{F}_3+\tilde{F}_4")
	lines!(t, SSAJL.reconstruct(F_ssa_L20,5:20), label=L"\tilde{F}_5 + \cdots + \tilde{F}_{20}")
	lines!(t, SSAJL.reconstruct(F_ssa_L20,4), label=L"\tilde{F}_4")
	
	axislegend(ax91,framevisible = false,position = :ct)
	
	fg9
	
end

# ╔═╡ e4460467-8dbd-49a9-bc33-e2f9c9329609
md"""
- For ``L = 20`` we begin to see the trend and periodic components start to take shape. The single component ``\tilde{F}_1`` looks like the parabolic trend, and the group ``\tilde{F}_2 + \tilde{F}_3 + \tilde{F}_4`` is a very handsome periodicity, corresponding to the sum of the two periodic components in the original definition of the toy time series.
    
    - The component ``\tilde{F}_4`` is indeed troublesome, looking like it contributes to both noise *and* periodicity. This suggests we need to increase the window length and see if we get an improved separation of noise and periodicity.
"""

# ╔═╡ 37ffd8ff-1342-4d05-a82a-de68648350ec
md"""
#### ``L = 40``

- Once again, we'll double the window length and inspect the w-correlation matrix first:
"""

# ╔═╡ 2004b09b-be66-44d5-a707-077ac0f670fc
begin
	F_ssa_L40 = SSAJL.SSA(F, 40)
	
	SSAJL.plot_wcorr(F_ssa_L40,ptitle=L"W-Correlation for Toy Time Series, $L=40$")
end

# ╔═╡ fdd667e0-224e-492d-83ce-d2e6a724d5ab
md"""
- The w-correlation matrix for ``L=40`` retains the two-block structure, with ``\tilde{F}_1, \ldots, \tilde{F}_6`` in the first block, and ``\tilde{F}_7, \ldots, \tilde{F}_{40}`` in the second. Let us group the components as follows:
```math
\begin{align*}
    \tilde{F}^{(1)} & = \tilde{F}_1 \\
    \tilde{F}^{(2)} & = \tilde{F}_2 + \tilde{F}_3 + \tilde{F}_4 \\
    \tilde{F}^{(3)} & = \tilde{F}_5 + \tilde{F}_6 \\
    \tilde{F}^{(4)} & = \tilde{F}_7 + \cdots + \tilde{F}_{40} \\
\end{align*}
```
- Once again, it can be argued that this grouping is not ideal, given the non-negligible w-correlation between, for example,  ``\tilde{F}_1`` and ``\tilde{F}_2``.
"""

# ╔═╡ acfc1bdd-032b-4579-803d-f2698d883179
begin
	fg10 = Figure()
	
	ax110 = Axis(fg10[1, 1],
	    xlabel = "t",
	    ylabel = L"\tilde{F}_i(t)",
	    title = "Component Groupings for Toy Time Series, L=40",
	)
	
	lines!(t, SSAJL.reconstruct(F_ssa_L40,1), label=L"\tilde{F}_1")
	lines!(t, SSAJL.reconstruct(F_ssa_L40,[2,3,4]), label=L"\tilde{F}_2+\tilde{F}_3+\tilde{F}_4")
	lines!(t, SSAJL.reconstruct(F_ssa_L40,[5,6]), label=L"\tilde{F}_5+\tilde{F}_6+")
	lines!(t, SSAJL.reconstruct(F_ssa_L40,7:40), label=L"\tilde{F}_7 + \cdots + \tilde{F}_{40}")
	
	axislegend(ax110,framevisible = false,position = :ct)
	
	fg10
end

# ╔═╡ 08b62fc5-5a03-4c2b-9229-ad4f5804544b
md"""
- Interestingly, the trend component ``\tilde{F}^{(0)}`` has started to deteriorate at ``L=40``, with notable 'kinks' in the time series.
- The periodicity in the toy series is now separated into two periodic components of differing amplitudes and frequencies, however with significant deteriorations near the beginning and end of the time series.
- At this stage our window length is 20% of the length of the time series, but the poor quality of the separated components suggests that we still need to increase the window size.
"""

# ╔═╡ b8a9faf2-c55c-4a9d-8bc1-61b23553e1e0
md"""
#### ``L = 60``

- We are now approaching our first choice of ``L=70``, so it's worth investigating how the decomposition converges to our original results. Inspect the w-correlation matrix first:
"""

# ╔═╡ cc03510a-f0a2-4f36-906d-2497cda3e9f4
begin
	F_ssa_L60 = SSAJL.SSA(F, 60)
	
	SSAJL.plot_wcorr(F_ssa_L60,ptitle=L"W-Correlation for Toy Time Series, $L=60$")
end

# ╔═╡ 69f0ccc4-c142-4beb-a7a1-c205ce20622b
md"""
- As with the original ``L=70`` result, the w-correlation matrix is now composed of two separate blocks: ``\tilde{F}_1`` to ``\tilde{F}_7``, and ``\tilde{F}_8`` to ``\tilde{F}_{60}``.
- From experience now, it is clear that ``\tilde{F}^{\text{(signal)}} = \sum_{i=1}^7 \tilde{F}_i`` will be the combined trend and periodic components ('signal'), and ``\tilde{F}^{\text{(noise)}} = \sum_{i=8}^{60} \tilde{F}_i`` will be the noise:
"""

# ╔═╡ 69477ae2-6b9f-4215-a65e-59409b96f6dc
begin
	fg11 = Figure()
	
	ax111 = Axis(fg11[1, 1],
	    xlabel = "t",
	    ylabel = L"\tilde{F}_i(t)",
	    title = "Signal and Noise Components of Toy Time Series, L = 60",
	)
	
	lines!(t, SSAJL.reconstruct(F_ssa_L60,1:7), label=L"\tilde{F}(signal)")
	lines!(t, SSAJL.reconstruct(F_ssa_L60,8:60), label=L"\tilde{F}(noise)")
	
	axislegend(ax111,framevisible = false,position = :ct)
	
	fg11
end

# ╔═╡ fa24cfcf-ea10-42c3-b4b2-63a61327aa2e
md"- It is also of interest viewing the w-correlation matrix for components 1–7:"

# ╔═╡ 04602bfd-4203-4834-b5ce-c575a4a58ca2
SSAJL.plot_wcorr(F_ssa_L60,max_idx=7,ptitle=L"W-Correlation for Toy Time Series, $L=60$ (zoomed)")

# ╔═╡ 07306231-459f-48f2-b872-230517523dad
md"- In order to understand why there is non-negligible w-correlation between most of the first seven components, it'll be prudent to plot all of them at once:"

# ╔═╡ 01d5f464-aa4f-4d94-81f4-e16fe84ed6b2
begin
	cpmm12 = SSAJL.components_to_matrix(F_ssa_L60,7)
	
	
	fg12 = Figure()
	
	ax112 = Axis(fg12[1, 1],
	    xlabel = "t",
	    ylabel = L"\tilde{F}_i(t)",
	    title = "The First 7 Components of the Toy Time Series, L=60",
	)
	
	for i in 1:7
	    lines!(t, cpmm12[:,i], label=L"F_{%$i}")
	end
	
	axislegend(ax112,framevisible = false,position = :ct)
	
	fg12
end

# ╔═╡ 90f36511-2f1f-422d-a8b2-a605c50a7482
md"""
- It is clear from the plot above that the trend component ``\tilde{F}_2`` contains oscillations that are the same frequency as the periodic components ``\tilde{F}_3`` and ``\tilde{F}_4``, resulting in substantial w-correlations.
- As we've already witnessed, increasing the window length to ``L = 70`` is enough to remove these oscillations almost entirely.
    - It is certainly possible to increase the window length further—up to ``L = N/2 = 100`` in this case—however, beyond ``L=70``, the time series components generated by the leading seven eigentriples are (relatively) insensitive to the window length.

- In SSA, there are no hard and fast rules to setting the perfect window length, beyond ``2 \le L \le N/2``. However, longer window lengths (up to 30–45% of the length of the time series) are sometimes required to adequately separate underlying periodicities from the overall trend. Some trial and error is needed, but it is often easy to start at a 'large-enough' window length, and work from there.
"""

# ╔═╡ 2d596c28-74bc-4ff8-a030-fbac18dbceb0
md"""
## Reference & Material

Material and papers related to the topics discussed in this lecture.


- [Hassani (2007) - "Singular Spectrum Analysis: Methodology and Comparison"](https://jds-online.org/journal/JDS/article/1027/info)
- [Rico et al. (2025) - "Singular spectrum analysis of Fermi-LAT blazar light curves: A systematic search for periodicity and trends in the time domain"](https://ui.adsabs.harvard.edu/abs/2025A%26A...697A..35R/abstract)
"""

# ╔═╡ 12f422fb-acc1-4434-af63-63b3b13ebcef
md"""
## Further Material

Papers for examining more closely some of the discussed topics.

- [Deng (2014) - "Time Series Decomposition Using Singular Spectrum Analysis "](https://dc.etsu.edu/etd/2352/)
- [Golyandina et al. (2001) - "Analysis of Time Series Structure: SSA and Related Techniques"](https://www.crcpress.com/Analysis-of-Time-Series-Structure-SSA-and-Related-Techniques/Golyandina-Nekrutkin-Zhigljavsky/p/book/9781584881940)
"""

# ╔═╡ f94e2f05-862a-4cd3-ac7e-7ca818c882dd
md"""
### Credits
***

This notebook is obtained and further elaborated from: [https://www.kaggle.com/code/jdarcy/introducing-ssa-for-time-series-decomposition](https://www.kaggle.com/code/jdarcy/introducing-ssa-for-time-series-decomposition).
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
    <td><a href="./open?path=Lectures/Lecture-SingularSpectrumAnalysis/Lecture-MotionSensors.jl">Science case about motion sensor data</a></td>
  </tr>
  <tr>
	<td>html</td>
    <td><a href="Lectures/Lecture-NonParametricPeriodograms/Lecture-NnParametricPeriodograms.html">Lecture about non-parametric periodograms</a></td>
<td><a href="Lectures/Lecture-SingularSpectrumAnalysis/Lecture-MotionSensors.html">Science case about motion sensor data</a></td>
  </tr>
</table>


"""

# ╔═╡ 206474b8-0811-4785-8a71-acdcfd20b76c
md"""
**Copyright**

This notebook is provided as [Open Educational Resource](https://en.wikipedia.org/wiki/Open_educational_resources). Feel free to use the notebook for your own purposes. The text is licensed under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/), the code of the examples, unless obtained from other properly quoted sources, under the [MIT license](https://opensource.org/licenses/MIT). Please attribute the work as follows: *Stefano Covino, Time Domain Astrophysics - Lecture notes featuring computational examples, 2026*.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
CommonMark = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CairoMakie = "~0.15.9"
CommonMark = "~0.8.15"
LaTeXStrings = "~1.4.0"
PlutoUI = "~0.7.61"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.5"
manifest_format = "2.0"
project_hash = "c82c12050b4356fb944cb2359f2d87d68691c7da"

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
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonMark]]
deps = ["Crayons", "PrecompileTools"]
git-tree-sha1 = "3faae67b8899797592335832fccf4b3c80bb04fa"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.15"

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
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

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

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

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
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

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

[[deps.MIMEs]]
git-tree-sha1 = "1833212fd6f580c20d4291da9c1b4e8a655b128e"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.0.0"

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
git-tree-sha1 = "135492b7e97fc86d9b132b96a54d2d3dd3e0c6a8"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.4.8+0"

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
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

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
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "7e71a55b87222942f0f9337be62e26b1f103d3e4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.61"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

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
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

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
git-tree-sha1 = "e2a7072fc0cdd7949528c1455a3e5da4122e1153"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.56+0"

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
# ╟─4d477519-c44f-434c-b7e0-8daaa5009358
# ╟─ec67de24-d88a-46fa-ae46-c0cd7b797adc
# ╠═6a1315d1-9a6d-4ce0-b1c0-3fe22beb1ec2
# ╟─3ddd0f61-79d0-473c-8fec-a0e0c3fc72bf
# ╟─5029a214-0841-40fb-b397-4a2e1047bfb7
# ╟─404060d3-23ec-400b-84cf-779e63b90293
# ╟─72cf6fcf-2e2f-4dee-994b-ce2bfef51802
# ╠═0646411b-bbbe-4b6d-bf94-33c4e9522ed5
# ╟─2f2430a4-ea1f-4373-b5a7-10db4756b29b
# ╟─e210c43a-1205-40b8-8989-7975ea85e648
# ╠═a759da24-db69-402a-bde3-25edbae59834
# ╠═37551c89-38d4-46e6-9109-52b3a065ed73
# ╟─99fb1b5c-2857-43a3-b69b-fd55a7e35e71
# ╟─b393afda-d52b-4b48-aee0-cf1a42054afe
# ╟─61850790-88da-4cf9-8ed3-f96282f04e2d
# ╟─35390351-c3d9-4218-aa8a-02cabda3c192
# ╟─86d2e08b-bfab-46a3-9db6-6ef8a22147b2
# ╠═490e3a2a-cf23-42e3-9530-3a6fad82bb30
# ╟─aaebbbcd-6bb6-46e0-8564-8029a4f2d1cf
# ╠═2b3f952d-db20-4790-8482-9c6e981f60e0
# ╠═ef9e9f9a-7598-4e22-9add-f1dd2a74c516
# ╟─3d354f7c-ea01-4465-9b0e-b911801d742b
# ╠═54273740-9fb0-46b6-a928-c74b9077eaa3
# ╟─018667ce-b84d-46f2-9fce-dd890fc9a557
# ╟─f05b56b9-801c-42e1-9062-4aa3234bc482
# ╟─911757a4-a582-469d-8a2a-9273065e8d83
# ╠═1d406a03-8158-4b3f-a734-d1b80c259847
# ╠═38fc4698-f890-4104-b37b-0d938a02f012
# ╟─d45f4fbd-8ef0-47db-8278-72cd1230ab5e
# ╠═5385afbd-33a9-4dba-af14-dfe98c841611
# ╟─18745681-9e0c-46f1-b9fe-d00b61fbafb8
# ╠═4c55c72d-3c4e-4fb9-8fa0-4dd6d6259637
# ╟─9993023e-588a-4ab8-999a-377702cc93a3
# ╠═7c912664-b822-46e6-9ff8-7743c6569447
# ╠═af222a9a-cf1f-4926-9ace-65723b03fc69
# ╠═8060ad99-5a09-4c77-8686-eab07dd08204
# ╟─ef8aa983-c67b-432a-84a6-1c1e45f92e22
# ╟─099bc6ae-ba2d-45f8-be0d-40fb83326b3f
# ╠═4c70cf9d-54f7-42ae-86bf-1df2230d986b
# ╟─5d5ff0d5-f24c-4260-8c86-b9f5162940b9
# ╠═72a030bd-518c-47e2-a020-d460f4062f6d
# ╟─fd32e9e5-546a-4be3-90f9-9d36738152b1
# ╠═7fcb3776-fa5c-46f7-abee-688d4ea28466
# ╟─919a747c-c657-45b1-b578-cf9f532e54bf
# ╟─097f54d6-2a9c-4eeb-8008-07f087d592fc
# ╠═cd347311-bb1e-4404-a93c-3ac7a37e04e2
# ╟─692bd29f-8949-42a6-885a-cada66630eac
# ╠═bd0ac540-1ca1-42f7-b0df-7c41b315b1f5
# ╟─128845f9-9def-4f79-8562-e8ec14120594
# ╠═88b4bbbe-90df-434e-8ab9-3c2783733db4
# ╟─d0ecf62d-4f12-4921-a43b-392fa41f039f
# ╠═ff15d461-5923-4bf0-9ce7-739f9040375c
# ╟─4177b9f8-957f-4b2e-9ae3-60b370ff5fc4
# ╠═0b935d48-7d19-49ad-a907-da757ecdec27
# ╟─e4460467-8dbd-49a9-bc33-e2f9c9329609
# ╟─37ffd8ff-1342-4d05-a82a-de68648350ec
# ╠═2004b09b-be66-44d5-a707-077ac0f670fc
# ╟─fdd667e0-224e-492d-83ce-d2e6a724d5ab
# ╠═acfc1bdd-032b-4579-803d-f2698d883179
# ╟─08b62fc5-5a03-4c2b-9229-ad4f5804544b
# ╟─b8a9faf2-c55c-4a9d-8bc1-61b23553e1e0
# ╠═cc03510a-f0a2-4f36-906d-2497cda3e9f4
# ╟─69f0ccc4-c142-4beb-a7a1-c205ce20622b
# ╠═69477ae2-6b9f-4215-a65e-59409b96f6dc
# ╟─fa24cfcf-ea10-42c3-b4b2-63a61327aa2e
# ╠═04602bfd-4203-4834-b5ce-c575a4a58ca2
# ╟─07306231-459f-48f2-b872-230517523dad
# ╠═01d5f464-aa4f-4d94-81f4-e16fe84ed6b2
# ╟─90f36511-2f1f-422d-a8b2-a605c50a7482
# ╟─2d596c28-74bc-4ff8-a030-fbac18dbceb0
# ╟─12f422fb-acc1-4434-af63-63b3b13ebcef
# ╟─f94e2f05-862a-4cd3-ac7e-7ca818c882dd
# ╠═b36fd613-95c8-44bf-876d-4eb345c26f08
# ╟─206474b8-0811-4785-8a71-acdcfd20b76c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
