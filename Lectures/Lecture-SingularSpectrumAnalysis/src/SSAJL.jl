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
