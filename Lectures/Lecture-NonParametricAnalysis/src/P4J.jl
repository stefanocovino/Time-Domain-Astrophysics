
using LinearAlgebra, Statistics, StatsBase

export Periodogram, fit!, finetune!, periodogram, best_period, best_periods

abstract type PeriodogramMethod end

struct QMI <: PeriodogramMethod
    variant::Symbol  # :euclidean or :cauchy_schwarz
end
struct PDM <: PeriodogramMethod end
struct LK <: PeriodogramMethod end
struct MHAoV <: PeriodogramMethod end

mutable struct Periodogram{M<:PeriodogramMethod}
    method::M
    frequencies::Vector{Float64}
    scores::Vector{Float64}
    best_frequency::Float64
    best_score::Float64
    t::Vector{Float64}
    y::Vector{Float64}
    dy::Vector{Float64}
    whiten::Bool
    coarse_step::Float64
    refined_frequencies::Vector{Float64}
    refined_scores::Vector{Float64}
end

# ==================== Constructor ====================
function Periodogram(method::String)
    m = lowercase(method)
    if m in ("qmi", "qmieu")
        meth = QMI(:euclidean)
    elseif m in ("qmics", "qmi")
        meth = QMI(:cauchy_schwarz)
    elseif m == "pdm"
        meth = PDM()
    elseif m in ("lk", "stringlength", "lksl")
        meth = LK()
    elseif m in ("mhaov", "aov")
        meth = MHAoV()
    else
        error("Metodo sconosciuto: $method. Supported: qmieu, qmics, pdm, lk, mhaov")
    end
    Periodogram(meth, Float64[], Float64[], NaN, NaN, Float64[], Float64[], Float64[], false, 0.0, Float64[], Float64[])
end


# ==================== QMI ====================
wrapped_cauchy(Δφ::Float64, hφ::Float64) = (1 - exp(-2hφ)) / (1 + exp(-2hφ)) / (1 + exp(-2hφ * cos(2π * Δφ)))

function _qmi_score(t, y, dy, f, variant, whiten)
    N = length(t)
    phase = mod.(t .* f, 1.0)
    w = @. 1.0 / (dy^2); w ./= sum(w)
    loc = sum(y .* w)
    scale = sqrt(sum(w .* (y .- loc).^2))
    y_std = whiten ? @.((y - loc) / scale) : y
    dy_std = whiten ? @.(dy / scale) : dy

    hm = 0.9 * min(sqrt(sum(w .* (y_std .- sum(y_std .* w)).^2)), 
                   weighted_iqr(y_std, w) / 1.349) * N^(-0.2)
    hp = 1.0

    IP_M = IP_Φ = IP_ΦM = 0.0
    @inbounds for i in 1:N
        yi, σi, ϕi = y_std[i], dy_std[i], phase[i]
        for j in 1:N
            Δm = yi - y_std[j]
            Δφ = ϕi - phase[j]
            denom = hm*hm + σi*σi + dy_std[j]*dy_std[j]
            g = exp(-0.5 * Δm*Δm / denom) / sqrt(2π * denom)
            wc = wrapped_cauchy(Δφ, hp)
            IP_M += g
            IP_Φ += wc
            IP_ΦM += g * wc
        end
    end
    IP_M /= N*N; IP_Φ /= N*N; IP_ΦM /= N*N

    IP_M_prod = IP_Φ_prod = 0.0
    @inbounds for i in 1:N, j in 1:N
        Δm = y_std[i] - y_std[j]
        Δφ = phase[i] - phase[j]
        denom = hm*hm + dy_std[i]*dy_std[i] + dy_std[j]*dy_std[j]
        g = exp(-0.5 * Δm*Δm / denom) / sqrt(2π * denom)
        wc = wrapped_cauchy(Δφ, hp)
        IP_M_prod += g
        IP_Φ_prod += wc
    end
    IP_M_prod /= N*N; IP_Φ_prod /= N*N

	if variant === :euclidean
        # Versione scalata per avvicinarsi al comportamento di P4J Python
        raw = IP_ΦM - 2 * (IP_M_prod * IP_Φ_prod) + IP_Φ * IP_M
        return raw
    else  # :cauchy_schwarz
        return log(IP_ΦM) - 2*log(IP_M_prod * IP_Φ_prod) + log(IP_Φ) + log(IP_M)
    end
end

function weighted_iqr(x, w)
    idx = sortperm(x)
    x, w = x[idx], w[idx] ./ sum(w)
    cumw = cumsum(w)
    q25 = x[searchsortedfirst(cumw, 0.25)]
    q75 = x[searchsortedfirst(cumw, 0.75)]
    q75 - q25
end

# ==================== PDM con pesi ====================
function _pdm_score(t::AbstractVector, y::AbstractVector, dy::AbstractVector, f::Real)
    phase = mod.(t .* f, 1.0)
    nbins = 10
    bin_edges = range(0.0, 1.0; length=nbins + 1)
    
    s2_within = 0.0
    w_total = 0.0
    
    @inbounds for b in 1:nbins
        mask = (bin_edges[b] .<= phase) .& (phase .< bin_edges[b+1])
        idx = findall(mask)
        n = length(idx)
        n < 2 && continue
        
        y_bin = @view y[idx]
        w_bin = @. 1.0 / (dy[idx]^2)
        w_sum = sum(w_bin)
        mean_bin = sum(y_bin .* w_bin) / w_sum
        bin_var = sum(w_bin .* (y_bin .- mean_bin).^2) / w_sum
        
        s2_within += bin_var * (n - 1)
        w_total += w_sum
    end
    
    w_total < 1e-8 && return 1.0
    
    w_all = @. 1.0 / (dy^2)
    w_sum_all = sum(w_all)
    mean_all = sum(y .* w_all) / w_sum_all
    s2_total = sum(w_all .* (y .- mean_all).^2) / w_sum_all
    
    theta = s2_within / (s2_total * (length(y) - 1))
    return theta
end

# ==================== LK allineato con P4J Python ====================
function _lk_score(t::AbstractVector, y::AbstractVector, dy::AbstractVector, f::Real)
    N = length(t)
    N < 2 && return 0.0

    phase = mod.(t .* f, 1.0)
    idx = sortperm(phase)

    sl = 0.0
    @inbounds for i in 1:N-1
        j = idx[i]
        k = idx[i+1]
        Δφ = phase[k] - phase[j]
        Δy = y[k] - y[j]
        sl += sqrt(Δφ^2 + Δy^2)
    end

    # Chiusura del ciclo (esattamente come in P4J Python)
    j = idx[end]
    k = idx[1]
    Δφ = phase[k] - phase[j]
    Δφ < 0 && (Δφ += 1.0)
    Δy = y[k] - y[j]
    sl += sqrt(Δφ^2 + Δy^2)

	# Normalizzazione che porta i valori vicini a quelli di Python P4J
    # (sl / N è la più comune; a volte si usa sl / (N * std(y)))
    
	#sl_normalized = sl / N                     # prova prima questa
 	sl_normalized = sl / (N * std(y))        # alternativa se serve

    return -sl_normalized      # segno negativo → più alto = meglio
end

# ==================== MHAoV ====================
function _mhaov_score(t::AbstractVector, y::AbstractVector, dy::AbstractVector, f::Real)
    N = length(t)
    N < 4 && return 0.0

    phase = 2π .* mod.(t .* f, 1.0)
    
    # Design matrix per 1 armonica (costante + sin + cos)
    X = [ones(N) sin.(phase) cos.(phase)]
    
    β = X \ y
    yhat = X * β
    
    rss = sum((y .- yhat).^2)          # Residual Sum of Squares
    tss = sum((y .- mean(y)).^2)       # Total Sum of Squares
    
    if rss < 1e-12
        return 1000.0                  # segnale praticamente perfetto
    end

    #power = (tss - rss) / rss          # versione base (F-like)

    # Alternativa ancora più "amplificata" se serve (prova prima quella sopra)
    power = (N - 3) * (tss - rss) / rss
	#power = (tss / rss) - 1 

    return power
end

# ==================== Score dispatch & refinement ====================
function _compute_score(pg::Periodogram, f::Float64)
    if pg.method isa QMI
        _qmi_score(pg.t, pg.y, pg.dy, f, pg.method.variant, pg.whiten)
    elseif pg.method isa PDM
        -_pdm_score(pg.t, pg.y, pg.dy, f)
    elseif pg.method isa LK
        _lk_score(pg.t, pg.y, pg.dy, f)
    elseif pg.method isa MHAoV
        _mhaov_score(pg.t, pg.y, pg.dy, f)
    else
        error("Metodo non supportato")
    end
end

function find_local_maxima(scores::AbstractVector, n::Int=10)
    length(scores) < 3 && return [argmax(scores)]
    local_idx = [i for i in 2:length(scores)-1 if scores[i] > scores[i-1] && scores[i] > scores[i+1]]
    isempty(local_idx) && return [argmax(scores)]
    top = local_idx[partialsortperm(scores[local_idx], 1:min(n, length(local_idx)), rev=true)]
    top
end

function finetune!(pg::Periodogram; fresolution::Float64=1e-5, n_local_optima::Int=10)
    length(pg.scores) == 0 && error("Esegui fit! prima")
    local_idx = find_local_maxima(pg.scores, n_local_optima)

    refined_f = Float64[]
    refined_s = Float64[]
    step = pg.coarse_step

    for i in local_idx
        fc = pg.frequencies[i]
        freqs_fine = collect(range(max(1e-6, fc-step), fc+step; step=fresolution))
        scores_fine = [_compute_score(pg, f) for f in freqs_fine]
        best_i = argmax(scores_fine)
        push!(refined_f, freqs_fine[best_i])
        push!(refined_s, scores_fine[best_i])
    end

    sort_idx = sortperm(refined_s, rev=true)
    pg.refined_frequencies = refined_f[sort_idx]
    pg.refined_scores = refined_s[sort_idx]

    if !isempty(refined_f)
        best_i = argmax(refined_s)
        pg.best_frequency = refined_f[best_i]
        pg.best_score = refined_s[best_i]
    end
    pg
end

# ==================== Fit! ====================
function fit!(pg::Periodogram, t::AbstractVector, y::AbstractVector;
              dy::Union{Nothing,AbstractVector}=nothing,
              fmin::Real=1e-4, fmax::Real=10.0, resolution::Real=1e-4,
              whiten::Bool=false, n_local_optima::Int=10)

    dy = isnothing(dy) ? fill(median(abs.(y .- median(y)))/0.6745, length(t)) : collect(dy)

    pg.t = collect(t)
    pg.y = collect(y)
    pg.dy = dy
    pg.whiten = whiten
    pg.coarse_step = resolution

    freqs = collect(range(fmin, fmax; step=resolution))
    scores = [_compute_score(pg, f) for f in freqs]   # usa dispatch

    pg.frequencies = freqs
    pg.scores = scores
    best_idx = argmax(scores)
    pg.best_frequency = freqs[best_idx]
    pg.best_score = scores[best_idx]

    if n_local_optima > 0
        finetune!(pg; fresolution=1e-5, n_local_optima=n_local_optima)
    else
        pg.refined_frequencies = Float64[]
        pg.refined_scores = Float64[]
    end

    pg
end

# ==================== Accessors ====================
periodogram(pg::Periodogram) = (periods = 1.0 ./ pg.frequencies, scores = pg.scores)

best_period(pg::Periodogram) = isempty(pg.refined_frequencies) ? 
    1.0 / pg.best_frequency : 1.0 / pg.refined_frequencies[1]

function best_periods(pg::Periodogram; n::Int=5)
    if !isempty(pg.refined_frequencies)
        n = min(n, length(pg.refined_frequencies))
        pg.refined_frequencies[1:n], pg.refined_scores[1:n]
    else
        idx = partialsortperm(pg.scores, 1:n, rev=true)
        1.0 ./ pg.frequencies[idx], pg.scores[idx]
    end
end
