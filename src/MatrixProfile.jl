module MatrixProfile
using Statistics, LinearAlgebra
using DSP
using LoopVectorization
using ProgressMeter
using RecipesBase

using Distances
export Euclidean

using SlidingDistancesBase
import SlidingDistancesBase: floattype, lastlength, distance_profile, distance_profile!

export matrix_profile, distance_profile, motifs, anomalies, mpdist, mpdist_profile, onsets, snippets, seqs, resample


struct Profile{TT,TP,QT}
    T::TT
    P::TP
    I::Vector{Int}
    m::Int
    Q::QT
end

"""
    profile = matrix_profile(T, m; showprogress=true)

Return the matrix profile and the profile indices of time series `T` with window length `m`. See fields `profile.P, profile.I`. You can also plot the profile. Uses the STOMP algorithm.

Reference: [Matrix profile II](https://www.cs.ucr.edu/~eamonn/STOMP_GPU_final_submission_camera_ready.pdf).
"""
function matrix_profile(T::AbstractVector{<:Number}, m::Int; showprogress=true)
    n   = lastlength(T)
    l   = n-m+1
    n > 2m+1 || throw(ArgumentError("Window length too long, maximum length is $((n+1)÷2)"))
    μ,σ  = running_mean_std(T, m)
    QT   = window_dot(getwindow(T, m, 1), T)
    QT₀  = copy(QT)
    D    = distance_profile(Euclidean(), QT, μ, σ, m)
    P    = copy(D)
    I    = ones(Int, l)
    prog = Progress((l - 1) ÷ 10, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
    @inbounds for i = 2:l
        for j = l:-1:2
            @fastmath QT[j] = QT[j-1] - T[j-1] * T[i-1] + T[j+m-1] * T[i+m-1] # Consider updating to eqs. 3-7 https://www.cs.ucr.edu/~eamonn/SCAMP-camera-ready-final1.pdf
            # QT[j] = muladd(-T[j-1], T[i-1],  muladd(T[j+m-1], T[i+m-1], QT[j-1]))
            # The expression with fastmath appears to be both more accurate and faster than both muladd and fma
        end
        QT[1] = QT₀[i]
        distance_profile!(D, Euclidean(), QT, μ, σ, m, i)
        update_min!(P, I, D, i)
        showprogress && i % 5 == 0 && next!(prog)
    end
    Profile(T, P, I, m, nothing)
end


function matrix_profile(A::AbstractVector{<:Number}, T::AbstractVector{<:Number}, m::Int; showprogress=true)
    n   = length(A)
    l   = n-m+1
    lT   = lastlength(T)-m+1
    # n > 2m+1 || throw(ArgumentError("Window length too long, maximum length is $((n+1)÷2)"))
    μT,σT  = running_mean_std(T, m)
    μA,σA  = running_mean_std(A, m)
    QT     = window_dot(getwindow(A, m, 1), T)
    @assert length(QT) == lT
    QT₀  = copy(QT)
    D    = distance_profile(Euclidean(), QT, μA, σA, μT, σT, m)
    P    = copy(D)
    I    = ones(Int, lT)
    prog = Progress((l - 1) ÷ 10, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
    @inbounds for i = 2:l
        for j = lT:-1:2
            @fastmath QT[j] = QT[j-1] - T[j-1] * A[i-1] + T[j+m-1] * A[i+m-1]
        end
        QT[1] = dot(getwindow(A, m, i), getwindow(T,m,1)) #QT₀[i]
        distance_profile!(D, Euclidean(), QT, μA, σA, μT, σT, m, i)
        update_min!(P, I, D, i, false)
        showprogress && i % 5 == 0 && next!(prog)
    end
    Profile(T, P, I, m, A)
end

matrix_profile(T::AbstractVector{<:Number}, m::Int, ::Euclidean; kwargs...) =
    matrix_profile(T, m; kwargs...)
matrix_profile(A::AbstractVector{<:Number}, T::AbstractVector{<:Number}, m::Int, ::Euclidean; kwargs...) =
        matrix_profile(A, T, m; kwargs...)


function distance_profile!(D::AbstractVector{S},::Euclidean, QT::AbstractVector{S}, μ, σ, m::Int, i::Int) where S <: Number
    @assert i <= length(D)
    @avx for j = eachindex(D)
        frac = (QT[j] - m*μ[i]*μ[j]) / (m*σ[i]*σ[j])
        D[j] = sqrt(max(2m*(1-frac), 0))
    end
    D[i] = typemax(eltype(D))
    D
end


function distance_profile!(D::AbstractVector{S},::Euclidean, QT::AbstractVector{S}, μA, σA, μT, σT, m::Int, i::Int) where S <: Number
    @assert i <= length(μA)
    @avx for j = eachindex(D,QT,μT,σT)
        frac = (QT[j] - m*μA[i]*μT[j]) / (m*σA[i]*σT[j])
        D[j] = sqrt(max(2m*(1-frac), 0))
    end
    D
end

distance_profile(::Euclidean,
    QT::AbstractVector{S},
    μ::AbstractVector{S},
    σ::AbstractVector{S},
    m::Int,
) where {S<:Number} = distance_profile!(similar(μ), Euclidean(), QT, μ, σ, m, 1)

distance_profile(::Euclidean, QT::AbstractVector{S}, μA, σA, μT, σT, m::Int) where {S<:Number} =
    distance_profile!(similar(μT), Euclidean(), QT, μA, σA, μT, σT, m, 1)


function update_min!(P, I, D, i, sym=true)
    n = lastlength(P)
    @assert n == lastlength(I) == lastlength(D) "Lengths are not consistent, $(eachindex(P,I,D))"
    @inbounds for j in 1:n
        sym && j == i && continue
        if D[j] < P[j]
            P[j] = D[j]
            I[j] = i
        end
    end
end

"""
Input: A query Q, and a user provided time series T
Output: The dot product between Q and all subsequences in T
"""
function window_dot(Q, T)
    n   = length(T)
    m   = length(Q)
    QT  = conv(reverse(Q), T)
    return QT[m:n]
end

function running_mean_std(x::AbstractArray{T}, m) where T
    @assert length(x) >= m
    n = length(x)-m+1
    s = ss = zero(T)
    μ = similar(x, n)
    σ = similar(x, n)
    @fastmath @inbounds for i = 1:m
        s  += x[i]
        ss += x[i]^2
    end
    μ[1] = s/m
    σ[1] = sqrt(ss/m - μ[1]^2)
    @fastmath @inbounds for i = 1:n-1 # fastmath making it more accurate here as well, but not faster
        s -= x[i]
        ss -= x[i]^2
        s += x[i+m]
        ss += x[i+m]^2
        μ[i+1] = s/m
        σ[i+1] = sqrt(ss/m - μ[i+1]^2)
    end
    μ,σ
end

function moving_mean(x::AbstractArray{T}, m) where T
    filtfilt(ones(m), [m], x)
end

function znorm(x::AbstractVector)
    x = x .- mean(x)
    x ./= std(x, mean=0, corrected=false)
end




"""
    distance_profile(::Euclidean, Q, T)

Compute the z-normalized Euclidean distance profile corresponding to sliding `Q` over `T`
"""
function distance_profile!(D::AbstractVector{S}, ::Euclidean, Q::AbstractVector{S}, T::AbstractVector{S}) where S <: Number
    m = length(Q)
    μ,σ  = running_mean_std(T, m)
    QT   = window_dot(znorm(Q), T)
    @avx for j = eachindex(D)
        frac = QT[j] / (m*σ[j])
        D[j] = sqrt(max(2m*(1-frac), 0))
    end
    D
end
distance_profile(::Euclidean, Q::AbstractArray{S}, T::AbstractArray{S}) where {S} =
    distance_profile!(similar(T, length(T) - length(Q) + 1), Euclidean(), Q, T)

## General data and distance ===================================================

function matrix_profile(T, m::Int, dist; showprogress=true)
    n   = lastlength(T)
    l   = n-m+1
    # n > 2m+1 || throw(ArgumentError("Window length too long, maximum length is $((n+1)÷2)"))
    P    = distance_profile(dist, getwindow(T,m,1), T)
    P[1] = typemax(floattype(P))
    D    = similar(P)
    I    = ones(Int, l)
    prog = Progress((l - 1) ÷ 10, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
    @inbounds for i = 2:l
        Ti = getwindow(T,m,i)
        distance_profile!(D, dist, Ti, T)
        update_min!(P, I, D, i)
        showprogress && i % 5 == 0 && next!(prog)
    end
    Profile(T, P, I, m, nothing)
end

function matrix_profile(A::AbstractArray{S}, B::AbstractArray{S}, m::Int, dist; showprogress=true) where S
    n   = lastlength(A)
    l   = n-m+1
    lT   = lastlength(B)-m+1
    n > 2m+1 || throw(ArgumentError("Window length too long, maximum length is $((n+1)÷2)"))
    P    = distance_profile(dist, getwindow(A,m,1), B)
    # P[1] = typemax(eltype(P))
    D    = similar(P)
    I    = ones(Int, lT)
    prog = Progress((l - 1) ÷ 10, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
    @inbounds for i = 2:l
        Ai = getwindow(A,m,i)
        distance_profile!(D, dist, Ai, B)
        update_min!(P, I, D, i, false)
        showprogress && i % 5 == 0 && next!(prog)
    end
    Profile(B, P, I, m, A)
end


include("motifs.jl")
include("plotting.jl")



"""
    mpdist(A, B, m, k = ((length(A) + length(B)) - 2m) ÷ 20)

The MP distance between `A` and `B` using window length `M` and returning the `k`th smallest value.
"""
function mpdist(A,B,m,d=Euclidean(),k=(length(A)+length(B)-2m)÷20)
    p1 = matrix_profile(A,B,m,d)
    p2 = matrix_profile(B,A,m,d)
    partialsort([p1.P; p2.P], k)
end

"""
    mpdist(MP::AbstractArray, th::Real, N::Int)

The MP distance given a precomputed matrix profile, calculating `k` as `th*N` where `N` is the length of the data used to create `MP`.
"""
function mpdist(MP::AbstractArray, th::Real, N::Int)
    k = ceil(Int, th*N)
    filter!(isfinite, MP)
    if isempty(MP)
        return typemax(eltype(MP))
    elseif length(MP) >= k
        return partialsort(MP, k)
    else
        return maximum(MP)
    end
end

"""
    mpdist_profile(T::AbstractVector, S::Int, m::Int)

All MP distance profiles between subsequences of length `S` in `T` using internal window length `m`.
"""
function mpdist_profile(T::AbstractVector,S::Int, m::Int, d=Euclidean())
    S >= m || throw(ArgumentError("S should be > m"))
    n = lastlength(T)
    pad = S * ceil(Int, n / S) - n
    T = [T;zeros(pad)]
    @showprogress 1 "MP dist profile" map(1:S:n-S) do i
        d = mpdist_profile(getwindow(T,S,i), T, m, d)
    end
end



"""
    mpdist_profile(A::AbstractVector, B::AbstractVector, m::Int)

MP distance profile between two time series using window length `m`.
"""
function mpdist_profile(A::AbstractVector, B::AbstractVector, m::Int, d=Euclidean())
    # % A the longer time series
    # % B the shorter time series
    th = 0.05
    l1 = lastlength(A)
    l2 = lastlength(B)
    if l1 < l2
        A, B = B, A
    end

    D               = mpdistmat(A, B, m, d)
    rows, cols      = size(D)
    moving_means    = similar(D)
    right_marginals = minimum(D, dims=2)
    for i = 1:cols
        moving_means[:,i] = moving_mean(D[:,i], cols)
    end

    l                = lastlength(A)-lastlength(B)+1
    N_right_marginal = lastlength(B)-m+1
    left_marginal    = zeros(N_right_marginal)

    map(1:l) do i
        right_marginal = right_marginals[i:N_right_marginal+i-1]
        left_marginal .= moving_means[i+cols÷2,:]
        m_profile = [left_marginal; right_marginal]
        mpdist(m_profile, th, 2length(B))
    end
end


function mpdistmat(A::AbstractVector, B::AbstractVector, m::Int, d)
    N = lastlength(B)-m +1
    D = similar(A, lastlength(A)-m+1, N)
    for i = 1:N
        distance_profile!(D[!,i], Euclidean(), getwindow(B, m, i), A)
    end
    D
end



"""
    snippets(T, k, S; m = max(S ÷ 10, 4))

Summarize time series `T` by extracting `k` snippets of length `S`
The parameter `m` controls the window length used internally.
"""
function snippets(T, k, S, d=Euclidean(); m = max(S÷10, 4))
    D = mpdist_profile(T,S,m,d)
    Q = fill(Inf, length(D[1]))
    n = length(T)
    minI = 0
    onsets = zeros(Int, k)
    snippet_profiles = similar(D, k)
    for j = 1:k
        minA = Inf
        for i = 1:length(D)
            A = sum(min.(D[i], Q))
            if A < minA
                minA = A
                minI = i
            end
        end
        snippet_profiles[j] = D[minI]
        Q = min.(D[minI], Q)
        onsets[j] = (minI-1)*S+1
    end
    tot_min = reduce((x,y)->min.(x,y), snippet_profiles, init=Inf*D[1])
    Cfracs = map(1:k) do j
        spj = snippet_profiles[j]
        f = count(spj .== tot_min)
        Cfrac = f/(n-S+1)
    end

    profile = Profile(T,tot_min,Int[],S,nothing)
    mot = Subsequence(profile, onsets, "Snippet")
    profile, mot, Cfracs
end




end
