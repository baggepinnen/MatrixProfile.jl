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
    prog = Progress((l - 1) ÷ 5, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
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
    prog = Progress((l - 1) ÷ 5, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
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

function moving_mean!(μ,x::AbstractArray{T}, m) where T
    # filt(ones(m), [m], x)
    @assert length(x) >= m
    n = length(x)-m+1
    s = zero(T)
    @fastmath @inbounds for i = 1:m
        s  += x[i]
    end
    μ[1] = s/m
    @fastmath @inbounds for i = 1:n-1 # fastmath making it more accurate here as well, but not faster
        s -= x[i]
        s += x[i+m]
        μ[i+1] = s/m
    end
    μ
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
    prog = Progress((l - 1) ÷ 5, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
    @inbounds for i = 2:l
        Ti = getwindow(T,m,i)
        distance_profile!(D, dist, Ti, T)
        update_min!(P, I, D, i)
        showprogress && i % 5 == 0 && next!(prog)
    end
    Profile(T, P, I, m, nothing)
end

function matrix_profile(A::AbstractArray{S}, B::AbstractArray{S}, m::Int, dist; showprogress=true) where S
    n  = lastlength(A)
    l  = n-m+1
    lT = lastlength(B)-m+1
    n > 2m+1 || throw(ArgumentError("Window length too long, maximum length is $((n+1)÷2)"))
    P    = distance_profile(dist, getwindow(A,m,1), B)
    # P[1] = typemax(eltype(P))
    D    = similar(P)
    I    = ones(Int, lT)
    prog = Progress((l - 1) ÷ 5, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
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
include("mpdist.jl")



end
