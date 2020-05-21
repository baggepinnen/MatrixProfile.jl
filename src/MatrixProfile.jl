module MatrixProfile
using Statistics, LinearAlgebra
using LoopVectorization
using ProgressMeter
using RecipesBase

using Distances

using SlidingDistancesBase
import SlidingDistancesBase: floattype, lastlength, distance_profile, distance_profile!, window_dot
export Euclidean, ZEuclidean

export matrix_profile, distance_profile, motifs, anomalies, mpdist, mpdist_profile, onsets, snippets, seqs, resample


struct Profile{TT,TP,QT}
    T::TT
    P::TP
    I::Vector{Int}
    m::Int
    Q::QT
end

"""
    profile = matrix_profile(T, m, [dist = ZEuclidean()]; showprogress=true)

Return the matrix profile and the profile indices of time series `T` with window length `m`. See fields `profile.P, profile.I`. You can also plot the profile. If `dist = ZEuclidean()` the STOMP algorithm will be used.

Reference: [Matrix profile II](https://www.cs.ucr.edu/~eamonn/STOMP_GPU_final_submission_camera_ready.pdf).
"""
function matrix_profile(T::AbstractVector{<:Number}, m::Int; showprogress=true)
    n   = lastlength(T)
    l   = n-m+1
    n > 2m+1 || throw(ArgumentError("Window length too long, maximum length is $((n+1)÷2)"))
    μ,σ  = sliding_mean_std(T, m)
    QT   = window_dot(getwindow(T, m, 1), T)
    QT₀   = copy(QT)
    D    = distance_profile(ZEuclidean(), QT, μ, σ, m)
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
        distance_profile!(D, ZEuclidean(), QT, μ, σ, m, i)
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
    μT,σT  = sliding_mean_std(T, m)
    μA,σA  = sliding_mean_std(A, m)
    QT     = window_dot(getwindow(A, m, 1), T)
    @assert length(QT) == lT
    QT₀  = copy(QT)
    D    = distance_profile(ZEuclidean(), QT, μA, σA, μT, σT, m)
    P    = copy(D)
    I    = ones(Int, lT)
    prog = Progress((l - 1) ÷ 5, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
    @inbounds for i = 2:l
        for j = lT:-1:2
            @fastmath QT[j] = QT[j-1] - T[j-1] * A[i-1] + T[j+m-1] * A[i+m-1]
        end
        QT[1] = dot(getwindow(A, m, i), getwindow(T,m,1)) #QT₀[i]
        distance_profile!(D, ZEuclidean(), QT, μA, σA, μT, σT, m, i)
        update_min!(P, I, D, i, false)
        showprogress && i % 5 == 0 && next!(prog)
    end
    Profile(T, P, I, m, A)
end

matrix_profile(T::AbstractVector{<:Number}, m::Int, ::ZEuclidean; kwargs...) =
    matrix_profile(T, m; kwargs...)
matrix_profile(A::AbstractVector{<:Number}, T::AbstractVector{<:Number}, m::Int, ::ZEuclidean; kwargs...) =
        matrix_profile(A, T, m; kwargs...)



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
include("mpdist.jl")
include("plotting.jl")



end
