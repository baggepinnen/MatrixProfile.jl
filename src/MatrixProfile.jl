module MatrixProfile
using Statistics, LinearAlgebra
using DSP
using LoopVectorization
using ProgressMeter
using RecipesBase

using SlidingDistancesBase
import SlidingDistancesBase: floattype, lastlength, distance_profile, distance_profile!

export matrix_profile, distance_profile, motifs, anomalies


struct Profile{TT,TP,BT}
    T::TT
    P::TP
    I::Vector{Int}
    m::Int
    B::BT
end

"""
    profile = matrix_profile(T, m; showprogress=true)

Return the matrix profile and the profile indices of time series `T` with window length `m`. See fields `profile.P, profile.I`. You can also plot the profile. Uses the STOMP algorithm.

Reference: [Matrix profile II](https://www.cs.ucr.edu/~eamonn/STOMP_GPU_final_submission_camera_ready.pdf).
"""
function matrix_profile(T::AbstractVector{<:Number}, m::Int; showprogress=true)
    n   = length(T)
    l   = n-m+1
    n > 2m+1 || throw(ArgumentError("Window length too long, maximum length is $(2m)"))
    μ,σ  = running_mean_std(T, m)
    QT   = window_dot(view(T, 1:m), T)
    QT₀  = copy(QT)
    D    = distance_profile(QT, μ, σ, m)
    P    = copy(D)
    I    = ones(Int, l)
    prog = Progress((l - 1) ÷ 10, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
    @inbounds for i = 2:l
        for j = l:-1:2
            QT[j] = QT[j-1] - T[j-1] * T[i-1] + T[j+m-1] * T[i+m-1]
        end
        QT[1] = QT₀[i]
        distance_profile!(D, QT, μ, σ, m, i)
        update_min!(P, I, D, i)
        showprogress && i % 10 == 0 && next!(prog)
    end
    Profile(T, P, I, m, nothing)
end


function matrix_profile(A::AbstractVector{<:Number}, T::AbstractVector{<:Number}, m::Int; showprogress=true)
    n   = length(A)
    l   = n-m+1
    lT   = length(T)-m+1
    n > 2m+1 || throw(ArgumentError("Window length too long, maximum length is $(2m)"))
    μT,σT  = running_mean_std(T, m)
    μA,σA  = running_mean_std(A, m)
    QT   = window_dot(view(A, 1:m), T)
    QT₀  = copy(QT)
    D    = distance_profile(QT, μA, σA, μT, σT, m)
    P    = copy(D)
    I    = ones(Int, lT)
    prog = Progress((l - 1) ÷ 10, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
    @inbounds for i = 2:l
        for j = lT:-1:2
            QT[j] = QT[j-1] - T[j-1] * A[i-1] + T[j+m-1] * A[i+m-1]
        end
        # QT[1] = QT₀[i]
        distance_profile!(D, QT, μA, σA, μT, σT, m, i)
        update_min!(P, I, D, i, false)
        showprogress && i % 10 == 0 && next!(prog)
    end
    Profile(T, P, I, m, A)
end


function distance_profile!(D::AbstractVector{S}, QT::AbstractVector{S}, μ, σ, m::Int, i::Int) where S <: Number
    @assert i <= length(D)
    @avx for j = eachindex(D)
        frac = (QT[j] - m*μ[i]*μ[j]) / (m*σ[i]*σ[j])
        D[j] = sqrt(max(2m*(1-frac), 0))
    end
    D[i] = typemax(eltype(D))
    D
end


function distance_profile!(D::AbstractVector{S}, QT::AbstractVector{S}, μA, σA, μT, σT, m::Int, i::Int) where S <: Number
    @assert i <= length(μA)
    @avx for j = eachindex(D,QT)
        frac = (QT[j] - m*μA[i]*μT[j]) / (m*σA[i]*σT[j])
        D[j] = sqrt(max(2m*(1-frac), 0))
    end
    D
end

distance_profile(
    QT::AbstractVector{S},
    μ::AbstractVector{S},
    σ::AbstractVector{S},
    m::Int,
) where {S<:Number} = distance_profile!(similar(μ), QT, μ, σ, m, 1)

distance_profile(QT::AbstractVector{S}, μA, σA, μT, σT, m::Int) where {S<:Number} =
    distance_profile!(similar(μT), QT, μA, σA, μT, σT, m, 1)


function update_min!(P, I, D, i, sym=true)
    @inbounds for j in eachindex(P,I,D)
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
    @assert length(x) > m
    n = length(x)-m+1
    s = ss = zero(T)
    μ = similar(x, n)
    σ = similar(x, n)
    @inbounds for i = 1:m
        s  += x[i]
        ss += x[i]^2
    end
    μ[1] = s/m
    σ[1] = sqrt(ss/m - μ[1]^2)
    @inbounds for i = 1:n-1
        s -= x[i]
        ss -= x[i]^2
        s += x[i+m]
        ss += x[i+m]^2
        μ[i+1] = s/m
        σ[i+1] = sqrt(ss/m - μ[i+1]^2)
    end
    μ,σ
end

function znorm(x::AbstractVector)
    x = x .- mean(x)
    x ./= std(x, mean=0, corrected=false)
end



## General data and distance ===================================================

"""
    distance_profile(Q, T)

Compute the z-normalized Euclidean distance profile corresponding to sliding `Q` over `T`
"""
function distance_profile!(D::AbstractVector{S}, Q::AbstractVector{S}, T::AbstractVector{S}) where S <: Number
    m = length(Q)
    μ,σ  = running_mean_std(T, m)
    QT   = window_dot(znorm(Q), T)
    @avx for j = eachindex(D)
        frac = QT[j] / (m*σ[j])
        D[j] = sqrt(max(2m*(1-frac), 0))
    end
    D
end
distance_profile(Q, T) = distance_profile!(similar(T, length(T)-length(Q)+1), Q, T)


function matrix_profile(T, m::Int, dist; showprogress=true)
    n   = length(T)
    l   = n-m+1
    n > 2m+1 || throw(ArgumentError("Window length too long, maximum length is $(2m)"))
    P    = distance_profile(dist, getwindow(T,m,1), T)
    P[1] = typemax(eltype(P))
    D    = similar(P)
    I    = ones(Int, l)
    prog = Progress((l - 1) ÷ 10, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
    @inbounds for i = 2:l
        Ti = getwindow(T,m,i)
        distance_profile!(D, dist, Ti, T)
        update_min!(P, I, D, i)
        showprogress && i % 10 == 0 && next!(prog)
    end
    Profile(T, P, I, m, nothing)
end

function matrix_profile(A::AbstractArray{S}, B::AbstractArray{S}, m::Int, dist; showprogress=true) where S
    n   = length(A)
    l   = n-m+1
    n > 2m+1 || throw(ArgumentError("Window length too long, maximum length is $(2m)"))
    P    = distance_profile(dist, getwindow(A,m,1), B)
    # P[1] = typemax(eltype(P))
    D    = similar(P)
    I    = ones(Int, l)
    prog = Progress((l - 1) ÷ 10, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
    @inbounds for i = 2:l
        Ai = getwindow(A,m,i)
        distance_profile!(D, dist, Ai, B)
        update_min!(P, I, D, i, false)
        showprogress && i % 10 == 0 && next!(prog)
    end
    Profile(A, P, I, m, B)
end


include("motifs.jl")
include("plotting.jl")


end
