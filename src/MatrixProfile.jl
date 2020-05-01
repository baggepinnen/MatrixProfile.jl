module MatrixProfile
using Statistics, LinearAlgebra
using DSP
using LoopVectorization
using ProgressMeter
using RecipesBase

export matrix_profile, distance_profile, motifs


struct Profile{TT,TP}
    T::TT
    P::TP
    I::Vector{Int}
    m::Int
end

"""
    profile = matrix_profile(T, m; showprogress=true)

Return the matrix profile and the profile indices of time series `T` with window length `m`. See fields `profile.P, profile.I`. You can also plot the profile. Uses the STOMP algorithm.

Reference: [Matrix profile II](https://www.cs.ucr.edu/~eamonn/STOMP_GPU_final_submission_camera_ready.pdf).
"""
function matrix_profile(T, m; showprogress=true)
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
            @fastmath QT[j] = QT[j-1]-T[j-1]*T[i-1]+T[j+m-1]*T[i+m-1]
        end
        QT[1] = QT₀[i]
        distance_profile!(D, QT, μ, σ, m, i)
        update_min!(P, I, D, i)
        showprogress && i % 10 == 0 && next!(prog)
    end
    Profile(T, P, I, m)
end

"""
    distance_profile(Q, T)

Compute the Euclidean distance profile corresponding to sliding `Q` over `T`
"""
function distance_profile!(D, Q, T)
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

function distance_profile!(D, QT, μ, σ, m, i)
    @assert i <= length(D)
    @avx for j = eachindex(D)
        frac = (QT[j] - m*μ[i]*μ[j]) / (m*σ[i]*σ[j])
        D[j] = sqrt(max(2m*(1-frac), 0))
    end
    D[i] = typemax(eltype(D))
    D
end

distance_profile(QT, μ, σ, m) = distance_profile!(similar(QT), QT, μ, σ, m, 1)


function update_min!(P, I, D, i)
    @inbounds for j in eachindex(P,I,D)
        j == i && continue
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
# function window_dot(Q, T::AbstractArray{S}) where S
#     n    = length(T)
#     m    = length(Q)
#     Ta   = [T; zeros(S,n)]
#     Qr   = reverse(Q)
#     Qra  = [Qr; zeros(S,2n-m)]
#     Qraf = rfft(Qra)
#     Taf  = rfft(Ta)
#     QT   = irfft(Qraf .* Taf, length(Qra))
#     return QT[m:n]
# end


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

function remove_found_motifs!(P, found_motifs, th)
    for mo in found_motifs
        for o in mo.onsets
            inds = max(1, o-th):min(length(P), o+th)
            P[inds] .= typemax(eltype(P))
        end
    end
end

function motifs(p::Profile, k, r, th = 0, found_motifs = Motif[])
    length(found_motifs) == k && return found_motifs
    m = p.m
    P = copy(p.P)
    remove_found_motifs!(P, found_motifs, th)
    d, i = findmin(P)
    inds = 0:m-1
    distance_profile!(P, view(p.T, inds .+ i), p.T)
    remove_found_motifs!(P, found_motifs, th)
    onsets = findall(<=(d * r), P)
    push!(found_motifs, Motif(p, onsets))
    motifs(p::Profile, k, r, th, found_motifs)
end

struct Motif
    motifs
    onsets
end

function Motif(p::Profile, onsets)
    m = p.m
    inds = 0:m-1
    n = length(onsets)
    vecs = map(1:n) do i
        p.T[inds .+ onsets[i]]
    end
    Motif(vecs, onsets)
end


function znorm(x)
    x = x .- mean(x)
    x ./= std(x, mean=0, corrected=false)
end

function dist(p::Profile, i, j)
    inds = 0:p.m-1
    sum(abs2, znorm(p.T[i .+ inds]) - znorm(p.T[j .+ inds]))
end


@recipe function plot(p::Profile)
    link --> :x
    layout --> (2,1)
    legend --> false
    @series begin
        title --> "T"
        label --> "T"
        subplot --> 1
        p.T
    end
    @series begin
        title --> "Matrix profile"
        subplot --> 2
        p.P
    end
end

_append_inf(x) = vec([x; fill(Inf, 1, size(x,2))])

@recipe function plot(p::Profile, motifs::Vector{Motif})

    @series begin
        p
    end

    inds = 0:p.m-1
    linewidth --> 2
    for (j,m) in enumerate(motifs)
        @series begin
            legend --> true
            group := j
            label --> j
            subplot := 1
            _append_inf(m.onsets' .+ inds), _append_inf(reduce(hcat, m.motifs))
        end
    end

end

end
