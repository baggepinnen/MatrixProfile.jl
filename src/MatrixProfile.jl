module MatrixProfile
using Statistics, LinearAlgebra
using LoopVectorization
using ProgressMeter
using RecipesBase

using Distances

using SlidingDistancesBase
import SlidingDistancesBase: floattype, lastlength, distance_profile, distance_profile!, window_dot
export Euclidean, ZEuclidean

export matrix_profile, distance_profile, motifs, anomalies, mpdist, mpdist_profile, onsets, snippets, seqs, resample, segment, segment_profile
export damp, mass


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

"""
    profile = matrix_profile(A, T, m, dist=ZEuclidean(); showprogress=true)

Mutual matrix profile between `A` and `T`.
"""
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
include("segment.jl")
include("plotting.jl")

using AbstractFFTs

"""
    mass(x::AbstractVector{T}, y::AbstractVector{T}, k)

# Arguments:
- `x`: Data
- `y`: query
- `k`: window size, must be at least `length(y)`
"""
function mass(x::AbstractVector{T}, y::AbstractVector{T}, k=length(y)) where T
    n = length(x)
    m = length(y)
    m <= k || throw(ArgumentError("k = $k is too small, must be at least length(y) = $m"))
    dist = T[]
    
    μy = mean(y)
    σy = std(y, corrected=false)
    μx, σx = sliding_mean_std(x, m-1, 0)
    @assert length(μx) == length(σx) == n
       
    y = reverse(y)
    y = [y; zeros(T, length(m+1:k))]
    
    j = 0
    P = @views plan_rfft(x[1:k])
    O = zero(T)
    local iP
    Y = P * y
    Z = similar(Y)
    for outer j = 1:k-m+1:n-k+1
        X = P * x[j:j+k-1]
        Z .= X.*Y
        if j == 1
            iP = @views plan_irfft(Z, length(y))
        end
        z = iP * Z
        
        d = @views @. sqrt(max(2*(m-(z[m:k]-m*μx[m+j-1:j+k-1]*μy)/(σx[m+j-1:j+k-1]*σy)), O))
        append!(dist, d)
    end
    
    j = j+k-m
    k = n-j
    if k >= m
        X = @views rfft(x[j+1:n])
        y = y[1:k]
        Y = rfft(y)
        Z = Y
        Z .= X.*Y
        z = irfft(Z, length(y))
        
        @views d = @. sqrt(max(2*(m-(z[m:k]-m*μx[j+m:n]*μy)/(σx[j+m:n]*σy)), O))
        append!(dist, d)
    end
    return dist
end




# x = randn(100000000);
# y = rand(400);
# k = 2^23;
# @time p = MatrixProfile.mass(x,y,k);

# x = randn(1000000)
# y = rand(16)
# k = 2^5
# @time p = MatrixProfile.mass(x,y,k);


"""
    damp(T, m, ind = length(T) ÷ 10)

DAMP algorithm from https://www.cs.ucr.edu/~eamonn/DAMP_long_version.pdf

# Arguments:
- `T`: Time series
- `m`: Subsequence length (choose as approximate period of repeating patterns)
- `ind`: Location of split point between training and test data, defaults to 10% of the data.

# Returns:
-  `aMP`: Left approximate Matrix Profile. Large values indicate anomalies.
"""
function damp(T, m, ind = length(T) ÷ 10)
    ind < length(T) || throw(ArgumentError("ind must be less than length(T)"))
    m < length(T) || throw(ArgumentError("m must be (quite a bit) less than length(T)"))
    PV = ones(Int, length(T)-m+1)
    aMP = zeros(length(T)-m+1)
    BSF = zero(eltype(T)) # The current best discord score
    # Scan all subsequences in the test data
    for i = ind:length(T) - m + 1
        if PV[i] != 1 # Skip the pruned subsequence
            aMP[i] = aMP[i-1]
        else
            aMP[i], BSF = _dampb(T, m, i, BSF)
            _dampf!(T, m, i, BSF, PV)
        end
    end
    Profile(T, aMP, PV, m, nothing)
end


"""
    aMP[i], BSF = dampb(T, m, i, BSF)

Backward processing of damp algorithm

# Arguments:
- `T`: Time series
- `m`: Subsequence
- `i`: Index of current query
- `BSF`: Highest discord score so far
- `aMPi`: Discord value at position i
"""
function _dampb(T,m,i,BSF)
    aMPi = Inf
    prefix = nextpow(2, m) # Initial length of prefix
    while aMPi ≥ BSF
        if i-prefix+1 <= 1 #the search reaches the beginning of the time series
            aMPi = minimum(mass(T[1:i],T[i:i+m-1]))
            if aMPi > BSF # Update the current best discord score
                BSF = aMPi
            end
            break
        else
            # @show length(T), prefix, i, m
            aMPi = minimum(mass(T[i-prefix+1:i],T[i:i+m-1]))
            if aMPi < BSF
                break # Stop searching
            else # Double the length of prefix
                prefix *= 2
            end
        end
    end
    return aMPi, BSF
end

"""
    _dampf!(T, m, i, BSF, PV)

Forward processing of damp algorithm

# Arguments:
- `T`: Time series
- `m`: Subsequence length
- `i`: Index of current query
- `BSF`: Highest discord score so far
- `PV`: Pruned Vector
"""
function _dampf!(T, m, i, BSF, PV)
    lookahead = nextpow(2, m)
    start = i + m
    last = min(start + lookahead - 1, length(T))

    if last < length(T)#  && i+m-1 <= length(T) #the search does not reach the end of the time series
        Di = mass(T[start:last], T[i:i+m-1])
        indices = findall(Di .< BSF)
        indices .+= (start - 1)
        PV[indices] .= false
    end
end

end

