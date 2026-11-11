module MatrixProfile
using Statistics, LinearAlgebra
using LoopVectorization
using ProgressMeter
using RecipesBase

using Distances

using SlidingDistancesBase
import SlidingDistancesBase: floattype, lastlength, distance_profile, distance_profile!, window_dot
export Euclidean, ZEuclidean

export matrix_profile, distance_profile, motifs, anomalies, mpdist, mpdist_profile, onsets, snippets, seqs, segment, segment_profile
export damp, mass


"""
    Profile

Result of [`matrix_profile`](@ref) and [`damp`](@ref).

# Fields
- `T`: the time series the profile is indexed by. For an AB-join `matrix_profile(A, B, m)`, this is `B`.
- `P`: the matrix profile. `P[j]` is the distance between the window starting at `j` and its nearest neighbor.
- `I`: the profile index. `I[j]` is the onset of the nearest neighbor of the window starting at `j`. For an AB-join, this is an index into `A`. For [`damp`](@ref), this field holds the pruning vector.
- `m`: the window length.
- `Q`: the query series `A` of an AB-join, `nothing` otherwise.

A profile can be plotted with `plot(profile)` and passed to [`motifs`](@ref), [`anomalies`](@ref), [`segment`](@ref) and [`segment_profile`](@ref).
"""
struct Profile{TT,TP,QT}
    T::TT
    P::TP
    I::Vector{Int}
    m::Int
    Q::QT
end

"""
    profile = matrix_profile(T, m, [dist = ZEuclidean()]; showprogress=true, exclusion_zone = 0)

Return the matrix profile and the profile indices of time series `T` with window length `m`. See fields `profile.P, profile.I`. You can also plot the profile. If `dist = ZEuclidean()` the STOMP algorithm will be used.

- `exclusion_zone` denotes an integer number of samples around the trivial match to avoid. The paper suggests using `exclusion_zone = m ÷ 4`. This is likely most beneficial for time-series dominated by low frequencies. Only supported for `dist = ZEuclidean()`.
- `showprogress` determines whether a progress bar is displayed.

The window length must satisfy `length(T) > 2m + 1`. Integer-valued series are converted to floating point.

Windows that contain non-finite values (`NaN`, `Inf`) are at infinite distance from all other windows. Two constant windows are at distance zero, and a constant and a non-constant window are at distance `√m`, the distance between a z-normalized window and the zero vector. These conventions agree with those of the Python package stumpy.

Reference: [Matrix profile II](https://www.cs.ucr.edu/~eamonn/STOMP_GPU_final_submission_camera_ready.pdf).
"""
function matrix_profile(T0::AbstractVector{<:Number}, m::Int; showprogress=true, exclusion_zone=0)
    n   = lastlength(T0)
    l   = n-m+1
    n > 2m+1 || throw(ArgumentError("Window length too long, maximum length is $((n-2)÷2)"))
    T, finite, constant = degenerate_windows(T0, m)
    special = findall(.!finite .| constant)
    degenerate = !isempty(special)
    μ,σ  = sliding_mean_std(T, m)
    QT   = window_dot(getwindow(T, m, 1), T)
    QT₀   = copy(QT)
    D    = distance_profile!(similar(μ), ZEuclidean(), QT, μ, σ, m, 1; exclusion_zone)
    degenerate && correct_degenerate!(D, finite[1], constant[1], finite, constant, special, m, 1-exclusion_zone:1+exclusion_zone)
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
        distance_profile!(D, ZEuclidean(), QT, μ, σ, m, i; exclusion_zone)
        degenerate && correct_degenerate!(D, finite[i], constant[i], finite, constant, special, m, i-exclusion_zone:i+exclusion_zone)
        update_min!(P, I, D, i)
        showprogress && i % 5 == 0 && next!(prog)
    end
    Profile(T0, P, I, m, nothing)
end

"""
    T, finite, constant = degenerate_windows(T0, m)

Return `T0` with non-finite values replaced by zero, so that they do not propagate through the recursive update of the sliding dot products, together with indicators of the windows of length `m` that contain only finite values and of the windows that are constant.
"""
function degenerate_windows(T0::AbstractVector, m::Int)
    n = length(T0)
    l = n-m+1
    nonfinite = cumsum(.!isfinite.(T0))
    finite    = [nonfinite[j+m-1] - (j > 1 ? nonfinite[j-1] : 0) == 0 for j in 1:l]
    runlength = ones(Int, n) # Number of consecutive equal values ending at each index
    for k in 2:n
        T0[k] == T0[k-1] && (runlength[k] = runlength[k-1] + 1)
    end
    constant = [runlength[j+m-1] >= m for j in 1:l]
    T = all(finite) ? T0 : map(x -> isfinite(x) ? x : zero(x), T0)
    T, finite, constant
end

"""
    correct_degenerate!(D, qfinite, qconstant, finite, constant, special, m, excluded)

Correct the distance profile `D` of a query window with properties `qfinite, qconstant` against windows with properties `finite, constant`, see [`degenerate_windows`](@ref). `special` holds the indices of the windows that are non-finite or constant; if the query window is neither, only these entries are corrected. Entries with indices in `excluded` (the exclusion zone of a self-join) are not modified.
"""
function correct_degenerate!(D, qfinite, qconstant, finite, constant, special, m, excluded)
    S = eltype(D)
    @inbounds for j in (qfinite && !qconstant ? special : eachindex(D))
        j ∈ excluded && continue
        if !(qfinite && finite[j])
            D[j] = typemax(S)
        elseif qconstant && constant[j]
            D[j] = zero(S)
        elseif qconstant || constant[j]
            D[j] = sqrt(S(m))
        end
    end
    D
end

"""
    profile = matrix_profile(A, T, m, [dist = ZEuclidean()]; showprogress=true)

Matrix profile of an AB-join between `A` and `T`. Every window of `T` is compared to all windows of `A`; `profile.P[j]` is the distance between the window of `T` starting at `j` and its nearest neighbor in `A`, and `profile.I[j]` is the onset of that neighbor in `A`. The result is in general not symmetric, `matrix_profile(A, T, m) != matrix_profile(T, A, m)`.

Both series must have a length of at least `m`. Integer-valued series are converted to floating point. Windows with non-finite values and constant windows are treated as described in the self-join method above.
"""
function matrix_profile(A0::AbstractVector{S}, T0::AbstractVector{S}, m::Int; showprogress=true) where S <: Number
    n   = length(A0)
    l   = n-m+1
    n >= m || throw(ArgumentError("A is shorter than the window length m = $m"))
    lastlength(T0) >= m || throw(ArgumentError("T is shorter than the window length m = $m"))
    A, finiteA, constantA = degenerate_windows(A0, m)
    T, finiteT, constantT = degenerate_windows(T0, m)
    specialT = findall(.!finiteT .| constantT)
    degenerate = !isempty(specialT) || !all(finiteA) || any(constantA)
    lT   = lastlength(T)-m+1
    μT,σT  = sliding_mean_std(T, m)
    μA,σA  = sliding_mean_std(A, m)
    QT     = window_dot(getwindow(A, m, 1), T)
    @assert length(QT) == lT
    QT₀  = copy(QT)
    D    = distance_profile(ZEuclidean(), QT, μA, σA, μT, σT, m)
    degenerate && correct_degenerate!(D, finiteA[1], constantA[1], finiteT, constantT, specialT, m, 1:0)
    P    = copy(D)
    I    = ones(Int, lT)
    prog = Progress((l - 1) ÷ 5, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
    @inbounds for i = 2:l
        for j = lT:-1:2
            @fastmath QT[j] = QT[j-1] - T[j-1] * A[i-1] + T[j+m-1] * A[i+m-1]
        end
        QT[1] = dot(getwindow(A, m, i), getwindow(T,m,1)) #QT₀[i]
        distance_profile!(D, ZEuclidean(), QT, μA, σA, μT, σT, m, i)
        degenerate && correct_degenerate!(D, finiteA[i], constantA[i], finiteT, constantT, specialT, m, 1:0)
        update_min!(P, I, D, i, false)
        showprogress && i % 5 == 0 && next!(prog)
    end
    Profile(T0, P, I, m, A0)
end

matrix_profile(T::AbstractVector{<:Integer}, m::Int; kwargs...) =
    matrix_profile(float.(T), m; kwargs...)
function matrix_profile(A::AbstractVector{<:Number}, T::AbstractVector{<:Number}, m::Int; kwargs...)
    S = promote_type(float(eltype(A)), float(eltype(T)))
    matrix_profile(S.(A), S.(T), m; kwargs...)
end
matrix_profile(A::AbstractVector{S}, T::AbstractVector{S}, m::Int; kwargs...) where S <: Integer =
    matrix_profile(float.(A), float.(T), m; kwargs...)

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
    # The 2m+1 lower bound applies only to the self-join; for an AB-join A
    # is independent of B and can be as short as m. See the ZEuclidean
    # AB-join above for the same reason.
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
    mass(Q::Vector, T::Vector)

Compute the Z-normalized [`distance_profile`](@ref) using FFT. Equivalent to
```
d = distance_profile(ZEuclidean(), Q, T)
```

# Arguments:
- `Q`: query
- `T`: Data
"""
mass(Q::AbstractVector, T::AbstractVector) = distance_profile(ZEuclidean(), Q, T)

# The implementation below seems to perform much worse than the one above
# function mass(y::AbstractVector{T}, x::AbstractVector{T}, k=length(y)) where T
#     n = length(x)
#     m = length(y)
#     m <= k || throw(ArgumentError("k = $k is too small, must be at least length(y) = $m"))
#     dist = T[]
    
#     μy = mean(y)
#     σy = std(y, corrected=false)
#     μx, σx = sliding_mean_std(x, m-1, 0)
#     @assert length(μx) == length(σx) == n
       
#     y = reverse(y)
#     y = [y; zeros(T, length(m+1:k))]
    
#     j = 0
#     P = @views plan_rfft(x[1:k])
#     O = zero(T)
#     local iP
#     Y = P * y
#     Z = similar(Y)
#     for outer j = 1:k-m+1:n-k+1
#         X = P * x[j:j+k-1]
#         Z .= X.*Y
#         if j == 1
#             iP = @views plan_irfft(Z, length(y))
#         end
#         z = iP * Z
        
#         d = @views @. sqrt(max(2*(m-(z[m:k]-m*μx[m+j-1:j+k-1]*μy)/(σx[m+j-1:j+k-1]*σy)), O))
#         append!(dist, d)
#     end
    
#     j = j+k-m
#     k = n-j
#     if k >= m
#         error()
#         X = @views rfft(x[j+1:n])
#         y = y[1:k]
#         Y = rfft(y)
#         Z = Y
#         Z .= X.*Y
#         z = irfft(Z, length(y))
        
#         @views d = @. sqrt(max(2*(m-(z[m:k]-m*μx[j+m:n]*μy)/(σx[j+m:n]*σy)), O))
#         append!(dist, d)
#     end
#     return dist
# end


"""
    profile = damp(T, m, ind = max(length(T) ÷ 10, m))

DAMP algorithm from https://www.cs.ucr.edu/~eamonn/DAMP_long_version.pdf

# Arguments:
- `T`: Time series
- `m`: Subsequence length (choose as approximate period of repeating patterns)
- `ind`: Location of split point between training and test data, defaults to 10% of the data. Must satisfy `m ≤ ind ≤ length(T) - m + 1`.

# Returns:
A [`Profile`](@ref) where `profile.P` is the left approximate matrix profile. Large values indicate anomalies. Entries before `ind` belong to the training data and are zero. `profile.I` holds the pruning vector, an entry of zero indicates that the corresponding subsequence was pruned.
"""
function damp(T, m, ind = max(length(T) ÷ 10, m))
    m < length(T) || throw(ArgumentError("m must be (quite a bit) less than length(T)"))
    m <= ind <= length(T) - m + 1 || throw(ArgumentError("ind must satisfy m ≤ ind ≤ length(T) - m + 1"))
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
    aMPi, BSF = _dampb(T, m, i, BSF)

Backward processing of damp algorithm

# Arguments:
- `T`: Time series
- `m`: Subsequence length
- `i`: Index of current query
- `BSF`: Highest discord score so far

# Returns:
- `aMPi`: Discord value at position i
- `BSF`: Possibly updated best-so-far discord score
"""
function _dampb(T,m,i,BSF)
    aMPi = Inf
    prefix = nextpow(2, m) # Initial length of prefix
    @views while aMPi ≥ BSF
        if i-prefix+1 <= 1 #the search reaches the beginning of the time series
            aMPi = minimum(mass(T[i:i+m-1], T[1:i]))
            if aMPi > BSF # Update the current best discord score
                BSF = aMPi
            end
            break
        else
            # @show length(T), prefix, i, m
            aMPi = minimum(mass(T[i:i+m-1], T[i-prefix+1:i]))
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

    if last - start + 1 >= m # at least one full subsequence in the lookahead
        @views Di = mass(T[i:i+m-1], T[start:last])
        indices = findall(Di .< BSF)
        indices .+= (start - 1)
        PV[indices] .= false
    end
end

end

