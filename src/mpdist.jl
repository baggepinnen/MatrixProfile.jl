

"""
    mpdist(A, B, m, d = ZEuclidean(), k = max(1, (length(A) + length(B) - 2m) ÷ 20))

The MP distance between `A` and `B` using window length `m`. The two AB-join matrix profiles `matrix_profile(A, B, m, d)` and `matrix_profile(B, A, m, d)` are concatenated and their `k`th smallest value is returned.

Reference: [Matrix profile XII](https://www.cs.ucr.edu/~eamonn/MPdist_Expanded.pdf).
"""
function mpdist(A,B,m,d=ZEuclidean(),k=max(1, (length(A)+length(B)-2m)÷20))
    p1 = matrix_profile(A,B,m,d)
    p2 = matrix_profile(B,A,m,d)
    partialsort([p1.P; p2.P], k)
end

"""
    mpdist(MP::AbstractArray, th::Real, N::Int)

The MP distance given a precomputed (concatenated) matrix profile, calculating `k` as `ceil(th*N)` where `N` is the length of the data used to create `MP`. If `MP` has fewer than `k` entries, the largest finite entry is returned.
"""
function mpdist(MP::AbstractArray, th::Real, N::Int)
    k = ceil(Int, th*N)
    if !any(isfinite, MP)
        return typemax(eltype(MP))
    elseif length(MP) >= k
        return partialsort(MP, k)
    else
        return maximum(filter(isfinite, MP)) # Do not inplace this
    end
end

"""
    mpdist_profile(T::AbstractVector, S::Int, m::Int, [d = ZEuclidean()])

All MP distance profiles between subsequences of length `S` in `T` using internal window length `m`.

`T` is divided into `length(T) ÷ S` non-overlapping subsequences of length `S`; trailing samples that do not fill a complete subsequence are not used. The result is a vector with one entry per subsequence, each entry is the [`mpdist_profile`](@ref) of that subsequence against `T`.
"""
function mpdist_profile(T::AbstractVector{TT},S::Int, m::Int, d::DT = ZEuclidean()) where {TT,DT}
    S >= m || throw(ArgumentError("S should be ≥ m"))
    prev_fftw_threads = SlidingDistancesBase.DSP.FFTW.get_num_threads()
    SlidingDistancesBase.DSP.FFTW.set_num_threads(1)
    try
        n = length(T)
        # Imani et al., Matrix Profile XIII, Table II: floor(n/S) non-overlapping
        # windows T[(i-1)S+1 : iS]; trailing n mod S samples are not used.
        n_windows = n ÷ S
        last_start = (n_windows - 1) * S + 1
        N = S-m+1
        D = Matrix{float(TT)}(undef, n-m+1, N)
        sliding_means = similar(D)
        m_profile = similar(D, 2N)

        prog = Progress(n_windows, dt=1, desc="MP dist profile")
        map(1:S:last_start) do i
            dp = mpdist_profile(getwindow(T,S,i), T, m, d, D, sliding_means, m_profile)
            next!(prog)
            dp
        end
    finally
        SlidingDistancesBase.DSP.FFTW.set_num_threads(prev_fftw_threads)
    end
end



"""
    mpdist_profile(Q::AbstractVector, T::AbstractVector, m::Int, [d = ZEuclidean()])

MP distance profile of the query `Q` with respect to the time series `T`, using internal window length `m`. Entry `i` of the result is the MP distance between `Q` and the subsequence of `T` of length `length(Q)` starting at `i`. `Q` must not be longer than `T`.
"""
function mpdist_profile(
    B::AbstractVector,
    A::AbstractVector,
    m::Int,
    d = ZEuclidean(),
    D = Matrix{float(eltype(A))}(undef, lastlength(A)-m+1, lastlength(B)-m+1),
    sliding_means = similar(D),
    m_profile = similar(D, 2size(D,2))
)
    lA = lastlength(A)
    lB = lastlength(B)
    lA >= lB ||
    throw(ArgumentError("The first argument should be shorter than the second."))
    mpdistmat!(D,A,B,m,d)
    th = floattype(B)(0.05)

    rows, cols = size(D)
    right_margs = vec(minimum(D, dims = 2))
    for i = 1:cols
        sliding_mean!(sliding_means[!, i], D[!, i], cols)
    end

    l = lA - lB + 1
    map(1:l) do i
        m_profile[1:cols] .= @view sliding_means[i+cols÷2, :] # left marg
        copyto!(m_profile, cols+1, right_margs, i, cols)
        mpdist(m_profile, th, 2 * lastlength(B))
    end
end

function mpdistmat!(D, A::AbstractVector, B::AbstractVector, m::Int, d)
    N = lastlength(B)-m +1
    Threads.@threads for i = 1:N # Not thread safe for more than one FFTW thread
        distance_profile!(D[!,i], d, getwindow(B, m, i), A)
    end
    D
end



"""
    snippets(T, k, S, d = ZEuclidean(); m = max(S ÷ 10, 4), th = Inf)

Summarize time series `T` by extracting `k` snippets of length `S` using
the Snippet-Finder algorithm of Imani et al., *Matrix Profile XIII*.

Returns a [`Snippets`](@ref) object that can be plotted with `plot(snips)`.
Candidate snippets are the `length(T) ÷ S` non-overlapping subsequences of
length `S`, so `T` must contain at least `k` of them.

# Keyword arguments
- `m`: window length used internally by `mpdist_profile`.
- `th`: **non-standard extension** to the paper. When `th < Inf`, each
  newly selected snippet onset is required to be more than `th` samples
  away from every previously selected onset. The default `Inf` disables
  the constraint and reproduces the paper's algorithm verbatim (Table I,
  line 7). Use a finite `th` only if you specifically want to force a
  minimum temporal spacing between snippets — it is not part of
  Snippet-Finder as published.
"""
function snippets(T, k, S, d=ZEuclidean(); m = max(S÷10, 4), th=Inf)
    lastlength(T) >= k*S || throw(ArgumentError("T must contain at least k = $k subsequences of length S = $S"))
    D      = mpdist_profile(T,S,m,d)
    INF    = typemax(floattype(T))
    Q      = fill(INF, length(D[1]))
    n      = lastlength(T)
    minI   = 0
    onsets = zeros(Int, k)
    snippet_profiles = similar(D, k)
    constrain_spacing = th < Inf
    for j = 1:k
        minA = INF
        for i = 1:length(D)
            A = sum(min.(D[i], Q))
            if A < minA && (!constrain_spacing || all(abs.(((i-1)*S+1) .- onsets[1:j-1]) .> th))
                minA = A
                minI = i
            end
        end
        snippet_profiles[j] = D[minI]
        Q .= min.(D[minI], Q)
        onsets[j] = (minI-1)*S+1
    end
    tot_min = reduce((x,y)->min.(x,y), snippet_profiles, init=snippet_profiles[1])
    Cfracs = map(1:k) do j
        spj = snippet_profiles[j]
        f = count(spj .== tot_min)
        Cfrac = f/length(spj)
    end

    profile = Profile(T,tot_min,Int[],S,nothing)
    snips = Subsequence(profile, onsets, "Snippet")
    Snippets(T,tot_min, snips, Cfracs, snippet_profiles, S)
end


"""
    Snippets

Result of [`snippets`](@ref).

# Fields
- `T`: the time series.
- `minprofile`: the elementwise minimum of the MP distance profiles of the selected snippets.
- `snippets`: vector of `Subsequence`s holding the selected snippets and their onsets.
- `Cfracs`: the fraction of `T` represented by each snippet.
- `profiles`: the MP distance profile of each snippet.
- `S`: the snippet length.

`segment(snips)` returns, for every index of `minprofile`, the number of the snippet that represents `T` at that index.
"""
struct Snippets
    T
    minprofile
    snippets
    Cfracs
    profiles
    S
end

"""
    segment(s::Snippets)

Return a vector holding, for every index of `s.minprofile`, the number of the snippet that represents the time series at that index.
"""
segment(s::Snippets) = [findfirst(==(s.minprofile[i]), getindex.(s.profiles, i)) for i in eachindex(s.minprofile)]
