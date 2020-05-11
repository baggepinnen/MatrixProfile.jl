

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
    if isempty(MP)
        return typemax(eltype(MP))
    elseif length(MP) >= k
        return partialsort(MP, k)
    else
        MP = filter(isfinite, MP) # Do not inplace this
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
    N = S-m+1
    D = similar(T, lastlength(T)-m+1, N)
    moving_means = similar(D)
    m_profile = similar(T, 2N)

    @showprogress 1 "MP dist profile" map(1:S:n-S) do i
        d = mpdist_profile(getwindow(T,S,i), T, m, d; D=D, moving_means=moving_means, m_profile=m_profile)
    end
end



"""
    mpdist_profile(A::AbstractVector, B::AbstractVector, m::Int)

MP distance profile between two time series using window length `m`.
"""
function mpdist_profile(
    B::AbstractVector,
    A::AbstractVector,
    m::Int,
    d = Euclidean();
    D = similar(A, lastlength(A)-m+1, lastlength(B)-m+1),
    moving_means = similar(D),
    m_profile = similar(D, 2size(D,2))
)
    mpdistmat!(D,A,B,m,d)
    th = 0.05
    lA = lastlength(A)
    lB = lastlength(B)
    lA >= lB ||
    throw(ArgumentError("The first argument should be shorter than the second."))

    rows, cols = size(D)
    right_margs = vec(minimum(D, dims = 2))
    for i = 1:cols
        moving_mean!(moving_means[!, i], D[!, i], cols)
    end

    l = lA - lB + 1
    map(1:l) do i
        m_profile[1:cols] .= @view moving_means[i+cols÷2, :] # left marg
        copyto!(m_profile, cols+1, right_margs, i, cols)
        mpdist(m_profile, th, 2 * lastlength(B))
    end
end

function mpdistmat(A::AbstractVector, B::AbstractVector, m::Int, d)
    N = lastlength(B)-m +1
    D = similar(A, lastlength(A)-m+1, N)
    mpdistmat!(D, A::AbstractVector, B::AbstractVector, m::Int, d)
end
function mpdistmat!(D, A::AbstractVector, B::AbstractVector, m::Int, d)
    N = lastlength(B)-m +1
    for i = 1:N # Not thread safe
        distance_profile!(D[!,i], Euclidean(), getwindow(B, m, i), A)
    end
    D
end



"""
    snippets(T, k, S; m = max(S ÷ 10, 4))

Summarize time series `T` by extracting `k` snippets of length `S`
The parameter `m` controls the window length used internally.
"""
function snippets(T, k, S, d=Euclidean(); m = max(S÷10, 4), th=S)
    D = mpdist_profile(T,S,m,d)
    INF = typemax(floattype(T))
    Q = fill(INF, length(D[1]))
    n = lastlength(T)
    minI = 0
    onsets = zeros(Int, k)
    snippet_profiles = similar(D, k)
    for j = 1:k
        minA = INF
        for i = 1:length(D)
            A = sum(min.(D[i], Q))
            if A < minA #&& all(abs.(((i-1)*S+1) .- onsets[1:j-1]) .> th)
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
    mot = Subsequence(profile, onsets, "Snippet")
    profile, mot, Cfracs
end
