

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
    right_marginals = vec(minimum(D, dims=2))
    for i = 1:cols
        moving_mean!(moving_means[!,i], D[!,i], cols)
    end

    l                = lastlength(A)-lastlength(B)+1
    N_right_marginal = lastlength(B)-m+1
    left_marginal    = similar(D, N_right_marginal)

    map(1:l) do i
        right_marginal = getwindow(right_marginals, N_right_marginal, i)
        left_marginal .= @view moving_means[i+cols÷2,:]
        m_profile = [left_marginal; right_marginal]
        mpdist(m_profile, th, 2length(B))
    end
end


function mpdistmat(A::AbstractVector, B::AbstractVector, m::Int, d)
    N = lastlength(B)-m +1
    D = similar(A, lastlength(A)-m+1, N)
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
        Q .= min.(D[minI], Q)
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
