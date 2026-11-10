abstract type AbstractSubSequence end

struct Subsequence{VT <: AbstractVector} <: AbstractSubSequence
    seq::VT
    onset::Int
    type::String
end

struct Motif <: AbstractSubSequence
    seqs::Vector{Subsequence}
end

const SubSeqType = Union{AbstractSubSequence, AbstractVector{<:AbstractSubSequence}}

"""
    onsets(m)

Return the onsets (start indices in the time series) of the occurrences of motif `m`, or of the subsequence(s) `m` returned by [`anomalies`](@ref) or [`snippets`](@ref).
"""
onsets(m::Motif) = getfield.(m.seqs, :onset)
"""
    seqs(m)

Return the subsequences belonging to motif `m`, or to the subsequence(s) `m` returned by [`anomalies`](@ref) or [`snippets`](@ref), as a vector of vectors.
"""
seqs(m::Motif) = getfield.(m.seqs, :seq)
seqlength(m::Motif) = length(m.seqs[1].seq)
subseqtype(::Union{Motif, Vector{<:Motif}}) = "Motif"

onsets(m::Subsequence) = m.onset
seqs(m::Subsequence) = [m.seq]
seqlength(m::Subsequence) = length(m.seq)
subseqtype(m::Subsequence) = m.type

onsets(m::Array{<:Subsequence}) = onsets.(m)
seqs(m::Array{<:Subsequence}) = getfield.(m, :seq)
seqlength(m::Array{<:Subsequence}) = length(m[1].seq)
subseqtype(m::Array{<:Subsequence}) = subseqtype(m[1])

"""
    motifs(p::Profile, k, found_motifs = Motif[]; r=2, th = p.m, dist = ZEuclidean())

- `k` is the number of motifs to extract
- `r` controls how similar two windows must be to belong to the same motif. A higher value leads to more windows being grouped together.
- `th` is a threshold on how nearby in time two motifs are allowed to be.
- `dist` is the distance used to find the occurrences of a motif.

Returns a vector of `Motif`s. Use [`onsets`](@ref) and [`seqs`](@ref) to obtain the onsets and subsequences of each motif, and `plot(p, motifs)` to visualize them.

Also see the function `anomalies(profile)` to find anomalies in the data, sometimes called *discords*.
"""
function motifs(p::Profile, k, found_motifs = Motif[]; r=2, th = p.m, dist = ZEuclidean())
    length(found_motifs) == k && return found_motifs
    m = p.m
    P = copy(p.P)
    remove_found_motifs!(P, found_motifs, th)
    d, i = findmin(P)
    distance_profile!(P, dist, getwindow(p.T, m, i), p.T)
    remove_found_motifs!(P, found_motifs, th)
    perm = sortperm(P)
    onsets = [perm[1]]
    j = 2
    while j <= length(perm) && P[perm[j]] < (d + 1e-5)*r
        all(abs.(perm[j] .- onsets) .> th) && push!(onsets, perm[j])
        j += 1
    end
    sort!(onsets)
    push!(found_motifs, Motif(Subsequence(p, onsets, "Motif")))
    motifs(p::Profile, k, found_motifs; r=r, th=th, dist=dist)
end


function Subsequence(p::Profile, onsets, type::String = "Subsequence")
    m    = p.m
    inds = 0:m-1
    vecs = map(eachindex(onsets)) do i
        p.T[inds .+ onsets[i]] # Don't use getwindow
    end
    Subsequence.(vecs, onsets, type)
end

function remove_found_motifs!(P, found_motifs, th)
    for mo in found_motifs
        for o in onsets(mo)
            inds = max(1, o-th):min(length(P), o+th)
            P[inds] .= typemax(eltype(P))
        end
    end
end


"""
    anomalies(p::Profile, k)

Return the `k` subsequences with the largest matrix-profile values, i.e., the subsequences that are most dissimilar to their nearest neighbors, sometimes called *discords*. The result is a vector of `Subsequence`s, use [`onsets`](@ref) and [`seqs`](@ref) to access their onsets and values and `plot(p, anomalies(p, k))` to visualize them.

Overlapping windows are not excluded; consecutive windows around a single anomaly may all be returned.
"""
function anomalies(p::Profile, k)
    onsets = partialsortperm(p.P, 1:k, rev=true)
    Subsequence(p, copy(onsets), "Anomaly")
end
