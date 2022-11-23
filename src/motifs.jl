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

onsets(m::Motif) = getfield.(m.seqs, :onset)
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
    while P[perm[j]] < (d + 1e-5)*r && j < length(perm)
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


function anomalies(p::Profile, k)
    onsets = partialsortperm(p.P, 1:k, rev=true)
    Subsequence(p, copy(onsets), "Anomaly")
end
