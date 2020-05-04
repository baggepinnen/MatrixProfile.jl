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
subseqtype(::Motif) = "Motif"

onsets(m::Subsequence) = m.onset
seqs(m::Subsequence) = [m.seq]
seqlength(m::Subsequence) = length(m.seq)
subseqtype(m::Subsequence) = m.type

onsets(m::Array{<:Subsequence}) = onsets.(m)
seqs(m::Array{<:Subsequence}) = getfield.(m, :seq)
seqlength(m::Array{<:Subsequence}) = seqlength(length(m[1]))
subseqtype(m::Array{<:Subsequence}) = subseqtype(m[1])

function motifs(p::Profile, k, r, th = 0, found_motifs = Motif[])
    length(found_motifs) == k && return found_motifs
    m = p.m
    P = copy(p.P)
    remove_found_motifs!(P, found_motifs, th)
    d, i = findmin(P)
    distance_profile!(P, getwindow(p.T, m, i), p.T)
    remove_found_motifs!(P, found_motifs, th)
    onsets = findall(<=((d + 1e-5)*r), P)
    i ∈ onsets || push!(onsets, i)
    sort!(onsets)
    push!(found_motifs, Motif(Subsequence(p, onsets, "Motif")))
    motifs(p::Profile, k, r, th, found_motifs)
end

function Subsequence(p::Profile, onsets, type::String = "Subsequence")
    m    = p.m
    inds = 0:m-1
    n    = length(onsets)
    vecs = map(1:n) do i
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
