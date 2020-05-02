struct Motif
    motifs
    onsets
end

function motifs(p::Profile, k, r, th = 0, found_motifs = Motif[])
    length(found_motifs) == k && return found_motifs
    m = p.m
    P = copy(p.P)
    remove_found_motifs!(P, found_motifs, th)
    d, i = findmin(P)
    distance_profile!(P, getwindow(p.T, m, i), p.T)
    remove_found_motifs!(P, found_motifs, th)
    onsets = findall(<=((d + sqrt(eps()))*r), P)
    onsets = sort!(unique(push!(onsets, i)))
    push!(found_motifs, Motif(p, onsets))
    motifs(p::Profile, k, r, th, found_motifs)
end

function Motif(p::Profile, onsets)
    m    = p.m
    inds = 0:m-1
    n    = length(onsets)
    vecs = map(1:n) do i
        p.T[inds .+ onsets[i]] # Don't use getwindow
    end
    Motif(vecs, onsets)
end

function remove_found_motifs!(P, found_motifs, th)
    for mo in found_motifs
        for o in mo.onsets
            inds = max(1, o-th):min(length(P), o+th)
            P[inds] .= typemax(eltype(P))
        end
    end
end


function anomalies(p::Profile, k)
    onsets = partialsortperm(p.P, 1:k, rev=true)
    Motif(p, copy(onsets))
end
