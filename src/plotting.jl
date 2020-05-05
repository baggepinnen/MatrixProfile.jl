_append_inf(x) = vec([x; fill(Inf, 1, size(x,2))])

@recipe function plot(p::Profile)
    link --> :x
    layout --> (2,1)
    legend --> false
    @series begin
        title --> "T"
        label --> "T"
        subplot --> 1
        p.T isa AbstractVector{<:Number} ? p.T : reduce(hcat, p.T)'
    end
    @series begin
        title --> "Matrix profile"
        subplot --> 2
        p.P
    end
end




@recipe function plot(p::Profile, motifs::SubSeqType)

    motifs isa Vector || (motifs = [motifs])

    @series begin
        p
    end

    inds = 0:p.m-1
    linewidth --> 2
    for (j,m) in enumerate(motifs)
        @series begin
            legend --> true
            group := j
            label --> string(subseqtype(m), " ", j)
            subplot := 1
            _append_inf(onsets(m)' .+ inds), _append_inf(reduce(hcat, seqs(m)))
        end
    end
end


@recipe function plot(sequences::SubSeqType)

    sequences isa Vector || (sequences = [sequences])
    layout --> length(sequences)
    for (j,s) in enumerate(sequences)
        @series begin
            legend --> false
            group := j
            label --> string(subseqtype(s), " ", j)
            title --> string(subseqtype(s), " ", j)
            subplot := j
            inds = 0:seqlength(s)-1
            _append_inf(onsets(s)' .+ inds), _append_inf(reduce(hcat, seqs(s)))
        end
    end
end

@recipe function plot(motifs::Vector{Motif})

    motifs isa Vector || (motifs = [motifs])
    layout --> length(motifs)
    for (j,m) in enumerate(motifs)
        @series begin
            legend --> false
            group := j
            label --> string(subseqtype(m), " ", j)
            title --> string(subseqtype(m), " ", j)
            subplot := j
            inds = 0:seqlength(m)-1
            reduce(hcat, seqs(m))
        end
    end
end
