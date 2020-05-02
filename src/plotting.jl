_append_inf(x) = vec([x; fill(Inf, 1, size(x,2))])

@recipe function plot(p::Profile)
    link --> :x
    layout --> (2,1)
    legend --> false
    @series begin
        title --> "T"
        label --> "T"
        subplot --> 1
        p.T
    end
    @series begin
        title --> "Matrix profile"
        subplot --> 2
        p.P
    end
end


@recipe function plot(p::Profile, motifs::Vector{Motif})

    @series begin
        p
    end

    inds = 0:p.m-1
    linewidth --> 2
    for (j,m) in enumerate(motifs)
        @series begin
            legend --> true
            group := j
            label --> j
            subplot := 1
            _append_inf(m.onsets' .+ inds), _append_inf(reduce(hcat, m.motifs))
        end
    end
end
