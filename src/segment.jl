"expected_arc(i,n) this polynomial indicates the expected number of NN-arcs that would pass over index `i` for a series of length`n` if the series was completely random. This is used to normalize the `segment_profile` to account for bias towards the edges."
expected_arc(i,n) = 2*(-i^2*n/(n-1)^2 + i*n/(n-1)) # polynomial coefficients solved for by using constraints e(0) = 0, e(n-1) = 0, e(n/2) = n/2

"""
    segment_profile(p::Profile)

Calculate the MP-index of a matrix profile. This index has the same length as the profile, and tells you how many nearest-neighbor arcs passes over index `i`. It's normalized by the expected number of arcs for a completely random time series, so that a value of 1 indicates a poor segmentation point, and a value close to 0 indicates a likely segmentation point.
"""
function segment_profile(p::Profile)
    I = p.I
    n = length(I)
    mpi = zeros(eltype(p.P), n)
    for i in 1:n
        Ii = I[i]
        if Ii > i
            for j = i:Ii-1
                mpi[j] += 1
            end
        else
            for j = Ii+1:i
                mpi[j] += 1
            end
        end
    end
    mpi .= min.(mpi ./ expected_arc.(0:n-1, n), 1)
    mpi
end

"""
    i = segment(p::Profile)

Returns an index `i` indicating the most likely segmentation point of profile `p`, i.e., the point which the fewst nearest-neighbor arcs passes over.

Ref: Matrix Profile VIII: Domain Agnostic Online Semantic
Segmentation at Superhuman Performance Levels
"""
function segment(p::Profile)
    argmin(segment_profile(p))
end

