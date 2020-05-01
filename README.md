# MatrixProfile

[![Build Status](https://github.com/baggepinnen/MatrixProfile.jl/workflows/CI/badge.svg)](https://github.com/baggepinnen/MatrixProfile.jl/actions)
[![Coverage](https://codecov.io/gh/baggepinnen/MatrixProfile.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/baggepinnen/MatrixProfile.jl)

Time-series analysis using the matrix profile. The matrix profile `P` tells you which sub-sequences of a time series `T` are similar to each other, and which are most dissimilar from all other. This will allow you to find repeated patterns, or *motifs*, as well as finding outliers.

The function `stomp` returns the matrix profile and profile indices. Here's an example where we insert a repeated pattern in an otherwise random time series.
```julia
t   = range(0, stop=1, step=1/10)
y0  = sin.(2pi .* t)
T   = [randn(20); y0; randn(20); y0; randn(20)]
window_length = length(y0)
P,I = stomp(T, window_length)
plot(T, layout=(2,1), title="T", legend=false, link=:x)
plot!(P, sp=2, title="Matrix profile", legend=false) # Should have minima at 21 and 52
```
![matrix_profile](mp.svg)

The matrix profile have two sharp minima at the onsets of the repeated pattern. The parameter `window_length` determines how long pattern to search for.

## Runtime
`stomp` benefits greatly in speed from the use of `Flaot32` instead of `Float64`. The computational time scales as the square of the length of `T`, but is invariant to the window length. Calculating the matrix profile of `2^17 ≈ 100k` points takes about a minute on a laptop.

## Clustering
This is not implemented yet, but a common technique is to find the two time windows that are closest to each other, and group them together with all points within a small factor, say 2, of the distance between them. This group the makes up a motif. This can be done recursively after having removed the first motif in order to find more of them.

## References
The `stomp` algorithm is detailed in the paper [Matrix profile II](https://www.cs.ucr.edu/~eamonn/STOMP_GPU_final_submission_camera_ready.pdf).
