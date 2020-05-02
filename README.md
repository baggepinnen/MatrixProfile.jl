# MatrixProfile

[![Build Status](https://github.com/baggepinnen/MatrixProfile.jl/workflows/CI/badge.svg)](https://github.com/baggepinnen/MatrixProfile.jl/actions)
[![Coverage](https://codecov.io/gh/baggepinnen/MatrixProfile.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/baggepinnen/MatrixProfile.jl)

Time-series analysis using the matrix profile. The matrix profile `P` tells you which sub-sequences of a time series `T` are similar to each other, and which are most dissimilar from all other. This will allow you to find repeated patterns, or *motifs*, as well as finding outliers.

The function `matrix_profile` returns the matrix profile and profile indices. Here's an example where we insert a repeated pattern in an otherwise random time series.
```julia
using MatrixProfile, Plots
t  = range(0, stop=1, step=1/10)
y0 = sin.(2pi .* t)
T  = [randn(20); y0; randn(20); y0; randn(20)]
window_length = length(y0)
profile = matrix_profile(T, window_length)
plot(profile) # Should have minima at 21 and 52
```
![matrix_profile](figures/mp.svg)

The matrix profile have two sharp minima at the onsets of the repeated pattern. The parameter `window_length` determines how long pattern to search for.

### Runtime
`matrix_profile` benefits greatly in speed from the use of `Float32` instead of `Float64`. The computational time scales as the square of the length of `T`, but is invariant to the window length. Calculating the matrix profile of `2^17 ≈ 100k` points takes about a minute on a laptop.

## Motif grouping
Using the fake data from the example above, we can do
```julia
k = 2; r = 2; th = 5;
mot = motifs(profile, k, r, th)
plot(profile, mot)
```
- `k` is the number of motifs to extract
- `r` controls how similar two windows must be to belong to the same motif. A higher value leads to more windows being grouped together.
- `th` is a threshold on how far nearby in time two motifs are allowed to be.
![motif_plot](figures/motifs.svg)


## Arbitrary metrics and input types
The matrix profile can be computed for any object that has a "time axis" as a notion of distance. The examples so far have dealt with one-dimensional arrays of real numbers with the Euclidean metric, for which the matrix profile is particularly efficient to compute. We do not have to limit ourselves to this setting, though, and `matrix_profile` accepts any array-like object and any distance function on the form `dist(x,y)`. The interface looks like this
```julia
profile = matrix_profile(T, m, dist)
```
If `T` is a high-dimensional array, time is considered to be the last axis. `T` can also be a vector of any arbitrary julia objects for which the function `dist(x,y)` is defined. Note that if `T` has a long time dimensions, the matrix profile will be expensive to compute, 𝒪(n²log(n)). This method does not make use of the STOMP algorithm, since this is limited to one-dimensional data under the Euclidean metric.



## References
The STOMP algorithm used in `matrix_profile` is detailed in the paper [Matrix profile II](https://www.cs.ucr.edu/~eamonn/STOMP_GPU_final_submission_camera_ready.pdf).
