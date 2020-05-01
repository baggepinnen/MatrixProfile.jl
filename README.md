# MatrixProfile

[![Build Status](https://github.com/baggepinnen/MatrixProfile.jl/workflows/CI/badge.svg)](https://github.com/baggepinnen/MatrixProfile.jl/actions)
[![Coverage](https://codecov.io/gh/baggepinnen/MatrixProfile.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/baggepinnen/MatrixProfile.jl)


The function `stomp` returns the matrix profile and profile indices.
```julia
t   = range(0, stop=1, step=1/10)
y0  = sin.(2pi .* t)
T   = [randn(50); y0; randn(50); y0; randn(50)]
window_length = length(y0)
P,I = stomp(T, window_length)
plot(T, layout=2)
plot!(P, sp=2) # Should have minima at 51 and 112
```

`stomp` benefits greatly in speed from the use of `Flaot32` instead of `Float64`.

Reference: [Matrix profile II](https://www.cs.ucr.edu/~eamonn/STOMP_GPU_final_submission_camera_ready.pdf).
