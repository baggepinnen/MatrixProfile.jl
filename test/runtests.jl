using MatrixProfile
using Test, Statistics, LinearAlgebra, Plots

using MatrixProfile: znorm

@testset "MatrixProfile.jl" begin

   x = randn(30)
   for w = 2:10
       m,s = MatrixProfile.running_mean_std(x, w)
       @test length(m) == length(s) == 30-w+1
       @test m[1] ≈ mean(x[1:w])
       @test s[1] ≈ std(x[1:w], corrected=false)

       @test m[2] ≈ mean(x[2:w+1])
       @test s[2] ≈ std(x[2:w+1], corrected=false)
   end


   t = range(0, stop=1, step=1/10)
   y0 = sin.(2pi .* t)

   T = [randn(50); y0; randn(50); y0; randn(50)]

   profile = matrix_profile(T, length(y0))
   P,I = profile.P, profile.I
   @test_nowarn plot(profile)
   # plot(T, layout=2)
   # plot!(P, sp=2)

   m = findmin(P)
   @test m[1] < 1e-6
   @test m[2] == 51 || m[2] == 112

   Q = randn(5)
   T = [randn(5); Q; randn(5)]
   D = MatrixProfile.distance_profile(Q,T)
   @test D[6] < 1e-6
   @test D[1] ≈ norm(znorm(Q) - znorm(T[1:5]))

   @time matrix_profile(randn(Float32, 2^15), 256)



   @testset "Motifs" begin
       @info "Testing Motifs"
       @time mot = motifs(profile, 2, 2, 5)
       @test mot[1].onsets == [51, 112]
       @test_nowarn plot(profile, mot)
   end

 end






# Benchmark
# @btime matrix_profile($randn(Float32, 10_000), 20)
# @profiler matrix_profile(randn(10_000), 20)
# Q = randn(5)
# T = randn(10)
# d1 = DynamicAxisWarping.window_dot(Q,T)
# d2 = DynamicAxisWarping.window_dot2(Q,T)
# d3 = DynamicAxisWarping.window_dot3(Q,T)
# [d1 d2 d3]
