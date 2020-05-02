using MatrixProfile
using Test, Statistics, LinearAlgebra, Plots

using MatrixProfile: znorm

normdist = (x,y)->norm(znorm(x)-znorm(y))

@testset "MatrixProfile.jl" begin


    @testset "running stats" begin
        @info "Testing running stats"

        x = randn(30)
        for w = 2:10
            m,s = MatrixProfile.running_mean_std(x, w)
            @test length(m) == length(s) == 30-w+1
            @test m[1] ≈ mean(x[1:w])
            @test s[1] ≈ std(x[1:w], corrected=false)

            @test m[2] ≈ mean(x[2:w+1])
            @test s[2] ≈ std(x[2:w+1], corrected=false)
        end

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


   # Test Euclidean between two series
   profile3 = matrix_profile(T, T, length(y0))
   @test all(profile3.P .< 1e-6)


   profile5 = matrix_profile(T[1:end÷2], T, length(y0))
   @test sum(profile5.P[1:length(T)÷2-length(y0)]) < 1e-5
   profile6 = matrix_profile(T, T[1:end÷2], length(y0))
   @test sum(profile6.P) < 1e-5


   @testset "Generic matrix profile" begin
       @info "Testing Generic matrix profile"

       # Test generic one serie
       profile2 = matrix_profile(T, length(y0), normdist)
       @test profile2.P ≈  profile.P
       @test all(eachindex(profile.I)) do i
           profile2.I[i] ∈ (51,112) || profile2.I[i] == profile.I[i]
       end



       # Test generic between two series
       profile4 = matrix_profile(T, T, length(y0), normdist)
       @test all(profile4.P .< 1e-6)


   end



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

       an = MatrixProfile.anomalies(profile, 3)
       @test length(an.motifs) == length(an.onsets) == 3
       @test an.onsets[1] == argmax(profile.P)
   end






 end






# Benchmark
# @btime matrix_profile($(randn(Float32, 10_000)), 20)
#
# @time matrix_profile((randn(Float32, 100_000)),  20)
#
# @time matrix_profile((randn(Float32, 10_000)), (randn(Float32, 1_000)), 20, (x,y)->norm(x-y))
#
# @profiler matrix_profile(randn(10_000), 20)
# Q = randn(5)
# T = randn(10)
# d1 = DynamicAxisWarping.window_dot(Q,T)
# d2 = DynamicAxisWarping.window_dot2(Q,T)
# d3 = DynamicAxisWarping.window_dot3(Q,T)
# [d1 d2 d3]


# Eval precision
# a = randn(Float64, 200_000)
# p1 = matrix_profile(a,  20)
# p2 = matrix_profile(Float32.(a),  20)
# plot(abs.(p1.P - p2.P)./p1.P)
