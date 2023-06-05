using MatrixProfile
using Test, Statistics, LinearAlgebra, Plots

using SlidingDistancesBase
using MatrixProfile: Profile

normdist = (x,y)->norm(znorm(x)-znorm(y))

function naive_matrix_profile(A,B,m)
    res = map(1:length(B)-m+1) do i
        findmin(distance_profile(ZEuclidean(), getwindow(B,m,i), A))
    end
    Profile(B,getindex.(res, 1), getindex.(res, 2), m, A)
end

function profeq(p1,p2)
    cond = p1.P ≈ p2.P
    cond && return true
    @show norm(p1.P - p2.P)
    false
end


@testset "MatrixProfile.jl" begin

   t = range(0, stop=1, step=1/10)
   y0 = sin.(2pi .* t)
   T = [randn(50); y0; randn(50); y0; randn(50)]
   A = sign.([randn(50); y0; randn(50); y0; randn(50)])[1:end÷2] .+ 0.01.*randn.()

   profile = @inferred matrix_profile(T, length(y0))
   P,I = profile.P, profile.I
   @test_nowarn plot(profile)
   # plot(T, layout=2)
   # plot!(P, sp=2)

   m = findmin(P)
   @test m[1] < 1e-6
   @test m[2] == 51 || m[2] == 112


   # Test Euclidean between two series
   profile3 = @inferred matrix_profile(T, T, length(y0))
   @test all(profile3.P .< 1e-6)


   profile5 = matrix_profile(T[1:end÷2], T, length(y0))
   @test sum(profile5.P[1:length(T)÷2-length(y0)]) < 1e-5

   profile6 = matrix_profile(T, T[1:end÷2], length(y0))
   @test sum(profile6.P) < 1e-5


   profile7n = naive_matrix_profile(T, A, length(y0))
   profile7n2 = @inferred matrix_profile(T, A, length(y0), normdist)
   profile7 = matrix_profile(T, A, length(y0))
   @test profeq(profile7, profile7n)
   @test profeq(profile7n, profile7n2)

   profile8n = naive_matrix_profile(A, T, length(y0))
   profile8n2 = @inferred matrix_profile(A, T, length(y0), normdist)
   profile8 = @inferred matrix_profile(A, T, length(y0))
   @test profeq(profile8, profile8n)
   @test profeq(profile8n, profile8n2)


   @testset "Generic matrix profile" begin
       @info "Testing Generic matrix profile"

       # Test generic one serie
       profile2 = @inferred matrix_profile(T, length(y0), normdist)
       @test profile2.P ≈  profile.P
       @test all(eachindex(profile.I)) do i
           profile2.I[i] ∈ (51,112) || profile2.I[i] == profile.I[i]
       end

       # Test generic between two series
       profile4 = @inferred matrix_profile(T, T, length(y0), normdist)
       @test all(profile4.P .< 1e-6)
   end





   @testset "Motifs" begin
       @info "Testing Motifs"
       t = range(0, stop=1, step=1/10)
       y0 = sin.(2pi .* t)
       T = [randn(50); y0; randn(50); y0; randn(50)]
       A = sign.([randn(50); y0; randn(50); y0; randn(50)])[1:end÷2] .+ 0.01.*randn.()
       profile = matrix_profile(T, length(y0))

       @time mot = motifs(profile, 2, r=2, th=5)
       @test MatrixProfile.subseqtype(mot) == "Motif"
       @test onsets(mot[1]) == [51, 112]
       @test_nowarn plot(mot)
       @test_nowarn plot(profile, mot)

       an = anomalies(profile, 3)
       @test length(an) == 3
       @test length(seqs(an)) == 3
       @test MatrixProfile.subseqtype(an) == "Anomaly"
       @test onsets(an)[1] == argmax(profile.P)
       @test_nowarn plot(an)
       @test MatrixProfile.seqlength(an) == length(y0)
   end

   @testset "damp" begin
        @info "Testing damp"
        t = 0:0.1:100
        T = sin.(t)
        T[500:550] .*= 0
        aMP = damp(T, 100, 200)
        @test argmax(aMP.P) ≈ 500 atol=1
   end



   @time matrix_profile(randn(Float32, 2^15), 256)

   @testset "mpdist" begin
       @info "Testing mpdist"

       t = 1:0.01:3
       A = @. sin(2pi * t) + 0.1 * randn()
       B = @. sign(sin(2pi * t)) + 0.1 * randn()
       m = 4
       # p1 = matrix_profile(A,B,m)
       # p2 = matrix_profile(B,A,m)
       # partialsort([p1.P; p2.P], 20)
       # plot(p1)
       # plot!(p2)
       @test @inferred(mpdist(A,B, m)) > 0
       @test mpdist(A,A, m) < 10sqrt(eps())
       @test mpdist(B,B, m) < 10sqrt(eps())
       T = [A; B]
       p = mpdist_profile(T, 50, 5)
       @test_nowarn plot(p)

       snips = snippets(T, 2, 50, m=5)
       @test_nowarn plot(snips)

   end


   @testset "segment_profile and segmentation" begin
       @info "Testing segment_profile and segmentation"

        @test MatrixProfile.expected_arc(0,10) == 0
        @test MatrixProfile.expected_arc(9,10) == 0
        @test MatrixProfile.expected_arc(9/2,10) ≈ 5

        x = randn(1000)
        p = matrix_profile(x, 20)
        s = segment_profile(p)
        @test minimum(s) > 0.5

        x = [sin.(1:0.1:100); sign.(sin.(1:0.1:100))]
        x .+= 0.01 .* randn.()
        p = p = matrix_profile(x, 20)
        s = segment_profile(p)
        val, ind = findmin(s)
        @test abs(ind-length(s)÷2) < 10
        @test val < 0.02

        @test segment(p) == ind
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
# a = randn(Float64, 50_000)
# p1 = matrix_profile(a,  20)
# p2 = matrix_profile(Float32.(a),  20)
# @show mean(abs.(p1.P - p2.P)./p1.P)
# plot(abs.(p1.P - p2.P)./p1.P)

# path = "/tmp/MixedBag/01911m_02019m_III_7680_200.txt"
# T = parse.(Int, split(join(Char.(read(path))), ','))
# profile, snips, Cfracs = snippets(T, 3, 20, 10)
# plot(plot(snips, layout = (1, 3), size = (800, 200)), plot(profile, snips, legend = false), layout=(2,1))
