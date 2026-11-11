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

   profileez = @inferred matrix_profile(T, length(y0), exclusion_zone=10)
   Pez,Iez = profileez.P, profileez.I
   @test_nowarn plot(profileez)
   # plot(T, layout=2)
   # plot!(P, sp=2)

   m = findmin(Pez)
   @test m[1] < 1e-6
   @test m[2] == 51 || m[2] == 112
   # Distances near zero carry rounding errors of order sqrt(eps) due to cancellation in sqrt(2m(1-ρ))
   @test all(Pez .>= P .- 1e-6)


   # Test Euclidean between two series
   profile3 = @inferred matrix_profile(T, T, length(y0))
   @test all(profile3.P .< 1e-6)


   profile5 = matrix_profile(T[1:end÷2], T, length(y0))
   @test all(profile5.P[1:length(T)÷2-length(y0)] .< 1e-6)

   profile6 = matrix_profile(T, T[1:end÷2], length(y0))
   @test all(profile6.P .< 1e-6)


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

       # Regression: generic AB-join must accept a short A (n <= 2m+1).
       # A short query is a perfectly valid mutual matrix-profile use case.
       m = length(y0)
       short_A = T[51:51+m+1]                # length m+2, well below 2m+1
       @test_nowarn matrix_profile(short_A, T, m, normdist)
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


   @testset "mass" begin
       @info "Testing mass"
       T = randn(100)
       Q = randn(10)
       d = distance_profile(ZEuclidean(), Q, T)
       dmass = mass(Q, T)
       @test d ≈ dmass rtol=1e-3
   end

   @testset "Input validation and edge cases" begin
       @info "Testing input validation and edge cases"
       # A unique, smooth segment at the start of the series: without the exclusion
       # zone, its neighboring windows are trivial matches of the first window.
       T = randn(400); T[1:60] .= range(0, 1, length=60).^2
       ez = 10
       p = matrix_profile(T, 30; exclusion_zone=ez, showprogress=false)
       @test all(abs.(p.I .- (1:length(p.I))) .> ez)

       @test_throws ArgumentError matrix_profile(randn(20), 10, showprogress=false)
       @test !isdefined(MatrixProfile, :resample)

       Ti = rand(1:10, 100)
       @test matrix_profile(Ti, 10, showprogress=false).P ≈ matrix_profile(float.(Ti), 10, showprogress=false).P
       Ai = rand(1:10, 50)
       @test matrix_profile(Ai, Ti, 10, showprogress=false).P ≈ matrix_profile(float.(Ai), float.(Ti), 10, showprogress=false).P
       A32 = randn(Float32, 50)
       @test matrix_profile(A32, float.(Ti), 10, showprogress=false).P ≈ matrix_profile(Float64.(A32), float.(Ti), 10, showprogress=false).P
       @test_throws ArgumentError matrix_profile(randn(5), randn(100), 10, showprogress=false)

       @test mpdist([Inf, Inf], 0.9, 10) == Inf
       @test_throws ArgumentError snippets(randn(100), 3, 40)

       Td = sin.(0:0.1:30)
       @test length(damp(Td, 50).P) == length(Td) - 50 + 1 # default split point is smaller than m
       @test_throws ArgumentError damp(Td, 50, 10)
   end

   @testset "Comparison with stumpy" begin
       @info "Testing comparison with stumpy"
       # Reference values computed with stumpy 1.14.1:
       # stumpy.stump(T, m)[:, 0] and stumpy.stump(T, m, A, ignore_trivial=False)[:, 0],
       # with stumpy's default exclusion zone ceil(m/4).
       series(n, c) = [sin(0.37k) + 0.3cos(c*k + 0.01k^2) for k in 1:n]
       T = series(80, 2.1); T[20:35] .= 0.5; T[60] = NaN
       A = series(50, 1.3); A[10:24] .= -1.0; A[40] = NaN
       m = 8
       P_self = [1.8880019391642575, 1.676186000949532, 1.0505714316064334, 1.0543076185474054, 0.8317304101322976, 0.8601796945227194, 1.1961136501503076, 1.3860469598992893, 1.901916286543613, 1.6487601148702706, 1.3920978692459827, 0.9277662248964594, 0.9846045080499953, 1.188804612410343, 1.2064305258865833, 1.4454297746372962, 2.319330800086516, 2.8284271247461903, 1.6386344576710505, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.7176807181964557, 2.03574817513157, 1.3149340077313285, 1.754380757418532, 1.5512275742846233, 1.8128322003117296, 1.6656655638531843, 0.9710613806966734, 0.798248511465681, 0.7410147026138494, 0.5463374522282325, 0.5386779339412084, 0.9224161722355307, 1.544354886810089, 2.1144325258151193, 1.2372837365158695, 0.9920004622394213, 0.7002644482924371, 0.6250618364191982, 0.6227779082541649, 0.8120539224795632, 0.9956760356144023, 1.4598703443652248, 1.6656655638531843, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 1.2372837365158695, 0.9920004622394213, 0.7002644482924371, 0.6250618364191982, 0.6227779082541649, 0.8120539224795632, 1.0696181081304736, 1.9385208451044904, 1.742344721977506, 1.1032711188399653, 0.9034883311152503, 0.5463374522282325, 0.5386779339412084]
       P_AB = [1.2738419655665218, 2.0040084147477777, 1.2496392134759629, 0.8711595215974859, 0.78444520152703, 0.4923332905577403, 1.048280316840652, 1.197651296676646, 2.0627119203319504, 1.222030062745129, 1.01524702211665, 1.1122611210573603, 0.9915521078055421, 0.9989065782300426, 1.80413660737282, 1.9820418050223374, 1.9853151954484802, 2.0474207852910657, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3546485456500324, 1.4899995178679013, 1.6639552940161033, 1.771575443240779, 2.5241372884722266, 2.373186281184167, 1.7548051200822896, 1.4435865929651799, 1.0547382086360406, 0.9862486314494404, 0.9253107163288791, 1.1412980946948, 1.3752804704931398, 1.3134899700757183, 1.5334055392916475, 1.1446144148661193, 1.0132821624352624, 0.81810669447043, 0.9318155907020974, 0.9891230314664651, 1.2364709705003178, 2.0289516122786853, 1.8184819735037847, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 1.7893226205960924, 1.577766854043758, 1.2722891341616482, 1.2087964676976362, 1.2087487192108168, 1.2980713610980792, 1.6442934936975597, 2.162142749960906, 2.3410673206079027, 1.7055358583103377, 1.398604158914797, 1.1167980795813304, 1.0940402820904795]
       p = matrix_profile(T, m; exclusion_zone=cld(m, 4), showprogress=false)
       pab = matrix_profile(A, T, m; showprogress=false)
       for (P, R) in ((p.P, P_self), (pab.P, P_AB))
           @test all(isinf.(P) .== isinf.(R))
           @test P[isfinite.(R)] ≈ R[isfinite.(R)] atol=1e-6
       end
       # A non-finite value affects only the windows that contain it
       @test findall(isinf, p.P) == 53:60
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
       # The self-distance is zero up to rounding errors of order sqrt(eps)
       @test mpdist(A,A, m) < 1e-6
       @test mpdist(B,B, m) < 1e-6

       # Regression: short inputs send default k to 0; mpdist must not error.
       short = randn(10)
       @test_nowarn mpdist(short, short, 4)
       T = [A; B]
       S_in = 50
       p = mpdist_profile(T, S_in, 5)
       @test_nowarn plot(p)
       # Imani et al., Matrix Profile XIII, Table II: floor(n/S) windows.
       @test length(p) == length(T) ÷ S_in

       snips = snippets(T, 2, 50, m=5)
       @test_nowarn plot(snips)

       # Regression test: mpdist_profile / snippets must not mutate the caller's T.
       T_in = [A; B]
       T_ref = copy(T_in)
       mpdist_profile(T_in, 50, 5)
       @test T_in == T_ref
       snippets(T_in, 2, 50, m=5)
       @test T_in == T_ref
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
        p = matrix_profile(x, 20)
        s = segment_profile(p)
        val, ind = findmin(s)
        @test abs(ind-length(s)÷2) < 10
        @test val < 0.02

        @test segment(p) == ind

        # Regression: when no arcs cover the first/last index (mpi[1] or
        # mpi[end] == 0) the previous implementation produced NaN via
        # 0/0 division at the endpoints. Construct a Profile where each
        # window is its own NN to force this exact case.
        let n_idx = 8
            I_self = collect(1:n_idx)
            P_zero = zeros(Float64, n_idx)
            p_self = Profile(randn(n_idx + 4), P_zero, I_self, 5, nothing)
            s_self = segment_profile(p_self)
            @test !any(isnan, s_self)
        end
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



# T = randn(100000)
# Q = randn(10000)
# d = @btime distance_profile(ZEuclidean(), $Q, $T);
# dmass = @btime mass($Q, $T);
# @test d ≈ dmass rtol=1e-3

