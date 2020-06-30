using MatrixProfile
# Download data from https://sites.google.com/site/snippetfinderinfo/home/MixedBag.zip
datapath = "/home/fredrikb/Downloads/MixedBag/"
cd(datapath)

## Some data
T = parse.(Float64, split(join(Char.(read("01911m_02019m_III_7680_200.txt"))), ','))
snips = snippets(T, 3, 200, m=50)
plot(snips)


## Power demand
T = parse.(Float64, readlines("Powerdemand_12_4500_200.txt"))
snips = snippets(T, 4, 24, m=24)
plot(snips)

profile = matrix_profile(T, 24)
mot = motifs(profile, 4, r=3, th=48)
plot(profile, mot, legend=false)
plot(mot)


## Walking styles
T = parse.(Float64, readlines("PAMAP_Subject4_NormalWalking_NordicWalking_17001_500.txt"))
T = T[1:2:end]
snips = snippets(T, 4, 200)
plot(snips)

profile = matrix_profile(T, 50)
mot = motifs(profile, 4, r=3)
plot(mot)
plot(profile, mot, legend=false)


## Robot dog
T = parse.(Float64, readlines("RoboticDogActivityY_64_4000_400.txt"))
snips = snippets(T, 4, 100)
plot(snips)

profile = matrix_profile(T, 50)
mot = motifs(profile, 4, r=2, th=200)
plot(profile, mot, legend=false)
plot(mot)

using DynamicAxisWarping

dist = DTW(3)
normalizer = ZNormalizer
profile2 = matrix_profile(T, 50, dist, normalizer=normalizer)


using Distances
struct ZWrapper{D} <: Distances.Metric
    dist::D
end
function Distances.evaluate(d::ZWrapper, x, y)
    evaluate(d.dist, znorm(x), znorm(y))
end
(d::ZWrapper)(x,y) = evaluate(d,x,y)

mot2 = motifs(profile2, 4, r=3, th=200, dist=ZWrapper(dist))
plot(profile2, mot2, legend=false)
plot(mot2)

snips = snippets(T, 4, 100, dist)
plot(snips)

@btime snippets($(T[1:2000]), 4, 100)
@profiler snippets((T[1:2000]), 4, 100)
