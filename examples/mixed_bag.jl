using MatrixProfile
# Download data from https://sites.google.com/site/snippetfinderinfo/home/MixedBag.zip
datapath = "/home/fredrikb/Downloads/MixedBag/"
cd(datapath)

## Some data
T = parse.(Int, split(join(Char.(read("01911m_02019m_III_7680_200.txt"))), ','))
profile, snips, Cfracs = snippets(T, 3, 100, m=50)
plot(plot(snips, layout = (1, 3), size = (800, 200)), plot(profile, snips, legend = false))

## Power demand
T = parse.(Float32, readlines("Powerdemand_12_4500_200.txt"))
profile, snips, Cfracs = snippets(T, 4, 24, m=12)
plot(plot(snips, size = (800, 200), xrotation=45), plot(profile, snips, legend = false, link=:none))

profile = matrix_profile(T, 24)
mot = motifs(profile, 4, r=3, th=48)
plot(profile, mot, legend=false)
plot(mot)


## Walking styles
T = parse.(Float32, readlines("PAMAP_Subject4_NormalWalking_NordicWalking_17001_500.txt"))
T = T[1:2:end]
profile, snips, Cfracs = snippets(T, 4, 200)
plot(plot(snips, size = (800, 200), xrotation=45), plot(profile, snips, legend = false, link=:none))

profile = matrix_profile(T, 50)
mot = motifs(profile, 4, r=3)
plot(mot)
plot(profile, mot, legend=false)


## Robot dog
T = parse.(Float32, readlines("RoboticDogActivityY_64_4000_400.txt"))
profile, snips, Cfracs = snippets(T, 4, 100)
plot(plot(snips, size = (800, 200), xrotation=45), plot(profile, snips, legend = false, link=:none))

profile = matrix_profile(T, 50)
mot = motifs(profile, 4, r=2, th=200)
plot(profile, mot, legend=false)
plot(mot)

using DynamicAxisWarping
using MatrixProfile: znorm

dist = DTWDistance(DTW(5))
normalizer = Val(ZNormalizer)
profile2 = matrix_profile(T, 50, dist, normalizer=normalizer)

mot2 = motifs(profile2, 4, r=3, th=200, dist=(x,y)->dist(znorm(x),znorm(y)))
plot(profile2, mot2, legend=false)
plot(mot2)
