module MatrixProfile
using Statistics
using DSP
using LoopVectorization

export stomp


"""
    P,I = stomp(T, m)

Return the matrix profile and the profile indices of time series `T` with window length `m`.

Reference: [Matrix profile II](https://www.cs.ucr.edu/~eamonn/STOMP_GPU_final_submission_camera_ready.pdf).
"""
function stomp(T, m)
    n   = length(T)
    l   = n-m+1
    n > 2m+1 || throw(ArgumentError("Window length too long, maximum length is $(2m)"))
    μ,σ = running_mean_std(T, m)
    QT  = window_dot(view(T, 1:m), T)
    QT₀ = copy(QT)
    D   = distance_profile(QT, μ, σ, m)
    P   = copy(D)
    I   = ones(Int, l)
    @inbounds for i = 2:l
        for j = l:-1:2
            @fastmath QT[j] = QT[j-1]-T[j-1]*T[i-1]+T[j+m-1]*T[i+m-1]
        end
        QT[1] = QT₀[i]
        distance_profile!(D, QT, μ, σ, m, i)
        update_min!(P, I, D, i)
    end
    P, I
end

function distance_profile!(D, QT, μ, σ, m, i)
    @assert i <= length(D)
    @avx for j = eachindex(D)
        frac = (QT[j] - m*μ[i]*μ[j]) / (m*σ[i]*σ[j])
        D[j] = sqrt(max(2m*(1-frac), 0))
    end
    D[i] = typemax(eltype(D))
    D
end

distance_profile(QT, μ, σ, m) = distance_profile!(similar(QT), QT, μ, σ, m, 1)


function update_min!(P, I, D, i)
    @inbounds for j in eachindex(P,I,D)
        j == i && continue
        if D[j] < P[j]
            P[j] = D[j]
            I[j] = i
        end
    end
end

"""
Input: A query Q, and a user provided time series T
Output: The dot product between Q and all subsequences in T
"""
function window_dot(Q, T::AbstractArray{S}) where S
    n   = length(T)
    m   = length(Q)
    Qr  = reverse(Q)
    Qra = Qr
    QT  = conv(Qr, T)
    return QT[m:n]
end
# function window_dot(Q, T::AbstractArray{S}) where S
#     n    = length(T)
#     m    = length(Q)
#     Ta   = [T; zeros(S,n)]
#     Qr   = reverse(Q)
#     Qra  = [Qr; zeros(S,2n-m)]
#     Qraf = rfft(Qra)
#     Taf  = rfft(Ta)
#     QT   = irfft(Qraf .* Taf, length(Qra))
#     return QT[m:n]
# end


function running_mean_std(x::AbstractArray{T}, m) where T
    @assert length(x) > m
    n = length(x)-m+1
    s = ss = zero(T)
    μ = similar(x, n)
    σ = similar(x, n)
    @inbounds for i = 1:m
        s += x[i]
        ss += x[i]^2
    end
    μ[1] = s/m
    σ[1] = sqrt(ss/m - μ[1]^2)
    @inbounds for i = 1:n-1
        s -= x[i]
        ss -= x[i]^2
        s += x[i+m]
        ss += x[i+m]^2
        μ[i+1] = s/m
        σ[i+1] = sqrt(ss/m - μ[i+1]^2)
    end
    μ,σ
end


end
