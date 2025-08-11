"""
aic(Rxx, K)

Estimates number of sources using the Akaike information criterion (AIC).

arguments:
----------
    Rxx: covariance matrix of the array which is used for estimation
    K: sample size

References:
-----------
M. Wax and T. Kailath, ‘Detection of signals by information theoretic criteria’, IEEE Trans. Acoust., vol. 33, no. 2, pp. 387–392, Apr. 1985.

H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function aic(Rxx, K)
    N = size(Rxx, 1)
    λ = sort(real(eigen(Rxx).values), rev=true)
    λ = max.(λ, 0)

    aic = zeros(N)

    for d in 0:N-1
        λ_n = λ[d+1:end]
        L = ((N-d)K)*log(mean(λ_n)/StatsBase.geomean(λ_n))
        aic[d+1] = L + d*(2N-d)
    end

    # filter out invalid (Inf, -Inf, NaN) results
    orders = (0:N-1)[map(!, (isnan.(aic) .|| isinf.(aic)))]
    aic = aic[map(!, (isnan.(aic) .|| isinf.(aic)))]
    return orders[argmin(aic)]
end