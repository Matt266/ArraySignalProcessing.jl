"""
mdl(Rxx, K)

Estimates number of sources using the Minimum Description Length (MDL) model selection criteria.

arguments:
----------
    Rxx: covariance matrix of the array which is used for estimation
    K: sample size

References:
-----------
M. Wax and T. Kailath, ‘Detection of signals by information theoretic criteria’, IEEE Trans. Acoust., vol. 33, no. 2, pp. 387–392, Apr. 1985.

H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function mdl(Rxx, K)
    N = size(Rxx, 1)
    λ = real(eigvals(Rxx, sortby= λ -> -abs(λ)))
    λ = max.(λ, 0)
    mdl = zeros(N)

    for d in 0:N-1
        λ_n = λ[d+1:end]
        L = ((N-d)K)*log(mean(λ_n)/StatsBase.geomean(λ_n))
        mdl[d+1] = L + 0.5*(d*(2N-d)+1)log(K)
    end

    # filter out invalid (Inf, -Inf, NaN) results
    orders = (0:N-1)[map(!, (isnan.(mdl) .|| isinf.(mdl)))]
    mdl = mdl[map(!, (isnan.(mdl) .|| isinf.(mdl)))]
    return orders[argmin(mdl)]
end