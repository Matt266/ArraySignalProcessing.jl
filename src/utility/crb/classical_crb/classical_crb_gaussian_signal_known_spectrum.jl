"""
classical_crb_gaussian_signal_known_spectrum(Rss, nvar, am, θ, K=1)
Calculates the unconditional (complex gaussian zero-mean random signal model) classical Cramér-Rao Bound (CRB) for a complex Gaussian random observation vector x
and KNOWN source spectrum Rss (and known noise variance).

arguments:
---------
    Rss: DxD covariance matrix of the source signals
    nvar: noise variance
    am: Array Manifold Matrix am(θ) as function of the parameter vector θ
    θ: parameter vector of the Array Manifold Matrix to evaluate the CRB
    K: number of sample snapshots (default is 1)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002. (8.85), (8.200), (8.201), (8.202), (8.204), (8.205)
Check text paragraph between (8.89) and (8.90): Sx for Kx, also signal s is a sample function
from a zero-mean stationary complex Gaussian random process so m = 0
"""
function classical_crb_gaussian_signal_known_spectrum(Rss, nvar, am, θ::AbstractVector, K::Number=1)
    N = size(am(θ), 1)

    function Kx(θ_full)
        A = am(θ_full)
        return A*Rss*A'+ nvar * collect(I(N))
    end
    
    # signal and noise are zero-mean random processes
    m = fill!(similar(Rss, N, K), 0)

    return classical_crb(θ, Kx, m, K)
end

"""
classical_crb_gaussian_signal_known_spectrum(Rss, nvar, am, θw, θu, K=1) 
Calculates the unconditional (complex gaussian zero-mean random signal model) classical Cramér-Rao Bound (CRB) for a complex Gaussian random observation vector x
and KNOWN source spectrum Rss (and known noise variance).

arguments:
---------
    Rss: DxD covariance matrix of the source signals
    nvar: noise variance
    am: Array Manifold Matrix am(θw, θu) as function of the parameter vectors θw and θu
    θw: vector of wanted parameters 
    θu: vector of unwanted parameters
    K: number of sample snapshots (default is 1)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
Check text paragraph between (8.89) and (8.90): Sx for Kx, also signal s is a sample function
from a zero-mean stationary complex Gaussian random process so m = 0
"""
function classical_crb_gaussian_signal_known_spectrum(Rss, nvar, am, θw::AbstractVector, θu::AbstractVector, K::Number=1) 
    N = size(am(θw, θu), 1)
    
    function Kx(θ_full)
        θw_val = θ_full[1:length(θw)]
        θu_val = θ_full[(length(θw)+1):end]

        A = am(θw_val, θu_val)
        return A*Rss*A'+ nvar * collect(I(N))
    end
    
    # signal and noise are zero-mean random processes
    m = fill!(similar(Rss, N, K), 0)

    return classical_crb(θw, θu, Kx, m; K=K)
end