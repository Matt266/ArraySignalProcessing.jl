"""
unconditional_signals(Rss, N; norm=true)

Generates signals corresponding to the unconditional signal model (stochastic signals) but with given spatial correlation

arguments:
----------
    Rss: Covariance matrix of the signals
    N: number of snapshots to generate
    norm: normalize so source power adds up to unit power

References:
------------
A. J. Barabell, J. Capon, D. F. DeLong, K. D. Senne, and J. R. Johnson, ‘Performance Comparison of Superresolution Array Processing Algorithms’, MIT Lincoln Laboratory, Lexington, MA, May 1984.
"""
function unconditional_signals(Rss, N; norm=true)
    # number signals
    d = size(Rss, 1)

    # generate random signals
    w = (randn(d, N) + 1im.*randn(d, N))/sqrt(2)

    # normalize so source power adds up to unit power
    if norm
        Rss = Rss/tr(Rss)
    end

    # correlate sources 
    s = cholesky(Rss).L * w
    return s
end

# assuming unit power signal, converts SNR in dB to noise variance σₙ²
function snr2nvar(SNR)
    return 10^-(SNR/10)
end

function unconditional_signals(pa::AbstractPhasedArray, Rss, N, SNR, angles, f, c=c_0; norm=true, kwargs...)
    M = length(pa)

    A = steer(pa, angles, f, c; kwargs...)
    s = unconditional_signals(Rss, N; norm=norm)
    s = convert(typeof(A), s)

    # if sum of signals is not normed to unit power 
    # amplify noise by signal power to reach desired SNR  
    nvar = snr2nvar(SNR) * (norm ? 1 : tr(Rss))
    n = sqrt(nvar/2)*(randn(M, N) + 1im*randn(M, N))
    n = convert(typeof(A), n)
    return A*s+n
end