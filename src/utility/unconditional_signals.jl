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
    w = (randn(d, N) + 1im*randn(d, N))/sqrt(2)

    # normalize so source power adds up to unit power
    if norm
        Rss = Rss/tr(Rss)
    end

    # correlate sources 
    s = cholesky(Rss).L * w
    return s
end