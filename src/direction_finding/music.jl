"""
music(pa::AbstractPhasedArray, Rxx, d, f, angles; c=c_0)

Calculates the MUSIC spectrum for direction of arrival (DoA) estimation.

arguments:
----------
    pa: Array for which the MUSIC spectrum is computed
    Rxx: Covariance matrix of the received signals
    d: Number of signal sources (model order)
    f: Center/operating frequency
    angles: 1xD vector or 2xD matrix of azimuth and elevation angles. For 1xD input, elevation is assumed zero.
    c: Propagation speed of the wave (default: c_0)

References:
-----------
H. Krim and M. Viberg, ‘Two decades of array signal processing research: the parametric approach’, IEEE Signal Process. Mag., vol. 13, no. 4, pp. 67–94, Jul. 1996.
"""
function music(pa::AbstractPhasedArray, Rxx, d, angles, f, c=c_0; W=I, kwargs...)
    eigs = eigen(Rxx)
    U = eigs.vectors[:, sortperm(eigs.values; rev=true)]
    Un = U[:, d+1:size(U)[2]]

    A = steer(pa, angles, f, c; kwargs...)

    # weighted music
    #P = a'*a/(a'*Un*W*Un'*a)
    B = Un'*A
    P = vec(sum(abs2, A; dims=1) ./ sum(conj(B) .* (W * B); dims=1))
    return real(P)
end