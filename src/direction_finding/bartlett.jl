"""
bartlett(pa::AbstractPhasedArray, Rxx, angles; w=nothing, c=c_0)

Calculates the bartlett spectrum for direction of arrival estimation.
This is identically to steering a bartlett (delay-and-sum) beamformer and measuring the 
output power. Directly outputing the power spectrum for a given angle 
and data is just more convenient for DoA estimation.

arguments:
----------
    pa: AbstractPhasedArray to calculate the estimator for
    Rxx: covariance matrix of the array which is used for estimation
    f: center/operating frequency
    angles: 1xD or 2xD matrix of steering directions.
        - If 1xD: azimuth angles in radians, elevation assumed zero.
        - If 2xD: [azimuth; elevation] for D directions in radians.
    w: vector of taper weights for the array (e.g., chebyshev window) 
    c: propagation speed of the wave

References:
-----------
H. Krim and M. Viberg, ‘Two decades of array signal processing research: the parametric approach’, IEEE Signal Process. Mag., vol. 13, no. 4, pp. 67–94, Jul. 1996.
"""
function bartlett(pa::AbstractPhasedArray, Rxx, f, angles; w=nothing, c=c_0, coords=:azel)
    
    if isnothing(w)
        W = I
    else 
        assert(length(w) == length(pa))
        W = diagm(vec(w))
    end

    A = steer(pa, f, angles; c=c, coords=coords)

    # a'*W*Rxx*W'*a
    P = vec(sum(conj(A) .* (W*Rxx*W' * A), dims=1))
    return real(P)
end