"""
capon(pa::AbstractPhasedArray, Rxx, f, ϕ, θ=0; fs=nothing, c=c_0)

Calculates the capon spectrum for direction of arrival estimation.
This is identically to steering a capon beamformer and measuring the 
output power. Directly outputing the power spectrum for a given angle 
and data is just more convenient for DoA estimation.  

arguments:
----------
    pa: Array to evaluate for
    Rxx: Covariance matrix (MxM)
    angles: 1xD or 2xD matrix of [azimuth; elevation] angles in radians
    c: Propagation speed (default: c_0)

returns:
--------
    P: Power spectrum evaluated at each direction

References:
-----------
H. Krim and M. Viberg, ‘Two decades of array signal processing research: the parametric approach’, IEEE Signal Process. Mag., vol. 13, no. 4, pp. 67–94, Jul. 1996.
"""
function capon(pa::AbstractPhasedArray, Rxx, f, angles; c=c_0, coords=:azel)
    A = steer(pa, f, angles; c=c, coords=coords)
    #P = 1/(a'*inv(Rxx)*a)
    P = vec(1 ./ sum(conj(A) .* (inv(Rxx) * A), dims=1))
    return real(P)
end