"""
capon(am::AbstractArrayManifold, Rxx, angles, f, c=c_0; kwargs...)

Calculates the capon spectrum for direction of arrival estimation.
This is identically to steering a capon beamformer and measuring the 
output power. Directly outputing the power spectrum for a given angle 
and data is just more convenient for DoA estimation.  

arguments:
----------
    am: Array (Manifold) to evaluate for
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
function capon(am::AbstractArrayManifold, Rxx, angles, f, c=c_0; kwargs...)
    A = am(angles, f, c; kwargs...)
    #P = 1/(a'*inv(Rxx)*a)
    P = vec(1 ./ sum(conj(A) .* (inv(Rxx) * A), dims=1))
    return real(P)
end