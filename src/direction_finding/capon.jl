"""
capon(A::AbstractMatrix, Rxx)

Calculates the capon spectrum for direction of arrival estimation.
This is identically to steering a capon beamformer and measuring the 
output power. Directly outputing the power spectrum for a given angle 
and data is just more convenient for DoA estimation.  

arguments:
----------
    A: Array Manifold Matrix to calculate the estimator for
    Rxx: Covariance matrix (MxM)

returns:
--------
    P: Power spectrum evaluated at each direction

References:
-----------
H. Krim and M. Viberg, ‘Two decades of array signal processing research: the parametric approach’, IEEE Signal Process. Mag., vol. 13, no. 4, pp. 67–94, Jul. 1996.
"""
function capon(A::AbstractMatrix, Rxx)
    #P = 1/(a'*inv(Rxx)*a)
    P = vec(1 ./ sum(conj(A) .* (inv(Rxx) * A), dims=1))
    return real(P)
end