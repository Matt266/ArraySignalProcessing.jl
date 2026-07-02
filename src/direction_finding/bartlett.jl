"""
bartlett(A::AbstractMatrix, Rxx; W=I)

Calculates the bartlett spectrum for direction of arrival estimation.
This is identically to steering a bartlett (delay-and-sum) beamformer and measuring the 
output power. Directly outputing the power spectrum for a given angle 
and data is just more convenient for DoA estimation.

arguments:
----------
    A: Array Manifold Matrix to calculate the estimator for
    Rxx: covariance matrix of the array which is used for estimation
    W: diagonal matrix of taper weights for the array (e.g., chebyshev window) 

References:
-----------
H. Krim and M. Viberg, ‘Two decades of array signal processing research: the parametric approach’, IEEE Signal Process. Mag., vol. 13, no. 4, pp. 67–94, Jul. 1996.
"""
function bartlett(A::AbstractMatrix, Rxx; W=I)
    # a'*W*Rxx*W'*a
    P = vec(sum(conj(A) .* (W*Rxx*W' * A), dims=1))
    return real(P)
end