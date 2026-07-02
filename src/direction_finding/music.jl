"""
music(A::AbstractMatrix, Rxx, d; W=I)

Calculates the MUSIC spectrum for direction of arrival (DoA) estimation.

arguments:
----------
    A: Array Manifold Matrix to calculate the estimator for
    Rxx: Covariance matrix of the received signals
    d: Number of signal sources (model order)
    W: Weighting matrix for the MUSIC algorithm (default: I)
    c: Propagation speed of the wave (default: c_0)

References:
-----------
H. Krim and M. Viberg, ‘Two decades of array signal processing research: the parametric approach’, IEEE Signal Process. Mag., vol. 13, no. 4, pp. 67–94, Jul. 1996.
"""
function music(A::AbstractMatrix, Rxx, d; W=I)
    eigs = eigen(Rxx)
    U = eigs.vectors[:, sortperm(abs.(eigs.values); rev=true)]
    Un = U[:, d+1:size(U)[2]]

    # weighted music
    #P = a'*a/(a'*Un*W*Un'*a)
    B = Un'*A
    P = vec(sum(abs2, A; dims=1) ./ sum(conj(B) .* (W * B); dims=1))
    return real(P)
end