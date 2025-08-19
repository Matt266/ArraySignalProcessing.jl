
"""

Arguments:
-----------
C: Constraint matrix (e.g., each column a steering vector to set a gain in that direction)
G: Gain matrix (each column is a gain vector with desired responses for each set of angels)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function lcmv_weights(Rnn, C, G)
    # w' = g' * [C' * Rnn^-1 * C]^-1 * C'* Rnn^-1
    # (Rnn and inv(Rnn) are hermitian)
    return inv(Rnn)*C*inv(C'*inv(Rnn)*C)*G
end

lcmp_weights(Rxx, C, G) = lcmv_weights(Rxx, C, G)

"""
Calculate quiescent beamformer weights WQ, Blocking Matrix B and an initial set of adaptive weights (WA)
for usage with a Generalized Sidelobe Canceller (GSC)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function lcmv_gsc(Rxx, C, G)
    N, M = size(C)
    WQ = pinv(C')*G

    PC = C*pinv(C)
    U, _, _ = svd(I-PC)
    B = U[:, 1:(N - M)]

    WA = (B'*Rxx*B) \ (B'*Rxx*WQ)

    # Total weights:
    # W = WQ - B * WA

    return WQ, B, WA
end

"""
Basic utility power spectrum for LCMV. 

E.g., "Sweep" C as wideband distortionless constraint over multiple angles for 
wideband DoA estimation. 

As this does not directly perform the sweep over multiple angles and one may want to 
set additional constraints (e.g., derivative), this is not directly a detection method but 
a utility that can be used to perform a LMCV-based detection.

References:
-----------

W. Liu and S. Weiss, Wideband Beamforming. Nashville, TN: John Wiley & Sons, 2010.
"""
function lcmv(Rxx, C, G)
    # g'*inv(C'*inv(Rxx)*C)*g for multiple g as Matrix G
    P = vec(sum(conj(G) .* (inv(C'*inv(Rxx)*C) * G), dims=1))
    return real(P)
end