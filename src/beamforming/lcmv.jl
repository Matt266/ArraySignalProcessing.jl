
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