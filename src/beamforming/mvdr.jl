
"""
A: Array Manifold Matrix

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function mvdr_weights(A::AbstractMatrix, Rnn)
    #(inv(Rnn)*a)/(a'*inv(Rnn)*a) 
    #   -> adapted to steering matrix where each column is a steering vector
    Rnn_inv = inv(Rnn)
    numerator = Rnn_inv * A
    denominators = sum(conj(A) .* numerator, dims=1)
    return numerator ./ denominators
end

mpdr_weights(A::AbstractMatrix, Rxx) = mvdr_weights(A, Rxx)

const capon_weights = mpdr_weights