
"""
References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function mvdr_weights(am::AbstractArrayManifold, Rnn, angles, f, c=c_0; kwargs...)
    v = am(angles, f, c; kwargs...)
    #(inv(Rnn)*v)/(v'*inv(Rnn)*v) 
    #   -> adapted to steering matrix where each column is a steering vector
    Rnn_inv = inv(Rnn)
    numerator = Rnn_inv * v
    denominators = sum(conj(v) .* numerator, dims=1)
    return numerator ./ denominators
end

mpdr_weights(am::AbstractArrayManifold, Rxx, angles, f, c=c_0; kwargs...) = mvdr_weights(am, Rxx, angles, f, c; kwargs...)

const capon_weights = mpdr_weights