
"""
References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function mvdr_weights(pa::AbstractPhasedArray, Rnn, angles, f, c=c_0; kwargs...)
    v = steer(pa, angles, f, c; kwargs...)
    #(inv(Rnn)*v)/(v'*inv(Rnn)*v) 
    #   -> adapted to steering matrix where each column is a steering vector
    Rnn_inv = inv(Rnn)
    numerator = Rnn_inv * v
    denominators = sum(conj(v) .* numerator, dims=1)
    return numerator ./ denominators
end

mpdr_weights(pa::AbstractPhasedArray, Rxx, angles, f, c=c_0; kwargs...) = mvdr_weights(pa, Rxx, angles, f, c; kwargs...)

const capon_weights = mpdr_weights