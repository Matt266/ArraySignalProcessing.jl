
"""
References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function mvdr_weights(pa::AbstractPhasedArray, Rnn, f, angles; c=c_0, coords=:azel)
    v = steer(pa, f, angles; c=c, coords=coords)
    return (inv(Rnn)*v)/(v'*inv(Rnn)*v)
end

mpdr_weights(pa::AbstractPhasedArray, Rxx, f, angles; c=c_0, coords=:azel) = mvdr_weights(pa, Rxx, f, angles; c=c, coords=coords)

const capon_weights = mpdr_weights