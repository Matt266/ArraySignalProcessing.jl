"""
References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function dsb_weights(pa::AbstractPhasedArray, f, angles; c=c_0, coords=:azel)
    return steer(pa, f, angles; c=c, coords=coords)/length(pa)
end