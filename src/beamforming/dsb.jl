"""
References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function dsb_weights(pa::AbstractPhasedArray, angles, f, c=c_0; kwargs...)
    return steer(pa, angles, f, c; kwargs...)/length(pa)
end