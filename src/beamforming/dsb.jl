"""
References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function dsb_weights(am::AbstractArrayManifold, angles, f, c=c_0; kwargs...)
    return am(angles, f, c; kwargs...)/length(am)
end