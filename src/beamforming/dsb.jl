"""
A: Array Manifold Matrix

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function dsb_weights(A::AbstractMatrix)
    return A/size(A, 1)
end