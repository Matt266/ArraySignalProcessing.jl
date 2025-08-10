function whitenoise(pa::AbstractPhasedArray, σ²)
    return σ²*I(length(pa))
end

"""
References:
-----------
W. Herbordt, Sound capture for human / machine interfaces, 2005th ed. Berlin, Germany: Springer, 2005.
"""
function diffnoise(pa::AbstractPhasedArray, σ², f, c=c_0)
    ω = 2π*f
    k = ω/c
    p(x, i) = [e for e in x.manifold.r[:,i]] 
    si(x) = sinc(x/π)
    Γ(x, n, m, k) = si(k*norm(p(x,m)-p(x,n)))
    n = 1:length(pa)
    return σ²*Γ.(Ref(pa), n, n', Ref(k))
end