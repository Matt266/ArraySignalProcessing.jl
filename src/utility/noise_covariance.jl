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

    # si(x) = sinc(x/π)
    si(x) = sinc.(x ./ π)

    #p(x, i) = [e for e in x.manifold.r[:,i]] 
    #Γ(x, n, m, k) = si(k*norm(p(x,m)-p(x,n)))
    function Γ(n, m)
        d = @view(pa.manifold.r[:, m]) .- @view(pa.manifold.r[:, n])
        si(k * norm(d))
    end

    n = 1:length(pa)

    #σ²*Γ.(Ref(pa), n, n', Ref(k))
    return σ² .* [Γ(i, j) for i in n, j in n] 
end