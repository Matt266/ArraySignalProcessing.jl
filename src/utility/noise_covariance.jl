function whitenoise(am::AbstractArrayManifold, σ²)
    return σ²*I(length(am))
end

"""
References:
-----------
W. Herbordt, Sound capture for human / machine interfaces, 2005th ed. Berlin, Germany: Springer, 2005.
"""
function diffnoise(am::AbstractArrayManifold, σ², f, c=c_0)
    ω = 2π*f
    k = ω/c

    # si(x) = sinc(x/π)
    si(x) = sinc.(x ./ π)

    #p(x, i) = [e for e in x.manifold.r[:,i]] 
    #Γ(x, n, m, k) = si(k*norm(p(x,m)-p(x,n)))
    function Γ(n, m)
        d = @view(am.r[:, m]) .- @view(am.r[:, n])
        si(k * norm(d))
    end

    n = 1:length(am)

    #σ²*Γ.(Ref(am), n, n', Ref(k))
    return σ² .* [Γ(i, j) for i in n, j in n] 
end

# TODO: this is a quick fix to have diffnoise work with TappedDelayLines
#       References are needed to check how diffnoise is modelled
#       for taps correctly!
function diffnoise(tdl::TappedDelayLineManifold, σ², f, c=c_0)
    J = tdl.num_taps
    Ts = 1 / tdl.fs
    ω = 2π*f

    # diffuse spatial noise correlation
    R_spatial = diffnoise(tdl.sub_manifold, σ², f, c)

    # temporal correlation between taps  
    R_time = [exp(-1im * ω * (i-j) * Ts) for i=0:J-1, j=0:J-1]

    # extend to TDL: Kronecker product
    return kron(R_time, R_spatial)
end

# TODO: this is a quick fix to have whitenoise work with TappedDelayLines
#       References are needed to check how whitenoise is modelled
#       for taps correctly!
function whitenoise(tdl::TappedDelayLineManifold, σ², f)
    J  = tdl.num_taps
    Ts = 1 / tdl.fs
    ω  = 2π * f

    M = length(tdl.sub_manifold)

    # temporal correlation between taps  
    R_time = [exp(-1im * ω * (i-j) * Ts) for i = 0:J-1, j = 0:J-1]

    # spatially white between sensors
    R_spatial_white = I(M) 

    # extend to TDL: Kronecker product
    return σ² * kron(R_time, R_spatial_white)
end