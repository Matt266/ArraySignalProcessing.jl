"""
esprit(Rzz, Δ, d, f; c=c_0, TLS = true, side = :left)

Calculates the TLS esprit estimator for the direction of arrival.
Returns a vector of tuples with each tuple holding the two ambigues 
DOAs corresponding to a source. 

arguments:
----------
    Rzz: covariance matrix of total array (both subarrays vertically concatenated),
    Δ: displacement vector between both subarrays
    d: number of sources
    f: center/operating frequency
    c: propagation speed of the wave
    TLS: calculates total least squares solution if 'true' (default),
        least squares if 'false'
    side: choose angles on the left (':left'), right (':right'), or both (':both') sides
        of the displacement vector to decide between the two ambigues angles per source.

References:
-----------
R. Roy and T. Kailath, ‘ESPRIT-estimation of signal parameters via rotational invariance techniques’, IEEE Trans. Acoust., vol. 37, no. 7, pp. 984–995, Jul. 1989.

H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function esprit(Rzz, Δ, d, f, c=c_0; TLS = true, side = :left)
    # number of sensors in the array (p)
    # and the subarrays (ps)
    p = size(Rzz, 1)
    ps = Int(p/2)

    eigs = eigen(Rzz)
    U = eigs.vectors[:, sortperm(abs.(eigs.values); rev=true)]

    Es = U[:,1:d]
    Ex = Es[(1:ps),:]
    Ey = Es[(1:ps).+(ps),:]

    # estimate Φ by exploiting the array symmetry
    if(TLS)
        # TLS-ESPRIT
        eigs = eigen(Hermitian([Ex Ey]'*[Ex Ey]))
        E = eigs.vectors[:, sortperm(abs.(eigs.values); rev=true)]
        
        E12 = E[1:d, (1:d).+d]
        E22 = E[(1:d).+d, (1:d).+d]
        Ψ = -E12*inv(E22)
    else
        # LS-ESPRIT
        Ψ = Ex \ Ey
    end

    Φ = eigvals(Array(Ψ), sortby= λ -> -abs(λ))

    # calculate the directions of arrival (DoAs) from Φ
    k = (2π*f)/c

    # orientation of displacement vector
    ϕ0 = atan(Δ[2], Δ[1])  

    # angle estimates to the left of Δ
    Θ1 = ϕ0 .+ acos.(min.(max.(angle.(Φ) ./ (k * norm(Δ)), -1), 1))
    Θ1 =  mod.(Θ1 .+ π, 2π) .- π

    # angle estimates to the right of Δ
    Θ2 = ϕ0 .- acos.(min.(max.(angle.(Φ) ./ (k * norm(Δ)), -1), 1))
    Θ2 = mod.(Θ2 .+ π, 2π) .- π

    # select ambigues angles left or right
    # respective to displacement vector Δ
    if side == :left
        return Θ1
    elseif side == :right
        return Θ2
    elseif side == :both
        return Θ1, Θ2
    else
       error("Invalid symbol for side: ':$(side)'. Valid options are: ':left', ':right', ':both'")
    end
    # for displacement along y-axis:
    # asin.(angle.(Φ)/(k*norm(Δ)))
end