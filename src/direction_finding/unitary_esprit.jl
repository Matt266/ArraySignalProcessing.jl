"""
unitary_esprit(X, J1, Δ, d, f; c=c_0, TLS = true, side = :left)

Calculates the DoAs using the unitary esprit.
Requires a centrosymmetric array geometry.

arguments:
----------
    X: data matrix of the array (NO CONCATENATION of subarrays)
    J1: selection matrix for the first subarray
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
M. Haardt and J. A. Nossek, ‘Unitary ESPRIT: how to obtain increased estimation accuracy with a reduced computational burden’, IEEE Trans. Signal Process., vol. 43, no. 5, pp. 1232–1242, May 1995.
"""
function unitary_esprit(X, J1, Δ, d, f; c=c_0, TLS = true, side = :left)
    # NxN exchange matrix
    II(N) = begin
        return rotl90(Matrix(I,N,N))
    end
    
    # unitary matrices
    Q(N) = begin
        @assert N>=2
        n = Int(floor(N/2))
        if N%2 == 0
            # N even
            [I(n) 1im*I(n); II(n) -1im*II(n)]/sqrt(2)
        else
            # N odd
            [I(n) zeros(n,1) 1im*I(n); zeros(1,n) sqrt(2) zeros(1, n); II(n) zeros(n,1) -1im*II(n)]/sqrt(2)
        end
    end

    m = size(J1)[1] # number elements in subarray
    M = size(J1)[2] # number elements in array
    K1 = Q(m)'*(J1+II(m)*J1*II(M))*Q(M)
    K2 = Q(m)'*1im*(J1-II(m)*J1*II(M))*Q(M)

    Y = Q(M)'*X
    U, _ = svd([real(Y) imag(Y)])
    Es = U[:,1:d]

    C1 = K1*Es
    C2 = K2*Es

    if(TLS)
        # TLS solution
        E, _ = svd([C1 C2]'*[C1 C2])
        E12 = E[1:d, (1:d).+d]
        E22 = E[(1:d).+d, (1:d).+d]
        Ψ = -E12*inv(E22)
    else
        # LS solution
        Ψ = C1 \ C2
    end

    Φ = eigvals(Ψ, sortby= λ -> -abs(λ))

    # calculate the directions of arrival (DoAs) from Φ
    Μ = 2atan.(real(Φ))
    k = (2π*f)/c

    # orientation of displacement vector
    ϕ0 = atan(Δ[2], Δ[1])  

    # angle estimates to the left of Δ
    Θ1 = ϕ0 .+ acos.(min.(max.(Μ ./ (k * norm(Δ)), -1), 1))
    Θ1 =  mod.(Θ1 .+ π, 2π) .- π

    # angle estimates to the right of Δ
    Θ2 = ϕ0 .- acos.(min.(max.(Μ ./ (k * norm(Δ)), -1), 1))
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
    # Θ = asin.(Μ/(k*Δ))
end