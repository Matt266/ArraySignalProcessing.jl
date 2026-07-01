
"""
classical_log_likelihood(θ::AbstractVector, Kx, m, X)

Calculates the classical log-likelihood function.

arguments:
---------
    θ: parameter vector to evaluate the log-likelihood at
    Kx: function that returns the covariance matrix of x given θ, or a static matrix
    m: function that returns the mean of x given θ, or a static vector
    X: matrix where each column is a snapshot of the complex Gaussian random observation vector x

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002. (8.9) & (8.25)
"""
classical_log_likelihood(θ::AbstractVector, Kx, m, X) = begin
    num_snapshots = size(X, 2)
    Kx_val = applicable(Kx, θ) ? Kx(θ) : Kx
    m_val = applicable(m, θ) ? m(θ) : m
    X_c = X .- m_val

    # sum over k for (x_k^H-m^H) * Kx^-1 * (x_k-m) in (8.9)
    # is equal to tr((X-M)^H Kx^-1 \ (X-M)) with M being m repeated to the number of snapshots
    # this is equal to dot(X_c, Kx_val \ X_c) since tr(A^H B) = sum(conj(A).*B) = dot(A,B)
    return -num_snapshots*real(logdet(π*Kx_val)) - real(dot(X_c, Kx_val \ X_c))
end

"""
classical_fim(θ, Kx, m; K=1)

Calculates the classical Fisher Information Matrix (FIM) for a complex Gaussian random observation vector x.

arguments:
---------
    θ: parameter vector to evaluate the FIM at
    Kx: function that returns the covariance matrix given θ, or a static matrix
    m: function that returns the mean vector given θ, or a static vector
    K: number of sample snapshots (default is 1)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002. (8.34)
"""
function classical_fim(θ::AbstractVector, Kx, m; K=1)
    P = length(θ)

    K_val = applicable(Kx, θ) ? Hermitian(Kx(θ)) : Hermitian(Kx)
    m_val = applicable(m, θ) ? m(θ) : m
    
    N = size(K_val, 1)

    J_K_raw = applicable(Kx, θ) ? Zygote.jacobian(Kx, θ)[1] : nothing
    J_K = isnothing(J_K_raw) ? fill!(similar(K_val, N^2, P), 0) : J_K_raw
    J_m_raw = applicable(m, θ) ? Zygote.jacobian(m, θ)[1] : nothing
    J_m = isnothing(J_m_raw) ? fill!(similar(m_val, length(m_val), P), 0) : J_m_raw

    # Reshape N^2 x P Jacobian to an N x (N*P) matrix
    ∂K = reshape(J_K, N, N * P)
    Kinv_∂K = reshape(K_val \ ∂K , N, N, P)
    
    # compute tr(A * B) as vec(A^T)^T * vec(B) for all pairwise combinations (i, j) of P:
    Term1 = real.(
            transpose(reshape(permutedims(Kinv_∂K, (2, 1, 3)), N^2, P)) * reshape(Kinv_∂K, N^2, P)
            )
    Term2 = 2 .* real.(J_m' * (K_val \ J_m))
    return (Term1 .+ Term2) .* K
end

"""
classical_fim(θ, Kx, m; K=1)

Calculates the classical Fisher Information Matrix (FIM) for a complex Gaussian random observation vector x.
Returns Tuple of FIM partitions for wanted and unwanted parameters: (FIM_ww, FIM_wu, FIM_uw, FIM_uu)

arguments:
---------
    θw: vector of wanted parameters to evaluate the FIM at
    θu: vector of unwanted parameters to evaluate the FIM at
    Kx: function that returns the covariance matrix given θ, or a static matrix
    m: function that returns the mean vector given θ, or a static vector
    K: number of sample snapshots (default is 1)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002. (8.40)
"""
function classical_fim(θw::AbstractVector, θu::AbstractVector, Kx, m; K=1)
    FIM = classical_fim(vcat(θw,θu), Kx, m; K=K)
    # Partition to:
    # [ FIM_ww FIM_wu
    #   FIM_uw FIM_uu ]
    w = length(θw)
    return FIM[1:w, 1:w], FIM[1:w, w+1:end], FIM[w+1:end, 1:w], FIM[w+1:end, w+1:end]
end

"""
classical_crb(θ::AbstractVector, Kx, m; K=1)
Calculates the classical Cramér-Rao Bound (CRB) for a complex Gaussian random observation vector x.

arguments:
---------
    θ: parameter vector to evaluate the CRB at
    Kx: function that returns the covariance matrix given θ, or a static matrix
    m: function that returns the mean vector given θ, or a static vector
    K: number of sample snapshots (default is 1)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002. (8.25)
"""
function classical_crb(θ::AbstractVector, Kx, m; K=1) 
    return classical_fim(θ, Kx, m; K=K) \ I
end

"""
classical_crb(θw::AbstractVector, θu::AbstractVector, Kx, m; K=1) 
Calculates the classical Cramér-Rao Bound (CRB) for a complex Gaussian random observation vector x.

arguments:
---------
    θw: vector of wanted parameters to evaluate the CRB at
    θu: vector of unwanted parameters to evaluate the CRB at
    Kx: function that returns the covariance matrix given θ, or a static matrix
    m: function that returns the mean vector given θ, or a static vector
    K: number of sample snapshots (default is 1)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002. (8.43)
"""
function classical_crb(θw::AbstractVector, θu::AbstractVector, Kx, m; K=1) 
    FIM_ww, FIM_wu, FIM_uw, FIM_uu = classical_fim(θw, θu, Kx, m; K=K)
    return (FIM_ww - FIM_wu * (FIM_uu \ FIM_uw)) \ I
end