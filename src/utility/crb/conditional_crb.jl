"""
conditional_crb(s, nvar, am, θ; s_unwanted=true, nvar_unwanted=true)
Calculates the conditional (nonrandom signal model) classical Cramér-Rao Bound (CRB) for a complex Gaussian random observation vector x.

arguments:
---------
    s: DxK source signal matrix for K sample snapshots of D source signals
    nvar: Noise variance σₙ²
    am: Array Manifold Matrix am(θ) as function of the parameter vector θ
    θ: parameter vector of the Array Manifold Matrix to evaluate the CRB
    K: number of sample snapshots (default is 1)
    s_unwanted: treat source signals s as unwanted parameter (default: true)
    nvar_unwanted: treat noise variance (nvar) as unwanted parameter (default: true)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002. (8.85), (8.200), (8.201), (8.202), (8.204), (8.205)
"""
function conditional_crb(s, nvar, am, θ; s_unwanted=true, nvar_unwanted=true)
    # number of sensors
    N = size(am(θ), 1)
    D, K = size(s)
    
    # Van Trees (8.202): F = [Re[F(1)]^T, Im[F(1)]^T, ...]^T
    s_vec = vec(vcat(real(s), imag(s)))
    nvar_vec = [nvar]
    
    L_θ = length(θ)
    L_s = length(s_vec)
    L_n = 1

    if s_unwanted && nvar_unwanted
        θw = copy(θ)
        θu = vcat(s_vec, nvar_vec)
        
        idx_θ = 1:L_θ
        idx_s = L_θ+1 : L_θ+L_s
        idx_n = L_θ+L_s+1 : L_θ+L_s+L_n
        
    elseif !s_unwanted && !nvar_unwanted
        θw = vcat(θ, s_vec, nvar_vec)
        θu = nothing
        
        idx_θ = 1:L_θ
        idx_s = L_θ+1 : L_θ+L_s
        idx_n = L_θ+L_s+1 : L_θ+L_s+L_n
        
    elseif !s_unwanted && nvar_unwanted
        θw = vcat(θ, s_vec)
        θu = nvar_vec
        
        idx_θ = 1:L_θ
        idx_s = L_θ+1 : L_θ+L_s
        idx_n = L_θ+L_s+1 : L_θ+L_s+L_n
        
    else # s_unwanted && !nvar_unwanted
        θw = vcat(θ, nvar_vec)
        θu = s_vec
        
        idx_θ = 1:L_θ
        idx_n = L_θ+1 : L_θ+L_n
        idx_s = L_θ+L_n+1 : L_θ+L_n+L_s
    end

    function Kx(θ_full)
        nvar = θ_full[idx_n][1]
        return nvar * I(N)
    end
    
    function m(θ_full)
        θ = θ_full[idx_θ]
        s_vec = θ_full[idx_s]
        s_ri = reshape(s_vec, 2D, K)
        s = s_ri[1:D, :] .+ 1im .* s_ri[D+1:end, :]
        return am(θ) * s
    end

    if isnothing(θu)
        return classical_crb(θw, Kx, m; K=1)
    else
        return classical_crb(θw, θu, Kx, m; K=1)
    end
end

"""
conditional_crb(θw::AbstractVector, θu::AbstractVector, Kx, m; K=1) 
Calculates the conditional (nonrandom signal model) classical Cramér-Rao Bound (CRB) for a complex Gaussian random observation vector x.

arguments:
---------
    s: DxK source signal matrix for K sample snapshots of D source signals
    nvar: Noise variance σₙ²
    am: Array Manifold Matrix am(θw, θu) as function of the parameter vectors θw and θu
    θw: vector of wanted parameters 
    θu: vector of unwanted parameters
    s_unwanted: treat source signals s as unwanted parameter (default: true)
    nvar_unwanted: treat noise variance (nvar) as unwanted parameter (default: true)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function conditional_crb(s, nvar, am, θw, θu; s_unwanted=true, nvar_unwanted=true) 
end
