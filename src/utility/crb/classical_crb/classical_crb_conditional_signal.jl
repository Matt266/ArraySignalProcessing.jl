"""
classical_crb_conditional_signal(s, nvar, am, θ; s_unwanted=true, nvar_unwanted=true)
Calculates the conditional (unkown nonrandom source signal) classical Cramér-Rao Bound (CRB) for a complex Gaussian random observation vector x.

arguments:
---------
    s: DxK source signal matrix for K sample snapshots of D source signals
    nvar: noise variance (nvar)
    am: Array Manifold Matrix am(θ) as function of the parameter vector θ
    θ: parameter vector of the Array Manifold Matrix to evaluate the CRB
    s_unwanted: treat source signals s as unwanted parameter (default: true)
    nvar_unwanted: treat noise variance (nvar) as unwanted parameter (default: true)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002. (8.85), (8.200), (8.201), (8.202), (8.204), (8.205)
"""
function classical_crb_conditional_signal(s, nvar, am, θ::AbstractVector; s_unwanted=true, nvar_unwanted=true)
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
        nvar_val = θ_full[idx_n][1]
        return nvar_val * collect(I(N))
    end
    
    function m(θ_full)
        θ_val = θ_full[idx_θ]
        s_vec_val = θ_full[idx_s]
        s_ri = reshape(s_vec_val, 2D, K)
        s_mat = s_ri[1:D, :] .+ 1im .* s_ri[D+1:end, :]
        return am(θ_val) * s_mat
    end

    if isnothing(θu)
        return classical_crb(θw, Kx, m, K)
    else
        return classical_crb(θw, θu, Kx, m, K)
    end
end

"""
classical_crb_conditional_signal(s, nvar, am, θw, θu; s_unwanted=true, nvar_unwanted=true)
Calculates the conditional (unkown nonrandom source signal) classical Cramér-Rao Bound (CRB) for a complex Gaussian random observation vector x.

arguments:
---------
    s: DxK source signal matrix for K sample snapshots of D source signals
    nvar: noise variance (nvar)
    am: Array Manifold Matrix am(θw, θu) as function of the parameter vectors θw and θu
    θw: vector of wanted parameters 
    θu: vector of unwanted parameters
    s_unwanted: treat source signals s as unwanted parameter (default: true)
    nvar_unwanted: treat noise variance (nvar) as unwanted parameter (default: true)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function classical_crb_conditional_signal(s, nvar, am, θw::AbstractVector, θu::AbstractVector; s_unwanted=true, nvar_unwanted=true)
    N = size(am(θw, θu), 1)
    D, K = size(s)
    
    # Van Trees (8.202): F = [Re[F(1)]^T, Im[F(1)]^T, ...]^T
    s_vec = vec(vcat(real(s), imag(s)))
    nvar_vec = [nvar]
    
    L_θw = length(θw)
    L_θu = length(θu)
    L_s  = length(s_vec)
    L_n  = 1

    if s_unwanted && nvar_unwanted
        θw_combined = copy(θw)
        θu_combined = vcat(θu, s_vec, nvar_vec)
        
        L_wc = length(θw_combined)
        
        idx_θw = 1:L_θw
        idx_θu = L_wc+1 : L_wc+L_θu
        idx_s  = L_wc+L_θu+1 : L_wc+L_θu+L_s
        idx_n  = L_wc+L_θu+L_s+1 : L_wc+L_θu+L_s+L_n
        
    elseif !s_unwanted && !nvar_unwanted
        θw_combined = vcat(θw, s_vec, nvar_vec)
        θu_combined = copy(θu)
        
        L_wc = length(θw_combined)
        
        idx_θw = 1:L_θw
        idx_s  = L_θw+1 : L_θw+L_s
        idx_n  = L_θw+L_s+1 : L_θw+L_s+L_n
        
        idx_θu = L_wc+1 : L_wc+L_θu
        
    elseif !s_unwanted && nvar_unwanted
        θw_combined = vcat(θw, s_vec)
        θu_combined = vcat(θu, nvar_vec)
        
        L_wc = length(θw_combined)
        
        idx_θw = 1:L_θw
        idx_s  = L_θw+1 : L_θw+L_s
        
        idx_θu = L_wc+1 : L_wc+L_θu
        idx_n  = L_wc+L_θu+1 : L_wc+L_θu+L_n
        
    else # s_unwanted && !nvar_unwanted
        θw_combined = vcat(θw, nvar_vec)
        θu_combined = vcat(θu, s_vec)
        
        L_wc = length(θw_combined)
        
        idx_θw = 1:L_θw
        idx_n  = L_θw+1 : L_θw+L_n
        
        idx_θu = L_wc+1 : L_wc+L_θu
        idx_s  = L_wc+L_θu+1 : L_wc+L_θu+L_s
    end

    function Kx(θ_full)
        nvar_val = θ_full[idx_n][1]
        return nvar_val * collect(I(N))
    end
    
    function m(θ_full)
        θw_val = θ_full[idx_θw]
        θu_val = θ_full[idx_θu]
        s_vec_val = θ_full[idx_s]
        s_ri = reshape(s_vec_val, 2D, K)
        s_mat = s_ri[1:D, :] .+ 1im .* s_ri[D+1:end, :]
        return am(θw_val, θu_val) * s_mat
    end

    return classical_crb(θw_combined, θu_combined, Kx, m, K)
end