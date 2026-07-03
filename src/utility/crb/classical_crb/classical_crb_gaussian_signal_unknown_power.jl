"""
classical_crb_gaussian_signal_unknown_power(svar, nvar, am, θ, K=1; svar_unwanted=true, nvar_unwanted=true)
Calculates the unconditional (complex gaussian zero-mean random signal model) classical Cramér-Rao Bound (CRB) for a complex Gaussian random observation vector x
and uncorrelated source signals with unknown power.

arguments:
---------
    svar: size D vector of source signal variances
    nvar: noise variance
    am: Array Manifold Matrix am(θ) as function of the parameter vector θ
    θ: parameter vector of the Array Manifold Matrix to evaluate the CRB
    K: number of sample snapshots (default is 1)
    svar_unwanted: treat source signals s as unwanted parameter (default: true)
    nvar_unwanted: treat noise variance (nvar) as unwanted parameter (default: true)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002. (8.85), (8.200), (8.201), (8.202), (8.204), (8.205)
Check text paragraph between (8.89) and (8.90): Sx for Kx, also signal s is a sample function
from a zero-mean stationary complex Gaussian random process so m = 0
"""
function classical_crb_gaussian_signal_unknown_power(svar, nvar, am, θ::AbstractVector, K::Number=1; svar_unwanted=true, nvar_unwanted=true)
    N = size(am(θ), 1)
    
    svar_vec = isa(svar, Number) ? [svar] : svar
    nvar_vec = [nvar]
    
    L_θ = length(θ)
    L_svar = length(svar_vec)
    L_n = 1

    if svar_unwanted && nvar_unwanted
        θw = copy(θ)
        θu = vcat(svar_vec, nvar_vec)
        
        idx_θ = 1:L_θ
        idx_svar = L_θ+1 : L_θ+L_svar
        idx_n = L_θ+L_svar+1 : L_θ+L_svar+L_n
        
    elseif !svar_unwanted && !nvar_unwanted
        θw = vcat(θ, svar_vec, nvar_vec)
        θu = nothing
        
        idx_θ = 1:L_θ
        idx_svar = L_θ+1 : L_θ+L_svar
        idx_n = L_θ+L_svar+1 : L_θ+L_svar+L_n
        
    elseif !svar_unwanted && nvar_unwanted
        θw = vcat(θ, svar_vec)
        θu = nvar_vec
        
        idx_θ = 1:L_θ
        idx_svar = L_θ+1 : L_θ+L_svar
        idx_n = L_θ+L_svar+1 : L_θ+L_svar+L_n
        
    else # s_unwanted && !nvar_unwanted
        θw = vcat(θ, nvar_vec)
        θu = svar_vec
        
        idx_θ = 1:L_θ
        idx_n = L_θ+1 : L_θ+L_n
        idx_svar = L_θ+L_n+1 : L_θ+L_n+L_svar
    end

    function Kx(θ_full)
        θ_val = θ_full[idx_θ]
        svar_vec_val = θ_full[idx_svar]
        nvar_val = θ_full[idx_n][1]

        # (8.160)
        A = am(θ_val)
        A_scaled = svar_vec_val' .* A
        return A_scaled*A'+ nvar_val * collect(I(N))
    end
    
    # signal and noise are zero-mean random processes
    m = fill!(similar(svar_vec, N, K), 0)

    if isnothing(θu)
        return classical_crb(θw, Kx, m, K)
    else
        return classical_crb(θw, θu, Kx, m, K)
    end
end

"""
classical_crb_gaussian_signal_unknown_power(svar, nvar, am, θw, θu, K=1; svar_unwanted=true, nvar_unwanted=true) 
Calculates the unconditional (complex gaussian zero-mean random signal model) classical Cramér-Rao Bound (CRB) for a complex Gaussian random observation vector x
and uncorrelated source signals with unknown power.

arguments:
---------
    svar: size D vector of source signal variances
    nvar: noise variance
    am: Array Manifold Matrix am(θw, θu) as function of the parameter vectors θw and θu
    θw: vector of wanted parameters 
    θu: vector of unwanted parameters
    K: number of sample snapshots (default is 1)
    svar_unwanted: treat source signals s as unwanted parameter (default: true)
    nvar_unwanted: treat noise variance (nvar) as unwanted parameter (default: true)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
Check text paragraph between (8.89) and (8.90): Sx for Kx, also signal s is a sample function
from a zero-mean stationary complex Gaussian random process so m = 0
"""
function classical_crb_gaussian_signal_unknown_power(svar, nvar, am, θw::AbstractVector, θu::AbstractVector, K::Number=1; svar_unwanted=true, nvar_unwanted=true) 
    N = size(am(θw, θu), 1)

    svar_vec = isa(svar, Number) ? [svar] : svar
    nvar_vec = [nvar]
    
    L_θw = length(θw)
    L_θu = length(θu)
    L_svar  = length(svar_vec)
    L_n  = 1

    if svar_unwanted && nvar_unwanted
        θw_combined = copy(θw)
        θu_combined = vcat(θu, svar_vec, nvar_vec)
        
        L_wc = length(θw_combined)
        
        idx_θw = 1:L_θw
        idx_θu = L_wc+1 : L_wc+L_θu
        idx_svar  = L_wc+L_θu+1 : L_wc+L_θu+L_svar
        idx_n  = L_wc+L_θu+L_svar+1 : L_wc+L_θu+L_svar+L_n
        
    elseif !svar_unwanted && !nvar_unwanted
        θw_combined = vcat(θw, svar_vec, nvar_vec)
        θu_combined = copy(θu)
        
        L_wc = length(θw_combined)
        
        idx_θw = 1:L_θw
        idx_svar  = L_θw+1 : L_θw+L_svar
        idx_n  = L_θw+L_svar+1 : L_θw+L_svar+L_n
        
        idx_θu = L_wc+1 : L_wc+L_θu
        
    elseif !svar_unwanted && nvar_unwanted
        θw_combined = vcat(θw, svar_vec)
        θu_combined = vcat(θu, nvar_vec)
        
        L_wc = length(θw_combined)
        
        idx_θw = 1:L_θw
        idx_svar  = L_θw+1 : L_θw+L_svar
        
        idx_θu = L_wc+1 : L_wc+L_θu
        idx_n  = L_wc+L_θu+1 : L_wc+L_θu+L_n
        
    else # s_unwanted && !nvar_unwanted
        θw_combined = vcat(θw, nvar_vec)
        θu_combined = vcat(θu, svar_vec)
        
        L_wc = length(θw_combined)
        
        idx_θw = 1:L_θw
        idx_n  = L_θw+1 : L_θw+L_n
        
        idx_θu = L_wc+1 : L_wc+L_θu
        idx_svar  = L_wc+L_θu+1 : L_wc+L_θu+L_svar
    end

    function Kx(θ_full)
        θw_val = θ_full[idx_θw]
        θu_val = θ_full[idx_θu]
        svar_vec_val = θ_full[idx_svar]
        nvar_val = θ_full[idx_n][1]

        # (8.160)
        A = am(θw_val, θu_val)
        A_scaled = svar_vec_val' .* A
        return A_scaled*A'+ nvar_val * collect(I(N))
    end
    
    # signal and noise are zero-mean random processes
    m = fill!(similar(svar_vec, N, K), 0)

    return classical_crb(θw_combined, θu_combined, Kx, m, K)
end