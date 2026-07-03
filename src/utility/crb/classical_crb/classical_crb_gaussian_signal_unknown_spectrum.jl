"""
classical_crb_gaussian_signal_unknown_spectrum(Rss, nvar, am, θ, K=1; Rss_unwanted=true, nvar_unwanted=true)
Calculates the unconditional (complex gaussian zero-mean random signal model) classical Cramér-Rao Bound (CRB) for a complex Gaussian random observation vector x
and unknown source spectrum Rss.

arguments:
---------
    Rss: DxD covariance matrix of the source signals
    nvar: noise variance
    am: Array Manifold Matrix am(θ) as function of the parameter vector θ
    θ: parameter vector of the Array Manifold Matrix to evaluate the CRB
    K: number of sample snapshots (default is 1)
    Rss_unwanted: treat source signals s as unwanted parameter (default: true)
    nvar_unwanted: treat noise variance (nvar) as unwanted parameter (default: true)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002. (8.85), (8.200), (8.201), (8.202), (8.204), (8.205)
Check text paragraph between (8.89) and (8.90): Sx for Kx, also signal s is a sample function
from a zero-mean stationary complex Gaussian random process so m = 0
"""
function classical_crb_gaussian_signal_unknown_spectrum(Rss, nvar, am, θ::AbstractVector, K::Number=1; Rss_unwanted=true, nvar_unwanted=true)
    N = size(am(θ), 1)
    D = size(Rss, 1)
    
    # (8.87)
    Rss_vec = vec(vcat(real(Rss), imag(Rss)))
    nvar_vec = [nvar]
    
    L_θ = length(θ)
    L_Rss = length(Rss_vec)
    L_n = 1

    if Rss_unwanted && nvar_unwanted
        θw = copy(θ)
        θu = vcat(Rss_vec, nvar_vec)
        
        idx_θ = 1:L_θ
        idx_Rss = L_θ+1 : L_θ+L_Rss
        idx_n = L_θ+L_Rss+1 : L_θ+L_Rss+L_n
        
    elseif !Rss_unwanted && !nvar_unwanted
        θw = vcat(θ, Rss_vec, nvar_vec)
        θu = nothing
        
        idx_θ = 1:L_θ
        idx_Rss = L_θ+1 : L_θ+L_Rss
        idx_n = L_θ+L_Rss+1 : L_θ+L_Rss+L_n
        
    elseif !Rss_unwanted && nvar_unwanted
        θw = vcat(θ, Rss_vec)
        θu = nvar_vec
        
        idx_θ = 1:L_θ
        idx_Rss = L_θ+1 : L_θ+L_Rss
        idx_n = L_θ+L_Rss+1 : L_θ+L_Rss+L_n
        
    else # s_unwanted && !nvar_unwanted
        θw = vcat(θ, nvar_vec)
        θu = Rss_vec
        
        idx_θ = 1:L_θ
        idx_n = L_θ+1 : L_θ+L_n
        idx_Rss = L_θ+L_n+1 : L_θ+L_n+L_Rss
    end

    function Kx(θ_full)
        θ_val = θ_full[idx_θ]
        Rss_vec_val = θ_full[idx_Rss]
        Rss_ri = reshape(Rss_vec_val, 2D, D)
        Rss_mat = Rss_ri[1:D, :] .+ 1im .* Rss_ri[D+1:end, :]
        nvar_val = θ_full[idx_n][1]

        A = am(θ_val)
        return A*Rss_mat*A'+ nvar_val * collect(I(N))
    end
    
    # signal and noise are zero-mean random processes
    m = fill!(similar(Rss, N, K), 0)

    if isnothing(θu)
        return classical_crb(θw, Kx, m, K)
    else
        return classical_crb(θw, θu, Kx, m, K)
    end
end

"""
classical_crb_gaussian_signal_unknown_spectrum(Rss, nvar, am, θw, θu, K=1; Rss_unwanted=true, nvar_unwanted=true) 
Calculates the unconditional (complex gaussian zero-mean random signal model) classical Cramér-Rao Bound (CRB) for a complex Gaussian random observation vector x
and unknown source spectrum Rss.

arguments:
---------
    Rss: DxD covariance matrix of the source signals
    nvar: noise variance
    am: Array Manifold Matrix am(θw, θu) as function of the parameter vectors θw and θu
    θw: vector of wanted parameters 
    θu: vector of unwanted parameters
    K: number of sample snapshots (default is 1)
    Rss_unwanted: treat source signals s as unwanted parameter (default: true)
    nvar_unwanted: treat noise variance (nvar) as unwanted parameter (default: true)

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
Check text paragraph between (8.89) and (8.90): Sx for Kx, also signal s is a sample function
from a zero-mean stationary complex Gaussian random process so m = 0
"""
function classical_crb_gaussian_signal_unknown_spectrum(Rss, nvar, am, θw::AbstractVector, θu::AbstractVector, K::Number=1; Rss_unwanted=true, nvar_unwanted=true) 
    N = size(am(θw, θu), 1)
    D = size(Rss, 1)

    Rss_vec = vec(vcat(real(Rss), imag(Rss)))
    nvar_vec = [nvar]
    
    L_θw = length(θw)
    L_θu = length(θu)
    L_Rss  = length(Rss_vec)
    L_n  = 1

    if Rss_unwanted && nvar_unwanted
        θw_combined = copy(θw)
        θu_combined = vcat(θu, Rss_vec, nvar_vec)
        
        L_wc = length(θw_combined)
        
        idx_θw = 1:L_θw
        idx_θu = L_wc+1 : L_wc+L_θu
        idx_Rss  = L_wc+L_θu+1 : L_wc+L_θu+L_Rss
        idx_n  = L_wc+L_θu+L_Rss+1 : L_wc+L_θu+L_Rss+L_n
        
    elseif !Rss_unwanted && !nvar_unwanted
        θw_combined = vcat(θw, Rss_vec, nvar_vec)
        θu_combined = copy(θu)
        
        L_wc = length(θw_combined)
        
        idx_θw = 1:L_θw
        idx_Rss  = L_θw+1 : L_θw+L_Rss
        idx_n  = L_θw+L_Rss+1 : L_θw+L_Rss+L_n
        
        idx_θu = L_wc+1 : L_wc+L_θu
        
    elseif !Rss_unwanted && nvar_unwanted
        θw_combined = vcat(θw, Rss_vec)
        θu_combined = vcat(θu, nvar_vec)
        
        L_wc = length(θw_combined)
        
        idx_θw = 1:L_θw
        idx_Rss  = L_θw+1 : L_θw+L_Rss
        
        idx_θu = L_wc+1 : L_wc+L_θu
        idx_n  = L_wc+L_θu+1 : L_wc+L_θu+L_n
        
    else # s_unwanted && !nvar_unwanted
        θw_combined = vcat(θw, nvar_vec)
        θu_combined = vcat(θu, Rss_vec)
        
        L_wc = length(θw_combined)
        
        idx_θw = 1:L_θw
        idx_n  = L_θw+1 : L_θw+L_n
        
        idx_θu = L_wc+1 : L_wc+L_θu
        idx_Rss  = L_wc+L_θu+1 : L_wc+L_θu+L_Rss
    end

    function Kx(θ_full)
        θw_val = θ_full[idx_θw]
        θu_val = θ_full[idx_θu]
        Rss_vec_val = θ_full[idx_Rss]
        Rss_ri = reshape(Rss_vec_val, 2D, D)
        Rss_mat = Rss_ri[1:D, :] .+ 1im .* Rss_ri[D+1:end, :]
        nvar_val = θ_full[idx_n][1]

        A = am(θw_val, θu_val)
        return A*Rss_mat*A'+ nvar_val * collect(I(N))
    end
    
    # signal and noise are zero-mean random processes
    m = fill!(similar(Rss, N, K), 0)

    return classical_crb(θw_combined, θu_combined, Kx, m, K)
end