"""
ols(Y, A, d)

Orthogonal Least Squares (OLS) DOA estimation. Returns a vector representing the estimated, on-grid, sparse, spatial power spectrum of the signals. Estimated 
DOAs are the angles corresponding to indices of the non-zero values of the output spectrum.

arguments:
----------
    Y: Data matrix of the array
    A: Dictionary matrix of array response vectors from the angle grid 
    d: number of sources

References:
-----------
A. Hashemi and H. Vikalo, "Sparse linear regression via generalized orthogonal least-squares," 2016 IEEE Global Conference on Signal and Information Processing (GlobalSIP), Washington, DC, USA, 2016, pp. 1305-1309, doi: 10.1109/GlobalSIP.2016.7906052.

S. F. Cotter, B. D. Rao, Kjersti Engan and K. Kreutz-Delgado, "Sparse solutions to linear inverse problems with multiple measurement vectors," in IEEE Transactions on Signal Processing, vol. 53, no. 7, pp. 2477-2488, July 2005, doi: 10.1109/TSP.2005.849172.
"""
function ols(Y, A, d)
    r = copy(Y)
    Λ = Int[]

    cost(i) = begin
        Ψ = A[:, [Λ; i]]
        PΨ = Ψ*inv(Ψ'*Ψ)*Ψ'
        return norm((I-PΨ)*Y)^2
    end

    for _ in 1:d
        candidates = setdiff(1:size(A, 2), Λ)
        
        # select N indices that give the smallest cost and add them to the support
        idx = candidates[argmin(cost.(candidates))]
        push!(Λ, idx...)

        Ψ = A[:, Λ]
        X = Ψ \ Y
        r = Y - Ψ*X
    end
    
    Ψ = A[:, Λ]
    s_eltype = real(promote_type(eltype(A), eltype(Y)))
    s = fill!(similar(A, s_eltype, size(A,2)), zero(s_eltype))
    s[Λ] .= vec(sum(abs2, Ψ \ Y, dims=2))
    return s
end