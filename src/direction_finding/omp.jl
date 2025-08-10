"""
omp(Y, A, d)

Orthogonal Matching Pursuit (OMP) DOA estimation. Returns a vector representing the estimated, on-grid, sparse, spatial power spectrum of the signals. Estimated 
DOAs are the angles corresponding to indices of the non-zero values of the output spectrum.

arguments:
----------
    Y: Data matrix of the array
    A: Dictionary matrix of array response vectors from the angle grid 
    d: number of sources

References:
-----------
V. N. Xuan, K. Hartmann, W. Weihs, and O. Loffeld, ‘Modified orthogonal matching pursuit for multiple measurement vector with joint sparsity in super-resolution compressed sensing’, in 2017 51st Asilomar Conference on Signals, Systems, and Computers, Pacific Grove, CA, USA, 2017.

J. Chen and X. Huo, ‘Sparse representations for multiple measurement vectors (MMV) in an over-complete dictionary’, in Proceedings. (ICASSP ’05). IEEE International Conference on Acoustics, Speech, and Signal Processing, 2005, Philadelphia, Pennsylvania, USA, 2006.
"""
function omp(Y, A, d)
    r = copy(Y)
    Λ = Int[]

    for _ in 1:d
        corr = A'*r
        idx = argmax(sqrt.(vec(sum(abs2, corr, dims=2))))
        push!(Λ, idx)

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