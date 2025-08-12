# Copy of NormL21 from ProximalOperators.jl:
# https://github.com/JuliaFirstOrder/ProximalOperators.jl/blob/master/src/functions/normL21.jl
#
# But with minor changes to work with CUDA
import ProximalOperators: prox!
struct MyNormL21{R, I}
    lambda::R
    dim::I
    function MyNormL21{R,I}(lambda::R, dim::I) where {R, I}
            new(lambda, dim)
    end
end

is_convex(f::Type{<:MyNormL21}) = true

MyNormL21(lambda::R=1, dim::I=1) where {R, I} = MyNormL21{R, I}(lambda, dim)

function (f::MyNormL21)(X)
    return f.lambda * sum(sqrt, sum(abs2, X, dims=f.dim))
end

function prox!(Y, f::MyNormL21, X, gamma)
    gl = gamma * f.lambda
    shrink(norm) = max(1 - gl / norm, zero(real(eltype(X))))
    slice_norms = sqrt.(sum(abs2, X, dims=f.dim))
    scal = shrink.(slice_norms)

    @. Y = scal * X 
    return f.lambda * sum(scal .* slice_norms)
end

function prox_naive(f::MyNormL21, X, gamma)
    gl = gamma * f.lambda
    shrink(norm) = max(1 - gl / norm, zero(real(eltype(X))))
    slice_norms = sqrt.(sum(abs2, X, dims=f.dim))
    scal = shrink.(slice_norms)

    Y = scal .* X 
    return Y, f.lambda * f(Y)
end

"""
lasso(Y, A, λ=1e-2; max_iter=300, tol=1e-6)

LASSO DOA estimation. Returns a vector representing the estimated, on-grid, spatial power spectrum of the signals. Estimated 
DOAs are the grid positions for which the spectrum crosses a certain threshold, as shown in the 'LASSO.ipynb' example.    

arguments:
----------
    Y: Data matrix of the array (Group LASSO with L21-penalty, when Y has >1 columns)
    A: Dictionary matrix of array response vectors from the angle grid 
    λ: Regularization parameter for the LASSO problem
    maxit: maximum iterations for optimization
    tol: tolerance for optimization 

References:
-----------
Z. Yang, J. Li, P. Stoica, and L. Xie, ‘Sparse methods for direction-of-arrival estimation’, arXiv [cs.IT], 30-Sep-2016.
"""
function lasso(Y, A, λ=1e-2; kwargs...)
    X0 = fill!(similar(Y, eltype(Y), size(A,2), size(Y,2)), zero(eltype(Y)))
    f = X -> sum(abs2, A * X - Y) # least squares
    g = MyNormL21(λ, 2)
    ffb = ProximalAlgorithms.FastForwardBackward(; kwargs...)
    solution, _ = ffb(x0=X0, f=f, g=g)
    return vec(sum(abs2, solution, dims=2))
end

"""
M: ULA length
N: Number Snapshots
SNR: SNR in dB 

References:
-----------
Sparse Methods for Direction-of-Arrival Estimation. Zai Yang, Jian Li, Petre Stoica, and Lihua Xie. January 10, 2017
Optimal λ for ULA (p.45; below eq. (190))
Where source [85] is referenced for that:
    Y. Li, Y. Chi, Off-the-grid line spectrum denoising and estimation with multiple measurement vectors,
    IEEE Transactions on Signal Processing 64 (5) (2016) 1257–1269.
"""
function λ_opt(M, N, SNR)
    return sqrt(M*(N+log(M)+sqrt(2N*log(M)))*sqrt(snr2nvar(SNR)))
end