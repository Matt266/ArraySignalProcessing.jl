# for solving lasso with LeastSquares from ProximalOperators.jl
function value_and_gradient(f::ProximalOperators.LeastSquares, X)
    val = f(X)
    grad = similar(X)
    gradient!(grad, f, X)
    return val, grad
end

"""
lasso(Y, A, λ=1e-2; max_iter=300, tol=1e-6)

LASSO DOA estimation. Returns a vector representing the estimated, on-grid, spatial power spectrum of the signals. Estimated 
DOAs are the grid positions for which the spectrum crosses a certain threshold, as shown in the 'LASSO.ipynb' example.    

arguments:
----------
    Y: Data matrix of the array
    A: Dictionary matrix of array response vectors from the angle grid 
    λ: Regularization parameter for the LASSO problem
    maxit: maximum iterations for optimization
    tol: tolerance for optimization 

References:
-----------
Z. Yang, J. Li, P. Stoica, and L. Xie, ‘Sparse methods for direction-of-arrival estimation’, arXiv [cs.IT], 30-Sep-2016.
"""
function lasso(Y, A, λ=1e-2; maxit=100, tol=1e-6, kwargs...)
    f = LeastSquares(A, Y)
    g = NormL21(λ, 2)

    X0_eltype = promote_type(eltype(A), eltype(Y))
    X0 = fill!(similar(A, X0_eltype, size(A,2), size(Y,2)), zero(X0_eltype))
    ffb = ProximalAlgorithms.FastForwardBackward(;maxit=maxit, tol=tol, kwargs...)
    solution, _ = ffb(x0=X0, f=f, g=g)
    return vec(sum(abs2, solution, dims=2))
end