"""
bpdn(Y, A, η=1e-2)

Basis Pursuit Denoising (BPDN) DOA estimation. Returns a vector representing the estimated, on-grid, spatial power spectrum of the signals. Estimated 
DOAs are the grid positions for which the spectrum crosses a certain threshold, as shown in the 'BPDN.ipynb' example.    

arguments:
----------
    Y: Data matrix of the array
    A: Dictionary matrix of array response vectors from the angle grid 
    η: Constraint and upper bound on the noise energy (η≥||E||_2). 

References:
-----------
Z. Yang, J. Li, P. Stoica, and L. Xie, ‘Sparse methods for direction-of-arrival estimation’, arXiv [cs.IT], 30-Sep-2016.
"""
function bpdn(Y, A, η=1e-2)
    X = ComplexVariable(size(A)[2], size(Y)[2])
    p = minimize(sum([norm(X[i, :], 2) for i in axes(X, 1)]), norm(A*X-Y, 2) <= η)
    Convex.solve!(p, SCS.Optimizer)
    return norm.(eachrow(evaluate(X)),2).^2
end