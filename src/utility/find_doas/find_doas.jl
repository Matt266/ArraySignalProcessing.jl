"Finds peaks in the 1D gridded spectrum"
function find_doas(grid::AbstractVector, spectrum::AbstractVector, d)
    # all peaks
    peaks = Peaks.findmaxima(spectrum)
    pidx = peaks.indices
    pval = peaks.heights 

    # all peaks -> up to d highest peaks -> corresponding grid values
    return grid[pidx[sortperm(pval, rev=true)[1:min(d, end)]]]
end

"Finds peaks in the 1D power spectrum (spectrum = pfunc(angles)) via initial grid search
 and a following optimization based refinment"
function find_doas(grid::AbstractVector, pfunc::Function, d;
                    optimizer = Fminbox(LBFGS()), problem_kwargs=(), solve_kwargs=())
    # spectrum over grid
    spectrum = pfunc(transpose(grid))

    # all grid peaks
    peaks = Peaks.findmaxima(spectrum)
    pidx = peaks.indices
    pval = peaks.heights 

    # all peak indices -> indices of up to d highest peaks
    rough_indices = pidx[sortperm(pval, rev=true)[1:min(d, end)]]

    # invert to search maxima with cost minimization
    # and turn one-element output array (spectrum of single angle)
    # into scalar cost output for optimization
    cost_func(angle, args...) = -pfunc(angle)[1]

    refined_doas = []
    for idx in rough_indices
        f = OptimizationFunction(cost_func, AutoZygote())
        # search for peak only within the region to the adjacent grid points
        p = OptimizationProblem(f, [grid[idx]], (); lb=[grid[max(idx-1, 1)]],ub=[grid[min(idx+1, end)]], problem_kwargs...)
        s = solve(p, optimizer; solve_kwargs...)
        append!(refined_doas, s.u)
    end
    return refined_doas
end