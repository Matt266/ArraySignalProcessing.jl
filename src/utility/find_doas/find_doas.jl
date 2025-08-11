"Finds peaks in the 1D gridded spectrum"
function find_doas(grid, spectrum, d)
    # all peaks
    peaks = Peaks.findmaxima(spectrum)
    pidx = peaks.indices
    pval = peaks.heights 

    # all peaks -> up to d highest peaks -> corresponding grid values
    return grid[pidx[sortperm(pval, rev=true)[1:min(d, end)]]]
end