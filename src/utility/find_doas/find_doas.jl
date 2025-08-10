function filter_plateaus(indices)
    N = length(indices[1])
    filtered_indices = []
    indices = Set(indices)

    while !isempty(indices)
        start_idx = pop!(indices)

        # breadth-first search (BFS) to find all connected neighbors
        queue = [start_idx]
        plateau_group = [start_idx]
        min_norm_idx = start_idx
        min_norm = norm(Tuple(start_idx))
        while !isempty(queue)
            current_idx = popfirst!(queue)
            for neighbor_offset in CartesianIndices(ntuple(_ -> -1:1, N))
                # Skip the current index itself
                if all(Tuple(neighbor_offset) .== 0)
                    continue
                end

                # Check if the neighbor is in the set of unvisited indices
                neighbor_idx = current_idx + neighbor_offset
                if neighbor_idx in indices
                    push!(plateau_group, neighbor_idx)
                    push!(queue, neighbor_idx)
                    delete!(indices, neighbor_idx)

                    # Update the index with the minimum norm
                    current_norm = norm(Tuple(neighbor_idx))
                    if current_norm < min_norm
                        min_norm = current_norm
                        min_norm_idx = neighbor_idx
                    end
                end
            end
        end
        # Add the minimum-norm index from the plateau to the filtered list
        push!(filtered_indices, min_norm_idx)
    end

    return filtered_indices
end

function gridsearch(func::Function, grids::Tuple)
    grid_indices = CartesianIndices(map(grid -> 1:length(grid), grids))
    grid_data = [func(map( (i, g) -> g[i], Tuple(idx), grids)...) for idx in grid_indices]

    raw_peak_indices = []
    for idx in grid_indices
        is_local_max = true
        current_value = grid_data[idx]

        for neighbor_offset in CartesianIndices(ntuple(_ -> -1:1, length(grids)))
            # Skip the current index itself
            if all(Tuple(neighbor_offset) .== 0)
                continue
            end

            neighbor_idx = idx + neighbor_offset

            # Ensure neighbor is within grid bounds
            if checkbounds(Bool, grid_data, neighbor_idx)
                if grid_data[neighbor_idx] > current_value
                    is_local_max = false
                    break
                end
            end
        end

        if is_local_max
            push!(raw_peak_indices, idx)
        end
    end

    # Filter the raw indices to handle plateaus
    filtered_indices = filter_plateaus(raw_peak_indices)

    # Convert filtered indices back to coordinates and return
    peak_coords = [map((i, g) -> g[i], Tuple(idx), grids) for idx in filtered_indices]

    return filtered_indices, peak_coords
end

function refine_peaks(func::Function, peak_indices, grids)

    refined_peak_coords = []
    obj_func(x) = -func(x...)

    for current_index in peak_indices
        current_coords = map((idx, grid) -> grid[idx], Tuple(current_index), grids)

        lower_bounds = map(
            (j, grid) -> grid[max(1, current_index[j] - 1)], 
            1:length(grids), grids
        )
        upper_bounds = map(
            (j, grid) -> grid[min(length(grid), current_index[j] + 1)], 
            1:length(grids), grids
        )

        is_refinable = all(map((l, u) -> l != u, lower_bounds, upper_bounds))

        if is_refinable
            initial_guess = collect(current_coords)
            rhobeg = minimum(upper_bounds - lower_bounds) / 4.0
            result, _ = prima(obj_func, initial_guess; xl=lower_bounds, xu=upper_bounds, rhobeg=rhobeg)
            push!(refined_peak_coords, result)
        else 
            push!(refined_peak_coords, collect(current_coords))
        end
    end

    return refined_peak_coords
end

function merge_close_peaks(func::Function, peak_coords, merge_distance)
    merged_peak_coords = []

    for peak_coord in peak_coords
        peak_val = func(peak_coord...)

        is_merged = false
        for (i, (unique_peak, unique_val)) in enumerate(merged_peak_coords)
            # Check if the distance is within the tolerance
            if norm(peak_coord - unique_peak) < merge_distance
                is_merged = true
                # Keep the peak with the higher function value
                if peak_val > unique_val
                    merged_peak_coords[i] = (peak_coord, peak_val)
                end
                break
            end
        end

        if !is_merged
            push!(merged_peak_coords, (peak_coord, peak_val))
        end
    end

    return [p[1] for p in merged_peak_coords]
end

"""
Automates spectral search for maxima corresponding to DoAs like encountered in MUSIC.
Define a function func(coords...) that returns the power spectrum (or any real valued spectrum)
for the specified coordinates. The d highest maxima will be returned. You must give 
an initial grid for each coordinate axes that is fine enough to resolve each peak. 
"""
function find_doas(func::Function, d::Int, grids...; merge_distance = 0.01)




    function select_peaks(func::Function, d::Int, peak_coords)
        # Store peak values and their coordinates as tuples
        peak_info = [(func(p...), p) for p in peak_coords]

        # Sort the list in descending order based on the peak value (the first element of the tuple)
        sort!(peak_info, by = first, rev = true)
        return [p[2] for p in peak_info[1:min(d, end)]]
    end

    initial_peak_indices, _ = gridsearch(func, grids)

    if isempty(initial_peak_indices)
        return Matrix{Float64}(undef, length(grids), 0)
    end

    refined_peak_coords = refine_peaks(func, initial_peak_indices, grids)

    merged_peaks = merge_close_peaks(func, refined_peak_coords, merge_distance)

    top_peaks = select_peaks(func, d, merged_peaks)
    sort!(top_peaks, by = p -> Tuple(p))
    return reduce(hcat, top_peaks)
end