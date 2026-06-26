"""
Creates an ArrayManifold from recorded samples.
'responses' array holding the complex amplitudes should have following dimension order:
    propagation_speeds - frequencies - spatial coordinates - array elements

E.g., for AzEl this would be (check implemented Wavefront types for supported formats):
    propagation_speeds - frequencies - azimuths - elevations - array elements

Note that each dimension of size 1 (or less) can be dropped. So, if only azimuths for a
single frequency are recorded, responses reduces to:
    azimuths - array elements

You must sample over at least one dimension. If other dependencies (e.g., temperature) are 
be needed, either create an array of SampledArrayManifolds or implement a custom ArrayManifold subtype.
"""
struct SampledArrayManifold{NT, CT, KA, MT, PT} <: AbstractArrayManifold
    num_elements::NT
    coords_type::CT
    keep_axes::KA
    itp_mag::MT
    itp_phase::PT
end

Adapt.@adapt_structure SampledArrayManifold

function SampledArrayManifold(responses::AbstractArray; c_grid=[0], f_grid=[0], coords_grid=AzEl(0))
    num_elements = size(responses)[end]

    # sort all axes
    coords_perms = sortperm.(unique.(eachrow(coords_grid.coords)))
    coords_axes = map(i -> unique(eachrow(coords_grid.coords)[i])[coords_perms[i]], 1:length(coords_perms))

    c_perms = sortperm(c_grid)
    c_grid = c_grid[c_perms]

    f_perms = sortperm(f_grid)
    f_grid = f_grid[f_perms]

    element_perms = Vector(1:num_elements)
    element_axes = Vector(1:num_elements)

    # tuple with all axes of length >1
    grid_perms = [c_perms, f_perms, coords_perms..., element_perms]
    grid_axes = [c_grid, f_grid, coords_axes..., element_axes]
    keep_axes = length.(grid_axes) .> 1
    grid_perms = (grid_perms[keep_axes]...,)
    grid_axes = (grid_axes[keep_axes]...,)

    # responses with all axes of length >1
    flattened_responses = dropdims(responses; dims=Tuple(findall(size(responses) .<= 1)))

    # permute responsens to match sorted axes again
    flattened_responses = flattened_responses[grid_perms...]
    
    itp_mag = extrapolate(interpolate(grid_axes, abs.(flattened_responses), Gridded(Interpolations.Linear())),Interpolations.Flat())
    itp_phase = extrapolate(interpolate(grid_axes, angle.(flattened_responses), Gridded(Interpolations.Linear())),Interpolations.Flat())
    return SampledArrayManifold(num_elements, typeof(coords_grid), keep_axes, itp_mag, itp_phase)
end

# TODO: implement Wavefront type conversion logic.
# Currently will fail if called with other format than the sampled one,
# which for now is the intended behaviour 
# TODO: probably in each Wavefront constructer: constrain values (e.g., azimuth to -π...+π and elevation to 0...π).
function (a::SampledArrayManifold)(angles::Wavefront, f::Number, c::Number=c_0)
    angles_conv = convert(a.coords_type, angles)

    f_vec = f isa Number ? (f,) : f
    c_vec = c isa Number ? (c,) : c

    combinations = Iterators.product(eachcol(angles_conv.coords), f_vec, c_vec)

    # matrix of tuples: each row an element index e, each col an column from the angles 
    # all axes with size <= 1 will be dropped from the points
    A_matrix = [ 
        begin
            query_point = Tuple([c_val, f_val, col..., e][a.keep_axes])
            a.itp_mag(query_point...) * exp(1im * a.itp_phase(query_point...))
        end
        for e in 1:a.num_elements, (col, f_val, c_val) in combinations 
    ]
    
    return A_matrix
end

function Base.length(a::SampledArrayManifold)
    return a.num_elements
end
