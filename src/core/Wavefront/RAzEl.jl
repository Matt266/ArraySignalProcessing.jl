struct RAzEl{T, A<:AbstractMatrix{T}} <: SphericalWave
    coords::A
end

Adapt.@adapt_structure RAzEl

function RAzEl(coords::AbstractMatrix)
    M, D = size(coords)
    M <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one row."))
    M >= 4 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at most three rows:  [range...; azimuth...; elevation...]"))
    D <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one column"))

    padded_coords = similar(coords, eltype(coords), 3, D)
    padded_coords[1:M, :] = coords

    if M < 3
        padded_coords[M+1:3, :] .= 0
    end

    return RAzEl{eltype(padded_coords), typeof(padded_coords)}(padded_coords)
end

function RAzEl(r_vec::AbstractVector)
    return RAzEl(reshape(r_vec, 1, :))
end

function RAzEl(r...)
    return RAzEl(collect(promote(r...)))
end