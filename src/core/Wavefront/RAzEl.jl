struct RAzEl{T, A<:AbstractMatrix{T}} <: SphericalWave
    coords::A
    RAzEl{T, A}(coords::A) where {T, A<:AbstractMatrix{T}} = new{T, A}(coords)
end

Adapt.@adapt_structure RAzEl

function RAzEl(coords::AbstractMatrix)
    M, D = size(coords)
    M <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one row."))
    M >= 4 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at most three rows:  [range...; azimuth...; elevation...]"))
    D <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one column"))

    if M == 1
        padded_coords = vcat(coords, Zeros{eltype(coords)}(2, D))
    elseif M == 2
        padded_coords = vcat(coords, Zeros{eltype(coords)}(1, D))
    else
        padded_coords = coords
    end

    return RAzEl{eltype(padded_coords), typeof(padded_coords)}(padded_coords)
end

function RAzEl(r_vec::AbstractVector)
    return RAzEl(reshape(r_vec, 1, :))
end

function RAzEl(r...)
    return RAzEl(collect(promote(r...)))
end