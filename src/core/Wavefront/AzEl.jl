struct AzEl{T, A<:AbstractMatrix{T}} <: PlaneWave
    coords::A
    AzEl{T, A}(coords::A) where {T, A<:AbstractMatrix{T}} = new{T, A}(coords)
end

Adapt.@adapt_structure AzEl

function AzEl(coords::AbstractMatrix)
    M, D = size(coords)

    M <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one row."))
    M >= 3 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at most two rows: [azimuths...; elevations...]"))
    D <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one column"))

    padded_coords = similar(coords, eltype(coords), 2, D)
    padded_coords[1:M, :] = coords
    
    if M == 1
        padded_coords[2:2, :] .= 0
    end

    return AzEl{eltype(padded_coords), typeof(padded_coords)}(padded_coords)
end

function AzEl(azimuths::AbstractVector)
    return AzEl(reshape(azimuths, 1, :))
end

function AzEl(azimuths...)
    return AzEl(collect(promote(azimuths...)))
end