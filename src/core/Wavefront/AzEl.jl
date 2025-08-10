struct AzEl{T<:AbstractMatrix} <: PlaneWave
    coords::T
    function AzEl(coords::AbstractMatrix)
        # coords = [azimuths...; elevations...] 2xD matrix
        M, D = size(coords)

        M <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one row."))
        M >= 3 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at most two rows: [azimuths...; elevations...]"))
        D <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one column"))

        padded_coords = similar(coords, 2, D)
        padded_coords[1:M, :] = coords

        if M == 1
            fill!(@view(padded_coords[2, :]), zero(eltype(padded_coords)))
        end

        return new{typeof(padded_coords)}(padded_coords)
    end
end
# list of azimuths
function AzEl(azimuths::AbstractVector)
    return AzEl(reshape(azimuths, 1, :))
end

# list of azimuths
function AzEl(azimuths...)
    return AzEl(collect(promote(azimuths...)))
end