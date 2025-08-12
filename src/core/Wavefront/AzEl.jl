struct AzEl{T<:AbstractMatrix} <: PlaneWave
    coords::T
    function AzEl(coords::AbstractMatrix)
        # coords = [azimuths...; elevations...] 2xD matrix
        M, D = size(coords)

        M <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one row."))
        M >= 3 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at most two rows: [azimuths...; elevations...]"))
        D <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one column"))

        if M == 1
            padded_coords = vcat(coords, zeros(eltype(coords), 1, D))
        else
            padded_coords = vcat(coords, zeros(eltype(coords), 0, D))
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