struct AzEl{T} <: PlaneWave where T <: AbstractMatrix
    coords::T
    function AzEl(coords::T) where T <: AbstractMatrix
        # coords = [azimuths...; elevations...] 2xD matrix
        M, D = size(coords)

        M <= 0 && throw(DomainError(M, "'coords' must have at least one row."))
        M >= 3 && throw(DomainError(M, "'coords' must have at most two rows: [azimuths...; elevations...]"))
        D <= 0 && throw(DomainError(D, "'coords' must have at least one column"))

        if M == 1
            padding = similar(coords, 1, D)
            fill!(padding, zero(eltype(coords)))
            coords = vcat(coords, padding)
        end

        return new{typeof(coords)}(coords)
    end
end

# list of azimtuhs
function AzEl(azimuths::T) where T <: AbstractVector
    return AzEl(reshape(azimuths, 1, :))
end

# list of azimtuhs
function AzEl(azimuths::T...) where T <: Number
    return AzEl(collect(azimuths))
end