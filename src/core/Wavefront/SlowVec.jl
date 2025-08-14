struct SlowVec{T<:AbstractMatrix} <: PlaneWave
    coords::T
    function SlowVec(coords::AbstractMatrix)
        # coords = [sx...; sy...; sz...] 3xD matrix

        M, D = size(coords)
        M <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one row."))
        M >= 4 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at most three rows: [sx...; sy...; sz...]"))
        D <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one column"))

        padded_coords = similar(coords, 3, D)
        padded_coords[1:M, :] = coords

        if M == 1
            padded_coords = vcat(coords, zeros(eltype(coords), 2, D))
        elseif M==2
            padded_coords = vcat(coords, zeros(eltype(coords), 1, D))
        else
            padded_coords = vcat(coords, zeros(eltype(coords), 0, D))
        end

        return new{typeof(padded_coords)}(padded_coords)
    end
end

# list of sx
function SlowVec(sx_vec::AbstractVector)
    return SlowVec(reshape(sx_vec, 1, :))
end

# list of sx
function SlowVec(sx...)
    return SlowVec(collect(promote(sx...)))
end