struct SlowVec{T, A<:AbstractMatrix{T}} <: PlaneWave
    coords::A
end

Adapt.@adapt_structure SlowVec

function SlowVec(coords::AbstractMatrix)
    M, D = size(coords)
    M <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one row."))
    M >= 4 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at most three rows: [sx...; sy...; sz...]"))
    D <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one column"))

    padded_coords = similar(coords, eltype(coords), 3, D)
    padded_coords[1:M, :] = coords

    if M < 3
        padded_coords[M+1:3, :] .= 0
    end

    return SlowVec{eltype(padded_coords), typeof(padded_coords)}(padded_coords)
end

function SlowVec(sx_vec::AbstractVector)
    return SlowVec(reshape(sx_vec, 1, :))
end

function SlowVec(sx...)
    return SlowVec(collect(promote(sx...)))
end