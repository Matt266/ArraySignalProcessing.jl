struct SlowVec{T, A<:AbstractMatrix{T}} <: PlaneWave
    coords::A
    SlowVec{T, A}(coords::A) where {T, A<:AbstractMatrix{T}} = new{T, A}(coords)
end

Adapt.@adapt_structure SlowVec

function SlowVec(coords::AbstractMatrix)
    M, D = size(coords)
    M <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one row."))
    M >= 4 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at most three rows: [sx...; sy...; sz...]"))
    D <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one column"))

    if M == 1
        padded_coords = vcat(coords, Zeros{eltype(coords)}(2, D))
    elseif M == 2
        padded_coords = vcat(coords, Zeros{eltype(coords)}(1, D))
    else
        padded_coords = coords
    end

    return SlowVec{eltype(padded_coords), typeof(padded_coords)}(padded_coords)
end

function SlowVec(sx_vec::AbstractVector)
    return SlowVec(reshape(sx_vec, 1, :))
end

function SlowVec(sx...)
    return SlowVec(collect(promote(sx...)))
end