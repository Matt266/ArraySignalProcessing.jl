struct WaveVec{T, A<:AbstractMatrix{T}} <: PlaneWave
    coords::A
    WaveVec{T, A}(coords::A) where {T, A<:AbstractMatrix{T}} = new{T, A}(coords)
end

Adapt.@adapt_structure WaveVec

function WaveVec(coords::AbstractMatrix)
    M, D = size(coords)
    M <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one row."))
    M >= 4 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at most three rows: [kx...; ky...; kz...]"))
    D <= 0 && throw(DimensionMismatch("'coords' has size $(size(coords)) but must have at least one column"))

    padded_coords = similar(coords, eltype(coords), 3, D)
    padded_coords[1:M, :] = coords

    if M < 3
        padded_coords[M+1:3, :] .= 0
    end

    return WaveVec{eltype(padded_coords), typeof(padded_coords)}(padded_coords)
end

function WaveVec(kx_vec::AbstractVector)
    return WaveVec(reshape(kx_vec, 1, :))
end

function WaveVec(kx...)
    return WaveVec(collect(promote(kx...)))
end