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

    if M == 1
        padded_coords = vcat(coords, Zeros{eltype(coords)}(2, D))
    elseif M == 2
        padded_coords = vcat(coords, Zeros{eltype(coords)}(1, D))
    else
        padded_coords = coords
    end

    return WaveVec{eltype(padded_coords), typeof(padded_coords)}(padded_coords)
end

function WaveVec(kx_vec::AbstractVector)
    return WaveVec(reshape(kx_vec, 1, :))
end

function WaveVec(kx...)
    return WaveVec(collect(promote(kx...)))
end