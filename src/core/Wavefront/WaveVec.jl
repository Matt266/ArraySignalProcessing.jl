struct WaveVec{T} <: PlaneWave where T <: AbstractMatrix
    coords::T
    function WaveVec(coords::T) where T <: AbstractMatrix
        # coords = [kx; ky; kz] 3xD matrix

        M, D = size(coords)
        M <= 0 && throw(DomainError(M, "'coords' must have at least one row."))
        M >= 4 && throw(DomainError(M, "'coords' must have at most three rows: [kx...; ky...; kz...]"))
        D <= 0 && throw(DomainError(D, "'coords' must have at least one column"))

        if M == 1
            padding = similar(coords, 2, D)
            fill!(padding, zero(eltype(coords)))
            coords = vcat(coords, padding)
        elseif M==2
            padding = similar(coords, 1, D)
            fill!(padding, zero(eltype(coords)))
            coords = vcat(coords, padding)
        end

        return new{typeof(coords)}(coords)
    end
end

# list of kx
function WaveVec(kx_vec::T) where T <: AbstractVector
    return WaveVec(reshape(kx_vec, 1, :))
end

# list of kx
function WaveVec(kx::T...) where T <: Number
    return WaveVec(collect(kx))
end