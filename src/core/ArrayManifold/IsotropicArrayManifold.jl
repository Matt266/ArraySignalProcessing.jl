struct IsotropicArrayManifold{T<:AbstractMatrix} <: AbstractArrayManifold
    r::T
    function IsotropicArrayManifold(r::AbstractMatrix)
        M, D = size(r)

        M <= 0 && throw(DimensionMismatch("'r' has size $(size(r)) but must have at least one row."))
        M >= 4 && throw(DimensionMismatch("'r' has size $(size(r)) but must have at most three rows: [x...; y...; z...]"))
        D <= 0 && throw(DimensionMismatch("'r' has size $(size(r)) but must have at least one column"))

        padded_r = similar(r, 3, D)
        padded_r[1:M, :] = r

        if M == 1
            padded_r = vcat(r, zeros(eltype(r), 2, D))
        elseif M == 2
            padded_r = vcat(r, zeros(eltype(r), 1, D))
        else
            padded_r = vcat(r, zeros(eltype(r), 0, D))
        end

        return new{typeof(padded_r)}(padded_r) 
    end
end

function IsotropicArrayManifold(r::AbstractVector)
    return IsotropicArrayManifold(reshape(r, 1, :))
end

function IsotropicArrayManifold(elements...)
    return IsotropicArrayManifold(collect(promote(elements...)))
end

function Base.length(a::IsotropicArrayManifold)
    return size(a.r)[2]
end


"""
References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""

function (a::IsotropicArrayManifold)(angles::AzEl, f, c=c_0)
    k = convert(eltype(angles.coords), 2π * f / c)

    az = transpose(angles.coords[1, :])
    el = transpose(angles.coords[2, :])

    ζ = [cos.(el) .* cos.(az);
         cos.(el) .* sin.(az);
         sin.(el)]

    φ = -(k .* (a.r' * ζ))

    return exp.(-1im .* φ)
end

"""
References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function (a::IsotropicArrayManifold)(angles::WaveVec, f, c=c_0)
    #k = 2π * f / c  
    #ζ = angles.coords ./ k
    #φ = k .* (a.r' * ζ)
    φ = a.r' * angles.coords
    return exp.(-1im .* φ)
end

"""
References:
-----------
D. H. Johnson and D. E. Dudgeon, Array Signal Processing. Philadelphia, PA: Prentice Hall, 1993.
"""
function (a::IsotropicArrayManifold)(angles::SlowVec, f, c=c_0)
    ω = convert(eltype(angles.coords), 2π * f)
    φ = a.r' * (ω * angles.coords)
    return exp.(-1im .* φ)
end

"""
References:
-----------
Z. Ebadi, A. M. Molaei, M. A. B. Abbasi, S. Cotton, A. Tukmanov and O. Yurduseven, "Near-Field Localization with an Exact Propagation Model in Presence of Mutual Coupling," 2024 IEEE 99th Vehicular Technology Conference (VTC2024-Spring), Singapore, Singapore, 2024, pp. 1-5, doi: 10.1109/VTC2024-Spring62846.2024.10683010.
"""

function (a::IsotropicArrayManifold)(angles::RAzEl, f, c=c_0)
    k = convert(eltype(angles.coords), 2π * f / c)

    r  = transpose(angles.coords[1, :])
    az = transpose(angles.coords[2, :])
    el = transpose(angles.coords[3, :])

    ζ = r .* [cos.(el) .* cos.(az);
              cos.(el) .* sin.(az);
              sin.(el)]

    #(3xD → 3x1xD)
    src_pos = reshape(ζ, 3, 1, size(ζ,2))

    #(3xM → 3xMx1)
    sens_pos = reshape(a.r, 3, size(a.r,2), 1)

    # distances from each sensor to each source (M×D)
    R = dropdims(sqrt.(sum(abs2, sens_pos .- src_pos; dims=1)), dims=1)

    amp = (r ./ R)
    φ = -(k .* (R .- r))

    return amp .* exp.(-1im .* φ)
end