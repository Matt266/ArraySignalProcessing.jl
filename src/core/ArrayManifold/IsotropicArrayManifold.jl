struct IsotropicArrayManifold{T, A<:AbstractMatrix{T}} <: AbstractArrayManifold
    r::A
    IsotropicArrayManifold{T, A}(r::A) where {T, A<:AbstractMatrix{T}} = new{T, A}(r)
end

Adapt.@adapt_structure IsotropicArrayManifold

function IsotropicArrayManifold(r::AbstractMatrix)
    M, D = size(r)

    M <= 0 && throw(DimensionMismatch("'r' has size $(size(r)) but must have at least one row."))
    M >= 4 && throw(DimensionMismatch("'r' has size $(size(r)) but must have at most three rows: [x...; y...; z...]"))
    D <= 0 && throw(DimensionMismatch("'r' has size $(size(r)) but must have at least one column"))

    if M == 1
        padded_r = vcat(r, Zeros{eltype(r)}(2, D))
    elseif M == 2
        padded_r = vcat(r, Zeros{eltype(r)}(1, D))
    else
        padded_r = r
    end

    return IsotropicArrayManifold{eltype(padded_r), typeof(padded_r)}(padded_r)
end

function IsotropicArrayManifold(r::AbstractVector)
    return IsotropicArrayManifold(reshape(r, 1, :))
end

function IsotropicArrayManifold(elements...)
    return IsotropicArrayManifold(collect(promote(elements...)))
end

function Base.length(a::IsotropicArrayManifold)
    return size(a.r, 2)
end


"""
References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""

function (a::IsotropicArrayManifold)(angles::AzEl, f, c=c_0)
    f_res = f isa Number ? f : reshape(f, 1, 1, :, 1)
    c_res = c isa Number ? c : reshape(c, 1, 1, 1, :)

    k = convert.(promote_type(eltype(angles.coords), Float32), 2π .* f_res ./ c_res)

    az = angles.coords[1:1, :]
    el = angles.coords[2:2, :]

    ζ = vcat(cos.(el) .* cos.(az),
            cos.(el) .* sin.(az),
            sin.(el))

    φ = -(k .* (a.r' * ζ))
    A_tensor = exp.(-1im .* φ)

    return reshape(A_tensor, size(A_tensor, 1), :)
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
    A_base = exp.(-1im .* φ)
    
    # pad dimensions to stay consistent with wavefront types on return
    T = promote_type(eltype(angles.coords), Float32)
    f_zeros = f isa Number ? 0 : fill!(similar(angles.coords, T, 1, 1, length(f), 1), 0)
    c_zeros = c isa Number ? 0 : fill!(similar(angles.coords, T, 1, 1, 1, length(c)), 0)
    A_tensor = A_base .+ f_zeros .+ c_zeros
    return reshape(A_tensor, size(A_tensor, 1), :)
end

"""
References:
-----------
D. H. Johnson and D. E. Dudgeon, Array Signal Processing. Philadelphia, PA: Prentice Hall, 1993.
"""
function (a::IsotropicArrayManifold)(angles::SlowVec, f, c=c_0)
    f_res = f isa Number ? f : reshape(f, 1, 1, :, 1)

    ω = convert.(promote_type(eltype(angles.coords), Float32), 2π .* f_res)
    φ = a.r' * (ω * angles.coords)
    A_base = exp.(-1im .* φ)

    # pad dimensions to stay consistent with wavefront types on return
    T = promote_type(eltype(angles.coords), Float32)
    c_zeros = c isa Number ? 0 : fill!(similar(angles.coords, T, 1, 1, 1, length(c)), 0)
    A_tensor = A_base .+ c_zeros
    return reshape(A_tensor, size(A_tensor, 1), :)
end

"""
References:
-----------
Z. Ebadi, A. M. Molaei, M. A. B. Abbasi, S. Cotton, A. Tukmanov and O. Yurduseven, "Near-Field Localization with an Exact Propagation Model in Presence of Mutual Coupling," 2024 IEEE 99th Vehicular Technology Conference (VTC2024-Spring), Singapore, Singapore, 2024, pp. 1-5, doi: 10.1109/VTC2024-Spring62846.2024.10683010.
"""

function (a::IsotropicArrayManifold)(angles::RAzEl, f, c=c_0)
    f_res = f isa Number ? f : reshape(f, 1, 1, :, 1)
    c_res = c isa Number ? c : reshape(c, 1, 1, 1, :)

    k = convert.(promote_type(eltype(angles.coords), Float32), 2π .* f_res ./ c_res)

    r  = angles.coords[1:1, :]
    az = angles.coords[2:2, :]
    el = angles.coords[3:3, :]

    ζ = r .* vcat(cos.(el) .* cos.(az),
              cos.(el) .* sin.(az),
              sin.(el))

    #(3xD → 3x1xD)
    src_pos = reshape(ζ, 3, 1, size(ζ,2))

    #(3xM → 3xMx1)
    sens_pos = reshape(a.r, 3, size(a.r,2), 1)

    # distances from each sensor to each source (M×D)
    R = dropdims(sqrt.(sum(abs2, sens_pos .- src_pos; dims=1)), dims=1)

    amp = (r ./ R)
    φ = (k .* (R .- r))

    A_tensor = amp .* exp.(-1im .* φ)
    return reshape(A_tensor, size(A_tensor, 1), :)
end