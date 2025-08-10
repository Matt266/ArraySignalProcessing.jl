struct IsotropicArrayManifold <: AbstractArrayManifold
    r::AbstractMatrix
    function IsotropicArrayManifold(r::AbstractMatrix)
        @assert 1 <= size(r)[1] <= 3 "Found invalid shape of input matrix: $(size(r)[1])x$(size(r)[2]). Each column are the 1D, 2D, or 3D coordinates of an element. Input matrix must be of size 1xN, 2xN, or 3xN respectively"

        if size(r)[1] == 1
            return new([r; zeros(2, size(r)[2])])
        elseif size(r)[1] == 2
            return new([r; zeros(1, size(r)[2])])
        else
            return new(r)
        end 
    end
end

function IsotropicArrayManifold(r::AbstractVector)
    return IsotropicArrayManifold(reshape(r, 1, length(r)))
end

IsotropicArrayManifold(elements::Number...) = IsotropicArrayManifold([e for e in elements])

function Base.length(a::IsotropicArrayManifold)
    return size(a.r)[2]
end


"""
References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""

function (a::IsotropicArrayManifold)(f, angles::AzEl; c=c_0)
    k = 2π * f / c

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
function (a::IsotropicArrayManifold)(f, angles::WaveVec; c=c_0)
    #k = 2π * f / c  
    #ζ = angles.coords ./ k
    #φ = k .* (a.r' * ζ)
    φ = a.r' * angles.coords
    return exp.(-1im .* φ)
end