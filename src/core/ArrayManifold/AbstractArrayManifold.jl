abstract type AbstractArrayManifold end

"""
steer(x::IsotropicArray, angles, f, c=c_0; coords=:azel)

Calculates the steering vector/matrix for given azimuth and elevation angles. 

arguments:
----------
    angles: Matrix with angles in the format specified in 'coords'
    f: Frequency of the signal in Hz.
    c: Propagation speed of the wave (default: c_0).
    coords: Coordinate system of 'angles':
         - ':azel' (default): interpret angles as azimuth/elevation pairs [az; el]
         - ':k': interpret as wavevector [kx; ky; kz]
returns:
--------
    Complex steering matrix of size MxD, where M is the number of array elements and D is the number of angle pairs.

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function (a::AbstractArrayManifold)(angles, f, c=c_0; coords=:AzEl)
    if coords == :AzEl
        return a(AzEl(angles), f, c)
    elseif coords == :WaveVec
        return a(WaveVec(angles), f, c)
    elseif  coords == :SlowVec
        return a(SlowVec(angles), f, c)
    elseif  coords == :RAzEl
        return a(RAzEl(angles), f, c)
    else
        throw(DomainError("'coords' must be ':AzEl', ':WaveVec', ':SlowVec', or ''RAzEl'; got: '$(coords)'"))
    end
end

# fallback so default value c=c_0 is used when possible
function (a::AbstractArrayManifold)(angles::Wavefront, f, c=c_0)
    return a(angles, f, c)
end

# returns M x (A*F*C) matrix which can be reshaped to M x A x F x C if required
function (a::AbstractArrayManifold)(angles::Wavefront, f::AbstractVector, c::AbstractVector)
    #return mapreduce(((fi, ci),) -> a(angles, fi, ci), hcat, Iterators.product(f, c))
    a_list = map(((fi, ci),) -> a(angles, fi, ci), Iterators.product(f, c))
    return hcat(a_list...)
end

# returns M x (A*F) matrix which can be reshaped to M x A x F if required
function (a::AbstractArrayManifold)(angles::Wavefront, f::AbstractVector, c::Number)
    #return mapreduce(fi -> a(angles, fi, c), hcat, f)
    a_list = map(fi -> a(angles, fi, c), f)
    return hcat(a_list...)
end

# returns M x (A*C) matrix which can be reshaped to M x A x C if required
function (a::AbstractArrayManifold)(angles::Wavefront, f::Number, c::AbstractVector)
    #return mapreduce(ci -> a(angles, f, ci), hcat, c)
    a_list = map(ci -> a(angles, f, ci), c)
    return hcat(a_list...)
end