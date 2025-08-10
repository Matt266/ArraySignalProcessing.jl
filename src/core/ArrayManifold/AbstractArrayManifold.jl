abstract type AbstractArrayManifold end

"""
steer(x::IsotropicArray, f, angles; c=c_0, coords=:azel)

Calculates the steering vector/matrix for given azimuth and elevation angles. 

arguments:
----------
    f: Frequency of the signal in Hz.
    angles: Matrix with angles in the format specified in 'coords'
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
function (a::AbstractArrayManifold)(f, angles; c=c_0, coords=:azel)
    if coords == :azel
        return a(f, AzEl(angles); c=c)
    elseif coords == :k
        return a(f, WaveVec(angles); c=c)
    else
        throw(ArgumentError("coords must be ':azel' or ':k'; got: '$(mode)'"))
    end
end