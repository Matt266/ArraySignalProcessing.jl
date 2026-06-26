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
function (a::AbstractArrayManifold)(angles::AbstractMatrix, f, c=c_0; coords=AzEl)
        return a(coords(angles), f, c)
end

function (a::AbstractArrayManifold)(angles::Wavefront, f, c=c_0)
    return a(angles, f, c)
end