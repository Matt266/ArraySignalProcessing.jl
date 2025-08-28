struct TappedDelayLineManifold{MT, NT, FT} <: AbstractArrayManifold where {MT<:AbstractArrayManifold, NT<:Int, FT<:Number}
    sub_manifold::MT
    num_taps::NT
    fs::FT
end

function Base.length(a::TappedDelayLineManifold)
    return length(a.sub_manifold)*a.num_taps
end

"""
Nested structure: M elements and J taps. Each tap has it's own weight vector. Those 
are stacked on top of each other to one large weight vector for the whole TappedDelayLine.

Columns than correspond to angles. 

Result is a large matrix: (M*J) x (A).

Basically does: taps_phaseshifts_as_vector ⊗ steering_vector
But is adjusted for multiple angles.

References:
-----------

W. Liu and S. Weiss, Wideband Beamforming. Nashville, TN: John Wiley & Sons, 2010.
"""
function (a::TappedDelayLineManifold)(angles::Wavefront, f::Number, c::Number=c_0)
    # array manifold matrix: M x A
    v_array = a.sub_manifold(angles, f, c)

    M = size(v_array, 1)
    A = length(angles)

    # number of taps and sample period
    J = a.num_taps
    Ts = convert(promote_type(eltype(f), Float32), 1 / a.fs)
    
    # J x 1
    v_taps = similar(v_array, complex(promote_type(eltype(f), Float32)), J)
    copyto!(v_taps, hcat(exp.(-1im .* 2π .* f .* (0:J-1) .* Ts)))

    # reshape to broadcast
    v_array_reshaped = reshape(v_array, M, 1, A)   # M x 1 x A
    v_taps_reshaped = reshape(v_taps, 1, J, 1)     # 1 x J x 1

    
    # broadcast outer product over taps
    v_full = v_array_reshaped .* v_taps_reshaped   # M × J × A

    # reshape to: (M*J) × A
    v = reshape(v_full, M*J, A)

    return v
end