abstract type AbstractPhasedArray end

struct PhasedArray{T<:AbstractArrayManifold} <: AbstractPhasedArray
    manifold::T
end

struct NestedArray{E<:AbstractPhasedArray, S<:AbstractVector{<:AbstractPhasedArray}} <: AbstractPhasedArray
    elements::E
    subarrays::S
    function NestedArray(elements::E, subarrays::S) where {E<:AbstractPhasedArray, S<:AbstractVector{<:AbstractPhasedArray}}
        (length(elements) == length(subarrays)) || throw(DimensionMismatch("Length of elements ($(length(elements))) and number of subarrays ($(length(subarrays))) do't match)"))
        return new{E, S}(elements, subarrays)
    end
end

function Base.length(pa::PhasedArray)
    return length(pa.manifold)
end

function Base.length(pa::NestedArray)
    return sum(length.(pa.subarrays))
end

function steer(pa::PhasedArray, angles, f, c=c_0; kwargs...)
    return pa.manifold(angles, f, c; kwargs...)
end

function steer(pa::NestedArray, angles, f, c=c_0; kwargs...)
    v_super = steer(pa.elements, angles, f, c; kwargs...)
    v_sub = map(sa -> steer(sa, angles, f, c; kwargs...), pa.subarrays)

    v_sub_stacked = vcat(v_sub...)
    indices = reduce(vcat, [fill(i, l) for (i, l) in enumerate(length.(pa.subarrays))])
    v_super_expanded = v_super[indices, :]

    return v_sub_stacked .* v_super_expanded
end

IsotropicArray(args...; kwargs...) = PhasedArray(IsotropicArrayManifold(args...; kwargs...))  

struct TappedDelayLine{MT, NT, FT} <: AbstractPhasedArray where {MT<:AbstractArrayManifold, NT<:Int, FT<:Number}
    manifold::MT
    num_taps::NT
    fs::FT
end

function Base.length(pa::TappedDelayLine)
    return length(pa.manifold)*pa.num_taps
end

"""
Nested structure: M elements and J taps. Each tap has it's own weight vector. Those 
are stacked on top of each other to one large weight vector for the whole TappedDelayLine.

Columns than correspond to angles, frequencies, and propagation speeds. 

Result is a large matrix: (M*J) x (A*F*C). Which can be reshaped into a (M*J)xAxFxC or 
even MxJxAxFxC tensor. 

Basically does: taps_phaseshifts_as_vector ⊗ steering_vector
But is adjusted for multiple angles, frequencies, and propagation speeds.

References:
-----------

W. Liu and S. Weiss, Wideband Beamforming. Nashville, TN: John Wiley & Sons, 2010.
"""
function steer(pa::TappedDelayLine, angles, f, c=c_0; kwargs...)
    # array manifold vector: M x (A*F*C)
    v_array = pa.manifold(angles, f, c; kwargs...)

    M = size(v_array, 1)
    A = length(angles)
    F = length(f)
    C = length(c)

    # number of taps and sample period
    J = pa.num_taps
    Ts = convert(eltype(f), 1 / pa.fs)

    # taps 'steering vectors': (J x F)
    ωs = 2π .* f

    # J x F
    v_taps = similar(v_array, complex(eltype(f)), J, length(f))
    copyto!(v_taps, hcat([exp.(-1im .* ω .* (0:J-1) .* Ts) for ω in ωs]...))  # vector of J-vectors

    # reshape to broadcast
    v_array_reshaped = reshape(v_array, M, 1, A, F, C)   # M x 1 x A x F x C
    v_taps_reshaped = reshape(v_taps, 1, J, 1, F, 1)     # 1 x J x 1 x F x 1

    
    # broadcast outer product over taps
    v_full = v_array_reshaped .* v_taps_reshaped         # M × J × A × F × C

    # reshape to: (M*J) × (A*F*C)
    v = reshape(v_full, M*J, A*F*C)

    return v
end
