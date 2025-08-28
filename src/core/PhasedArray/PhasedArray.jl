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
TappedDelayLine(manifold::AbstractArrayManifold, num_taps::Int, fs::Number) = PhasedArray(TappedDelayLineManifold(manifold, num_taps, fs))