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

function steer(pa::PhasedArray, args...; kwargs...)
    return pa.manifold(args...; kwargs...)
end

#TODO: rework the large reduce() line
function steer(pa::NestedArray, args...; kwargs...)
    v_super = steer(pa.elements, args...; kwargs...)
    v_sub = map(sa -> steer(sa, args...; kwargs...), pa.subarrays)

    v_sub_stacked = vcat(v_sub...)
    indices = reduce(vcat, [fill(i, l) for (i, l) in enumerate(length.(pa.subarrays))])
    v_super_expanded = v_super[indices, :]

    return v_sub_stacked .* v_super_expanded
end

IsotropicArray(args...; kwargs...) = PhasedArray(IsotropicArrayManifold(args...; kwargs...))  