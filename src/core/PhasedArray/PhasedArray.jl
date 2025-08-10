abstract type AbstractPhasedArray end

struct PhasedArray{T<:AbstractArrayManifold} <: AbstractPhasedArray
    manifold::T
end

struct NestedArray{E<:AbstractPhasedArray, S<:AbstractVector{<:AbstractPhasedArray}} <: AbstractPhasedArray
    elements::E
    subarrays::S
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
    return reduce(vcat, map((sup_row, sub_mat) -> sub_mat .* reshape(sup_row, 1, :), eachrow(v_super), v_sub))
end

IsotropicArray(args...; kwargs...) = PhasedArray(IsotropicArrayManifold(args...; kwargs...))  