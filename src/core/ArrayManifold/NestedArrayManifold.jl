struct NestedArrayManifold{
    E<:AbstractArrayManifold, 
    S<:Tuple{Vararg{<:AbstractArrayManifold}}, 
    I<:AbstractVector{<:Integer}
} <: AbstractArrayManifold
    super_manifold::E
    sub_manifolds::S
    indices::I
end

Adapt.@adapt_structure NestedArrayManifold

function NestedArrayManifold(super_manifold::AbstractArrayManifold, sub_manifolds)
    sub_tuple = Tuple(sub_manifolds)
    
    length(super_manifold) == length(sub_tuple) || 
        throw(DimensionMismatch("Length of super_manifold ($(length(super_manifold))) and number of sub_manifolds ($(length(sub_tuple))) do not match."))

    # Precalculate and store indices
    lengths = length.(sub_tuple)
    indices = reduce(vcat, [fill(i, l) for (i, l) in enumerate(lengths)])

    return NestedArrayManifold(super_manifold, sub_tuple, indices)
end

function Base.length(a::NestedArrayManifold)
    return length(a.indices)
end

function (a::NestedArrayManifold)(angles::Wavefront, f, c=c_0)
    v_super = a.super_manifold(angles, f, c)
    v_super_expanded = v_super[a.indices, :]
    v_sub_tuple = map(sa -> sa(angles, f, c), a.sub_manifolds)
    v_sub_stacked = vcat(v_sub_tuple...)
    return v_sub_stacked .* v_super_expanded
end