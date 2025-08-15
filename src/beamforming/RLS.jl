struct RLS{μT<:Real, PT<:AbstractMatrix} <: AbstractUpdateMethod
    μ::μT
    P::PT
end