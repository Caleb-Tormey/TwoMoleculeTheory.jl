using StaticArrays
using LinearAlgebra

export RadialGrid, SiteType, MolecularTopology, ThermodynamicState
export count_site_types, count_monomers

"""
    RadialGrid{T<:AbstractFloat}
Stores real-space (r) and reciprocal-space (k) grids.
"""
struct RadialGrid{T<:AbstractFloat}
    N::Int
    dr::T
    dk::T
    r::Vector{T}
    k::Vector{T}

    function RadialGrid(N::Int, dr::T) where {T<:AbstractFloat}
        dk = T(π) / (N * dr)
        r = collect(T(dr):T(dr):T(N * dr))
        k = collect(T(dk):T(dk):T(N * dk))
        new{T}(N, dr, dk, r, k)
    end
end

"""
    SiteType{T<:AbstractFloat}
Physical parameters of a single interaction site.
"""
struct SiteType{T<:AbstractFloat}
    id::Symbol
    σ::T
    ϵ::T
    mass::T
end

"""
    MolecularTopology{T}
Maps the physical chain to the theoretical matrix indices.
"""
struct MolecularTopology{T}
    sites::Vector{SiteType{T}}
    chain_sequence::Vector{Int}
    bonds::Vector{SVector{2, Int}}
    
    function MolecularTopology(sites::Vector{SiteType{T}}, sequence::Vector{Int}, bonds::Vector{SVector{2, Int}}) where {T}
        max_id = maximum(sequence)
        if max_id > length(sites)
            error("Topology sequence references site index $max_id, but only $(length(sites)) site types defined.")
        end
        new{T}(sites, sequence, bonds)
    end
end

"""
    ThermodynamicState{T<:AbstractFloat}
Macroscopic state (T, rho, beta).
"""
struct ThermodynamicState{T<:AbstractFloat}
    temp::T
    ρ::T
    β::T
    k_B::T

    function ThermodynamicState(temp::T, ρ::T; k_B::T=0.001987204) where {T<:AbstractFloat}
        beta = 1 / (k_B * temp)
        new{T}(temp, ρ, beta, k_B)
    end
end

# Helper functions
count_site_types(topo::MolecularTopology) = length(topo.sites)
count_monomers(topo::MolecularTopology) = length(topo.chain_sequence)