using Unitful
import AtomsBase

function _cell_to_atomic_system(
    lattice::Lattice{T},
    species::AbstractArray,
    positions::AbstractArray,
    arrays::Dict;
    cell_unit=u"Å",
    kwargs...,
) where {T}
    cm = cellmat(lattice)
    ndim = size(cm, 1)
    box = SVector{ndim,SVector{ndim,T}}(SVector{ndim,T}(col) for col in eachcol(cm))
    # Unpack :arrays
    kws = [Dict{Symbol,Any}() for _ = 1:length(species)]
    for (key, value) in arrays
        for (iat, prop) in enumerate(eachslice(value; dims=ndims(value)))
            kws[iat][key] = prop
        end
    end

    atoms = [
        AtomsBase.Atom(symbol, pos .* cell_unit; kw...) for
        (symbol, pos, kw) in zip(species, eachcol(positions), kws)
    ]
    AtomsBase.periodic_system(atoms, box .* cell_unit; kwargs...)
end

"""
    atomic_system(cell::Cell;cell_unit=u"Å")

Construct a `Cell` object into a `FlexibleSystem`.
The length unit for the `Cell` object passed is assuemd to be `cell_unit`.
"""
AtomsBase.atomic_system(cell::Cell; cell_unit=u"Å") = _cell_to_atomic_system(
    lattice(cell),
    species(cell),
    positions(cell),
    cell.arrays;
    cell_unit,
    metadata(cell)...,
)

"""
    Cell(system::AbstractSystem;cell_unit=u"Å")

Construct a `Cell` object from AbstractSystem.
The length unit for the `Cell` object returned well be in the `cell_unit`.
"""
function Cell(system::AtomsBase.AbstractSystem; cell_unit=u"Å")
    pos = hcat(map(x -> collect(ustrip.(cell_unit, x)), AtomsBase.position(system, :))...)
    cm = hcat(map(x -> collect(ustrip.(cell_unit, x)), AtomsBase.bounding_box(system))...)
    @assert all(AtomsBase.periodicity(system))
    out = Cell(Lattice(cm), map(Symbol, AtomsBase.atomic_symbol(system, :)), pos)
    # Store sys.data in metadata
    for key in keys(system)
        key in (:periodicity, :bounding_box) && continue
        out.metadata[key] = system[key]
    end

    for key in AtomsBase.atomkeys(system)
        key in (:position, :specie, :mass) && continue
        arr = system[:, key]
        # Convert vector of arrays/vectors back to matrix format
        if !isempty(arr) && first(arr) isa AbstractArray
            out.arrays[key] = hcat(arr...)
        else
            out.arrays[key] = arr
        end
    end
    out
end


# This is for compatibility with AtomsBase.jl 0.4
function AB.position(cell::Cell{T,N}, idx) where {T,N}
    [SVector{N,T}(x) .* u"Å" for x in eachcol(@view positions(cell)[:, idx])]
end

function AB.position(cell::Cell{T,N}, idx::Int) where {T,N}
    SVector{N,T}(positions(cell)[:, idx]) .* u"Å"
end

"""
    species(structure::Cell)

Return a Vector of species names.
"""
AB.species(structure::Cell) = AB.ChemicalSpecies.(structure.symbols)
AB.species(structure::Cell, idx) = AB.ChemicalSpecies.(structure.symbols[idx])


function AB.cell(cell::Cell{T}) where {T}
    AB.PeriodicCell(
        cell_vectors=NTuple{3}(SVector{3,T}(x) .* u"Å" for x in eachcol(cellmat(cell))),
        periodicity=(true, true, true),
    )
end

function Base.getindex(cell::Cell, idx::Int)
    pos = positions(cell)[:, idx] .* u"Å"
    AB.Atom(species(cell)[idx], pos)
end

# System property access
function Base.getindex(system::Cell, x::Symbol)
    if x === :bounding_box
        AtomsBase.bounding_box(system)
    elseif x === :periodicity
        AtomsBase.periodicity(system)
    else
        getindex(system.metadata, x)
    end
end

Base.getindex(system::Cell, ::Colon, x::Symbol) = getindex(system.arrays, x)

function Base.setindex!(system::Cell, value, key::Symbol)
    system.metadata[key] = value
end

function Base.setindex!(system::Cell, value::AbstractArray, ::Colon, key::Symbol)
    @assert size(value, ndims(value)) == length(system)
    system.arrays[key] = value
end

function Base.haskey(system::Cell, x::Symbol)
    x in (:bounding_box, :periodicity) || haskey(system.metadata, x)
end

function AB.bounding_box(cell::Cell{T}) where {T}
    NTuple{3}(SVector{3,T}(x) .* u"Å" for x in eachcol(cellmat(cell)))
end

Base.keys(cell::Cell) = (:bounding_box, :periodicity, keys(cell.metadata)...)
AB.periodicity(cell::Cell) = (true, true, true)
AB.n_dimensions(cell::Cell{T,N}) where {T,N} = N

const n_dimensions = AB.n_dimensions

function AB.mass(cell::Cell, idx)
    getproperty.(AB.element.(AB.species(cell, idx)), :atomic_mass)
end

# Distance squared between functions interface
distance_squared_between(s1::AtomsBase.Atom, s2::AtomsBase.Atom) = sum((s1.position .- s2.position) .^ 2)
distance_squared_between(s1::AtomsBase.Atom, s2::AtomsBase.Atom, shift) =
            sum((s1.position .- s2.position .- shift) .^ 2)
