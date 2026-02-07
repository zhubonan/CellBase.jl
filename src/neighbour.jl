#=
Code for building neighbour lists and extended point array
=#
import Base
using Base.Threads
export ExtendedPointArray,
    NeighbourList, eachneighbour, nions_extended, nions_orig, num_neighbours
export rebuild!, update!

"Maximum number of shift vectors"
const MAX_SHIFTS = 100000


"""
Represent an array of points after expansion by periodic boundary
"""
struct ExtendedPointArray{T,D}
    "Original point indices"
    indices::Vector{Int}
    "Index of the shift"
    shiftidx::Vector{Int}
    "Shift vectors"
    shiftvecs::Vector{SVector{D,T}}
    "Positions"
    positions::Vector{SVector{D,T}}
    "original Positions"
    orig_positions::Vector{SVector{D,T}}
    "Index of of the all-zero shift vector"
    rcut::T
    lattice::Matrix{T}
end

function Base.show(io::IO, s::ExtendedPointArray)
    print(
        io,
        "ExtendedPointArray of $(length(s.indices)) points from $(length(s.orig_positions)) points",
    )
end

"""
    ExtendedPointArray(cell::Cell, rcut)

Constructed an ExtendedPointArray from a given structure.
Implicitly, the positions are wrapped inside the unit cell, even if the actual 
in the original Cell is outside the unit cell. This ensures the correct neighbour list
begin constructed.
"""
function ExtendedPointArray(cell::Cell{T,D}, rcut) where {T,D}
    rcut = convert(Float64, rcut)
    ni = nions(cell)
    shifts = CellBase.shift_vectors(cellmat(lattice(cell)), rcut; safe=false)
    # Pad the shift vectors with zeros
    if D > 3
        shifts =
            [vcat(value, SVector{D - 3,T}([zero(T) for i = 1:D-3])) for value in shifts]
    end
    indices = zeros(Int, ni * length(shifts))
    shiftidx = zeros(Int, ni * length(shifts))
    pos_extended =
        zeros(eltype(positions(cell)), size(positions(cell), 1), ni * length(shifts))

    # Use the wrapped positions for building the point array
    wrapped = wrapped_spos(cell)
    i = 1
    for (idx, pos_orig) in enumerate(wrapped)   # Each original positions
        for (ishift, shiftvec) in enumerate(shifts)   # Each shift positions
            pos_extended[:, i] = pos_orig .+ shiftvec
            indices[i] = idx
            shiftidx[i] = ishift
            i += 1
        end
    end
    ExtendedPointArray{T,D}(
        indices,
        shiftidx,
        shifts,
        [SVector{size(pos_extended, 1)}(x) for x in eachcol(pos_extended)],
        wrapped,
        rcut,
        copy(cellmat(lattice(cell))),
    )
end

"""
    rebuild!(ea::ExtendedPointArray, cell)

Rebuild the ExtendedPointArray for an existing cell
"""
function rebuild!(ea::ExtendedPointArray, cell)
    lattice_change = !all(ea.lattice .== cellmat(cell))
    # Rebuild shift vectors from scratch is the lattice vectors have been changed
    if lattice_change
        _update_ea_with_lattice_change(ea, cell)
    else
        # no change in lattice shift all extended points with the displacements of the original ones
        _update_ea_no_lattice_change(ea, cell)
    end
    ea
end

"""
    _update_ea_with_lattice_change(ea, cell)

Update extended points with lattice shifts - need to rebuild from scratch
"""
function _update_ea_with_lattice_change(ea::ExtendedPointArray{T,D}, cell) where {T,D}
    newshifts = CellBase.shift_vectors(cellmat(lattice(cell)), ea.rcut; safe=false)
    if D > 3
        newshifts =
            [vcat(value, SVector{D - 3,T}([zero(T) for i = 1:D-3])) for value in newshifts]
    end
    # Number of vectors change
    nnew = length(newshifts)
    nold = length(ea.shiftvecs)
    dl = nnew - nold
    @assert nnew < MAX_SHIFTS "Too many lattice shifts request ($nnew) - possible ill shaped cell?"
    if dl != 0
        resize!(ea.shiftvecs, nnew)
        # Also resize the positions array and indices
        ntot = nnew * length(ea.orig_positions)
        resize!(ea.positions, ntot)
        resize!(ea.indices, ntot)
        resize!(ea.shiftidx, ntot)
    end
    ea.shiftvecs .= newshifts
    ea.lattice .= cellmat(lattice(cell))

    # Rebuild the extended arrays
    ea.orig_positions .= wrapped_spos(cell)
    nshift = length(ea.shiftvecs)
    Threads.@threads for idx = 1:length(ea.orig_positions)   # Each original positions
        pos_orig = ea.orig_positions[idx]
        for (ishift, shiftvec) in enumerate(ea.shiftvecs)   # Each shift positions
            i = (idx - 1) * nshift + ishift
            ea.positions[i] = pos_orig .+ shiftvec
            ea.indices[i] = idx
            ea.shiftidx[i] = ishift
        end
    end
end

"""
    _update_ea_no_lattice_change(ea, cell)

Update extended points without lattice shifts - apply displacements to all image points
"""
function _update_ea_no_lattice_change(ea, cell)
    spos = wrapped_spos(cell)
    ns = length(ea.shiftvecs)
    Threads.@threads for idx = 1:length(spos)
        pos = spos[idx]
        disp = pos - ea.orig_positions[idx]
        # Displace all images by this amount
        for j = 1:length(ea.shiftvecs)
            m = (idx - 1) * ns + j
            ea.positions[m] = ea.positions[m] + disp
        end
    end
    ea.orig_positions .= spos
end

"""
    NeighbourList{T,N}

Type for representing a neighbour list
"""
mutable struct NeighbourList{T,D}
    ea::ExtendedPointArray{T,D}
    "Extended indices of the neighbours"
    extended_indices::Matrix{Int}
    "Original indices of the neighbours"
    orig_indices::Matrix{Int}
    "Distance to the neighbours"
    distance::Matrix{T}
    "Vector displacement to the neighbours"
    vectors::Array{SVector{D,T},2}
    "Number of neighbours"
    nneigh::Vector{Int}
    "Maximum number of neighbours that can be stored"
    nmax::Int
    "Contains vector displacements or not"
    has_vectors::Bool
    rcut::T
    nmax_limit::Int
    last_rebuild_positions::Vector{SVector{D,T}}
    skin::T
end

"""Check if a static vector only contains zeros"""
allzeros(svec::SVector{1}) = (svec[1] == 0)
allzeros(svec::SVector{2}) = (svec[1] == 0) && (svec[2] == 0)
allzeros(svec::SVector{3}) = (svec[1] == 0) && (svec[2] == 0) && (svec[3] == 0)
allzeros(svec::SVector{4}) =
    (svec[1] == 0) && (svec[2] == 0) && (svec[3] == 0) && (svec[4] == 0)
allzeros(svec::SVector) = !any(x -> x != 0, svec)

"Number of ions in the original cell"
nions_orig(n::ExtendedPointArray) = length(n.orig_positions)
"Number of ions in the original cell"
nions_orig(n::NeighbourList) = nions_orig(n.ea)

"Number of ions in the extended cell"
nions_extended(n::ExtendedPointArray) = length(n.positions)
"Number of ions in the extended cell"
nions_extended(n::NeighbourList) = nions_extended(n.ea)


function Base.show(io::IO, n::NeighbourList)
    println(io, "NeighbourList of $(nions_orig(n)) ($(nions_extended(n)) extended)")
    print(io, "Current max neighbours $(maximum(n.nneigh))\n")
    print(io, "Cut off radius: $(n.rcut)")
end

"""
    NeighbourList(ea::ExtendedPointArray, rcut, nmax=1000; savevec=false)

Construct a NeighbourList from an extended point array for the points in the original cell
"""
function NeighbourList(
    ea::ExtendedPointArray{T,D},
    rcut,
    nmax=1000;
    savevec=false,
    ndim=length(ea.positions[1]),
    nmax_limit=5000,
    skin=-1.0,
) where {T,D}
    rcut = convert(Float64, rcut)

    # If using skin, update the rcut
    if skin > 0
        rcut = rcut + skin
    end

    norig = length(ea.orig_positions)
    extended_indices = zeros(Int, nmax, norig)
    orig_indices = zeros(Int, nmax, norig)
    distance = fill(-1.0, nmax, norig)
    nneigh = zeros(Int, norig)

    # Ensure the nmax_limit is at least four times that of nmax otherwise there is not much point
    if nmax_limit < nmax * 4
        nmax_limit = nmax * 4
    end

    @assert rcut >= ea.rcut "Cut off radius is large that that of the periodic image cut off."

    # Save vectors or not
    base = @SVector fill(-1.0, ndim)
    savevec ? vectors = fill(base, nmax, norig) : vectors = fill(base, 1, 1)
    nl = NeighbourList{T,D}(
        ea,
        extended_indices,
        orig_indices,
        distance,
        vectors,
        nneigh,
        nmax,
        savevec,
        rcut,
        nmax_limit,
        T[],
        skin,
    )
    rebuild!(nl, ea)
    nl
end

"""
    increase_nmax!(nl::NeighbourList, nmax)

Increase the maximum number of neighbours storable in the neighbour list
"""
function update_nmax!(nl::NeighbourList, nmax)
    norig = length(nl.ea.orig_positions)
    nl.extended_indices = zeros(Int, nmax, norig)
    nl.orig_indices = zeros(Int, nmax, norig)
    nl.distance = fill(-1.0, nmax, norig)
    nl.nmax = nmax
    if nl.has_vectors
        base = nl.vectors[1]
        nl.vectors = fill(base, nmax, norig)
    end
    rebuild!(nl, nl.ea)
end

"""
    _need_rebuild(nl::NeighbourList)

Check if a full rebuild is needed for the NeighbourList if the `skin` is used and
if any atom has move more than the skin.
"""
function _need_rebuild(nl::NeighbourList)
    if nl.skin < 0
        return true
    end
    # Initial build?
    isempty(nl.last_rebuild_positions) && return true
    # Compute the displacement
    for (orig, moved) in zip(nl.ea.orig_positions, nl.last_rebuild_positions)
        (norm(orig - moved) >= nl.skin) && return true
    end
    return false
end

"""
    rebuild!(nl::NeighbourList, ea::ExtendedPointArray)

Perform a full rebuild of the neighbour list from scratch for a given ExtendedPointArray.
Extended the neighbour storage space if necessary.
"""
function rebuild!(nl::NeighbourList, ea::ExtendedPointArray)

    # Check if we actually need rebuilding or not....
    if !_need_rebuild(nl)
        return update!(nl)
    end

    (; extended_indices, orig_indices, distance, nneigh, nmax, has_vectors, vectors, rcut) =
        nl

    # Reset
    fill!(vectors, vectors[1] .* 0.0)
    fill!(distance, -1.0)
    fill!(orig_indices, 0)
    fill!(extended_indices, 0)
    fill!(nneigh, 0)

    ineigh_max = Atomic{Int}(0)
    Threads.@threads for iorig = 1:length(ea.orig_positions)
        posi = ea.orig_positions[iorig]
        ineigh = 0
        for (j, posj) in enumerate(ea.positions)
            # Skip if it is the same point 
            (ea.indices[j] == iorig) && (allzeros(ea.shiftvecs[ea.shiftidx[j]])) && continue
            dist = distance_between(posj, posi)
            # Store the information if the distance is smaller than the cut off
            if dist < rcut
                ineigh += 1
                if ineigh <= nmax
                    distance[ineigh, iorig] = dist
                    # Store the extended index
                    extended_indices[ineigh, iorig] = j
                    # Store the index of the point in the original cell
                    orig_indices[ineigh, iorig] = ea.indices[j]
                    has_vectors && (vectors[ineigh, iorig] = posj .- posi)
                end
            end
        end
        # Store the total number of neighbours for this point
        if ineigh > nmax
            if nmax > ineigh_max[]
                ineigh_max[] = ineigh
            end

        end
        nneigh[iorig] = ineigh
    end
    imax = ineigh_max[]
    if imax > 0
        new_nmax = min(imax * 2, nl.nmax_limit)
        if new_nmax < imax
            throw(
                ErrorException(
                    "Cannot increase the reuqired beyong the `nmax_limit``: $(nl.nmax_limit)",
                ),
            )
        end
        # @warn "Too many neighbours - increasing the value of nmax to $(new_nmax) × 2 and rebuild"
        update_nmax!(nl, new_nmax)
        rebuild!(nl, ea)
    end
    # If using skin we need to copy the last positions
    if nl.skin > 0
        nl.last_rebuild_positions = collect(ea.orig_positions)
    end

    nl
end

"""
    rebuild(nl::NeighbourList, cell::Cell) 

Perform a full rebuild of the NeighbourList with the latest geometry of the cell
"""
function rebuild!(nl::NeighbourList, cell::Cell)
    rebuild!(nl.ea, cell)
    rebuild!(nl, nl.ea)
end

"""
Update the calculated distances/vectors but do not rebuild the whole neighbour list
"""
function update!(nl::NeighbourList)
    ea = nl.ea
    for (isite, nneight) in enumerate(nl.nneigh)
        inn = 1
        posi = ea.orig_positions[isite]
        while inn <= nneight
            jid_ext = nl.extended_indices[inn, isite]
            posj = ea.positions[jid_ext]
            dij = distance_between(posj, posi)
            # Update distance
            nl.distance[inn, isite] = dij
            # Update the vector
            nl.has_vectors && (nl.vectors[inn, isite] = posj - posi)
            inn += 1
        end
    end
end

"""
    update!(nl::NeighbourList, cell::Cell) 

Update the NeighbourList with the latest geometry of the Cell.
No rebuilding is performed.
"""
function update!(nl::NeighbourList, cell::Cell)
    rebuild!(nl.ea, cell)
    update!(nl)
end


NeighbourList(cell::Cell, rcut, nmax=100; savevec=false) =
    NeighbourList(ExtendedPointArray(cell, rcut), rcut, nmax; savevec)

"Number of neighbours for a point"
num_neighbours(nl::NeighbourList, iorig) = nl.nneigh[iorig]


struct Neighbour{T}
    index::Int
    index_extended::Int
    position::T
    shift_vector::T
    distance::Float64
end

"""
    get_neighbour(nl::NeighbourList, iorig::Int, idx::Int)

Return a single neighbour.
"""
function get_neighbour(nl::NeighbourList, iorig::Int, idx::Int)
    i = nl.orig_indices[idx, iorig]
    j = nl.extended_indices[idx, iorig]
    shift_vec = nl.ea.shiftvecs[(j-1)%length(nl.ea.shiftvecs)+1]
    pos = nl.ea.orig_positions[i]
    dist = nl.distance[idx, iorig]
    Neighbour(i, j, pos, shift_vec, dist)
end

"""
    get_neighbours(nl::NeighbourList, iorig::Int)

Return a Vector of Neighbour objects
"""
function get_neighbours(nl::NeighbourList, iorig::Int)
    [get_neighbour(nl, iorig, idx) for idx = 1:nl.nneigh[iorig]]
end


abstract type AbstractNLIterator end

"Iterator interface for going through all neighbours"
struct NLIterator{T,G} <: AbstractNLIterator
    nl::T
    iorig::Int
    unique::G
end

"Iterator interface for going through all neighbours"
struct NLIteratorWithVector{T} <: AbstractNLIterator
    nl::T
    iorig::Int
end


Base.length(nli::AbstractNLIterator) = num_neighbours(nli.nl, nli.iorig)

function Base.iterate(nli::NLIterator, state=1)
    nl = nli.nl
    iorig = nli.iorig
    if state > nl.nneigh[iorig]
        return nothing
    end
    # Handle request for unique primitive cell image
    if nli.unique
        while true
            if state > nl.nneigh[iorig]
                return nothing
            end
            oidx = nl.orig_indices[state, iorig]
            # For if we have seen this neighbour
            seen = any(x -> nl.orig_indices[x, iorig] == oidx, 1:state-1)
            if seen
                # Search the next neghbour
                state += 1
            else
                # Unique neighbour found - stop
                break
            end
        end
    end
    return (
        nl.orig_indices[state, iorig],
        nl.extended_indices[state, iorig],
        nl.distance[state, iorig],
    ),
    state + 1
end

function Base.iterate(nli::NLIteratorWithVector, state=1)
    nl = nli.nl
    iorig = nli.iorig
    if state > nl.nneigh[nli.iorig]
        return nothing
    end
    return (
        nl.orig_indices[state, iorig],
        nl.extended_indices[state, iorig],
        nl.distance[state, iorig],
        nl.vectors[state, iorig],
    ),
    state + 1
end

"""
Iterate the neighbours of a site in the original cell.
Returns a tuple of (original_index, extended_index, distance) for each iteration
"""
eachneighbour(nl::NeighbourList, iorig; unique=false) = NLIterator(nl, iorig, unique)


"""
Iterate the neighbours of a site in the original cell.
Returns a tuple of (original_index, extended_index, distance, vector) for each iteration
"""
function eachneighbourvector(nl::NeighbourList, iorig)
    @assert nl.has_vectors "NeighbourList is not build with distance vectors"
    NLIteratorWithVector(nl, iorig)
end
