using Printf
using PeriodicTable
using LinearAlgebra
import Base

export Cell,
    nions,
    positions,
    species,
    atomic_numbers,
    lattice,
    volume,
    get_cellmat,
    cellmat,
    cellvecs,
    wrap!,
    cellpar,
    natoms,
    sposarray
export set_scaled_positions!,
    get_scaled_positions,
    set_cellmat!,
    set_positions!,
    get_positions,
    get_lattice,
    get_cellmat

"""
A Cell represents a periodic structure in three-dimensional space.

Defined as:
```julia
mutable struct Cell{D, T}
    lattice::Lattice{T}                 # Lattice of the structure
    symbols::Vector{Symbol}
    positions::Matrix{T}
    arrays::Dict{Symbol, Any}        # Any additional arrays
    metadata::Dict{Symbol, Any}
end
```

"""
mutable struct Cell{T, D}  <: AB.AbstractSystem{D}
    lattice::Lattice{T}                 # Lattice of the structure
    symbols::Vector{Symbol}
    positions::Matrix{T}
    arrays::Dict{Symbol, Any}        # Any additional arrays
    metadata::Dict{Symbol,Any}
end

"""
    Cell(l::Lattice, symbols, positions) where T

Construct a Cell type from arrays
"""
function Cell(l::Lattice, symbols::Vector{Symbol}, positions::Matrix)
    arrays = Dict{Symbol, Any}()
    @assert length(symbols) == size(positions, 2)
    Cell{eltype(positions), size(positions, 1)}(l, symbols, positions, arrays, Dict{Symbol,Any}())
end

"""
    Cell(lat::Lattice, numbers::Vector{Int}, positions)

Constructure the Cell type from lattice, positions and numbers
"""
function Cell(lat::Lattice, numbers::Vector{T}, positions::Matrix) where {T<:Real}
    species = [Symbol(elements[i].symbol) for i in numbers]
    @assert length(numbers) == size(positions, 2)
    Cell(lat, species, positions)
end


"""
    Cell(lat::Lattice, numbers::Vector{Int}, positions::Vector)

Constructure the Cell type from lattice, positions and numbers
"""
function Cell(lat::Lattice, species_id, positions::Vector)
    @assert length(species_id) == length(positions)
    posmat = zeros(length(positions[1]), length(positions))
    for (i, vec) in enumerate(positions)
        posmat[:, i] = vec
    end
    Cell(lat, species_id, posmat)
end


"""
    clip(s::Cell, mask::AbstractVector)

Clip a structure with a given indexing array
"""
function clip(cell::Cell{T, N}, mask::AbstractVector) where {T, N}
    new_pos = positions(cell)[:, mask]
    new_symbols = species(cell)[mask]
    # Clip any additional arrays
    new_array = Dict{Symbol,Any}()
    for (key, array) in pairs(cell.arrays)
        new_array[key] = selectdim(array, ndims(array), mask)
    end
    Cell{T, N}(lattice(cell), new_symbols, new_pos, new_array, cell.metadata)
end

Base.getindex(cell::Cell, i::AbstractVector) = clip(cell, i)


"""
    Base.sort(cell::Cell)

Sort a `Cell` with the positions sorted by the species kinds.
"""
Base.sort(cell::Cell) = cell[sortperm(species(cell))]

"""
    Base.sort!(cell::Cell)

In-place sort a `Cell` with the positions sorted by the species kinds.

!!! note
    This will in-place-update all underlying Array - use with caution where the
    data is shared between multiple instances.

"""
function Base.sort!(cell::Cell)
    idx = sortperm(species(cell))
    species(cell) .= species(cell)[idx]
    cell.positions .= cell.positions[:, idx]
    for array in values(cell.arrays)
        array .= selectdim(array, ndims(array), idx)
    end
    cell
end


# Basic interface 
"""
    nions(cell::Cell)

Return number of atoms in a structure.
"""
nions(cell::Cell) = length(cell.symbols)

const natoms = nions
@doc """
    natoms(cell::Cell)

Return number of atoms in a structure.
"""
natoms

"""
    positions(cell::Cell)

Return positions (cartesian coordinates) of the atoms in a structure.
"""
positions(cell::Cell) = cell.positions


"""
    get_positions(cell::Cell)

Return a *copy* of the positions (cartesian coordinates) of the atoms in a structure.
"""
get_positions(cell::Cell) = copy(cell.positions)

"""
    sposarray(cell::Cell)

Return the positions as a Vector of static arrays.
The returned array can provide improved performance for certain type of operations.
"""
sposarray(structure::Cell{T, N}) where {T, N} =
    [SVector{N,T}(x) for x in eachcol(positions(structure))]

"""
    species(structure::Cell)

Return a Vector of species names.
"""
species(structure::Cell) = structure.symbols

"""
    get_species(structure::Cell)

Return a Vector (copy) of species names.
"""
get_species(structure::Cell) = copy(structure.symbols)

"""
    atomic_numbers(structure::Cell)

Return a Vector of the atomic numbers.
"""
atomic_numbers(structure::Cell) = Int[elements[x].number for x in species(structure)]

## Wrapper for the Lattice ###

"""
    lattice(structure::Cell)

Return the `Lattice` instance.
"""
lattice(structure::Cell) = structure.lattice

"""
    get_lattice(structure::Cell)

Return the `Lattice` instance (copy).
"""
get_lattice(structure::Cell) = deepcopy(structure.lattice)

"""
    cellpar(structure::Cell)

Return the lattice parameters.
"""
cellpar(structure::Cell) = cellpar(lattice(structure))

"""
    volume(structure::Cell)

Return the volume of the cell.
"""
volume(structure::Cell) = volume(lattice(structure))

"""
    cellmat(structure::Cell)

Return the matrix of lattice vectors.
"""
cellmat(structure::Cell) = cellmat(lattice(structure))

"""
    get_cellmat(structure::Cell)

Return the matrix of lattice vectors(copy).
"""
get_cellmat(structure::Cell) = get_cellmat(lattice(structure))

## END ##


"""
    array(structure::Cell, arrayname::Symbol)

Return the additional array stored in the `Cell` object.
"""
array(structure::Cell, arrayname::Symbol) = structure.arrays[arrayname]

"""
    arraynames(structure::Cell)

Return the names of additional arrays.
"""
arraynames(structure::Cell) = keys(structure.arrays)

"""
    metadata(structure::Cell)

Return the `metadata` dictionary. 
"""
metadata(structure::Cell) = structure.metadata

"""
    attachmetadata!(structure::Cell, metadata::Dict)

Replace `metadata` with an existing dictionary.
"""
attachmetadata!(structure::Cell, metadata::Dict) = structure.metadata = metadata

"""
    num_fu(structure::Cell)

Return the number of formula units.
"""
num_fu(structure::Cell) = formula_and_factor(structure)[2]

"""
    reduced_fu(structure::Cell)

Return the reduced formula.
"""
reduced_fu(structure::Cell) = formula_and_factor(structure)[1]

"""
    sorted_symbols(symbols)

Return sorted symbols by atomic numbers.
"""
function sorted_symbols(symbols)
    z_array = [elements[sp].number for sp in symbols]
    perm = sortperm(z_array)
    return symbols[perm]
end

"""
    formula_and_factor(structure::Cell)

Return reduced formula and associated factor.
"""
function formula_and_factor(structure::Cell)

    # Check if computed results already exists
    metadata_dict = metadata(structure)
    rformula = get(metadata_dict, :formula, :None)
    num_fu = get(metadata_dict, :num_fu, 0)

    if (rformula != :None) & (num_fu != 0)
        return rformula, num_fu
    end

    # Find the reduced formula units
    sp_array = species(structure)
    unique_sp = sorted_symbols(unique(sp_array))
    num_atoms = Array{Int}(undef, size(unique_sp))
    for i = 1:length(unique_sp)
        num_atoms[i] = count(x -> x == unique_sp[i], sp_array)
    end
    num_fu = gcd(num_atoms)
    num_atoms ./= num_fu
    args = Symbol[]
    for i = 1:length(unique_sp)
        push!(args, unique_sp[i])
        if num_atoms[i] > 1
            push!(args, Symbol(num_atoms[i]))
        end
    end
    rformula = Symbol(args...)
    metadata(structure)[:formula] = rformula
    metadata(structure)[:num_fu] = num_fu
    return rformula, num_fu
end

"""
    specindex(structure::Cell)

Return the unique species and integer based indices for each atom.
"""
function specindex(structure)
    # Mapping between the speices as symbols and as intgers
    unique_spec = unique(species(structure))
    spec_indices = [findfirst(x -> x == sym, unique_spec) for sym in species(structure)]
    return unique_spec, spec_indices
end

"""
    get_fraction_positions(cell::Cell)

Return fractional positions of the cell.
"""
function get_scaled_positions(cell::Cell)
    rec_cellmat(lattice(cell)) * positions(cell)
end

"""
    set_scaled_positions!(cell::Cell, scaled::Matrix)

Set scaled positions for a cell.
"""
function set_scaled_positions!(cell::Cell, scaled::Matrix)
    cell.positions .= cellmat(cell) * scaled
end

"""
    wrap!(cell::Cell)

Wrap an atom outside of the lattice back into the box defined by the lattice vectors.
"""
function wrap!(cell::Cell)
    scaled = get_scaled_positions(cell)
    scaled .-= floor.(scaled)
    set_scaled_positions!(cell, scaled)
end

"""
    wrap!(vec::AbstractVector, l::Lattice)

Wrap a vector back to the periodic box defined by the lattice vectors.
"""
function wrap!(vec::AbstractVector, l::Lattice)
    frac = l.rec * vec 
    frac .-= floor.(frac)
    vec .= l.matrix * frac
end

wrap!(vec::AbstractVector, c::Cell) = wrap!(vec, lattice(c))


"""
    wrapped_spos(cell)

Return a static array of wrapped positons.
"""
function wrapped_spos(cell::Cell{T, 3}) where {T}
    posarray = sposarray(cell)
    recmat = SMatrix{3,3}(rec_cellmat(lattice(cell)))
    cmat = SMatrix{3,3}(cellmat(cell))
    for i = 1:length(posarray)
        x = recmat * posarray[i]
        x -= floor.(x)
        x = cmat * x
        posarray[i] = x
    end
    posarray
end

"""
    wrapped_spos(cell)

Return a static array of wrapped positons.
"""
function wrapped_spos(cell::Cell)
    posarray = sposarray(cell)
    recmat = SMatrix{3,3}(rec_cellmat(lattice(cell)))
    cmat = SMatrix{3,3}(cellmat(cell))
    for i = 1:length(posarray)
        x = recmat * posarray[i][1:3]
        x -= floor.(x)
        x = cmat * x
        posarray[i] = vcat(x, posarray[i][4:end])
    end
    posarray
end

"""
    set_cellmat!(cell::Cell, mat;scale_positions=true)

Update the `Lattice` with a new matrix of lattice vectors. 
Scale of the existing postions if needed.
"""
function set_cellmat!(cell::Cell, mat; scale_positions=true)
    if scale_positions
        scaled_pos = get_scaled_positions(cell)
        set_cellmat!(lattice(cell), mat)
        set_scaled_positions!(cell, scaled_pos)
    else
        set_cellmat!(lattice(cell), mat)
    end
    cell
end

"""
    set_positions!(cell::Cell, pos)

Set the positions of the `Cell` with a new matrix.
"""
set_positions!(cell::Cell, pos) = cell.positions .= pos

"""
    rattle!(cell::Cell, amp)

Rattle the positions of the cell for a given maximum amplitude (uniform distribution).
"""
function rattle!(cell::Cell, amp)
    dev = rand(length(positions(cell)))
    for i in eachindex(cell.positions)
        cell.positions[i] += (rand() - 0.5) * 2amp
    end
end

function Base.show(io::IO, s::Cell)
    a, b, c, α, β, γ = cellpar(lattice(s))
    sym = join(map(string, species(s)))
    print(io, "Cell $(sym) with $(nions(s)) atoms, lattice parametrs: $a $b $c $α $β $γ")
end

function Base.show(io::IO, ::MIME"text/plain", s::Cell)
    println(io, "Cell with $(nions(s)) ions")
    println(io, "Lattice: ")
    cellmat = s.lattice.matrix
    posmat = positions(s)
    for i = 1:3
        println(
            io,
            @sprintf "%8.3f   %8.3f   %8.3f" cellmat[1, i] cellmat[2, i] cellmat[3, i]
        )
    end

    println(io, "Sites: ")
    sym = species(s)
    for i = 1:nions(s)
        symbol = sym[i]
        println(
            io,
            @sprintf "%4s  %10.5f  %10.5f  %10.5f" symbol posmat[1, i] posmat[2, i] posmat[
                3,
                i,
            ]
        )
    end
end

Base.length(cell::Cell) = natoms(cell)
#Base.getindex(cell::Cell, i::Int) = Site(@view(cell.positions[:, i]), i, cell.symbols[i])


function distance_matrix(cell::Cell; mic=true)
    if mic == true
        _distance_matrix_mic(cell)
    else
        _distance_matrix_no_mic(cell)
    end
end


"""
    _distance_matrix_no_mic(structure::Cell)

Compute the distnace matrix without minimum image convention, e.g. the PBC is not respected.
"""
function _distance_matrix_no_mic(cell::Cell)
    # Compute the naive pair-wise vectors
    nn = nions(cell)
    pos = sposarray(cell)
    dmat = zeros(nn, nn)
    for i = 1:nn
        for j = i+1:nn
            d = norm(pos[j] - pos[i])
            dmat[i, j] = d
            dmat[j, i] = d
        end
    end
    dmat
end


"""
    _distance_matrix_mic(cell::Cell)

Compute the distance matrix for the given structure, using the minimum image convention (or not).

Note the returned matrix does not expand the cell. The distance matrix cannot be safety used for obtaining the minimum separations.
For example, a structure with a single atom would be a distance matrix containing only zero.

"""
function _distance_matrix_mic(cell::Cell)

    # Compute the naive pair-wise vectors
    nn = nions(cell)
    vecs = zeros(3, nn * nn)
    pos = sposarray(cell)
    dmat = zeros(nn, nn)
    ivec = 0
    for i = 1:nn
        for j = i+1:nn
            ivec += 1
            vecs[:, ivec] .= pos[j] .- pos[i]
        end
    end
    # Apply minimum image conventions
    _, dmic = mic(lattice(cell), @view(vecs[:, 1:ivec]))
    # Unpack computed distances
    ivec = 0
    for i = 1:nn
        for j = i+1:nn
            ivec += 1
            dmat[i, j] = dmic[ivec]
            dmat[j, i] = dmic[ivec]
        end
    end
    dmat
end

"""
    distance_squared_between(posmat::Matrix, i, j, svec::Matrix, ishift)

Return the squared distance between two positions stored in a matrix and shift vector.
"""
function distance_squared_between(posmat::Matrix, i, j, svec::Matrix, ishift)
    d2 = 0.0
    for n = 1:size(posmat, 1)
        d = posmat[n, j] - posmat[n, i] + svec[n, ishift]
        d2 += d * d
    end
    d2
end


"""
    check_minsep(structure::Cell, minsep::Dict)

Check if the minimum separation constraints are satisfied. Minimum separations are supplied
as an dictionary where the global version under the :global key. To supply the minimum separations
between A and B, pass Dict((:A=>:B)=>1.0).
Return true or false.
"""
function check_minsep(structure::Cell, minsep::Dict{T,Float64}) where {T}
    _, spec_indices, minsep_matrix = compute_minsep_mat(structure, minsep)
    dist_mat = distance_matrix(structure)
    ni = nions(structure)
    for i = 1:ni
        for j = i+1:ni
            si = spec_indices[i]
            sj = spec_indices[j]
            if minsep_matrix[si, sj] > dist_mat[i, j]
                return false
            end
        end
    end
    return true
end


"""
    minsep_matrix(structure::Cell, minsep::Dict)

Initialise the minimum separation matrix and species mapping.
Returns the unique species, integer indexed species and the minimum separation matrix.
"""
function compute_minsep_mat(structure, minsep::Dict{T,Float64}) where {T}
    unique_spec, spec_indices = specindex(structure)
    nunique = length(unique_spec)
    # Minimum separation matrix
    minsep_mat = zeros(nunique, nunique)
    for i = 1:nunique
        for j = i:nunique
            sA = unique_spec[i]
            sB = unique_spec[j]
            if (sA => sB) in keys(minsep)
                value = minsep[sA=>sB]
            elseif (sB => sA) in keys(minsep)
                value = minsep[sB=>sA]
            else
                value = get(minsep, :all, 1.0)
            end
            minsep_mat[i, j] = value
            minsep_mat[j, i] = value
        end
    end
    unique_spec, spec_indices, minsep_mat
end

"""
    make_supercell(structure::Cell, a, b, c)

Make a supercell

Currently only work with diagonal transform matrices.
TODO: Write function for the general cases....
"""
function make_supercell(structure::Cell, a, b, c)
    tmat = [
        a 0 0
        0 b 0
        0 0 c
    ]
    old_cell = cellmat(lattice(structure))
    new_cell = old_cell * tmat

    # Compute the shift vectors
    svec = shift_vectors(old_cell, a - 1, b - 1, c - 1, 0, 0, 0)
    nshifts = length(svec)

    current_pos = positions(structure)
    ns = natoms(structure)
    # New positions
    new_pos = zeros(3, nshifts * ns)
    for (i, shift) in enumerate(svec)
        for j = 1:ns
            idx = j + (i - 1) * ns  # New index
            for n = axes(new_pos, 1)
                @inbounds new_pos[n, idx] = current_pos[n, j] + shift[n]
            end
        end
    end
    new_spec = repeat(species(structure), nshifts)
    Cell(Lattice(new_cell), new_spec, new_pos)
end


"""
    fingerprint(s::Cell; dmat=distance_matrix(s), weighted=true, cut_bl=3.0)

Computed the fingerprint vector based on simple sorted pair-wise distances.
NOTE: Does not work for single atom cell!!
"""
function fingerprint(s::Cell; dmat=distance_matrix(s), weighted=true, cut_bl=3.0)
    # Cut off distance based on minimum bond length
    cut_bl = minimum(d for d in dmat if d > 0.0) * cut_bl
    # Allocate workspace
    nn, _ = size(dmat)
    dist = zeros(nn * nn)

    # Weighting matrix
    num = atomic_numbers(s)

    c = 0
    norm_weights = 0.0
    # Consider only the lower triangle
    for j = 1:nn
        for i = 1+j:nn
            d = dmat[i, j]
            d > cut_bl && continue
            c += 1
            dist[c] = d
            if weighted
                dist[c] *= num[i] * num[j]
                norm_weights += num[i] * num[j]
            end
        end
    end

    dist_out = dist[1:c]

    # Normlise
    if weighted
        dist_out ./= (norm_weights / c)
    end
    sort!(dist_out)
end


"""
    fingerprint_distance(f1::AbstractVector, f2::AbstractVector;lim=Inf)

Compute the deviation between two finger print vectors

Comparison is truncated by the size of the shortest vector of the two, or by the `lim`
key word.
"""
function fingerprint_distance(f1::AbstractVector, f2::AbstractVector; lim=Inf)
    l1 = length(f1)
    l2 = length(f2)
    comp = min(l1, l2)
    d = 0.0
    ncomp = 0
    for i = 1:comp
        f1[i] > lim && break
        f2[i] > lim && break
        d += abs(f1[i] - f2[i])
        ncomp += 1
    end
    d / comp
end
