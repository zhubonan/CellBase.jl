using Printf
using PeriodicTable
using LinearAlgebra
using Statistics
import Base: sort, sort!, repeat

export Cell,
    Cell3D,
    nions,
    positions,
    periodicity,
    species,
    atomic_numbers,
    masses,
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
export make_supercell,
    repeat,
    rotate,
    rotate!
export reduced_fu,
    num_fu,
    formula_and_factor,
    array,
    arraynames,
    distance_matrix,
    rattle!

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
mutable struct Cell{T,D} <: AB.AbstractSystem{D}
    lattice::Lattice{T}                 # Lattice of the structure
    symbols::Vector{Symbol}
    positions::Matrix{T}
    arrays::Dict{Symbol,Any}        # Any additional arrays
    metadata::Dict{Symbol,Any}
    periodicity::NTuple{D,Bool}
end

_n_dimensions(::Cell{T,D}) where {T,D} = D

_default_periodicity(::Val{D}) where {D} = ntuple(i -> i <= 3, D)

function _normalize_periodicity(periodicity, ::Val{D}) where {D}
    normalized =
        isnothing(periodicity) ? _default_periodicity(Val(D)) : Tuple(Bool.(collect(periodicity)))
    length(normalized) == D ||
        throw(ArgumentError("Periodicity must have length $D, got $(length(normalized))"))
    normalized
end

function _validate_periodicity(periodicity::NTuple{D,Bool}) where {D}
    D >= 3 || throw(ArgumentError("Cell positions must have at least 3 Cartesian dimensions"))
    expected = _default_periodicity(Val(D))
    periodicity == expected ||
        throw(
            ArgumentError(
                "Current CellBase support requires the first 3 dimensions to be periodic and any extra dimensions to be non-periodic",
            ),
        )
end

function _validate_cell_shape(l::Lattice, positions::Matrix, periodicity)
    size(cellmat(l)) == (3, 3) ||
        throw(ArgumentError("Current CellBase support requires a 3x3 lattice matrix"))
    size(positions, 1) == length(periodicity) ||
        throw(
            ArgumentError(
                "Position dimension $(size(positions, 1)) does not match periodicity length $(length(periodicity))",
            ),
        )
end

function _require_three_cartesian_dimensions(cell::Cell, context::AbstractString)
    size(positions(cell), 1) == 3 ||
        throw(
            ArgumentError(
                "$context only supports cells with 3 Cartesian coordinates per atom",
            ),
        )
end

function _warn_dropped_auxiliary_dimensions(cell::Cell, context::AbstractString)
    if size(positions(cell), 1) > 3
        @warn "$context drops auxiliary dimensions beyond the first 3 Cartesian coordinates"
    end
end

"""
    Cell3D{T}

Type alias for 3D cells, providing backward compatibility with code that uses `Cell{T}` annotations.

# Examples
```julia
# Old code (still works)
cell3d::Cell3D{Float64} = bulk("Cu")

# Equivalent to new type
cell::Cell{Float64,3} = bulk("Cu")
```
"""
const Cell3D{T} = Cell{T,3}

"""
    Cell(l::Lattice, symbols, positions) where T

Construct a Cell type from arrays
"""
function Cell(
    l::Lattice,
    symbols::Vector{Symbol},
    positions::Matrix,
    arrays::Dict{Symbol,Any},
    metadata::Dict{Symbol,Any};
    periodicity=nothing,
)
    @assert length(symbols) == size(positions, 2)
    D = size(positions, 1)
    periodicity_tuple = _normalize_periodicity(periodicity, Val(D))
    _validate_periodicity(periodicity_tuple)
    _validate_cell_shape(l, positions, periodicity_tuple)
    Cell{eltype(positions),D}(l, symbols, positions, arrays, metadata, periodicity_tuple)
end

function Cell(l::Lattice, symbols::Vector{Symbol}, positions::Matrix; periodicity=nothing)
    Cell(
        l,
        symbols,
        positions,
        Dict{Symbol,Any}(),
        Dict{Symbol,Any}();
        periodicity,
    )
end

"""
    Cell(lat::Lattice, numbers::Vector{Int}, positions)

Constructure the Cell type from lattice, positions and numbers
"""
function Cell(lat::Lattice, numbers::Vector{T}, positions::Matrix; periodicity=nothing) where {T<:Real}
    species = [Symbol(elements[i].symbol) for i in numbers]
    @assert length(numbers) == size(positions, 2)
    Cell(lat, species, positions; periodicity)
end


"""
    Cell(lat::Lattice, numbers::Vector{Int}, positions::Vector)

Constructure the Cell type from lattice, positions and numbers
"""
function Cell(lat::Lattice, species_id, positions::Vector; periodicity=nothing)
    @assert length(species_id) == length(positions)
    posmat = zeros(length(positions[1]), length(positions))
    for (i, vec) in enumerate(positions)
        posmat[:, i] = vec
    end
    Cell(lat, species_id, posmat; periodicity)
end


"""
    clip(s::Cell, mask::AbstractVector)

Clip a structure with a given indexing array
"""
function clip(cell::Cell{T,N}, mask::AbstractVector) where {T,N}
    new_pos = positions(cell)[:, mask]
    new_symbols = species(cell)[mask]
    # Clip any additional arrays
    new_array = Dict{Symbol,Any}()
    for (key, array) in pairs(cell.arrays)
        new_array[key] = selectdim(array, ndims(array), mask)
    end
    Cell(lattice(cell), new_symbols, new_pos, new_array, cell.metadata; periodicity=periodicity(cell))
end

Base.getindex(cell::Cell, i::AbstractVector) = clip(cell, i)


"""
    _compute_sort_indices(cell::Cell, by::Symbol) -> Vector{Int}

Compute sort indices based on the sorting criterion.
"""
function _compute_sort_indices(cell::Cell, by::Symbol)
    if by in (:symbol, :species)
        # Sort by element symbol
        return sortperm(species(cell))
    elseif by in (:number, :z, :atomic_number)
        # Sort by atomic number
        z = atomic_numbers(cell)
        return sortperm(z)
    elseif by == :position
        # Sort by position (x, then y, then z)
        pos = positions(cell)
        return sortperm(1:natoms(cell), by=i -> (pos[1, i], pos[2, i], pos[3, i]))
    else
        throw(ArgumentError("Unknown sort criterion: $by. Use :symbol, :number, or :position"))
    end
end


"""
    Base.sort(cell::Cell; by::Symbol=:symbol) -> Cell
    Base.sort(cell::Cell, tags::AbstractVector) -> Cell

Return a new Cell with sorted atomic order.

# Arguments
- `cell`: Input Cell instance
- `by`: Sorting criterion (default: :symbol)
  - `:symbol` or `:species`: Sort by element symbol
  - `:number`, `:z`, or `:atomic_number`: Sort by atomic number
  - `:position`: Sort by position (x, then y, then z)
- `tags`: Custom vector of values for sorting (one per atom). Lower values come first.

# Returns
A new `Cell` instance with sorted atoms.

# Examples
```julia
# Sort by element symbol (default)
sort(cell)
sort(cell, by=:symbol)

# Sort by atomic number
sort(cell, by=:number)

# Sort by position
sort(cell, by=:position)

# Sort by custom tags
sort(cell, [2, 1, 4, 3])
```
"""
Base.sort(cell::Cell; by::Symbol=:symbol) = cell[_compute_sort_indices(cell, by)]

function Base.sort(cell::Cell, tags::AbstractVector)
    if length(tags) != natoms(cell)
        throw(ArgumentError("Tags vector length ($(length(tags))) must match number of atoms ($(natoms(cell)))"))
    end
    return cell[sortperm(tags)]
end


"""
    Base.sort!(cell::Cell; by::Symbol=:symbol) -> Cell
    Base.sort!(cell::Cell, tags::AbstractVector) -> Cell

In-place sort atoms in a Cell.

Updates positions, species, and all additional arrays in place.

!!! warning
    This will modify shared data. Use with caution.

# Arguments
- `cell`: Cell instance to sort (modified in place)
- `by`: Sorting criterion (default: :symbol)
  - `:symbol` or `:species`: Sort by element symbol
  - `:number`, `:z`, or `:atomic_number`: Sort by atomic number
  - `:position`: Sort by position (x, then y, then z)
- `tags`: Custom vector of values for sorting (one per atom)

# Returns
The modified `cell` instance.

# Examples
```julia
# Sort by element symbol
sort!(cell)
sort!(cell, by=:symbol)

# Sort by atomic number
sort!(cell, by=:number)

# Sort by custom tags
sort!(cell, [2, 1, 4, 3])
```
"""
function Base.sort!(cell::Cell; by::Symbol=:symbol)
    idx = _compute_sort_indices(cell, by)
    _sort_impl!(cell, idx)
end

function Base.sort!(cell::Cell, tags::AbstractVector)
    if length(tags) != natoms(cell)
        throw(ArgumentError("Tags vector length ($(length(tags))) must match number of atoms ($(natoms(cell)))"))
    end
    idx = sortperm(tags)
    _sort_impl!(cell, idx)
end


"""
    _sort_impl!(cell::Cell, idx::Vector{Int}) -> Cell

Internal implementation of in-place sorting with validation.
Validates that all arrays have consistent sizes before sorting.
"""
function _sort_impl!(cell::Cell, idx::Vector{Int})
    n = natoms(cell)

    # Validate positions array
    if size(cell.positions, 2) != n
        throw(ArgumentError("Positions array size $(size(cell.positions, 2)) doesn't match natoms($n)"))
    end

    # Validate additional arrays
    for (name, arr) in pairs(cell.arrays)
        if ndims(arr) >= 1 && size(arr, ndims(arr)) != n
            throw(ArgumentError("Array '$name' size $(size(arr, ndims(arr))) doesn't match natoms($n)"))
        end
    end

    # Perform the sort
    species(cell) .= species(cell)[idx]
    cell.positions .= cell.positions[:, idx]
    for array in values(cell.arrays)
        if ndims(array) >= 1
            array .= selectdim(array, ndims(array), idx)
        end
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

periodicity(cell::Cell) = cell.periodicity


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
sposarray(structure::Cell{T,N}) where {T,N} =
    [SVector{N,T}(x) for x in eachcol(positions(structure))]

"""
    species(structure::Cell)

Return a Vector of species names.
"""
species(structure::Cell) = structure.symbols

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
function array(structure::Cell, arrayname::Symbol)
    if !haskey(structure.arrays, arrayname)
        available = join(keys(structure.arrays), ", ")
        throw(ArgumentError("Array '$arrayname' not found in Cell. Available arrays: $available"))
    end
    return structure.arrays[arrayname]
end

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

Return fractional positions of the periodic subspace.

For ordinary 3D cells this returns a `3 x N` matrix of fractional coordinates.
For hyper cells (`D > 3`), this still returns `3 x N` and only covers the first
three periodic Cartesian coordinates.
"""
function get_scaled_positions(cell::Cell)
    rec_cellmat(lattice(cell)) * @view positions(cell)[1:3, :]
end

"""
    set_scaled_positions!(cell::Cell, scaled::Matrix)

Set fractional coordinates for the periodic subspace.

`scaled` must have shape `3 x N`. For hyper cells, only the first three
periodic Cartesian coordinates are updated; auxiliary coordinates are preserved.
"""
function set_scaled_positions!(cell::Cell, scaled::Matrix)
    size(scaled, 1) == 3 ||
        throw(ArgumentError("Scaled positions must have exactly 3 rows, got $(size(scaled, 1))"))
    size(scaled, 2) == natoms(cell) ||
        throw(
            ArgumentError(
                "Scaled positions must have $(natoms(cell)) columns, got $(size(scaled, 2))",
            ),
        )
    @view(cell.positions[1:3, :]) .= cellmat(cell) * @view(scaled[1:3, :])
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

function wrap!(vec::AbstractVector, c::Cell)
    length(vec) >= 3 || throw(ArgumentError("Vector must have at least 3 components"))
    frac = rec_cellmat(lattice(c)) * @view(vec[1:3])
    frac .-= floor.(frac)
    vec[1:3] .= cellmat(c) * frac
    vec
end


"""
    wrapped_spos(cell)

Return a static array of wrapped positons.
"""
function wrapped_spos(cell::Cell{T,3}) where {T}
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
    for i in eachindex(cell.positions)
        cell.positions[i] += (rand() - 0.5) * 2 * amp
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
        line =
            @sprintf "%4s  %10.5f  %10.5f  %10.5f" symbol posmat[1, i] posmat[2, i] posmat[
                3,
                i,
            ]
        if n_dimensions(s) > 3
            extra = join([@sprintf("%10.5f ", i) for i in posmat[4:end, i]], "")
            line = line * "  ($extra  )"
        end
        println(io, line)
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
    extra_vecs = max(size(positions(cell), 1) - 3, 0) > 0 ? zeros(size(positions(cell), 1) - 3, nn * nn) : nothing
    pos = sposarray(cell)
    dmat = zeros(nn, nn)
    ivec = 0
    for i = 1:nn
        for j = i+1:nn
            ivec += 1
            vecs[:, ivec] .= pos[j][1:3] .- pos[i][1:3]
            if !isnothing(extra_vecs)
                extra_vecs[:, ivec] .= pos[j][4:end] .- pos[i][4:end]
            end
        end
    end
    # Apply minimum image conventions
    vmic, dmic = mic(lattice(cell), @view(vecs[:, 1:ivec]))
    if !isnothing(extra_vecs)
        for idx = 1:ivec
            dmic[idx] = sqrt(dmic[idx]^2 + sum(abs2, @view extra_vecs[:, idx]))
        end
    end
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
    make_supercell(structure, tmat; wrap=false)
end


"""
    make_supercell(cell::Cell{T,D}, P::AbstractMatrix{<:Integer};
                   wrap=true, order="cell-major") where {T,D}

Generate a supercell by applying a general transformation matrix `P` to the input cell.

The transformation is described by a 3×3 integer matrix **P** (column-major format):
new_cell_matrix = old_cell_matrix * P

This is a generalization of the diagonal make_supercell that works with any
integer transformation matrix, including non-diagonal (shear) transformations.

# Arguments
- `cell`: Input Cell instance
- `P`: 3×3 integer transformation matrix (column-major, cell vectors are columns)
- `wrap`: Whether to wrap positions back into the new cell (default: true)
- `order`: Ordering of atoms - "cell-major" (default) or "atom-major"
  - "cell-major": [atom1_cell1, atom2_cell1, ..., atom1_cell2, atom2_cell2, ...]
  - "atom-major": [atom1_cell1, atom1_cell2, ..., atom2_cell1, atom2_cell2, ...]

# Returns
A new `Cell` instance representing the supercell.

# Examples
```julia
# 2x2x2 supercell via transformation matrix
P = [2 0 0; 0 2 0; 0 0 2]
supercell = make_supercell(cell, P)

# Create a 2x1x1 supercell with shear
P = [2 1 0; 0 1 0; 0 0 1]
supercell = make_supercell(cell, P)

# Different atom ordering
supercell = make_supercell(cell, P, order="atom-major")
```
"""
function make_supercell(cell::Cell{T,D}, P::AbstractMatrix{<:Real};
                        wrap=true, order="cell-major") where {T,D}
    # Validate P is 3×3 matrix
    if size(P) != (3, 3)
        throw(ArgumentError("Transformation matrix P must be 3×3, got $(size(P))"))
    end

    # Ensure integer matrix
    P_int = round.(Int, P)
    if any(P_int .!= P)
        throw(ArgumentError("Transformation matrix P must contain only integers"))
    end

    # Compute new cell: new_cell = old_cell * P (column-major)
    old_cell = cellmat(cell)
    new_cell = old_cell * P_int

    # Compute determinant to get expected number of atoms
    det_P = abs(det(P_int))
    n_target = round(Int, det_P * natoms(cell))

    # Generate all lattice points in the supercell
    # We need integer solutions to P * n = i where i are the lattice indices
    lattice_points_frac = _lattice_points_in_supercell(P_int)
    N = size(lattice_points_frac, 2)

    # Convert fractional shifts to Cartesian
    lattice_shifts = old_cell * P_int * lattice_points_frac  # 3 × N matrix

    # Get current positions and species
    current_pos = positions(cell)
    ns = natoms(cell)

    # Create new positions based on order parameter
    if order == "cell-major"
        # [atom1_shift1, atom2_shift1, ..., atom1_shift2, atom2_shift2, ...]
        new_pos = zeros(T, D, N * ns)
        for i = 1:N
            shift = lattice_shifts[:, i]
            for j = 1:ns
                idx = j + (i - 1) * ns
                new_pos[1:3, idx] .= current_pos[1:3, j] .+ shift
                if D > 3
                    new_pos[4:end, idx] .= current_pos[4:end, j]
                end
            end
        end
    elseif order == "atom-major"
        # [atom1_shift1, atom1_shift2, ..., atom2_shift1, atom2_shift2, ...]
        new_pos = zeros(T, D, N * ns)
        for j = 1:ns
            for i = 1:N
                idx = i + (j - 1) * N
                shift = lattice_shifts[:, i]
                new_pos[1:3, idx] .= current_pos[1:3, j] .+ shift
                if D > 3
                    new_pos[4:end, idx] .= current_pos[4:end, j]
                end
            end
        end
    else
        throw(ArgumentError("Invalid order '$order'. Use 'cell-major' or 'atom-major'."))
    end

    # Tile the species
    new_spec = repeat(species(cell), N)

    # Create new cell
    new_cell_obj = Cell(Lattice(new_cell), new_spec, new_pos)

    # Copy over additional arrays
    for (name, arr) in pairs(cell.arrays)
        if ndims(arr) >= 1
            n_dim = ndims(arr)
            if order == "cell-major"
                # Tile along the last dimension (atom dimension) - repeat each element N times
                # For cell-major: [a1, a2, a3] with N=2 -> [a1, a1, a2, a2, a3, a3]
                new_shape = ntuple(d -> d == n_dim ? size(arr, d) * N : size(arr, d), n_dim)
                new_arr = similar(arr, new_shape)
                for i = 1:ns
                    for j = 1:N
                        idx_out = i + (j - 1) * ns
                        selectdim(new_arr, n_dim, idx_out) .= selectdim(arr, n_dim, i)
                    end
                end
            else  # atom-major
                # For atom-major: [a1, a2, a3] with N=2 -> [a1, a2, a3, a1, a2, a3]
                shape = size(arr)
                new_shape = ntuple(d -> d == n_dim ? shape[d] * N : shape[d], n_dim)
                new_arr = similar(arr, new_shape)
                for j = 1:N
                    for i = 1:ns
                        idx_out = i + (j - 1) * ns
                        idx_in = i
                        selectdim(new_arr, n_dim, idx_out) .= selectdim(arr, n_dim, idx_in)
                    end
                end
            end
            new_cell_obj.arrays[name] = new_arr
        end
    end

    # Copy metadata
    new_cell_obj.metadata = copy(cell.metadata)

    # Verify atom count
    if n_target != natoms(new_cell_obj)
        throw(ErrorException("Number of atoms in supercell: $(natoms(new_cell_obj)), expected: $n_target"))
    end

    # Wrap positions if requested
    if wrap
        wrap!(new_cell_obj)
    end

    return new_cell_obj
end


"""
    _lattice_points_in_supercell(P::Matrix{Int}) -> Matrix{Float64}

Generate all lattice points (in fractional coordinates of the supercell) within
a supercell defined by transformation matrix P.

Returns a 3×N matrix where each column is a lattice point in fractional
coordinates of the new cell.

This uses a brute-force search over a bounding box determined by the transformation
matrix columns.
"""
function _lattice_points_in_supercell(P::Matrix{Int})
    # The number of lattice points equals the absolute determinant of P
    det_P = round(Int, abs(det(Float64.(P))))

    # Determine search bounds based on the transformation matrix
    # For a transformation P, lattice points n satisfy: P⁻¹n ∈ [0,1)³
    # We need bounds large enough to cover the parallelepiped defined by P
    # Use sum of absolute values of each row to get conservative bounds
    row_sums = [sum(abs.(P[i, :])) for i in 1:3]
    max_bound = maximum(row_sums)
    # Add some margin for skewed cells
    bounds = -max_bound:max_bound

    P_inv = inv(Float64.(P))
    tol = 1e-10

    # Pre-allocate with expected size (may be larger, will trim later)
    points = Vector{Float64}[]
    sizehint!(points, det_P * 2)

    for i in bounds
        for j in bounds
            for k in bounds
                n = Float64[i, j, k]
                frac = P_inv * n

                # Check if fractional coordinates are in [0, 1)
                if all(x -> -tol <= x < 1.0 - tol, frac)
                    # Normalize to [0, 1) range
                    frac_norm = frac .- floor.(frac .+ tol)
                    push!(points, frac_norm)
                end
            end
        end
    end

    # Remove duplicates (shouldn't happen but just in case)
    unique_points = unique(points)

    npoints = length(unique_points)

    if npoints != det_P
        # If we didn't find the right number of points, try with larger bounds
        if npoints < det_P
            # Recursive call with larger bounds
            return _lattice_points_in_supercell_extended(P, max_bound * 2)
        end
        @warn "Expected $det_P lattice points, found $npoints"
    end

    if npoints == 0
        return zeros(3, 0)
    end

    # Convert to matrix efficiently
    points_matrix = Matrix{Float64}(undef, 3, npoints)
    for i in 1:npoints
        points_matrix[:, i] .= unique_points[i]
    end

    # Sort for consistent ordering (by x, then y, then z)
    idx = sortperm(1:npoints, by=i -> (points_matrix[1, i], points_matrix[2, i], points_matrix[3, i]))

    return points_matrix[:, idx]
end


"""
    _lattice_points_in_supercell_extended(P::Matrix{Int}, bound::Int) -> Matrix{Float64}

Extended search for lattice points with larger bounds. Called recursively if the
initial search doesn't find enough points.
"""
function _lattice_points_in_supercell_extended(P::Matrix{Int}, bound::Int)
    det_P = round(Int, abs(det(Float64.(P))))
    bounds = -bound:bound

    P_inv = inv(Float64.(P))
    tol = 1e-10

    points = Vector{Float64}[]
    sizehint!(points, det_P * 2)

    for i in bounds
        for j in bounds
            for k in bounds
                n = Float64[i, j, k]
                frac = P_inv * n

                if all(x -> -tol <= x < 1.0 - tol, frac)
                    frac_norm = frac .- floor.(frac .+ tol)
                    push!(points, frac_norm)
                end
            end
        end
    end

    unique_points = unique(points)
    npoints = length(unique_points)

    if npoints != det_P && bound < 1000
        # Try again with even larger bounds
        return _lattice_points_in_supercell_extended(P, bound * 2)
    end

    if npoints == 0
        return zeros(3, 0)
    end

    points_matrix = Matrix{Float64}(undef, 3, npoints)
    for i in 1:npoints
        points_matrix[:, i] .= unique_points[i]
    end

    idx = sortperm(1:npoints, by=i -> (points_matrix[1, i], points_matrix[2, i], points_matrix[3, i]))
    return points_matrix[:, idx]
end


"""
    repeat(cell::Cell{T,D}, reps::NTuple{3,Int}; kwargs...) where {T,D}
    repeat(cell::Cell{T,D}, a::Int, b::Int, c::Int; kwargs...) where {T,D}
    repeat(cell::Cell{T,D}, n::Int; kwargs...) where {T,D}

Repeat the cell along the lattice vector directions.

This creates a supercell by repeating the original cell the specified number of
times along each lattice direction. It is equivalent to `make_supercell` with a
diagonal transformation matrix.

# Arguments
- `cell`: Input Cell instance
- `reps`: Tuple of (a, b, c) repeat factors along each lattice vector
- `a, b, c`: Individual repeat factors
- `n`: Single integer for uniform repeat (equivalent to (n, n, n))
- `kwargs`: Additional arguments passed to `make_supercell` (e.g., `order`, `wrap`)

# Returns
A new `Cell` instance representing the repeated supercell.

# Examples
```julia
# Repeat 2x2x2
repeat(cell, 2, 2, 2)
repeat(cell, (2, 2, 2))
repeat(cell, 2)  # equivalent to (2, 2, 2)

# Different ordering
repeat(cell, 2, 2, 2, order="atom-major")

# Don't wrap positions
repeat(cell, 2, 2, 2, wrap=false)
```
"""
function repeat(cell::Cell{T,D}, reps::NTuple{3,Int}; kwargs...) where {T,D}
    a, b, c = reps

    # Check for invalid repeats along undefined lattice vectors
    old_cell = cellmat(cell)
    for (rep, vec) in zip(reps, eachcol(old_cell))
        if rep != 1 && norm(vec) < 1e-10
            throw(ArgumentError("Cannot repeat along undefined lattice vector"))
        end
    end

    # Create diagonal transformation matrix
    P = diagm([a, b, c])

    return make_supercell(cell, P; kwargs...)
end

repeat(cell::Cell{T,D}, a::Int, b::Int, c::Int; kwargs...) where {T,D} =
    repeat(cell, (a, b, c); kwargs...)

repeat(cell::Cell{T,D}, n::Int; kwargs...) where {T,D} =
    repeat(cell, (n, n, n); kwargs...)


"""
    masses(cell::Cell)

Return a vector of atomic masses for all atoms in the cell (in atomic mass units).
"""
function masses(cell::Cell)
    return [elements[sym].atomic_mass for sym in species(cell)]
end


"""
    _get_rotation_center(cell::Cell{T,D}, center::Union{NTuple{3,<:Real},Symbol}) where {T,D}

Get the rotation center coordinates based on the center specification.

# Arguments
- `cell`: Cell instance
- `center`: Center specification - can be:
  - Tuple of coordinates: (0.0, 0.0, 0.0)
  - Symbol: :center_of_mass, :com, :center_of_positions, :cop, :center_of_cell, :coc

# Returns
Vector of 3 coordinates for the rotation center
"""
function _get_rotation_center(cell::Cell{T,D}, center::Union{NTuple{3,<:Real},Symbol,AbstractVector}) where {T,D}
    if isa(center, Symbol)
        if center in (:center_of_mass, :com)
            # Center of mass
            m = masses(cell)
            return vec(sum(@view(positions(cell)[1:3, :]) .* reshape(m, 1, :), dims=2) ./ sum(m))
        elseif center in (:center_of_positions, :cop)
            # Center of positions (geometric center)
            return vec(mean(@view(positions(cell)[1:3, :]), dims=2))
        elseif center in (:center_of_cell, :coc)
            # Center of the unit cell
            return vec(sum(cellmat(lattice(cell)), dims=2) ./ 2)
        else
            throw(ArgumentError("Unknown center specification: $center. Use :center_of_mass, :center_of_positions, or :center_of_cell"))
        end
    else
        length(center) >= 3 ||
            throw(ArgumentError("Rotation center must have at least 3 coordinates"))
        return collect(Float64, center[1:3])
    end
end


"""
    rotate(cell::Cell{T,D}, angle::Real, axis::Union{AbstractVector,Symbol,String};
           center=(0.0, 0.0, 0.0), rotate_cell::Bool=false) where {T,D}
    rotate(cell::Cell{T,D}, v1::AbstractVector, v2::AbstractVector;
           center=(0.0, 0.0, 0.0), rotate_cell::Bool=false) where {T,D}

Rotate atoms around an axis by an angle, or align one vector to another.

The rotation can optionally include the cell vectors, allowing for rotation of
the entire crystal structure including its periodic boundaries.

# Arguments
- `cell`: Input Cell instance
- `angle`: Rotation angle in degrees
- `axis`: Rotation axis as vector or symbol/string ("x", "y", "z", "-x", etc.)
- `v1`, `v2`: Vectors to rotate from and to (aligns v1 with v2)
- `center`: Fixed center for rotation. Can be:
  - Tuple/Vector of coordinates: (0.0, 0.0, 0.0)
  - `:center_of_mass` or `:com`
  - `:center_of_positions` or `:cop`
  - `:center_of_cell` or `:coc`
- `rotate_cell`: If true, also rotate the lattice vectors (default: false)

# Returns
A new `Cell` instance with rotated positions (and optionally cell).

# Examples
```julia
# Rotate 90 degrees around z-axis
rotate(cell, 90, "z")
rotate(cell, 90, [0, 0, 1])
rotate(cell, 90, :z)

# Align x-axis with y-axis
rotate(cell, [1, 0, 0], [0, 1, 0])

# Rotate around center of mass
rotate(cell, 45, "z", center=:center_of_mass)

# Rotate cell as well
rotate(cell, 90, "z", rotate_cell=true)

# Rotate around specific point
rotate(cell, 30, [1, 1, 1], center=(1.0, 2.0, 3.0))
```
"""
function rotate(cell::Cell{T,D}, angle::Real, axis::Union{AbstractVector,Symbol,String};
                center=(0.0, 0.0, 0.0), rotate_cell::Bool=false) where {T,D}
    # Convert string/symbol axis to vector
    ax = if isa(axis, AbstractVector)
        axis
    else
        axis_from_string(axis)
    end

    # Get rotation center
    c = _get_rotation_center(cell, center)

    # Compute rotation matrix
    R = rotation_matrix(angle, ax)

    # Create new cell
    new_cell = deepcopy(cell)

    # Rotate positions: p' = R * (p - c) + c
    pos = positions(new_cell)
    for i = 1:natoms(new_cell)
        pos[1:3, i] .= R * (pos[1:3, i] .- c) .+ c
    end

    # Optionally rotate lattice vectors
    if rotate_cell
        old_cell_mat = cellmat(lattice(new_cell))
        new_cell_mat = R * old_cell_mat
        set_cellmat!(lattice(new_cell), new_cell_mat)
    end

    return new_cell
end


function rotate(cell::Cell{T,D}, v1::AbstractVector, v2::AbstractVector;
                center=(0.0, 0.0, 0.0), rotate_cell::Bool=false) where {T,D}
    # Compute rotation matrix that aligns v1 with v2
    R = rotation_matrix_align(v1, v2)

    # Get rotation center
    c = _get_rotation_center(cell, center)

    # Create new cell
    new_cell = deepcopy(cell)

    # Rotate positions: p' = R * (p - c) + c
    pos = positions(new_cell)
    for i = 1:natoms(new_cell)
        pos[1:3, i] .= R * (pos[1:3, i] .- c) .+ c
    end

    # Optionally rotate lattice vectors
    if rotate_cell
        old_cell_mat = cellmat(lattice(new_cell))
        new_cell_mat = R * old_cell_mat
        set_cellmat!(lattice(new_cell), new_cell_mat)
    end

    return new_cell
end


"""
    rotate!(cell::Cell{T,D}, angle::Real, axis::Union{AbstractVector,Symbol,String};
            center=(0.0, 0.0, 0.0), rotate_cell::Bool=false) where {T,D}
    rotate!(cell::Cell{T,D}, v1::AbstractVector, v2::AbstractVector;
            center=(0.0, 0.0, 0.0), rotate_cell::Bool=false) where {T,D}

In-place version of rotate. Modifies the cell directly.

!!! warning
    This modifies the cell in place. Use with caution when data is shared.
"""
function rotate!(cell::Cell{T,D}, angle::Real, axis::Union{AbstractVector,Symbol,String};
                 center=(0.0, 0.0, 0.0), rotate_cell::Bool=false) where {T,D}
    # Convert string/symbol axis to vector
    ax = if isa(axis, AbstractVector)
        axis
    else
        axis_from_string(axis)
    end

    # Get rotation center
    c = _get_rotation_center(cell, center)

    # Compute rotation matrix
    R = rotation_matrix(angle, ax)

    # Rotate positions: p' = R * (p - c) + c
    pos = positions(cell)
    for i = 1:natoms(cell)
        pos[1:3, i] .= R * (pos[1:3, i] .- c) .+ c
    end

    # Optionally rotate lattice vectors
    if rotate_cell
        old_cell_mat = cellmat(lattice(cell))
        new_cell_mat = R * old_cell_mat
        set_cellmat!(lattice(cell), new_cell_mat)
    end

    return cell
end


function rotate!(cell::Cell{T,D}, v1::AbstractVector, v2::AbstractVector;
                 center=(0.0, 0.0, 0.0), rotate_cell::Bool=false) where {T,D}
    # Compute rotation matrix that aligns v1 with v2
    R = rotation_matrix_align(v1, v2)

    # Get rotation center
    c = _get_rotation_center(cell, center)

    # Rotate positions: p' = R * (p - c) + c
    pos = positions(cell)
    for i = 1:natoms(cell)
        pos[1:3, i] .= R * (pos[1:3, i] .- c) .+ c
    end

    # Optionally rotate lattice vectors
    if rotate_cell
        old_cell_mat = cellmat(lattice(cell))
        new_cell_mat = R * old_cell_mat
        set_cellmat!(lattice(cell), new_cell_mat)
    end

    return cell
end


"""
    fingerprint(s::Cell; dmat=distance_matrix(s), weighted=true, cut_bl=3.0)

Computed the fingerprint vector based on simple sorted pair-wise distances.
NOTE: Does not work for single atom cell!!
"""
function fingerprint(s::Cell; dmat=distance_matrix(s), weighted=true, cut_bl=3.0)
    # Validate input
    nn, _ = size(dmat)
    if nn < 2
        throw(ArgumentError("fingerprint requires at least 2 atoms, got $nn"))
    end

    # Check if there are any non-zero distances
    nonzero_dists = [d for d in dmat if d > 0.0]
    if isempty(nonzero_dists)
        throw(ArgumentError("fingerprint requires at least one non-zero distance"))
    end

    # Cut off distance based on minimum bond length
    cut_bl = minimum(nonzero_dists) * cut_bl
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


"""
    add_dimensions(cell::Cell, ndims, dscale=3.0)

Add additional dimensions to the cell
"""
function add_dimensions(cell::Cell, ndims, dscale=3.0)
    pos = positions(cell)
    hyperpos = vcat(pos, (rand(ndims, size(pos, 2)) .- 0.5) .* 2dscale)
    out = Cell(lattice(cell), species(cell), hyperpos)
    out.metadata = cell.metadata
    out
end


"""
    remove_dimensions(cell::Cell)

Remove the extra dimensions in the cell
"""
function remove_dimensions(cell::Cell)
    pos = positions(cell)[1:3, :]
    out = Cell(lattice(cell), species(cell), pos)
    out.metadata = cell.metadata
    out
end


"""
    versioninfo([io::IO=stdout])

Print version and dependency information for CellBase.jl, useful for debugging and bug reports.

# Examples
```julia
versioninfo()  # Print to stdout
versioninfo(stderr)  # Print to stderr
```
"""
function versioninfo(io::IO=stdout)
    println(io, "CellBase.jl Version Information")
    println(io, "================================")
    println(io, "Julia Version: ", VERSION)
    println(io, "CellBase.jl: ", _pkg_version(CellBase))
    println(io, "")
    println(io, "Key Dependencies:")
    println(io, "  AtomsBase: ", _pkg_version(AB))
    println(io, "  Spglib: ", _pkg_version(Spglib))
    println(io, "  PeriodicTable: ", _pkg_version(PeriodicTable))
    println(io, "  StaticArrays: ", _pkg_version(StaticArrays))
    println(io, "")
    println(io, "AtomsBase compatibility: ", AtomsBase_version_compat())
    println(io, "Hyperdimensional support: enabled (D parameter in Cell{T,D})")
    println(io, "Cell3D{T} alias: available for backward compatibility")
    nothing
end

function _pkg_version(mod::Module)
    try
        path = pathof(mod)
        if path === nothing
            return "unknown"
        end
        pkg_dir = dirname(dirname(path))
        proj_file = joinpath(pkg_dir, "Project.toml")
        if isfile(proj_file)
            for line in eachline(proj_file)
                if startswith(line, "version = ")
                    version_str = strip(line[11:end])
                    return replace(version_str, "\"" => "")
                end
            end
        end
        return "unknown"
    catch
        return "unknown"
    end
end

function AtomsBase_version_compat()
    try
        # Check if we're using AtomsBase 0.4 or later by checking for position method with idx parameter
        if hasmethod(AB.position, (Cell{Float64,3}, Int))
            return "0.4+ (indexable with idx)"
        else
            return "0.3 (legacy)"
        end
    catch
        return "not available"
    end
end
