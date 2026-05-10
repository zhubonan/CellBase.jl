#=
Code for reading/writing ABACUS STRU files
=#

using Printf
import Unitful: ustrip

const BOHR_TO_ANGSTROM = 0.529177210544 # https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
const ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM

"""
    read_stru(io::IO)

Read a `Cell` object from an ABACUS STRU file.
"""
function read_stru(io::IO)
    lines = readlines(io)
    read_stru(lines)
end

"""
    read_stru(fname::AbstractString)

Read a `Cell` object from an ABACUS STRU file.
"""
function read_stru(fname::AbstractString)
    open(fname) do io
        read_stru(io)
    end
end

"""
    read_stru(lines::Vector{<:AbstractString})

Read a `Cell` object from an array of STRU file lines.
"""
function read_stru(lines::Vector{<:AbstractString})
    lattice_constant = 1.0
    lattice_vectors = Matrix{Float64}(undef, 3, 3)
    has_lattice_vectors = false
    species = Symbol[]
    positions = Matrix{Float64}(undef, 3, 0)
    coord_type = "Direct"
    metadata = Dict{Symbol,Any}()

    i = 1
    while i <= length(lines)
        line = strip(lines[i])

        if isempty(line) || startswith(line, "//")
            i += 1
            continue
        end

        if startswith(line, "ATOMIC_SPECIES")
            i += 1
            while i <= length(lines)
                l = strip(lines[i])
                if isempty(l) ||
                   startswith(l, "//") ||
                   !in(l[1], 'A':'Z') && !in(l[1], 'a':'z')
                    break
                end
                i += 1
            end
        elseif startswith(line, "NUMERICAL_ORBITAL")
            i += 1
            while i <= length(lines)
                l = strip(lines[i])
                if isempty(l) ||
                   startswith(l, "//") ||
                   contains(l, "LATTICE") ||
                   contains(l, "ATOMIC")
                    break
                end
                i += 1
            end
        elseif startswith(line, "LATTICE_CONSTANT")
            i += 1
            while i <= length(lines)
                l = strip(lines[i])
                if isempty(l) || startswith(l, "//")
                    i += 1
                    continue
                end
                lattice_constant = parse(Float64, split(l)[1])
                metadata[:lattice_constant_bohr] = lattice_constant
                i += 1
                break
            end
        elseif startswith(line, "LATTICE_VECTORS")
            has_lattice_vectors = true
            i += 1
            for row = 1:3
                while i <= length(lines)
                    l = strip(lines[i])
                    if isempty(l) || startswith(l, "//")
                        i += 1
                        continue
                    end
                    lattice_vectors[row, :] = parse.(Float64, split(l))
                    i += 1
                    break
                end
            end
        elseif startswith(line, "LATTICE_PARAMETERS")
            i += 1
            while i <= length(lines)
                l = strip(lines[i])
                if isempty(l) || startswith(l, "//") || contains(l, "ATOMIC")
                    break
                end
                i += 1
            end
        elseif startswith(line, "ATOMIC_POSITIONS")
            i += 1
            while i <= length(lines)
                l = strip(lines[i])
                if isempty(l) || startswith(l, "//")
                    i += 1
                    continue
                end
                coord_type = l
                i += 1
                break
            end

            species = Symbol[]
            pos_list = Vector{Float64}[]

            while i <= length(lines)
                l = strip(lines[i])
                if isempty(l) || startswith(l, "//")
                    i += 1
                    continue
                end

                if contains(l, "LATTICE") ||
                   contains(l, "NUMERICAL") ||
                   contains(l, "ATOMIC_SPECIES")
                    break
                end

                elem = Symbol(split(l)[1])
                i += 1

                while i <= length(lines)
                    l = strip(lines[i])
                    if isempty(l) || startswith(l, "//")
                        i += 1
                        continue
                    end
                    break
                end

                i += 1

                while i <= length(lines)
                    l = strip(lines[i])
                    if isempty(l) || startswith(l, "//")
                        i += 1
                        continue
                    end
                    break
                end

                natoms = parse(Int, split(lines[i])[1])
                i += 1

                for j = 1:natoms
                    while i <= length(lines)
                        l = strip(lines[i])
                        if isempty(l) || startswith(l, "//")
                            i += 1
                            continue
                        end
                        break
                    end

                    tokens = split(lines[i])
                    pos = parse.(Float64, tokens[1:3])
                    push!(pos_list, pos)
                    push!(species, elem)
                    i += 1
                end
            end

            if !isempty(pos_list)
                positions = hcat(pos_list...)
            end
            break
        end
        i += 1
    end

    if !has_lattice_vectors
        throw(ErrorException("LATTICE_VECTORS section not found in STRU file"))
    end

    lattice_ang = lattice_vectors * lattice_constant * BOHR_TO_ANGSTROM

    pos_ang = if startswith(lowercase(coord_type), "d")
        lattice_ang * positions
    elseif startswith(lowercase(coord_type), "c") &&
           contains(lowercase(coord_type), "angstrom")
        positions
    else
        positions * lattice_constant * BOHR_TO_ANGSTROM
    end

    cell = Cell(Lattice(lattice_ang), species, pos_ang)
    cell.metadata = metadata
    cell
end

"""
    write_stru(io::IO, cell::Cell)

Write cell using the ABACUS STRU format.
"""
function write_stru(io::IO, cell::Cell)
    _warn_dropped_auxiliary_dimensions(cell, "write_stru")
    lat_ang = cellmat(cell)

    println(io, "ATOMIC_SPECIES")
    unique_species = unique(species(cell))
    for sp in unique_species
        println(io, "$(sp) $(ustrip(elements[sp].atomic_mass)) $(sp).upf upf")
    end

    println(io, "")
    println(io, "LATTICE_CONSTANT")
    @printf io "%.15f\n" ANGSTROM_TO_BOHR

    println(io, "")
    println(io, "LATTICE_VECTORS")
    for row = 1:3
        @printf io "%15.10f %15.10f %15.10f\n" lat_ang[1, row] lat_ang[2, row] lat_ang[
            3,
            row,
        ]
    end

    println(io, "")
    println(io, "ATOMIC_POSITIONS")
    println(io, "Direct")

    pos_frac = inv(lat_ang) * @view(positions(cell)[1:3, :])

    sorted_species = sorted_symbols(unique(species(cell)))

    for sp in sorted_species
        println(io, sp)
        println(io, "0.0")

        indices = findall(x -> x == sp, species(cell))
        println(io, length(indices))

        for idx in indices
            @printf io "%15.10f %15.10f %15.10f\n" pos_frac[1, idx] pos_frac[2, idx] pos_frac[
                3,
                idx,
            ]
        end
    end
end

"""
    write_stru(fname::AbstractString, cell::Cell)

Write `cell` as a STRU file.
"""
function write_stru(fname::AbstractString, cell::Cell)
    open(fname, "w") do io
        write_stru(io, cell)
    end
end

export read_stru, write_stru
