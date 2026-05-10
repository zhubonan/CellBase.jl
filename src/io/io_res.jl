#=
Functions for reading SHELX files
=#

using Printf
import GZip

"""
Read an array containing the lines of the SHELX file
"""
function read_res(lines::Vector{<:AbstractString})
    cellpar = Array{Float64}(undef, 6)
    title_items = Dict{Symbol,Any}()
    line_no = 1

    # Find the number of atoms in the structure
    counter = -1
    nat = -1
    for line in lines
        if contains(line, "SFAC")
            counter = 0
            continue
        end
        if contains(line, "END")
            nat = counter
            break
        end
        if counter >= 0
            counter += 1
        end
    end

    # Storage space for per-atom properties
    species = fill(:na, nat)
    scaled_pos = zeros(3, nat)
    spins = zeros(3, nat)

    iatom = 1
    while line_no < length(lines)
        line = lines[line_no]
        tokens = split(strip(line))
        if tokens[1] == "TITL"
            title_items = parse_titl(line)
        elseif (tokens[1] == "CELL") & (length(tokens) == 8)
            for i = 3:8
                cellpar[i-2] = parse(Float64, tokens[i])
            end
        elseif tokens[1] == "SFAC"
            for atom_line in @view lines[line_no+1:end]
                if strip(atom_line) == "END"
                    break
                end
                atokens = split(strip(atom_line))
                species[iatom] = Symbol(atokens[1])
                for i = 3:5
                    scaled_pos[i-2, iatom] = parse(Float64, atokens[i])
                end
                if length(atokens) == 7
                    spins[iatom] = parse(Float64, atokens[7])
                else
                    spins[iatom] = 0.0
                end
                iatom += 1
            end
        end
        line_no += 1
    end

    lattice = Lattice(cellpar)
    cell = Cell(lattice, species, cellmat(lattice) * scaled_pos)

    # Attach spin only if there are any non-zero ones...
    if any(x -> x != 0, spins)
        cell.arrays[:spins] = spins
    end

    CellBase.attachmetadata!(cell, title_items)
    cell
end


function read_res(s::AbstractString)
    open(s) do handle
        read_res(readlines(handle))
    end
end


"""
    read_res_many(s::IO)

Read many SHELX files from an IO stream
"""
function read_res_many(s::IO)
    lines = String[]
    cells = Cell{Float64}[]
    for line in eachline(s)
        if contains(line, "TITL") && !isempty(lines)
            push!(cells, read_res(lines))
            empty!(lines)
        end
        push!(lines, line)
    end
    isempty(lines) || push!(cells, read_res(lines))
    cells
end

"""
    read_res_many(s::AbstractString)

Read many SHELX files from a packed file.
"""
function read_res_many(s::AbstractString)
    if endswith(s, ".gz")
        op = GZip.open
    else
        op = open
    end
    op(s) do handle
        read_res_many(handle)
    end
end

function parse_titl(s::AbstractString)
    pfloat(x) = parse(Float64, x)
    tokens = split(strip(s))[2:end]
    Dict(
        :label => tokens[1],
        :pressure => pfloat(tokens[2]),
        :volume => pfloat(tokens[3]),
        :enthalpy => pfloat(tokens[4]),
        :spin => pfloat(tokens[5]),
        :abs_spin => pfloat(tokens[6]),
        :natoms => pfloat(tokens[7]),
        :symm => tokens[8],
        :flag1 => tokens[9],
        :flag2 => tokens[10],
        :flag3 => tokens[11],
    )
end

"""
    write_res(io::IO, structure::Cell)

Write out SHELX format data including additional information stored in  `structure.metadata`.

Supported keys: `:label`, `:pressure`, `:volume`, `:enthalpy`, `:spin`, `:abs_spin`, `:symm`, `:flag1`, `flag2`, `flag3`.

The following fields are also treated as the `REM` line:
- `info`: If a `Dict` is supplied the key-value pairs will be written.
- `comments`:  Any additional comments to be written. Must be supplied as a `Vector{<:AbstractString}`.

The composition of the structure will be written as `REM Composition: <E1> <N1> <E2> <N2>`.
"""
function write_res(io::IO, structure::Cell)
    _warn_dropped_auxiliary_dimensions(structure, "write_res")
    infodict = structure.metadata
    titl = (
        label=get(infodict, :label, "CellBase-in-out"),
        pressure=get(infodict, :pressure, 0.0),
        volume=volume(structure),
        enthalpy=get(infodict, :enthalpy, 0.0),
        spin=get(infodict, :spin, 0.0),
        abs_spin=get(infodict, :abs_spin, 0.0),
        natoms=nions(structure),
        symm=get(infodict, :symm, "(n/a)"),
        flag1=get(infodict, :flag1, "n"),
        flag2=get(infodict, :flag2, "-"),
        flag3=get(infodict, :flag3, "1"),
    )
    titl_line = @sprintf("TITL %s %.10f %.10f %.10f %.3f %.3f %d %s %s %s %s\n", titl...)
    write(io, titl_line)

    # Write comment lines
    if :comments in keys(infodict)
        for line in infodict[:comments]
            println(io, "REM $(line)")
        end
    end

    # Write key value pairs
    if :info in keys(infodict) && isa(infodict[:info], Dict)
        for (key, value) in pairs(infodict[:info])
            println(io, "REM $(key): $(value)")
        end
    end

    # Write composition
    comp = Composition(structure)
    print(io, "REM Composition:")
    for (elem, n) in comp
        print(io, " ", elem, " ", n)
    end
    print(io, "\n")

    cell_line = @sprintf(
        "CELL 1.54180 %.6f %.6f %.6f %.6f %.6f %.6f\n",
        cellpar(lattice(structure))...
    )
    write(io, cell_line)
    write(io, "LATT -1\n")
    write(io, "SFAC ", join(map(string, unique(species(structure))), " "), "\n")
    fposmat = CellBase.get_scaled_positions(structure)
    # Wrap
    fposmat .-= floor.(fposmat)

    count = 1
    last_symbol = structure.symbols[1]
    if :spins in keys(structure.arrays)
        spin_array = structure.arrays[:spins]
        for (i, symbol) in enumerate(structure.symbols)
            if symbol != last_symbol
                count += 1
            end
            write(
                io,
                @sprintf(
                    "%-7s %2s %15.12f %15.12f %15.12f 1.0 %.3f\n",
                    symbol,
                    count,
                    fposmat[1, i],
                    fposmat[2, i],
                    fposmat[3, i],
                    spin_array[i]
                )
            )
            last_symbol = symbol
        end
    else
        for (i, symbol) in enumerate(structure.symbols)
            if symbol != last_symbol
                count += 1
            end
            write(
                io,
                @sprintf(
                    "%-7s %2s %15.12f %15.12f %15.12f 1.0\n",
                    symbol,
                    count,
                    fposmat[1, i],
                    fposmat[2, i],
                    fposmat[3, i]
                )
            )
            last_symbol = symbol
        end
    end
    write(io, "END\n")
end

"""
    write_res(fname::AbstractString, structure::Cell, mode="w")

Write out SHELX format data to a file
"""
function write_res(fname::AbstractString, structure::Cell, mode="w")
    open(fname, mode) do fh
        write_res(fh, structure)
    end
end
