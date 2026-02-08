#=
Functions for writing XYZ files
=#

## XYZ files

"""
    Write snapshots to a xyz file
"""
function write_xyz(fname, structures::Vector{Cell{T,D}}) where {T,D}
    lines = String[]
    for structure in structures
        push_xyz!(lines, structure)
    end

    open(fname, "w") do handle
        for line in lines
            write(handle, line)
            write(handle, "\n")
        end
    end
end


"""xyz lines for a single frame"""
function push_xyz!(lines, structure::Cell)
    ns = nions(structure)
    cell = cellmat(lattice(structure))
    push!(lines, "$ns")
    ax, ay, az = cell[:, 1]
    bx, by, bz = cell[:, 2]
    cx, cy, cz = cell[:, 3]

    # Include extra properties
    info_lines = []
    for (key, value) in CellBase.metadata(structure)
        push!(info_lines, "$(key)=\"$(value)\"")
    end
    info_string = join(info_lines, " ")

    comment_line = "Lattice = \"$ax $ay $az $bx $by $bz $cx $cy $cz\" Properties=\"species:S:1:pos:R:3\" $(info_string)"
    push!(lines, comment_line)

    # Write the atoms
    sp = species(structure)
    pos = positions(structure)
    for i = 1:ns
        specie = sp[i]
        x, y, z = pos[:, i]
        line = "$specie $x $y $z"
        push!(lines, line)
    end
    return lines
end

"""
Read XYZ file

NOTE: This is a pure julia implementation which may not work with all cases.
In the future, better to use ExtXYZ.jl interface.
"""
function read_xyz(io::IO; extra_col_map=nothing)
    frames = []
    while true
        line = readline(io)
        if line == ""
            break
        end
        natoms = parse(Int, line)
        attr = readline(io)
        # Read the atoms
        species = Symbol[:NULL for _ = 1:natoms]
        positions = Matrix{Float64}(undef, 3, natoms)
        extra_cols = Any[]
        for i = 1:natoms
            content = split(readline(io))
            species[i] = Symbol(content[1])
            for j = 1:3
                positions[j, i] = parse(Float64, content[1+j])
            end
            if !isnothing(extra_col_map)
                push!(extra_cols, parse.(extra_col_map, content[5:end]))
            end
        end
        # Process attr
        attrdict = _process_xyz_attr(attr)
        lattice = parse.(Float64, split(attrdict[:Lattice]))
        if length(lattice) == 9
            cmat = reshape(lattice, 3, 3) |> transpose |> collect
            cell = Cell(Lattice(cmat), species, positions)
            if length(lattice) == 6
                cell = Cell(Lattice(lattice...), species, positions)
            end
        else
            throw(ErrorException())
        end
        for (a, b) in attrdict
            cell.metadata[a] = b
        end
        cell.metadata[:extra_cols] = extra_cols
        push!(frames, cell)
    end
    frames
end

function read_xyz(f::AbstractString; kwargs...)
    open(f) do fh
        read_xyz(fh; kwargs...)
    end
end

"""
    _process_xyz_attr(attr)

Process ExtXYZ attribute string.
"""
function _process_xyz_attr(attr)
    pattern = r"(\w+)=[^\"]([^ ]+)"
    out = Dict{Symbol,String}()
    for m in eachmatch(pattern, attr)
        out[Symbol(m.captures[1])] = m.captures[2]
    end
    pattern = r"(\w+)=\"([^\"]+)\""
    for m in eachmatch(pattern, attr)
        out[Symbol(m.captures[1])] = m.captures[2]
    end
    out
end
