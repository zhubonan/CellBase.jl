
module CellIO

using CellBase: cellpar2mat

"""
Read a cell file
"""
function read_cell(fname::AbstractString)
    lines = readlines(fname)
    read_cell(lines)
end

"Clean lines by removing new line symbols and comments"
function clean_lines(lines_in)
    lines_out = String[]
    for line in lines_in
        tmp = strip(line)
        if length(tmp) == 0
            continue
        end
        icomment = findnext('#', tmp, 1)
        if icomment === nothing
            push!(lines_out, string(tmp))
        elseif icomment == 1
            continue
        else
            push!(lines_out, string(strip(tmp[1:icomment-1])))
        end
    end
    lines_out
end

"""
Separate the unit from the actual content of the block
"""
function separate_unit(block)
    # Parse the unit
    unit = ""
    if length(split(block[1])) == 1
        unit = uppercase(strip(block[1]))
        popfirst!(block)
    end
    unit, block
end


"Read cell related sections"
function read_cellmat(lines)
    block = find_block(lines, "LATTICE_CART")

    if length(block) > 0
        unit, block = separate_unit(block)
        !isempty(unit) &&
            @assert startswith(unit, "ANG") "Only support Angstrom for lattice but $(unit) found"
        cellmat = copy(read_num_block(block, 3, column_major=true))

    else
        cell_par = Array{Float64}(undef, 6)
        block = find_block(lines, "LATTICE_ABC")
        unit, block = separate_unit(block)
        !isempty(unit) &&
            @assert startswith(unit, "ANG") "Only support Angstrom  for lattice but $(unit) found"

        @assert length(block) > 0 "Neither lattice_cart or lattice_abc are present"
        cell_par_raw = read_num_block(block, 3)
        cell_par[1:3] = cell_par_raw[1, :]
        cell_par[4:6] = cell_par_raw[2, :]
        cellmat = cellpar2mat(cell_par...)
    end
    cellmat
end

"Read positions related sections"
function read_positions(lines, cellmat)
    positions = Vector{Float64}[]
    is_abs = false
    block = find_block(lines, "POSITIONS_ABS")
    if length(block) > 0
        is_abs = true
        unit, block = separate_unit(block)
        !isempty(unit) &&
            @assert startswith(unit, "ANG") "Only support Angstrom for positions but found $(unit)"
    else
        block = find_block(lines, "POSITIONS_FRAC")
    end
    @assert length(block) > 0 "Positions block found!"
    # Parse the position block
    ion_names = Symbol[]
    for line in block
        sline = split(line)
        push!(positions, map(x -> parse(Float64, x), sline[2:4]))
        push!(ion_names, Symbol(sline[1]))
    end

    # The positions are stored as column vector matrix
    posmat = hcat(positions...)
    if is_abs == false
        posmat = cellmat * posmat
    end
    posmat, ion_names
end


"""
Read content of a cell file
"""
function read_cell(lines_in::Vector{T}) where {T<:AbstractString}

    lines = clean_lines(lines_in)
    cellmat = read_cellmat(lines)
    posmat, ion_names = read_positions(lines, cellmat)

    return cellmat, posmat, ion_names
end


"""
Read in cell 
"""
function read_seed(fname::AbstractString)
    lines = readlines(fname)
    read_seed(lines)
end


"""
Read seed with the label parsed. At this stage no explansions are done
"""
function read_seed(lines_in::Vector{T}) where {T<:AbstractString}

    cell_structure = read_cell(lines_in)
    # read positions

    is_abs = false
    block = find_block(lines_in, "POSITIONS_ABS")
    if !isempty(block)
        is_abs = true
    else
        block = find_block(lines_in, "POSITIONS_FRAC")
    end
    ion_tags = parse_taglines(block)

    # Read global tags in the format #KEY=<SOMETHING>
    global_tags = Dict{Symbol,Any}()
    for line in lines_in
        m = match(r"^#(\w+)=(.*)$", line)
        if !isnothing(m)
            global_tags[Symbol(m.captures[1])] = strip(m.captures[2])
        end
    end

    return cell_structure, ion_tags, global_tags
end

"""
Parse the tag line, return a dictionary of the
tagline is in the form of <ion_set_name> % KEY=VALUE KEY=VALUE
"""
function _parse_tagline(i, line)
    num_tags = ["POSAMP", "MINAMP", "ANGAMP", "XAMP", "ZAMP", "RAD"]
    flag_tags = ["FIX", "NOMOVE", "PERM", "ADATOM", "ATHOLE"]
    not_implemented = ["COORD"]

    tokens = split(line, "%", limit=2)
    has_settings = false
    if length(tokens) > 1
        ion_set_name = strip(tokens[1])
        tokens = split(tokens[2])
        has_settings = true
    else
        if isempty(tokens[1])
            ion_set_name = "Ion-" * string(i)
        else
            ion_set_name = strip(tokens[1])
        end
    end

    out_dict = Dict{Symbol,Any}(:ion_set_name => ion_set_name)

    # Read the tokens
    if has_settings
        for token in tokens
            m = match(r"(\w+)=(.+)", token)
            if !isnothing(m)
                key = m.captures[1]
                value = m.captures[2]
                out_dict[Symbol(key)] = parse(Float64, value)
            elseif token in flag_tags
                out_dict[Symbol(token)] = true
            else
                throw(ErrorException("Token \"$(token)\" not understood"))
            end
        end
    end
    return out_dict
end

"""
Parse tag lines, return a vector of dictionary containing parsed tags 
"""
function parse_taglines(lines::Vector{T}) where {T<:AbstractString}
    tags = []
    i = 0
    for line in lines
        line = strip(line)
        if startswith(line, "#")
            continue
        end
        tokens = split(line, "#")
        if length(tokens) <= 1
            continue
        else
            i += 1
            ltmp = tokens[2]
            push!(tags, _parse_tagline(i, ltmp))
        end
    end
    return tags
end

"""
Read a numerical block in the form of a vector of String
"""
function read_num_block(
    lines::Vector{T},
    span::Int;
    ntype::Type=Float64,
    column_major=false,
) where {T<:AbstractString}
    if column_major
        block = Array{ntype}(undef, (span, length(lines)))
        for (n, line) in enumerate(lines)
            block[:, n] = map(p -> parse(ntype, p), split(line))
        end
    else
        block = Array{ntype}(undef, (length(lines), span))
        for (n, line) in enumerate(lines)
            block[n, :] = map(p -> parse(ntype, p), split(line))
        end
    end
    return block
end


"""
Find the block with given name, return a Vector of the lines
Comments are skipped
"""
function find_block(lines, block_name)
    in_block = false
    block_lines = String[]
    b_start = "%BLOCK " * block_name
    bname = ""
    for line in lines
        # Skip any comment line
        if startswith(line, "#")
            continue
        elseif startswith(uppercase(line), b_start)
            current_block_name = uppercase(split(line)[2])
            in_block = true
        elseif in_block && startswith(uppercase(line), "%ENDBLOCK")
            bname = uppercase(split(line)[2])
            @assert bname == block_name "Block <$(block_name)> ends as <$(bname)>"
            break
        elseif in_block
            push!(block_lines, line)
        end
    end
    return block_lines
end

"""
Get rid of blocks 
"""
function filter_block(lines)
    out_lines = String[]
    in_block = false
    for line in lines
        # Skip any comment line
        if startswith(uppercase(line), "%BLOCK")
            in_block = true
        elseif in_block && startswith(uppercase(line), "%ENDBLOCK")
            in_block = false
        elseif !in_block
            push!(out_lines, line)
        end
    end
    return out_lines
end

"""
Filter away block names contained in only
"""
function filter_block(lines, only::Vector{T}) where {T<:AbstractString}
    out_lines = String[]
    in_block = false
    bname = ""
    only_upper = [uppercase(x) for x in only]
    for line in lines
        # Skip any comment line
        if startswith(uppercase(line), "%BLOCK")
            bname = split(line)[2]
            if bname in only
                in_block = true
            else
                push!(out_lines, line)
            end
        elseif in_block && startswith(uppercase(line), "%ENDBLOCK")
            @assert split(line)[2] == bname "Incomplete block <$bname>, ends with $line"
            in_block = false
        elseif !in_block
            push!(out_lines, line)
        end
    end
    return out_lines
end

# Writer
function write_cell(handle::IO, latt, pos, species)
    write(handle, "%BLOCK LATTICE_CART\n")
    for (a, b, c) in eachcol(latt)
        write(handle, " $(a) $(b) $(c)\n")
    end
    write(handle, "%ENDBLOCK LATTICE_CART\n")

    write(handle, "%BLOCK POSITIONS_ABS\n")
    for (sym, (a, b, c)) in zip(species, eachcol(pos))
        write(handle, "$(sym) $(a) $(b) $(c)\n")
    end
    write(handle, "%ENDBLOCK POSITIONS_ABS\n")
end

function write_cell(fname::AbstractString, latt, pos, species)
    open(fname, "w") do handle
        write_cell(handle, latt, pos, species)
    end
end

end # Module cell IO

using .CellIO

function read_cell(fname)
    cellmat, posmat, species = CellIO.read_cell(fname)
    lattice = Lattice(cellmat)
    return Cell(lattice, species, posmat)
end

"""
    write_cell(fname, cell::Cell)
"""
function write_cell(fname, cell::Cell)
    _warn_dropped_auxiliary_dimensions(cell, "write_cell")
    CellIO.write_cell(fname, cellmat(cell), @view(positions(cell)[1:3, :]), species(cell))
end
