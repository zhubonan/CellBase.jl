#=
Functions for writing XYZ files
=#

using ExtXYZ

const _EXTXYZ_RESERVED_INFO_KEYS = Set(["Lattice", "Properties", "pbc"])
const _EXTXYZ_RESERVED_ARRAY_KEYS = Set(["species", "pos"])

## XYZ files

"""
    write_xyz(fname, structures)

Write one or more structures to XYZ / ExtXYZ format using the `ExtXYZ.jl`
backend.

For hyper cells, the first 3 Cartesian coordinates are written as `pos` and any
auxiliary coordinates are written as per-atom ExtXYZ properties named
`extra_dim_1`, `extra_dim_2`, and so on.
"""
function write_xyz(fname, structures::Vector{Cell{T,D}}) where {T,D}
    ExtXYZ.write_frames(fname, _cell_to_extxyz_frame.(structures))
end

"""
    push_xyz!(fname::AbstractString, structure::Cell; append=true)

Append a single frame to an XYZ / ExtXYZ file.
"""
function push_xyz!(fname::AbstractString, structure::Cell; append=true)
    ExtXYZ.write_frames(fname, [_cell_to_extxyz_frame(structure)]; append=append)
    return fname
end

"""
    push_xyz!(io::IOStream, structure::Cell)

Write a single XYZ / ExtXYZ frame to an open stream.
"""
function push_xyz!(io::IOStream, structure::Cell)
    ExtXYZ.write_frames(io, [_cell_to_extxyz_frame(structure)])
    return io
end

"""
    push_xyz!(lines::AbstractVector{<:AbstractString}, structure::Cell)

Deprecated compatibility shim for the former line-buffer API.
"""
function push_xyz!(lines::AbstractVector{<:AbstractString}, structure::Cell)
    Base.depwarn(
        "push_xyz!(lines::AbstractVector{<:AbstractString}, structure) is deprecated; use push_xyz!(io::IOStream, structure) or push_xyz!(path::AbstractString, structure) instead.",
        :push_xyz!,
    )
    mktemp() do path, io
        close(io)
        push_xyz!(path, structure)
        append!(lines, readlines(path))
    end
    return lines
end

function _cell_to_extxyz_frame(structure::Cell)
    pos = positions(structure)
    arrays = Dict{String,Any}(
        "species" => string.(species(structure)),
        "pos" => copy(@view(pos[1:3, :])),
    )
    for dim = 4:size(pos, 1)
        arrays["extra_dim_$(dim - 3)"] = vec(copy(@view(pos[dim, :])))
    end
    for (name, array) in pairs(structure.arrays)
        key = string(name)
        key in _EXTXYZ_RESERVED_ARRAY_KEYS && continue
        occursin(r"^extra_dim_\d+$", key) && continue
        value = _extxyz_array_value(array, natoms(structure))
        isnothing(value) && continue
        arrays[key] = value
    end

    info = Dict{String,Any}()
    for (key, value) in pairs(metadata(structure))
        string_key = string(key)
        string_key in _EXTXYZ_RESERVED_INFO_KEYS && continue
        info[string_key] = _extxyz_info_value(value)
    end

    Dict(
        "N_atoms" => natoms(structure),
        "cell" => Matrix{Float64}(cellmat(structure)),
        "pbc" => collect(periodicity(structure)[1:3]),
        "info" => info,
        "arrays" => arrays,
    )
end

function _extxyz_array_value(array, natoms_expected::Integer)
    if array isa AbstractVector{<:Symbol}
        length(array) == natoms_expected || return nothing
        return string.(collect(array))
    elseif array isa AbstractVector{<:AbstractString}
        length(array) == natoms_expected || return nothing
        return collect(array)
    elseif array isa AbstractVector{<:Union{Integer,AbstractFloat,Bool}}
        length(array) == natoms_expected || return nothing
        return collect(array)
    elseif array isa AbstractMatrix{<:Union{Integer,AbstractFloat,Bool}}
        size(array, 2) == natoms_expected || return nothing
        return Matrix(array)
    end
    return nothing
end

function _extxyz_info_value(value)
    if value isa Symbol
        return string(value)
    elseif value isa AbstractVector{<:Symbol}
        return string.(collect(value))
    elseif value isa Union{Integer,AbstractFloat,Bool,AbstractString}
        return value
    elseif value isa AbstractVector{<:Union{Integer,AbstractFloat,Bool}}
        return collect(value)
    elseif value isa AbstractMatrix{<:Union{Integer,AbstractFloat,Bool}}
        return Matrix(value)
    end
    return string(value)
end

"""
    read_xyz(io::IO; extra_col_map=nothing)
    read_xyz(fname::AbstractString; extra_col_map=nothing)

Read one or more XYZ / ExtXYZ frames and return them as `Vector{Cell}`.

When `extra_dim_*` ExtXYZ properties are present, they are reconstructed into
auxiliary coordinate rows `positions(cell)[4:end, :]`.
"""
function read_xyz(io::IO; extra_col_map=nothing)
    _extxyz_frames_to_cells(ExtXYZ.read_frames(io); extra_col_map)
end

function read_xyz(f::AbstractString; kwargs...)
    _extxyz_frames_to_cells(ExtXYZ.read_frames(f); kwargs...)
end

function _extxyz_frames_to_cells(frames; extra_col_map=nothing)
    [_extxyz_frame_to_cell(frame; extra_col_map) for frame in frames]
end

function _extxyz_frame_to_cell(frame; extra_col_map=nothing)
    arrays = frame["arrays"]
    pos3 = _extxyz_pos_matrix(arrays)
    posmat = _extxyz_positions_with_extra_dims(pos3, arrays)
    nat = size(posmat, 2)
    lattice_matrix = Matrix{Float64}(frame["cell"])
    cell = Cell(Lattice(lattice_matrix), _extxyz_species(arrays, nat), posmat)

    for (key, value) in pairs(get(frame, "info", Dict{String,Any}()))
        cell.metadata[Symbol(key)] = value
    end
    haskey(frame, "pbc") && (cell.metadata[:pbc] = frame["pbc"])

    extra_cols = Any[]
    for (key, value) in pairs(arrays)
        key in _EXTXYZ_RESERVED_ARRAY_KEYS && continue
        occursin(r"^extra_dim_\d+$", key) && continue
        converted = _cell_array_value(value, nat)
        isnothing(converted) || (cell.arrays[Symbol(key)] = converted)
    end
    if !isnothing(extra_col_map)
        ordered_extra_keys = _ordered_extra_property_keys(arrays)
        for atom_idx = 1:nat
            raw_values = String[]
            for key in ordered_extra_keys
                append!(raw_values, _extxyz_atom_property_strings(arrays[key], atom_idx))
            end
            push!(extra_cols, parse.(extra_col_map, raw_values))
        end
    end
    cell.metadata[:extra_cols] = extra_cols
    cell
end

function _extxyz_pos_matrix(arrays::Dict{String,Any})
    haskey(arrays, "pos") || throw(ArgumentError("ExtXYZ frame is missing the required 'pos' property"))
    pos = Matrix{Float64}(arrays["pos"])
    size(pos, 1) == 3 ||
        throw(ArgumentError("ExtXYZ positions must have shape 3xN, got $(size(pos))"))
    pos
end

function _extxyz_species(arrays::Dict{String,Any}, natoms_expected::Integer)
    if haskey(arrays, "species")
        species_values = arrays["species"]
        length(species_values) == natoms_expected ||
            throw(ArgumentError("ExtXYZ 'species' property length $(length(species_values)) does not match natoms($natoms_expected)"))
        return Symbol.(species_values)
    elseif haskey(arrays, "Z")
        atomic_numbers = vec(Int.(arrays["Z"]))
        length(atomic_numbers) == natoms_expected ||
            throw(ArgumentError("ExtXYZ 'Z' property length $(length(atomic_numbers)) does not match natoms($natoms_expected)"))
        return [Symbol(elements[z].symbol) for z in atomic_numbers]
    end
    throw(ArgumentError("ExtXYZ frame must contain either 'species' or 'Z' properties"))
end

function _extxyz_positions_with_extra_dims(pos3::Matrix{Float64}, arrays::Dict{String,Any})
    extra_dims = Dict{Int,Vector{Float64}}()
    for key in keys(arrays)
        match_obj = match(r"^extra_dim_(\d+)$", key)
        isnothing(match_obj) && continue
        extra_idx = parse(Int, only(match_obj.captures))
        values = vec(Float64.(arrays[key]))
        length(values) == size(pos3, 2) ||
            throw(ArgumentError("ExtXYZ property '$key' length $(length(values)) does not match natoms($(size(pos3, 2)))"))
        extra_dims[extra_idx] = values
    end
    isempty(extra_dims) && return pos3

    posmat = zeros(Float64, 3 + maximum(keys(extra_dims)), size(pos3, 2))
    posmat[1:3, :] .= pos3
    for (extra_idx, values) in pairs(extra_dims)
        posmat[3 + extra_idx, :] .= values
    end
    posmat
end

function _cell_array_value(value, natoms_expected::Integer)
    if value isa AbstractVector
        length(value) == natoms_expected || return nothing
        return collect(value)
    elseif value isa AbstractMatrix
        size(value, 2) == natoms_expected || return nothing
        return Matrix(value)
    end
    return nothing
end

function _ordered_extra_property_keys(arrays::Dict{String,Any})
    hyper_keys = Tuple{Int,String}[]
    other_keys = String[]
    for key in keys(arrays)
        key in _EXTXYZ_RESERVED_ARRAY_KEYS && continue
        match_obj = match(r"^extra_dim_(\d+)$", key)
        if isnothing(match_obj)
            push!(other_keys, key)
        else
            push!(hyper_keys, (parse(Int, only(match_obj.captures)), key))
        end
    end
    sort!(hyper_keys, by=first)
    sort!(other_keys)
    vcat(last.(hyper_keys), other_keys)
end

function _extxyz_atom_property_strings(value, atom_idx::Integer)
    if value isa AbstractVector
        return [string(value[atom_idx])]
    elseif value isa AbstractMatrix
        return string.(vec(value[:, atom_idx]))
    end
    return String[]
end
