#=
Code for writing/parsing POSCAR
=#

"""
    write_poscar(io::IO, cell)

Write cell using the POSCAR format.
"""
function write_poscar(io::IO, cell::Cell)
    _warn_dropped_auxiliary_dimensions(cell, "write_poscar")
    println(io, formula(cell))
    println(io, "  1.000000000000")
    latt = cellmat(cell)
    @printf io "    %20.16f %20.16f %20.16f\n" latt[1, 1] latt[2, 1] latt[3, 1]
    @printf io "    %20.16f %20.16f %20.16f\n" latt[1, 2] latt[2, 2] latt[3, 2]
    @printf io "    %20.16f %20.16f %20.16f\n" latt[1, 3] latt[2, 3] latt[3, 3]

    # Write species
    # Here we need to sort the cell
    syms = species(cell)
    usyms = unique(syms)
    # Array of symbol index
    sindex = [findfirst(x -> x == sym, usyms) for sym in syms]

    sort_idx = sortperm(sindex)
    println(io, "   ", join(map(string, usyms), "  "))

    # Check the number of species
    specie_counts = [count(x -> x == sym, syms) for sym in usyms]
    println(io, "   ", join(map(string, specie_counts), "  "))
    println(io, "Direct")
    scaled_pos = get_scaled_positions(cell)
    for idx in sort_idx
        x, y, z = scaled_pos[:, idx]
        @printf io "  %20.16f  %20.16f  %20.16f\n" x y z
    end

end

"""
    write_poscar(fname::AbstractString, cell::Cell)
Write `cell` as a `POSCAR` file.
"""
function write_poscar(fname::AbstractString, cell::Cell)
    open(fname, "w") do io
        write_poscar(io, cell)
    end
end

"""
    read_poscar(io::IO)
Read a `Cell` object from a `POSCAR` file.
"""
function read_poscar(io::IO)
    comment = readline(io) # comment
    scale = parse(Float64, readline(io))
    v1 = parse.(Float64, split(readline(io)))
    v2 = parse.(Float64, split(readline(io)))
    v3 = parse.(Float64, split(readline(io)))
    mat = hcat(v1, v2, v3)
    mat .*= scale
    species = map(Symbol, split(readline(io)))
    counts = map(x -> parse(Int, x), split(readline(io)))
    symbols = vcat([repeat([x], c) for (x, c) in zip(species, counts)]...)

    posmat = zeros(Float64, 3, sum(counts))
    poskind = readline(io)

    for i in axes(posmat, 2)
        tokens = split(readline(io))
        for j = 1:3
            posmat[j, i] = parse(Float64, tokens[j])
        end
    end
    if startswith(lowercase(poskind), "d")
        posmat = mat * posmat
    elseif startswith(lowercase(poskind), "c")
        posmat .*= scale
    else
        throw(ErrorException("Unknown coordinate kind: $(poskind)"))
    end

    out = Cell(Lattice(mat), symbols, posmat)
    out.metadata[:comments] = [comment]
    out
end

"""
    read_poscar(fname::AbstractString)
Read a `Cell` object from a `POSCAR` file.
"""
function read_poscar(fname::AbstractString)
    open(fname) do io
        read_poscar(io)
    end
end

export write_poscar, read_poscar
