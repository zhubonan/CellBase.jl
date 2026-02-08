# [File I/O](@id file-io)

This guide covers reading and writing crystal structure files in various formats supported by CellBase.jl.

```@setup fileio
using CellBase
```

## Supported Formats

CellBase.jl supports the following file formats:

| Format | Extension | Read | Write | Notes |
|--------|-----------|------|-------|-------|
| SHELX RES | `.res` | Yes | Yes | AIRSS-style output |
| VASP POSCAR | `POSCAR` | No | Yes | Write-only |
| XYZ | `.xyz` | No | Yes | Write-only |
| CASTEP | `.castep` | Yes | No | Includes energies and forces |
| ABACUS STRU | `STRU` | Yes | Yes | ABACUS input format |
| CELL | `.cell` | Yes | No | CASTEP cell files |

## Writing Structures

### VASP POSCAR Format

```@example fileio
# Create a structure
cell = Cell(Lattice(5.0, 5.0, 5.0), [:Na, :Cl], [0.0 2.5; 0.0 2.5; 0.0 2.5])

println("Structure to save:")
println("  Atoms: ", natoms(cell))
for i in 1:natoms(cell)
    println("  Atom ", i, " (", species(cell)[i], "): ", positions(cell)[:, i])
end

# Write to POSCAR format
write_poscar("POSCAR_test", cell)

# Read it back to verify
content = read("POSCAR_test", String)
println("\nPOSCAR content:")
println(content)
```

### XYZ Format

```@example fileio
# Write to XYZ format - save multiple frames
frames = Cell{Float64, 3}[cell, cell]  # Two identical frames
write_xyz("structure.xyz", frames)

# Read it back
content = read("structure.xyz", String)
println("XYZ content:")
println(content)
```

### SHELX RES Format

```@example fileio
# Write to RES format
write_res("output.res", cell)

# Read it back
content = read("output.res", String)
println("RES content (first 500 chars):")
println(content[1:min(500, length(content))])
```

### ABACUS STRU Format

```@example fileio
# Write to STRU format
write_stru("STRU_test", cell)

# Read it back
content = read("STRU_test", String)
println("STRU content:")
println(content)
```

## Working with Additional Data

When working with files, additional data is stored in `arrays` and `metadata`:

```@example fileio
# Create a cell with additional data
cell = Cell(Lattice(5.0, 5.0, 5.0), [:Si, :Si], [0.0 2.5; 0.0 2.5; 0.0 2.5])

# Add forces
cell.arrays[:forces] = [0.1 0.2; 0.1 0.2; 0.1 0.2]

# Add metadata
cell.metadata[:energy] = -10.5
cell.metadata[:pressure] = 1.0

# Show positions and forces
println("Atom data:")
for i in 1:natoms(cell)
    println("  Atom ", i, " (", species(cell)[i], "):")
    println("    Position: ", positions(cell)[:, i])
    println("    Force: ", cell.arrays[:forces][:, i])
end

# Check what additional arrays are available
println("\nAvailable arrays: ", collect(arraynames(cell)))
println("Metadata keys: ", keys(cell.metadata))
```

## Complete Example: Creating and Saving

Here's a complete example of creating a structure and saving it to multiple formats:

```@example fileio
# Create a NaCl structure
lattice = Lattice(5.64, 5.64, 5.64)
symbols = [:Na, :Cl, :Na, :Cl, :Na, :Cl, :Na, :Cl]
pos_mat = [
    0.0  2.82  2.82  0.0   2.82  0.0   0.0   2.82;
    0.0  0.0   2.82  2.82  2.82  0.0   2.82  0.0;
    0.0  2.82  0.0   2.82  0.0   2.82  2.82  0.0
]
nacl = Cell(lattice, symbols, pos_mat)

# Inspect the structure
println("NaCl structure:")
println("  Formula: ", reduced_fu(nacl))
println("  Number of atoms: ", natoms(nacl))
println("  Volume: ", volume(nacl), " Å³")
println("\nAtom positions:")
for i in 1:min(4, natoms(nacl))
    println("  Atom ", i, " (", species(nacl)[i], "): ", positions(nacl)[:, i])
end

# Save to different formats
write_poscar("POSCAR_nacl", nacl)
# XYZ format for multiple frames
frames = Cell{Float64, 3}[nacl]
write_xyz("nacl.xyz", frames)
write_res("nacl.res", nacl)
write_stru("STRU_nacl", nacl)

println("\nSaved to POSCAR_nacl, nacl.xyz, nacl.res, STRU_nacl")
```

## Batch Processing Example

Example of how to process multiple structures:

```@example fileio
# Create several structures
structures = [
    ("NaCl", Cell(Lattice(5.64, 5.64, 5.64), [:Na, :Cl], [0.0 2.82; 0.0 2.82; 0.0 2.82])),
    ("LiF", Cell(Lattice(4.02, 4.02, 4.02), [:Li, :F], [0.0 2.01; 0.0 2.01; 0.0 2.01])),
    ("KCl", Cell(Lattice(6.29, 6.29, 6.29), [:K, :Cl], [0.0 3.15; 0.0 3.15; 0.0 3.15]))
]

# Show structures before saving
for (name, struc) in structures
    println(name, " structure:")
    println("  Atoms: ", natoms(struc))
    println("  Position: ", positions(struc)[:, 1])
end

# Save each to multiple formats
for (name, struc) in structures
    write_poscar("POSCAR_$(name)", struc)
    # XYZ format requires Vector{Cell}
    xyz_frames = Cell{Float64, 3}[struc]
    write_xyz("$(name).xyz", xyz_frames)
    println("\nSaved $name to POSCAR_$name and $(name).xyz")
end
```

## Best Practices

1. **Always check format support**: Some formats are read-only or write-only
2. **Handle missing data**: Files may not contain all optional information  
3. **Validate structures**: After creating, check that the structure is reasonable:

```@example fileio
# Validation example
cell = Cell(Lattice(5.0, 5.0, 5.0), [:H, :H], [0.0 2.0; 0.0 0.0; 0.0 0.0])

println("Structure validation:")
println("  Atom positions:")
for i in 1:natoms(cell)
    println("    Atom ", i, ": ", positions(cell)[:, i])
end
println("\n  Has atoms: ", natoms(cell) > 0)
println("  Volume is positive: ", volume(cell) > 0)
println("  All positions finite: ", all(isfinite, positions(cell)))
```
