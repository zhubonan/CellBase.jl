# STRU File I/O for CellBase.jl

## Overview

CellBase.jl now supports reading and writing ABACUS STRU files. The STRU format is used by the ABACUS density functional theory software package.

## Quick Start

```julia
using CellBase

# Read a STRU file
cell = read_stru("structure.stru")

# Write a cell to STRU format
write_stru("output.stru", cell)
```

## API

### Reading STRU Files

```julia
# Read from file
cell = read_stru("filename.stru")

# Read from IO stream
open("filename.stru") do io
    cell = read_stru(io)
end
```

### Writing STRU Files

```julia
# Write to file
write_stru("output.stru", cell)

# Write to IO stream
open("output.stru", "w") do io
    write_stru(io, cell)
end
```

## Features Supported

### Reading

- **Coordinate Systems:**
  - `Direct` - Fractional coordinates (most common)
  - `Cartesian` - Cartesian in units of LATTICE_CONSTANT (Bohr)
  - `Cartesian_au` - Cartesian in Bohr
  - `Cartesian_angstrom` - Cartesian in Angstrom

- **Sections:**
  - `ATOMIC_SPECIES` - Element definitions (skipped - labels in ATOMIC_POSITIONS)
  - `NUMERICAL_ORBITAL` - Orbital files (skipped)
  - `LATTICE_CONSTANT` - Lattice scaling factor (stored in metadata)
  - `LATTICE_VECTORS` - Lattice vectors (full support)
  - `LATTICE_PARAMETERS` - Bravais lattice parameters (skipped)
  - `ATOMIC_POSITIONS` - Atomic positions (full support)

- **Other Features:**
  - Comments with `//`
  - Blank lines
  - Multiple species
  - Per-species magnetism (read but not stored)

### Writing

- Always writes in `Direct` (fractional) coordinates
- Species sorted by atomic number
- `LATTICE_CONSTANT` = ANGSTROM_TO_BOHR (≈1.889726)
- `LATTICE_VECTORS` written directly in Angstrom
- Default magnetism = 0.0 for all species

## Unit Conversion

CellBase.jl uses Angstrom as the internal unit. All STRU files are converted appropriately:

### Reading (STRU → Cell, Angstrom)

| Coordinate Type | Conversion |
|---------------|-------------|
| Direct | `pos = lattice * positions` |
| Cartesian | `pos = positions * lc_bohr * BOHR_TO_ANGSTROM` |
| Cartesian_au | `pos = positions * BOHR_TO_ANGSTROM` |
| Cartesian_angstrom | `pos = positions` (no conversion) |

### Writing (Cell → STRU)

| Value | STRU Output |
|-------|-------------|
| LATTICE_CONSTANT | `ANGSTROM_TO_BOHR` (≈1.889726) |
| LATTICE_VECTORS | lattice (Angstrom, no conversion) |
| ATOMIC_POSITIONS | Direct coordinates |

## Examples

### Example 1: Create and Write a Structure

```julia
using CellBase

# Create a simple SiO₂ structure
lat = Lattice(9.0, 9.0, 9.0)
cell_species = [:Si, :O, :O, :O]
atom_positions = [
    [0.0, 1.0, 2.0],  # Si
    [2.0, 3.0, 4.0],  # O
    [4.0, 5.0, 6.0],  # O
    [6.0, 7.0, 8.0],  # O
]
pos = hcat(atom_positions...)
cell = Cell(lat, cell_species, pos)

# Write to STRU file
write_stru("sio2.stru", cell)
```

### Example 2: Read and Analyze a STRU File

```julia
using CellBase

# Read STRU file
cell = read_stru("structure.stru")

# Access properties
println("Species: ", species(cell))
println("Number of atoms: ", nions(cell))
println("Lattice vectors:")
println(cellmat(cell))
println("Volume: ", volume(cell), " Å³")

# Access metadata (if available)
if haskey(cell.metadata, :lattice_constant_bohr)
    lc = cell.metadata[:lattice_constant_bohr]
    println("Lattice constant: $lc Bohr")
end
```

### Example 3: Round-Trip Verification

```julia
using CellBase

# Read a STRU file
original_cell = read_stru("input.stru")

# Write to new file
write_stru("output.stru", original_cell)

# Read back and verify
read_cell = read_stru("output.stru")

# Verify structure
println("Lattice preserved: ", isapprox(cellmat(read_cell), cellmat(original_cell)))
println("Species preserved: ", Set(species(read_cell)) == Set(species(original_cell)))
println("Positions preserved: ", isapprox(positions(sort(read_cell)), positions(sort(original_cell))))
```

### Example 4: Structure Manipulation

```julia
using CellBase

# Read STRU file
cell = read_stru("structure.stru")

# Sort by species
sorted_cell = sort(cell)

# Wrap positions to unit cell
wrap!(sorted_cell)

# Write sorted structure
write_stru("sorted.stru", sorted_cell)
```

## STRU File Format Example

```
ATOMIC_SPECIES
Si 28.085 Si.upf upf
O 15.999 O.upf upf

LATTICE_CONSTANT
1.889726124565062

LATTICE_VECTORS
   5.4300000000    0.0000000000    0.0000000000
   0.0000000000    5.4300000000    0.0000000000
   0.0000000000    0.0000000000    5.4300000000

ATOMIC_POSITIONS
Direct
Si
0.0
2
   0.0000000000    0.0000000000    0.0000000000
   0.2500000000    0.2500000000    0.2500000000
O
0.0
4
   0.1250000000    0.1250000000    0.1250000000
   0.3750000000    0.3750000000    0.1250000000
   0.1250000000    0.3750000000    0.3750000000
   0.3750000000    0.1250000000    0.3750000000
```

## Limitations

- **Not supported:**
  - `LATTICE_PARAMETERS` section (Bravais lattice types)
  - Per-atom magnetism, velocities, or spin constraints
  - Pseudopotential and orbital file information (not preserved in round-trip)
  - Coordinate types other than `Direct` for writing

- **Notes:**
  - Species are sorted by atomic number when writing
  - Small numerical errors (<0.01%) may occur in round-trip conversions
  - All internal positions are stored in Cartesian Angstrom coordinates

## Testing

Run the STRU I/O tests:

```bash
julia --project -e 'using Test; include("test/io_tests/test_stru.jl")'
```

## Demo Script

A comprehensive demo is available:

```bash
julia --project demo_stru.jl
```

This demo covers:
- Creating structures and writing to STRU
- Reading STRU files with different coordinate systems
- Round-trip verification
- Accessing metadata
- Structure manipulation
- Multiple species support

## References

- ABACUS documentation: https://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html
- CellBase.jl: https://github.com/your-org/CellBase.jl

## Contributors

- Implementation based on ABACUS STRU file format specification
- Follows CellBase.jl coding conventions and patterns
