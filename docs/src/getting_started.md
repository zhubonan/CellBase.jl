# [Getting Started](@id getting-started)

Welcome to CellBase.jl! This guide will help you install and start using the package for working with periodic crystal structures.

## Installation

CellBase.jl can be installed using Julia's package manager. Open Julia and run:

```julia
using Pkg
Pkg.add("CellBase")
```

Or use the package mode (press `]` in the Julia REPL):

```julia
] add CellBase
```

## Quick Start

Here's a simple example to get you started:

```@setup quickstart
using CellBase
```

```@example quickstart
# Create a simple cubic lattice with a=5.0 Å
lattice = Lattice(5.0, 5.0, 5.0)

# Create a Cell with two atoms
symbols = [:Na, :Cl]
pos_mat = [0.0 2.5; 0.0 2.5; 0.0 2.5]  # 3×2 matrix (column vectors)
cell = Cell(lattice, symbols, pos_mat)

# Inspect the cell
println("Number of atoms: ", natoms(cell))
println("Volume: ", volume(cell), " Å³")
println("Species: ", species(cell))
```

## Core Concepts

### Lattice

The `Lattice` type represents the periodic cell vectors. Lattice vectors are stored as **column vectors** in a 3×3 matrix:

```@example quickstart
# Create from lattice parameters
lattice = Lattice(5.0, 5.0, 5.0)  # Cubic cell

# Create from a matrix (column vectors)
mat = [5.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 5.0]
lattice2 = Lattice(mat)

# Access lattice vectors
a, b, c = cellvecs(lattice)
println("Lattice vector a: ", a)
println("Lattice vector b: ", b)
println("Lattice vector c: ", c)
```

### Cell

The `Cell` type combines a `Lattice` with atomic positions and species:

```@example quickstart
# Create a Cell
lattice = Lattice(4.0, 4.0, 4.0)
symbols = [:Si, :Si]
pos_mat = [0.0 2.0; 0.0 2.0; 0.0 2.0]
cell2 = Cell(lattice, symbols, pos_mat)

# Get information
println("Number of atoms: ", natoms(cell2))
println("Species: ", species(cell2))
println("Atomic numbers: ", atomic_numbers(cell2))

# Show positions for each atom
pos_mat2 = positions(cell2)
for i in 1:natoms(cell2)
    println("Atom ", i, " (", species(cell2)[i], "): ", pos_mat2[:, i])
end
```

### Positions

Positions are stored as a 3×N matrix where each column is an atom's Cartesian coordinates:

```@example quickstart
# Show positions for each atom
pos_mat3 = positions(cell2)
for i in 1:natoms(cell2)
    println("Atom ", i, " (", species(cell2)[i], "): ", pos_mat3[:, i])
end
```

You can also work with scaled (fractional) positions:

```@example quickstart
# Get scaled positions
scaled = get_scaled_positions(cell2)
println("Scaled positions (fractional):")
for i in 1:natoms(cell2)
    println("  Atom ", i, ": ", scaled[:, i])
end

# Set scaled positions
set_scaled_positions!(cell2, [0.0 0.5; 0.0 0.5; 0.0 0.5])
println("\nNew positions after setting fractional coords:")
for i in 1:natoms(cell2)
    println("  Atom ", i, ": ", positions(cell2)[:, i])
end

# Wrap atoms back into the cell
wrap!(cell2)
println("\nWrapped positions:")
for i in 1:natoms(cell2)
    println("  Atom ", i, ": ", positions(cell2)[:, i])
end
```

## Next Steps

- Learn about [Working with Cells](@ref working-with-cells) for operations like supercells, rotations, and sorting
- See [File I/O](@ref file-io) for reading and writing structure files
- Check out [Building Structures](@ref building-structures) for creating common crystal structures
- Explore the API Reference below for detailed function documentation
