# CellBase.jl

CellBase.jl is a Julia package for building and manipulating periodic crystal structures. It provides a simple yet powerful interface for working with atomic structures, with support for file I/O, symmetry operations, and common crystallographic tasks.

## What This Package Does

CellBase.jl provides:

- **Lattice operations**: Create and manipulate periodic lattice vectors
- **Atomic structure handling**: Store and manipulate atomic positions and species
- **File I/O**: Read and write various structure formats (RES, POSCAR, XYZ, STRU, CASTEP)
- **Structure building**: Create common crystal structures (FCC, BCC, HCP, diamond, etc.)
- **Symmetry analysis**: Integration with Spglib for space group operations
- **Supercell generation**: Create supercells from primitive cells
- **Neighbor lists**: Compute neighbor lists and distances
- **Hyperdimensional cells**: Support auxiliary non-periodic coordinates for advanced workflows

The target application is for small and periodic cells, with emphasis on both ease of use and performance.

## Getting Started

New to CellBase.jl? Start with the [Getting Started](@ref getting-started) guide to learn the basics.

```julia
using CellBase

# Create a simple NaCl structure
cell = bulk("NaCl", "rocksalt", a=5.64)

# Create a 2×2×2 supercell
supercell = repeat(cell, 2, 2, 2)

# Save to file
write_poscar("POSCAR", supercell)
```

## Guides

### [Getting Started](@ref getting-started)
Installation, basic concepts, and your first steps with CellBase.jl.

### [Working with Cells](@ref working-with-cells)
Learn how to manipulate Cell objects: supercells, sorting atoms, rotations, and more.
This guide also covers the current hyper-cell model and its `3 x N` scaled-coordinate API.

### [File I/O](@ref file-io)
Reading and writing structure files in various formats (RES, POSCAR, XYZ/ExtXYZ, STRU, CASTEP).

### [Building Structures](@ref building-structures)
Create common crystal structures using the `bulk` function: FCC, BCC, HCP, diamond, rocksalt, and more.

### [Spacegroup Operations](@ref spacegroup-operations)
Symmetry analysis using Spglib integration: finding space groups, standardizing cells, and more.

## API Reference

### Core Types

```@autodocs
Modules = [CellBase]
Order = [:type]
Pages = ["lattice.jl", "cell.jl", "site.jl", "neighbour.jl", "composition.jl"]
```

### Lattice Functions

```@autodocs
Modules = [CellBase]
Order = [:function]
Pages = ["lattice.jl", "periodic.jl"]
```

### Cell Functions

```@autodocs
Modules = [CellBase]
Order = [:function]
Pages = ["cell.jl"]
```

### Site Functions

```@autodocs
Modules = [CellBase]
Order = [:function]
Pages = ["site.jl"]
```

### File I/O Types

```@autodocs
Modules = [CellBase.SheapIO, CellBase.DotCastep]
Order = [:type]
Pages = ["io/io_sheap.jl", "io/io_dotcastep.jl"]
```

### Building Structures

```@autodocs
Modules = [CellBase]
Order = [:function]
Pages = ["build.jl"]
```

### File I/O

```@docs
CellBase.read_cell
CellBase.write_cell
CellBase.read_res
CellBase.write_res
CellBase.read_poscar
CellBase.write_poscar
CellBase.read_xyz
CellBase.write_xyz
CellBase.push_xyz!
CellBase.read_stru
CellBase.write_stru
CellBase.DotCastep.read_castep
CellBase.SheapIO.run_sheap
```

### Space Group Operations

```@autodocs
Modules = [CellBase]
Order = [:function]
Pages = ["spg.jl"]
```

### Neighbor Lists

```@autodocs
Modules = [CellBase]
Order = [:function]
Pages = ["neighbour.jl"]
```

### Math Utilities

```@autodocs
Modules = [CellBase]
Order = [:function]
Pages = ["mathutils.jl", "minkowski.jl"]
```

### Composition

```@autodocs
Modules = [CellBase]
Order = [:function]
Pages = ["composition.jl"]
```

## Index

```@index
```
