# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CellBase.jl is a Julia package for building and manipulating periodic crystal structures. It provides lattice operations, atomic position handling, file I/O for various formats (VASP, CASTEP, XYZ, etc.), and integration with Spglib and AtomsBase.

## Build, Test, and Development Commands

### Running Tests

```bash
# Run all tests
julia --project -e 'using Pkg; Pkg.test("CellBase")'

# Run tests with coverage
julia --project -e 'using Pkg; Pkg.test("CellBase", coverage=true)'

# Run a single test file interactively
julia --project -e 'include("test/test_cell.jl")'

# Run specific testset
julia --project -e 'using Test; include("test/runtests.jl")' -- --testset="Cell"
```

### Environment Setup

```bash
# Instantiate project dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'

# Add a new dependency
julia --project -e 'using Pkg; Pkg.add("PackageName")'

# Update dependencies
julia --project -e 'using Pkg; Pkg.update()'
```

### Documentation

```bash
# Build documentation locally
julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs/ docs/make.jl
```

### Code Formatting

```bash
# Format code using JuliaFormatter
julia --project -e 'using JuliaFormatter; format("src")'
julia --project -e 'using JuliaFormatter; format("test")'
```

Formatter configuration is in `.JuliaFormatter.toml`:
- `ignore = ["scripts"]` - Skip scripts directory
- `whitespace_in_kwargs = false` - No whitespace around keyword argument equals

## Architecture Overview

### Core Data Structures

**Cell{T,D}** - The central data structure representing a periodic crystal structure:
- `T`: Element type (usually Float64)
- `D`: Number of dimensions (typically 3, but supports hyperdimensional)
- Contains: `Lattice`, `symbols`, `positions`, `arrays`, `metadata`
- Extends `AtomsBase.AbstractSystem{D}` for interoperability

**Lattice{T}** - Represents periodic lattice vectors:
- Stores both `matrix` (real space) and `rec` (reciprocal space, pre-computed)
- Reciprocal lattice must be updated via `update_rec!()` after modifying matrix
- Supports both 3D and N-dimensional lattices

### File Organization

```
src/
├── CellBase.jl          # Main module entry point
├── cell.jl              # Core Cell struct and operations
├── lattice.jl           # Lattice struct and MIC (Minimum Image Convention)
├── build.jl             # bulk() constructor for common crystal structures
├── composition.jl       # Chemical formula and composition utilities
├── neighbour.jl         # Neighbor list computation
├── spg.jl              # Spglib integration for symmetry
├── minkowski.jl        # Minkowski reduction for safe MIC
├── periodic.jl         # Periodic boundary operations
├── site.jl             # Site abstraction
├── mathutils.jl        # Mathematical utilities (rotations, etc.)
├── reference_data.jl   # Elemental reference data
├── io/                 # File I/O operations
│   ├── io.jl          # I/O module and exports
│   ├── io_cell.jl     # AIRSS .cell format
│   ├── io_res.jl      # AIRSS .res format
│   ├── io_xyz.jl      # XYZ format
│   ├── io_poscar.jl   # VASP POSCAR format
│   ├── io_dotcastep.jl # CASTEP .castep file format
│   └── io_stru.jl     # ABACUS STRU format
└── external/
    └── atomsbase.jl    # AtomsBase interface implementation
```

### Key Architectural Patterns

**Type Parameters:** Both `Cell` and `Lattice` are parametric types:
- Use `Cell{T,D}` where `T` is element type and `D` is dimensions
- Enables N-dimensional positions (hyperdimensional relaxation support)
- StaticArrays (`SVector`, `SMatrix`) used for performance-critical paths

**Minimum Image Convention (MIC):**
- `mic(lattice, vectors)` computes minimum-image representation
- Falls back to safe Minkowski-reduction based method for skewed cells
- Always use MIC when computing distances/forces in periodic systems

**Bulk Crystal Builder:**
- `bulk("Cu")` - Auto-detects structure from reference database
- `bulk("MgO", "rocksalt")` - Explicit structure type
- Supports: sc, fcc, bcc, hcp, diamond, zincblende, wurtzite, rocksalt, etc.
- Reference data from SMACT package in `reference_data.jl`

**AtomsBase Integration:**
- `Cell` implements `AtomsBase.AbstractSystem{D}`
- Convert: `atomic_system(cell)` and `Cell(system)`
- Enables interoperability with other AtomsBase-compatible packages

**Supercell Generation:**
- `make_supercell(cell, P)` - General transformation matrix
- `repeat(cell, a, b, c)` - Diagonal supercell (shorthand)
- Supports both "cell-major" and "atom-major" ordering
- Handles additional arrays correctly

## Code Style Guidelines

### Imports

Use explicit imports at the top of each file:

```julia
using Printf
using PeriodicTable
using LinearAlgebra
import Base  # For method extension
import AtomsBase
const AB = AtomsBase  # Short alias pattern
```

**Patterns:**
- `using` for functions/types used frequently
- `import` when extending methods (e.g., `import Base: show`)
- Create short aliases for long module names (e.g., `const AB = AtomsBase`)

### Naming Conventions

- **Modules**: PascalCase (e.g., `CellBase`)
- **Types/Structs**: PascalCase (e.g., `Cell`, `Lattice`, `NeighbourList`)
- **Functions**: lowercase with underscores (e.g., `get_positions`, `niggli_reduce`)
- **Constants**: UPPERCASE or regular variables (no strict const pattern)
- **Type parameters**: Short uppercase (e.g., `T`, `D`, `N`)
- **Private functions**: Prefixed with underscore (e.g., `_distance_matrix_mic`)

### Type Signatures

Always use explicit type signatures for public functions:

```julia
function Cell(l::Lattice, symbols::Vector{Symbol}, positions::Matrix)
function clip(cell::Cell{T,N}, mask::AbstractVector) where {T,N}
```

### Documentation

Use Julia docstrings with triple quotes:

```julia
"""
    function_name(arg1::Type1, arg2::Type2) -> ReturnType

Brief description of what the function does.

# Arguments
- `arg1`: Description
- `arg2`: Description

# Examples
```julia
function_name(value1, value2)
```
"""
```

### Error Handling

Use descriptive error messages with `throw`:

```julia
if !haskey(structure.arrays, arrayname)
    available = join(keys(structure.arrays), ", ")
    throw(ArgumentError("Array '$arrayname' not found in Cell. Available arrays: $available"))
end
```

Use `@assert` for internal invariants:

```julia
@assert length(symbols) == size(positions, 2)
```

### Struct Definitions

Use parametric types for flexibility:

```julia
mutable struct Cell{T,D} <: AB.AbstractSystem{D}
    lattice::Lattice{T}
    symbols::Vector{Symbol}
    positions::Matrix{T}
    arrays::Dict{Symbol,Any}
    metadata::Dict{Symbol,Any}
end
```

### Function Patterns

- Use multiple dispatch to handle different input types
- Destructuring for tuple returns: `rformula, num_fu = formula_and_factor(structure)`
- Broadcasting with dot syntax: `scaled .-= floor.(scaled)`
- In-place operations with `!` suffix: `sort!`, `wrap!`, `set_positions!`

### Testing

Test files use `@testset` blocks with descriptive names:

```julia
using Test
using CellBase

@testset "Feature Name" begin
    @testset "Sub-feature" begin
        @test condition
        @test_throws ErrorType expression
    end
end
```

### Git Workflow

- CI runs on pushes to `master` and `static` branches
- PRs should target `master`
- Tests must pass with coverage reporting to Codecov

## Key Dependencies

- `AtomsBase` - Abstract system interface
- `Spglib` - Space group operations
- `StaticArrays` - Performance-critical arrays
- `PeriodicTable` - Element data
- `LinearAlgebra` - Matrix operations
- `Unitful` - Physical units (in AtomsBase integration)
- `Parameters` - Type-based parameter handling

## Important Implementation Notes

### Hyperdimensional Support
The package supports N-dimensional positions (not just 3D) - this is used for hyperdimensional relaxation in structure prediction. When working with hyperdimensional cells, be aware that:
- `size(positions, 1)` gives the dimension count
- Some operations (like rotation) may only work correctly on 3D subsets
- Use `remove_dimensions()` to extract 3D structure from hyperdimensional

### Reciprocal Lattice Updates
When modifying `lattice.matrix` directly, always call `update_rec!(lattice)` to recompute the reciprocal lattice. The `set_cellmat!()` function handles this automatically.

### Array Data Storage
Additional per-atom data (forces, charges, etc.) is stored in `cell.arrays::Dict{Symbol,Any}`. When cloning/clipping cells, ensure these arrays are handled correctly - they should be sliced along their last dimension (the atom dimension).
