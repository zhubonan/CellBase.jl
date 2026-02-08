# AGENTS.md - Agentic Coding Guidelines for CellBase.jl

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
# Format code using JuliaFormatter (if installed)
julia --project -e 'using JuliaFormatter; format("src")'
julia --project -e 'using JuliaFormatter; format("test")'
```

Formatter configuration is in `.JuliaFormatter.toml`:
- `ignore = ["scripts"]` - Skip scripts directory
- `whitespace_in_kwargs = false` - No whitespace around keyword argument equals

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

### File Organization

- `src/CellBase.jl` - Main module with includes
- `src/<feature>.jl` - Individual feature files
- `src/io/` - File I/O operations
- `src/external/` - Third-party integrations
- `test/test_<feature>.jl` - Corresponding test files
- `test/runtests.jl` - Test entry point with includes

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

## Notes

- The package supports N-dimensional positions (not just 3D)
- Uses StaticArrays (`SVector`, `SMatrix`) for performance
- Follows AtomsBase interface for interoperability
- Maintains compatibility with Julia 1.6+
