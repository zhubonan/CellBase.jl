# [Spacegroup Operations](@id spacegroup-operations)

This guide covers symmetry operations and spacegroup analysis using CellBase.jl's integration with Spglib.jl.

```@setup spg
using CellBase
using Spglib
```

## Overview

CellBase.jl provides seamless integration with [Spglib.jl](https://github.com/singularitti/Spglib.jl), which wraps the [spglib](https://spglib.github.io/spglib/) library for symmetry analysis. All Spglib functions work directly with CellBase `Cell` objects.

## Quick Example

Here's a quick example of using spacegroup operations:

```@example spg
# Create a structure
cu = CellBase.bulk("Cu")

# Get space group
sym = get_international(cu)
println("Space group: ", sym)
```

## Getting Symmetry Information

### Space Group Information

Get the international space group symbol and number:

```@example spg
# Create a high-symmetry structure (NaCl)
nacl = CellBase.bulk("NaCl", "rocksalt", a=5.64)

# Get space group symbol
sym = get_international(nacl)
println("Space group: ", sym)

# Get full symmetry dataset
dataset = get_dataset(nacl)
println("Space group number: ", dataset.spacegroup_number)
println("International symbol: ", dataset.international_symbol)
println("Hall symbol: ", dataset.hall_symbol)
```

### Symmetry Operations

Get the symmetry operations (rotations and translations) for a structure:

```@example spg
# Create an FCC structure (high symmetry)
cu = CellBase.bulk("Cu")

# Get symmetry operations
rotations, translations = get_symmetry(cu)
println("Number of symmetry operations: ", size(rotations, 3))

# Show first few operations
for i in 1:min(3, size(rotations, 3))
    println("\nOperation ", i, ":")
    println("Rotation:")
    display(rotations[:, :, i])
    println("Translation: ", translations[:, i])
end
```

## Cell Standardization

### Standardize Cell

Convert a cell to its standardized form:

```@example spg
# Create a structure that might not be in standard setting
lattice = CellBase.Lattice(5.0, 5.0, 5.0)
symbols = [:Si, :Si]
pos_mat = [0.1 0.6; 0.2 0.7; 0.3 0.8]  # Non-standard positions
nonstd_cell = CellBase.Cell(lattice, symbols, pos_mat)

println("Original cell:")
for i in 1:CellBase.natoms(nonstd_cell)
    println("  Atom ", i, " (", CellBase.species(nonstd_cell)[i], "): ", CellBase.positions(nonstd_cell)[:, i])
end

# Standardize the cell
std_cell = standardize_cell(nonstd_cell)
println("\nStandardized cell:")
for i in 1:CellBase.natoms(std_cell)
    println("  Atom ", i, " (", CellBase.species(std_cell)[i], "): ", CellBase.positions(std_cell)[:, i])
end
println("Number of atoms: ", CellBase.natoms(std_cell))
```

### Find Primitive Cell

Extract the primitive cell from a conventional cell:

```@example spg
# Create conventional FCC cell (4 atoms)
cu_conv = CellBase.bulk("Cu", "fcc", a=3.61, cubic=true)
println("Conventional cell: ", CellBase.natoms(cu_conv), " atoms")
println("First atom position: ", CellBase.positions(cu_conv)[:, 1])

# Find primitive cell
prim_cell = find_primitive(cu_conv)
println("\nPrimitive cell: ", CellBase.natoms(prim_cell), " atoms")
println("First atom position: ", CellBase.positions(prim_cell)[:, 1])
```

### Refine Cell

Refine a cell to higher precision:

```@example spg
# Create a slightly distorted cell
lattice = CellBase.Lattice(5.0, 5.0, 5.0)
symbols = [:Si, :Si]
pos_mat = [0.0 0.26; 0.0 0.25; 0.0 0.25]  # Slightly off
distorted = CellBase.Cell(lattice, symbols, pos_mat)

println("Distorted cell:")
for i in 1:CellBase.natoms(distorted)
    println("  Atom ", i, ": ", CellBase.positions(distorted)[:, i])
end

# Refine the cell
refined = refine_cell(distorted)
println("\nRefined cell:")
for i in 1:CellBase.natoms(refined)
    println("  Atom ", i, ": ", CellBase.positions(refined)[:, i])
end
```

## Cell Reduction

### Niggli Reduction

Reduce a cell to its Niggli reduced form:

```@example spg
# Create a triclinic-like cell
mat = [5.0 1.0 0.5; 0.0 4.0 0.3; 0.0 0.0 3.0]
lattice = CellBase.Lattice(mat)
cell = CellBase.Cell(lattice, [:H], reshape([2.5, 2.0, 1.5], 3, 1))

println("Original lattice matrix:")
println(CellBase.cellmat(CellBase.lattice(cell)))

# Niggli reduce
reduced = niggli_reduce(cell)
println("\nReduced lattice matrix:")
println(CellBase.cellmat(CellBase.lattice(reduced)))
```

## Working Example: Symmetry Analysis Pipeline

Here's a complete workflow for analyzing crystal symmetry:

```@example spg
# Step 1: Create a structure
si = CellBase.bulk("Si", "diamond", a=5.43)
println("Created silicon diamond structure")
println("  Atoms: ", CellBase.natoms(si))
println("  First atom position: ", CellBase.positions(si)[:, 1])

# Step 2: Get symmetry information
dataset = get_dataset(si)
println("\nSymmetry analysis:")
println("  Space group: ", dataset.international_symbol)
println("  Space group number: ", dataset.spacegroup_number)

# Step 3: Check if primitive
prim_si = find_primitive(si)
if CellBase.natoms(prim_si) < CellBase.natoms(si)
    println("\nConventional cell detected")
    println("  Primitive cell has ", CellBase.natoms(prim_si), " atoms")
else
    println("\nAlready a primitive cell")
end

# Step 4: Standardize
std_si = standardize_cell(si)
println("\nStandardized cell:")
println("  Atoms: ", CellBase.natoms(std_si))
```

## Checking Symmetry

### Wyckoff Positions

Get Wyckoff position information:

```@example spg
# Create a high-symmetry structure
cu = CellBase.bulk("Cu")
dataset = get_dataset(cu)

println("Wyckoff positions:")
for i in 1:CellBase.natoms(cu)
    println("  Atom ", i, " (", CellBase.species(cu)[i], "): ", dataset.wyckoffs[i])
end

println("\nEquivalent atoms:")
println(dataset.equivalent_atoms)
```

### Crystal System

Determine the crystal system:

```@example spg
# Test different structures
structures = [
    ("Cu (FCC)", CellBase.bulk("Cu")),
    ("Fe (BCC)", CellBase.bulk("Fe", "bcc", a=2.87)),
    ("Mg (HCP)", CellBase.bulk("Mg")),
    ("Si (Diamond)", CellBase.bulk("Si", "diamond", a=5.43)),
]

for (name, struc) in structures
    dataset = get_dataset(struc)
    println(rpad(name, 15), " => ", dataset.international_symbol)
end
```

## Practical Applications

### Comparing Structures

Use symmetry to compare if two structures are equivalent:

```@example spg
# Create two FCC cells with different lattice constants
cu1 = CellBase.bulk("Cu", "fcc", a=3.61)
cu2 = CellBase.bulk("Cu", "fcc", a=3.65)

# Get their space groups
sg1 = get_international(cu1)
sg2 = get_international(cu2)

println("Cu1 space group: ", sg1)
println("Cu2 space group: ", sg2)
println("Same space group: ", sg1 == sg2)

# Compare volumes
println("\nVolume difference: ", abs(CellBase.volume(cu1) - CellBase.volume(cu2)), " Å³")
```

### Finding Unique Atoms

Find crystallographically unique atoms:

```@example spg
# Create NaCl structure
nacl = CellBase.bulk("NaCl", "rocksalt", a=5.64)
dataset = get_dataset(nacl)

println("NaCl has ", CellBase.natoms(nacl), " atoms")
println("Atom positions:")
for i in 1:min(4, CellBase.natoms(nacl))
    println("  Atom ", i, " (", CellBase.species(nacl)[i], "): ", CellBase.positions(nacl)[:, i])
end
println("\nCrystallographically unique atoms: ", length(unique(dataset.equivalent_atoms)))
println("Equivalent atom mapping: ", dataset.equivalent_atoms)
```

## Best Practices

1. **Always standardize before comparison**: Use `standardize_cell` to ensure consistent representations
2. **Check for primitive cells**: Conventional cells may have more atoms than necessary
3. **Handle symmetry tolerance**: Spglib uses a default tolerance; be aware of numerical precision
4. **Validate results**: Always check that symmetry operations make physical sense

```@example spg
# Best practice example: validating a structure
cell = CellBase.bulk("Si", "diamond", a=5.43)

println("Structure:")
for i in 1:min(2, CellBase.natoms(cell))
    println("  Atom ", i, " (", CellBase.species(cell)[i], "): ", CellBase.positions(cell)[:, i])
end

# Check space group
dataset = get_dataset(cell)
println("\nSpace group: ", dataset.international_symbol)

# Verify it's a known structure type
if dataset.spacegroup_number == 227
    println("✓ Confirmed: Diamond structure (Fd-3m)")
end

# Check primitive vs conventional
if CellBase.natoms(find_primitive(cell)) == CellBase.natoms(cell)
    println("✓ This is a primitive cell")
else
    println("✓ This is a conventional cell")
end
```
