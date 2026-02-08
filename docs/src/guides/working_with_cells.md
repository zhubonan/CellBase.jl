# [Working with Cells](@id working-with-cells)

This guide covers common operations on `Cell` objects, including creating supercells, sorting atoms, rotating structures, and more.

```@setup cells
using CellBase
```

## Creating Supercells

### Simple Repetition

The easiest way to create a supercell is using the `repeat` function:

```@example cells
# Create a simple cubic cell
lattice = Lattice(4.0, 4.0, 4.0)
symbols = [:Cu, :Cu, :Cu, :Cu]
pos_mat = [
    0.0  0.0  2.0  2.0;
    0.0  2.0  0.0  2.0;
    0.0  2.0  2.0  0.0
]
cell = Cell(lattice, symbols, pos_mat)
println("Original cell has ", natoms(cell), " atoms")
println("Atom positions:")
for i in 1:natoms(cell)
    println("  Atom ", i, " (", species(cell)[i], "): ", positions(cell)[:, i])
end

# Create a 2×2×2 supercell
supercell = repeat(cell, 2, 2, 2)
println("\nSupercell has ", natoms(supercell), " atoms")
```

### General Transformation Matrix

For more complex supercells (including shear transformations), use `make_supercell`:

```@example cells
# 2×2×2 supercell via transformation matrix
P = [2 0 0; 0 2 0; 0 0 2]
supercell = make_supercell(cell, P)
println("Transformation supercell has ", natoms(supercell), " atoms")
println("First 4 atom positions:")
for i in 1:4
    println("  Atom ", i, ": ", positions(supercell)[:, i])
end
```

## Sorting Atoms

Sort atoms by element symbol, atomic number, or position:

```@example cells
# Create a mixed cell
lattice = Lattice(5.0, 5.0, 5.0)
symbols = [:O, :H, :H, :O, :H, :H]  # Two water molecules
pos_mat = rand(3, 6) .* 5.0
cell = Cell(lattice, symbols, pos_mat)

println("Before sorting:")
println("Species: ", species(cell))
for i in 1:min(3, natoms(cell))
    println("  Atom ", i, " (", species(cell)[i], "): ", positions(cell)[:, i])
end

# Sort by element symbol (default)
sorted_cell = sort(cell)
println("\nAfter sorting by symbol:")
println("Species: ", species(sorted_cell))
for i in 1:min(3, natoms(sorted_cell))
    println("  Atom ", i, " (", species(sorted_cell)[i], "): ", positions(sorted_cell)[:, i])
end
```

## Rotating Structures

Rotate atoms around an axis or align vectors:

```@example cells
# Create a simple cell for rotation
cell = Cell(Lattice(5.0, 5.0, 5.0), [:H, :H], [0.0 1.0; 0.0 0.0; 0.0 0.0])
println("Original positions:")
for i in 1:natoms(cell)
    println("  Atom ", i, " (", species(cell)[i], "): ", positions(cell)[:, i])
end

# Rotate 90 degrees around z-axis
rotated = rotate(cell, 90, "z")
println("\nAfter 90° rotation around z:")
for i in 1:natoms(rotated)
    println("  Atom ", i, " (", species(rotated)[i], "): ", positions(rotated)[:, i])
end
```

## Wrapping Positions

Ensure all atoms are within the unit cell:

```@example cells
# Create a cell with atoms outside
lattice = Lattice(5.0, 5.0, 5.0)
symbols = [:H, :H]
pos_mat = [6.0 8.0; 6.0 7.0; 6.0 9.0]  # Outside the cell
cell = Cell(lattice, symbols, pos_mat)

println("Before wrapping:")
for i in 1:natoms(cell)
    println("  Atom ", i, ": ", positions(cell)[:, i])
end

# Wrap all atoms back into the cell
wrap!(cell)
println("\nAfter wrapping:")
for i in 1:natoms(cell)
    println("  Atom ", i, ": ", positions(cell)[:, i])
end
```

## Selecting Atoms

Extract a subset of atoms:

```@example cells
# Create a cell with 6 atoms
lattice = Lattice(10.0, 10.0, 10.0)
symbols = [:Si for _ in 1:6]
pos_mat = rand(3, 6) .* 10.0
cell = Cell(lattice, symbols, pos_mat)
println("Original cell has ", natoms(cell), " atoms")

# Show first 3 atoms
println("First 3 atom positions:")
for i in 1:3
    println("  Atom ", i, ": ", positions(cell)[:, i])
end

# Select first 3 atoms
small_cell = cell[1:3]
println("\nSelected cell has ", natoms(small_cell), " atoms")
println("Selected atom positions:")
for i in 1:natoms(small_cell)
    println("  Atom ", i, ": ", positions(small_cell)[:, i])
end
```

## Working with Additional Arrays

Store and retrieve additional per-atom data:

```@example cells
# Create a cell
cell = Cell(Lattice(5.0, 5.0, 5.0), [:H, :H], [0.0 2.5; 0.0 2.5; 0.0 2.5])

# Store forces (3×N matrix)
cell.arrays[:forces] = rand(3, 2)

# Store charges (N-vector)
cell.arrays[:charges] = [0.5, -0.5]

# Access arrays - show forces for each atom
println("Forces on atoms:")
for i in 1:natoms(cell)
    println("  Atom ", i, " (", species(cell)[i], "): ", cell.arrays[:forces][:, i])
end
println("\nCharges: ", array(cell, :charges))

# List all array names
println("\nAvailable arrays: ", collect(arraynames(cell)))

# Metadata storage
cell.metadata[:energy] = -10.5
cell.metadata[:pressure] = 1.0
println("Metadata keys: ", keys(cell.metadata))
```

## Computing Distances

Calculate distances between atoms:

```@example cells
# Create a cell with known positions
lattice3 = Lattice(5.0, 5.0, 5.0)
symbols3 = [:H, :H]
pos_mat3 = [0.0 1.0; 0.0 0.0; 0.0 0.0]  # 1 Å apart
cell3 = Cell(lattice3, symbols3, pos_mat3)

println("Atom positions:")
for i in 1:natoms(cell3)
    println("  Atom ", i, ": ", positions(cell3)[:, i])
end

# Compute distance matrix (with minimum image convention)
dmat = distance_matrix(cell3)
println("\nDistance between atoms: ", dmat[1, 2], " Å")
```

## Formula and Composition

Get chemical formula information:

```@example cells
# Create a compound cell
lattice2 = Lattice(5.0, 5.0, 5.0)
symbols2 = [:Na, :Cl, :Na, :Cl]
pos_mat2 = rand(3, 4) .* 5.0
cell2 = Cell(lattice2, symbols2, pos_mat2)

println("Atom positions:")
for i in 1:natoms(cell2)
    println("  Atom ", i, " (", species(cell2)[i], "): ", positions(cell2)[:, i])
end

# Reduced formula
formula = reduced_fu(cell2)
println("\nReduced formula: ", formula)

# Number of formula units
n = num_fu(cell2)
println("Number of formula units: ", n)

# Get atomic masses
masses_vec = masses(cell2)
println("Atomic masses: ", masses_vec)
```

## Rattle (Random Displacements)

Add random displacements to atomic positions:

```@example cells
using Random

# Create a cell
cell_rattle = Cell(Lattice(5.0, 5.0, 5.0), [:H, :H], [0.0 2.0; 0.0 0.0; 0.0 0.0])
println("Original positions:")
for i in 1:natoms(cell_rattle)
    println("  Atom ", i, ": ", positions(cell_rattle)[:, i])
end

# Random displacement up to 0.1 Å (set seed for reproducibility)
Random.seed!(42)
rattle!(cell_rattle, 0.1)
println("\nAfter rattling:")
for i in 1:natoms(cell_rattle)
    println("  Atom ", i, ": ", positions(cell_rattle)[:, i])
end
```

## Modifying Lattice

Change the cell dimensions:

```@example cells
# Create a cell
cell4 = Cell(Lattice(5.0, 5.0, 5.0), [:H], reshape([2.5, 2.5, 2.5], 3, 1))
println("Atom position: ", positions(cell4)[:, 1])
println("Original lattice parameters: ", cellpar(cell4))
println("Original volume: ", volume(cell4))

# Set new lattice matrix (scale positions to maintain fractional coords)
new_mat = [6.0 0.0 0.0; 0.0 6.0 0.0; 0.0 0.0 6.0]
set_cellmat!(cell4, new_mat)
println("\nAfter scaling:")
println("Atom position: ", positions(cell4)[:, 1])
println("New lattice parameters: ", cellpar(cell4))
println("New volume: ", volume(cell4))
```
