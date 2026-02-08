# [Building Structures](@id building-structures)

This guide covers creating common crystal structures using the `bulk` function, which provides an easy way to generate standard crystal structures.

```@setup building
using CellBase
import CellBase: bulk
```

## Basic Usage

The `bulk` function creates crystal structures from chemical formulas and structure types:

```@example building
using CellBase

# Create FCC copper (auto-detects from reference database)
cu = CellBase.bulk("Cu")
println("Copper structure:")
println("Formula: ", reduced_fu(cu))
println("Number of atoms: ", natoms(cu))
println("Volume: ", volume(cu), " Å³")
```

```@example building
# Create BCC iron
fe = bulk("Fe")
println("Iron structure:")
println("Formula: ", reduced_fu(fe))
println("Number of atoms: ", natoms(fe))
```

```@example building
# Create HCP magnesium
mg = bulk("Mg")
println("Magnesium structure:")
println("Formula: ", reduced_fu(mg))
println("Number of atoms: ", natoms(mg))
println("Lattice parameters: ", cellpar(mg))
```

## Supported Structures

### Single Element Structures

| Structure | Symbol | Atoms | Description |
|-----------|--------|-------|-------------|
| Simple Cubic | `:sc` | 1 | Primitive simple cubic |
| FCC | `:fcc` | 1 | Face-centered cubic |
| BCC | `:bcc` | 1 | Body-centered cubic |
| HCP | `:hcp` | 2 | Hexagonal close-packed |
| Diamond | `:diamond` | 2 | Diamond cubic |
| Tetragonal | `:tetragonal` | 1 | Simple tetragonal |
| BCT | `:bct` | 1 | Body-centered tetragonal |
| Rhombohedral | `:rhombohedral` | 1 | Rhombohedral |
| Orthorhombic | `:orthorhombic` | 1 | Orthorhombic |

### Compound Structures

| Structure | Symbol | Atoms | Description |
|-----------|--------|-------|-------------|
| Zincblende | `:zincblende` | 2 | ZnS structure |
| Rocksalt | `:rocksalt` | 2 | NaCl structure |
| Cesium Chloride | `:cesiumchloride` | 2 | CsCl structure |
| Fluorite | `:fluorite` | 3 | CaF₂ structure |
| Wurtzite | `:wurtzite` | 4 | ZnS wurtzite |

## Creating Specific Structures

### FCC, BCC, and HCP

```@example building
# FCC with explicit lattice constant
si = bulk("Si", "fcc", a=5.43)
println("Silicon FCC:")
println("  Atoms: ", natoms(si))
println("  Volume: ", volume(si))
```

```@example building
# BCC
fe = bulk("Fe", "bcc", a=2.87)
println("Iron BCC:")
println("  Atoms: ", natoms(fe))
println("  Volume: ", volume(fe))
```

```@example building
# HCP (uses ideal c/a ratio by default)
mg = bulk("Mg", "hcp", a=3.21)
println("Magnesium HCP:")
println("  Atoms: ", natoms(mg))
println("  Lattice parameters: ", cellpar(mg))
```

```@example building
# HCP with explicit c/a ratio
ti = bulk("Ti", "hcp", a=2.95, covera=1.587)
println("Titanium HCP with c/a ratio:")
a, b, c, α, β, γ = cellpar(ti)
println("  a = ", a, ", c = ", c)
println("  c/a = ", c/a)
```

### Diamond Structure

```@example building
# Diamond cubic
c = bulk("C", "diamond", a=3.57)
println("Carbon diamond:")
println("  Atoms: ", natoms(c))
println("  Volume: ", volume(c))
```

### Binary Compounds

```@example building
# Rocksalt (NaCl)
nacl = bulk("NaCl", "rocksalt", a=5.64)
println("NaCl rocksalt:")
println("  Formula: ", reduced_fu(nacl))
println("  Atoms: ", natoms(nacl))
println("  Volume: ", volume(nacl))
```

```@example building
# Zincblende (GaAs)
gaas = bulk("GaAs", "zincblende", a=5.65)
println("GaAs zincblende:")
println("  Formula: ", reduced_fu(gaas))
println("  Atoms: ", natoms(gaas))
```

```@example building
# Cesium chloride structure
mgte = bulk("MgTe", "cesiumchloride", a=3.27)
println("MgTe cesium chloride:")
println("  Formula: ", reduced_fu(mgte))
println("  Atoms: ", natoms(mgte))
```

```@example building
# Fluorite (CaF₂)
caf2 = bulk("CaFF", "fluorite", a=5.46)
println("CaF₂ fluorite:")
println("  Formula: ", reduced_fu(caf2))
println("  Atoms: ", natoms(caf2))
```

## Conventional vs Primitive Cells

Create conventional (cubic/orthorhombic) cells instead of primitive cells:

```@example building
# FCC primitive cell (1 atom)
primitive = bulk("Cu", "fcc", a=3.61)
println("FCC primitive: ", natoms(primitive), " atoms")

# FCC conventional cubic cell (4 atoms)
conventional = bulk("Cu", "fcc", a=3.61, cubic=true)
println("FCC conventional: ", natoms(conventional), " atoms")
```

```@example building
# BCC conventional cell (2 atoms)
bcc_conv = bulk("Fe", "bcc", a=2.87, cubic=true)
println("BCC conventional: ", natoms(bcc_conv), " atoms")
```

```@example building
# Orthorhombic cells
ortho_fcc = bulk("Cu", "fcc", a=3.61, orthorhombic=true)
println("FCC orthorhombic: ", natoms(ortho_fcc), " atoms")
```

## Working with Created Structures

Once created, you can manipulate structures as normal `Cell` objects:

```@example building
# Create a structure
si = bulk("Si", "diamond", a=5.43)
println("Original silicon:")
println("  Atoms: ", natoms(si))

# Create a supercell
supercell = repeat(si, 2, 2, 2)  # 2×2×2 supercell
println("Supercell:")
println("  Atoms: ", natoms(supercell))

# Get information
println("Formula: ", reduced_fu(si))
println("Lattice parameters: ", cellpar(si))
```

## Complete Examples

### Example 1: Creating a Surface Slab

```@example building
# Create bulk MgO
mgo = bulk("MgO", "rocksalt", a=4.21)
println("Bulk MgO: ", natoms(mgo), " atoms")

# Make a 2×2×4 supercell for surface slab
slab = repeat(mgo, 2, 2, 4)
println("Slab: ", natoms(slab), " atoms")
```

### Example 2: Creating a Supercell for Defects

```@example building
# Create bulk silicon
si = bulk("Si", "diamond", a=5.43)
println("Primitive Si: ", natoms(si), " atoms")

# 3×3×3 supercell (216 atoms)
super_si = repeat(si, 3, 3, 3)
println("Supercell: ", natoms(super_si), " atoms")

# Sort atoms
sort!(super_si, by=:symbol)
println("Sorted supercell ready for defect study")
```

### Example 3: Comparing Structures

```@example building
# Create FCC and BCC iron
fe_fcc = bulk("Fe", "fcc", a=3.57)
fe_bcc = bulk("Fe", "bcc", a=2.87)

# Compare volumes
println("FCC volume: ", round(volume(fe_fcc), digits=3), " Å³")
println("BCC volume: ", round(volume(fe_bcc), digits=3), " Å³")
println("Volume difference: ", round(volume(fe_fcc) - volume(fe_bcc), digits=3), " Å³")
```
