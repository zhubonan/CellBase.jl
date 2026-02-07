#!/usr/bin/env julia
#=
Simple STRU I/O Examples
=#

using CellBase

println("=" ^ 60)
println("Simple STRU I/O Examples")
println("=" ^ 60)
println()

# ============================================================================
# Example 1: Read a STRU file
# ============================================================================
println("Example 1: Read a STRU file")
println("-" ^ 60)

cell = read_stru("test/io_tests/Si.stru")

println("Structure information:")
println("  Species: ", unique(species(cell)))
println("  Number of atoms: ", nions(cell))
println("  Volume: ", round(volume(cell), digits=2), " Å³")
println()

# ============================================================================
# Example 2: Create and write a structure
# ============================================================================
println("Example 2: Create and write a structure")
println("-" ^ 60)

# Create a simple cubic cell with 4 H2O molecules
lat = Lattice(10.0, 10.0, 10.0)
h2o_species = vcat([:O for _ in 1:4], [:H for _ in 1:8])
h2o_positions = hcat([
    [0.0, 0.0, 0.0],    # O1
    [2.5, 2.5, 2.5],    # O2
    [5.0, 5.0, 5.0],    # O3
    [7.5, 7.5, 7.5],    # O4
    [0.5, 0.5, 0.5],    # H1
    [-0.5, -0.5, -0.5], # H2
    [3.0, 3.0, 3.0],    # H3
    [2.0, 2.0, 2.0],    # H4
    [5.5, 5.5, 5.5],    # H5
    [4.5, 4.5, 4.5],    # H6
    [8.0, 8.0, 8.0],    # H7
    [7.0, 7.0, 7.0],    # H8
]...)

h2o_cell = Cell(lat, h2o_species, h2o_positions)

write_stru("h2o_cluster.stru", h2o_cell)
println("✓ Written to: h2o_cluster.stru")
println()

# ============================================================================
# Example 3: Round-trip with verification
# ============================================================================
println("Example 3: Round-trip verification")
println("-" ^ 60)

# Read STRU file
original = read_stru("test/io_tests/Si.stru")

# Write to new file
write_stru("si_copy.stru", original)

# Read back
copy = read_stru("si_copy.stru")

# Verify
lat_match = isapprox(cellmat(original), cellmat(copy))
pos_match = isapprox(positions(sort(original)), positions(sort(copy)), rtol=1e-3)

println("Lattice preserved: ", lat_match ? "✓" : "✗")
println("Positions preserved: ", pos_match ? "✓" : "✗")
println()

# ============================================================================
# Example 4: Working with different coordinate systems
# ============================================================================
println("Example 4: Different coordinate systems")
println("-" ^ 60)

# Read STRU files with different coordinate systems
cell_direct = read_stru("test/io_tests/Si.stru")
cell_cart = read_stru("test/io_tests/Si_cart.stru")
cell_ang = read_stru("test/io_tests/Si_ang.stru")

println("Direct coordinates volume: ", round(volume(cell_direct), digits=2), " Å³")
println("Cartesian (Bohr) volume: ", round(volume(cell_cart), digits=2), " Å³")
println("Cartesian_angstrom volume: ", round(volume(cell_ang), digits=2), " Å³")
println()

# ============================================================================
# Example 5: Accessing metadata
# ============================================================================
println("Example 5: Accessing STRU metadata")
println("-" ^ 60)

cell = read_stru("test/io_tests/Si.stru")

if haskey(cell.metadata, :lattice_constant_bohr)
    lc_bohr = cell.metadata[:lattice_constant_bohr]
    println("Lattice constant: ", lc_bohr, " Bohr")
    println("                ", round(lc_bohr * 0.529177, digits=3), " Å")
end
println()

println("=" ^ 60)
println("Examples completed!")
println("=" ^ 60)
