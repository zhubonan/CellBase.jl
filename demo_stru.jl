#!/usr/bin/env julia
#=
Demo of STRU file I/O functionality in CellBase.jl
=#

using CellBase

const BOHR_TO_ANGSTROM = 0.52917721092

println("=" ^ 70)
println("CellBase.jl STRU File I/O Demo")
println("=" ^ 70)
println()

# ============================================================================
# Demo 1: Creating a simple structure and writing to STRU
# ============================================================================
println("Demo 1: Creating a simple SiO₂ structure and writing to STRU")
println("-" ^ 70)

# Create a simple SiO₂ structure
lat = Lattice(9.0, 9.0, 9.0)
cell_species = [:Si, :O, :O, :O]
# Positions as column vectors: [x y z]^T for each atom
atom_positions = [[0.0, 1.0, 2.0], [2.0, 3.0, 4.0], [4.0, 5.0, 6.0], [6.0, 7.0, 8.0]]
pos = hcat(atom_positions...)
sio2_cell = Cell(lat, cell_species, pos)

println("Created SiO₂ cell:")
println("  Lattice:", cellmat(sio2_cell))
println("  Species: ", species(sio2_cell))
println("  Positions:", positions(sio2_cell))
println()

# Write to file
stru_filename = "demo_sio2.stru"
write_stru(stru_filename, sio2_cell)
println("✓ Written to: $stru_filename")
println()

# Show the file content
println("File content:")
println("-" ^ 70)
println(read(stru_filename, String))
println("-" ^ 70)
println()

# ============================================================================
# Demo 2: Reading STRU file
# ============================================================================
println("Demo 2: Reading STRU file")
println("-" ^ 70)

read_cell = read_stru(stru_filename)

println("Read from $stru_filename:")
println("  Lattice:", cellmat(read_cell))
println("  Species: ", species(read_cell))
println("  Positions:", positions(read_cell))
println()

# Verify round-trip
println("Round-trip verification:")
println("  Lattice match: ", isapprox(cellmat(read_cell), cellmat(sio2_cell)))
println("  Species match: ", Set(species(read_cell)) == Set(cell_species))
println("  Atom count match: ", length(species(read_cell)) == length(cell_species))
println("  Positions match (after sorting): ", isapprox(positions(sort(read_cell)), positions(sort(sio2_cell)), rtol=1e-3))
println()

# ============================================================================
# Demo 3: Reading STRU files with different coordinate systems
# ============================================================================
println("Demo 3: Reading STRU with different coordinate systems")
println("-" ^ 70)

# Read test STRU files with different coordinate types
println("Reading Direct coordinates...")
cell_direct = read_stru("test/io_tests/Si.stru")
println("  ✓ Species: ", unique(species(cell_direct)))
println("  ✓ Volume: ", round(volume(cell_direct), digits=2), " Å³")
println()

println("Reading Cartesian coordinates (Bohr)...")
cell_cart = read_stru("test/io_tests/Si_cart.stru")
println("  ✓ Species: ", unique(species(cell_cart)))
println("  ✓ Volume: ", round(volume(cell_cart), digits=2), " Å³")
println()

println("Reading Cartesian_angstrom...")
cell_ang = read_stru("test/io_tests/Si_ang.stru")
println("  ✓ Species: ", unique(species(cell_ang)))
println("  ✓ Volume: ", round(volume(cell_ang), digits=2), " Å³")
println()

# ============================================================================
# Demo 4: Converting structures
# ============================================================================
println("Demo 4: Structure manipulation and STRU output")
println("-" ^ 70)

# Scale a structure
println("Scaling lattice by factor of 2...")
si_cell = read_stru("test/io_tests/Si.stru")
original_lat = cellmat(si_cell)
scaled_lat = original_lat * 2.0
scaled_cell = Cell(Lattice(scaled_lat), species(si_cell), positions(si_cell) * 2.0)
println("  Original volume: ", round(volume(si_cell), digits=2), " Å³")
println("  Scaled volume: ", round(volume(scaled_cell), digits=2), " Å³")
println()

# Write scaled structure
write_stru("scaled_si.stru", scaled_cell)
println("✓ Scaled structure written to: scaled_si.stru")
println()

# ============================================================================
# Demo 5: Working with metadata
# ============================================================================
println("Demo 5: Accessing STRU metadata")
println("-" ^ 70)

# Read a STRU file
test_cell = read_stru("test/io_tests/Si.stru")

println("Metadata from STRU file:")
if haskey(test_cell.metadata, :lattice_constant_bohr)
    println("  ✓ Lattice constant: ", test_cell.metadata[:lattice_constant_bohr], " Bohr")
    println("  ✓ In Angstrom: ", test_cell.metadata[:lattice_constant_bohr] * BOHR_TO_ANGSTROM, " Å")
end
println()

# ============================================================================
# Demo 6: Formatting output structure
# ============================================================================
println("Demo 6: Customizing STRU output")
println("-" ^ 70)

# Sort the structure by species
sorted_cell = sort(sio2_cell)
println("Original species order: ", species(sio2_cell))
println("Sorted species order: ", species(sorted_cell))
println()

# Write sorted structure
write_stru("sorted_sio2.stru", sorted_cell)
println("✓ Sorted structure written to: sorted_sio2.stru")
println()

# ============================================================================
# Demo 7: Multiple species in STRU
# ============================================================================
println("Demo 7: Structure with multiple species")
println("-" ^ 70)

# Create a mixed-species structure
lat2 = Lattice(10.0, 10.0, 10.0)
cell_species2 = [:Fe, :Ni, :O, :O, :O]
# Positions as column vectors: [x y z]^T for each atom
atom_positions = [[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [5.0, 0.0, 0.0], [7.5, 0.0, 0.0], [0.0, 2.5, 0.0]]
pos2 = hcat(atom_positions...)
mixed_cell = Cell(lat2, cell_species2, pos2)

println("Mixed species structure:")
println("  Species: ", species(mixed_cell))
println("  Composition: ", unique(species(mixed_cell)))
println("  Atom counts:")
for sp in unique(species(mixed_cell))
    count = length(findall(x -> x == sp, species(mixed_cell)))
    println("    $sp: $count")
end
println()

write_stru("mixed_species.stru", mixed_cell)
println("✓ Written to: mixed_species.stru")
println()

println("File content:")
println("-" ^ 70)
println(read("mixed_species.stru", String))
println("-" ^ 70)
println()

# ============================================================================
# Cleanup
# ============================================================================
println("=" ^ 70)
println("Demo completed!")
println("=" ^ 70)
println()
println("Generated files:")
println("  - demo_sio2.stru")
println("  - scaled_si.stru")
println("  - sorted_sio2.stru")
println("  - mixed_species.stru")
println()
