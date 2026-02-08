using Test
using CellBase

@testset "Build Functions" begin
    
    @testset "Auto-detection from reference database" begin
        # Test auto-lookup for various elements
        
        # FCC elements - test with cubic=true for easy cellpar check
        cu = bulk("Cu", cubic=true)
        @test natoms(cu) == 4
        cp = cellpar(cu)
        @test isapprox(cp[1], 3.61, rtol=0.01)  # Reference value
        
        al = bulk("Al", cubic=true)
        @test natoms(al) == 4
        cp_al = cellpar(al)
        @test isapprox(cp_al[1], 4.05, rtol=0.01)
        
        # BCC elements - test with cubic=true
        fe = bulk("Fe", cubic=true)
        @test natoms(fe) == 2
        cp_fe = cellpar(fe)
        @test isapprox(cp_fe[1], 2.87, rtol=0.01)
        
        # HCP elements
        mg = bulk("Mg")
        @test natoms(mg) == 2
        cp_mg = cellpar(mg)
        @test isapprox(cp_mg[1], 3.21, rtol=0.01)
        @test isapprox(cp_mg[3], 3.21 * 1.624, rtol=0.01)  # c = a * c/a
        
        # Diamond elements
        si = bulk("Si", cubic=true)
        @test natoms(si) == 8
        cp_si = cellpar(si)
        @test isapprox(cp_si[1], 5.43, rtol=0.01)
    end
    
    @testset "Explicit crystalstructure as positional argument" begin
        # Test explicit structure specification with cubic=true for easy verification
        cu_fcc = bulk("Cu", "fcc", a=3.65, cubic=true)
        @test natoms(cu_fcc) == 4
        cp = cellpar(cu_fcc)
        @test isapprox(cp[1], 3.65, rtol=0.01)
        
        fe_bcc = bulk("Fe", "bcc", a=2.90, cubic=true)
        @test natoms(fe_bcc) == 2
        
        # Test with Symbol
        mg_hcp = bulk("Mg", :hcp, a=3.25, covera=1.60)
        @test natoms(mg_hcp) == 2
    end
    
    @testset "Simple Cubic (sc)" begin
        cell = bulk("Fe", "sc", a=2.87)
        @test natoms(cell) == 1
        @test species(cell) == [:Fe]
        # Check lattice is cubic
        cp = cellpar(cell)
        @test cp[1] ≈ 2.87
        @test cp[2] ≈ 2.87
        @test cp[3] ≈ 2.87
        @test cp[4:6] ≈ [90.0, 90.0, 90.0]
    end
    
    @testset "BCC (primitive)" begin
        cell = bulk("Fe", "bcc", a=2.87)
        @test natoms(cell) == 1
        @test species(cell) == [:Fe]
    end
    
    @testset "BCC (orthorhombic)" begin
        cell = bulk("Fe", "bcc", a=2.87, orthorhombic=true)
        @test natoms(cell) == 2
        @test species(cell) == [:Fe, :Fe]
        cp = cellpar(cell)
        @test cp[1:3] ≈ [2.87, 2.87, 2.87]
    end
    
    @testset "FCC (primitive)" begin
        cell = bulk("Cu", "fcc", a=3.61)
        @test natoms(cell) == 1
        @test species(cell) == [:Cu]
    end
    
    @testset "FCC (cubic)" begin
        cell = bulk("Cu", "fcc", a=3.61, cubic=true)
        @test natoms(cell) == 4
        @test species(cell) == [:Cu, :Cu, :Cu, :Cu]
        cp = cellpar(cell)
        @test cp[1:3] ≈ [3.61, 3.61, 3.61]
        @test all(cp[4:6] .≈ 90.0)
    end
    
    @testset "HCP" begin
        cell = bulk("Mg", "hcp", a=3.21, covera=1.624)
        @test natoms(cell) == 2
        @test species(cell) == [:Mg, :Mg]
        cp = cellpar(cell)
        @test cp[1] ≈ 3.21
        @test cp[2] ≈ 3.21
        @test isapprox(cp[3], 3.21 * 1.624, rtol=1e-3)
        @test isapprox(cp[6], 120.0, atol=1e-3)
    end
    
    @testset "HCP (orthorhombic)" begin
        cell = bulk("Mg", "hcp", a=3.21, covera=1.624, orthorhombic=true)
        @test natoms(cell) == 4
        @test count(s -> s == :Mg, species(cell)) == 4
        cp = cellpar(cell)
        @test all(cp[4:6] .≈ 90.0)  # Orthorhombic has 90 degree angles
    end
    
    @testset "Diamond" begin
        cell = bulk("Si", "diamond", a=5.43)
        @test natoms(cell) == 2
        @test species(cell) == [:Si, :Si]
    end
    
    @testset "Diamond (cubic)" begin
        cell = bulk("Si", "diamond", a=5.43, cubic=true)
        @test natoms(cell) == 8
        @test count(s -> s == :Si, species(cell)) == 8
    end
    
    @testset "Zincblende" begin
        cell = bulk("ZnS", "zincblende", a=5.41)
        @test natoms(cell) == 2
        @test :Zn in species(cell)
        @test :S in species(cell)
    end
    
    @testset "Zincblende (cubic)" begin
        cell = bulk("ZnS", "zincblende", a=5.41, cubic=true)
        @test natoms(cell) == 8
        @test count(s -> s == :Zn, species(cell)) == 4
        @test count(s -> s == :S, species(cell)) == 4
    end
    
    @testset "Rocksalt (NaCl)" begin
        cell = bulk("NaCl", "rocksalt", a=5.64)
        @test natoms(cell) == 2
        @test :Na in species(cell)
        @test :Cl in species(cell)
    end
    
    @testset "Rocksalt (cubic)" begin
        cell = bulk("NaCl", "rocksalt", a=5.64, cubic=true)
        @test natoms(cell) == 8
        @test count(s -> s == :Na, species(cell)) == 4
        @test count(s -> s == :Cl, species(cell)) == 4
    end
    
    @testset "Cesium Chloride" begin
        cell = bulk("CsCl", "cesiumchloride", a=4.11)
        @test natoms(cell) == 2
        @test :Cs in species(cell)
        @test :Cl in species(cell)
    end
    
    @testset "Fluorite (CaF2)" begin
        cell = bulk("CaFF", "fluorite", a=5.46)
        @test natoms(cell) == 3
        @test count(s -> s == :Ca, species(cell)) == 1
        @test count(s -> s == :F, species(cell)) == 2
    end
    
    @testset "Wurtzite" begin
        cell = bulk("ZnO", "wurtzite", a=3.25, covera=1.60)
        @test natoms(cell) == 4
        @test count(s -> s == :Zn, species(cell)) == 2
        @test count(s -> s == :O, species(cell)) == 2
    end
    
    @testset "Tetragonal" begin
        cell = bulk("Sn", "tetragonal", a=5.83, c=3.18)
        @test natoms(cell) == 1
        @test species(cell) == [:Sn]
        cp = cellpar(cell)
        @test cp[1] ≈ 5.83
        @test cp[2] ≈ 5.83
        @test cp[3] ≈ 3.18
    end
    
    @testset "BCT (body-centered tetragonal)" begin
        cell = bulk("In", "bct", a=3.25, c=4.95)
        @test natoms(cell) == 1
        @test species(cell) == [:In]
        cp = cellpar(cell)
        @test cp[1] ≈ 3.25
        @test cp[2] ≈ 3.25
        @test cp[3] ≈ 4.95
    end
    
    @testset "Orthorhombic" begin
        cell = bulk("Ga", "orthorhombic", a=4.52, b=4.52, c=7.66)
        @test natoms(cell) == 1
        @test species(cell) == [:Ga]
        cp = cellpar(cell)
        @test cp[1:3] ≈ [4.52, 4.52, 7.66]
    end
    
    @testset "Rhombohedral" begin
        cell = bulk("Bi", "rhombohedral", a=4.75, alpha=57.14)
        @test natoms(cell) == 1
        @test species(cell) == [:Bi]
        cp = cellpar(cell)
        @test cp[1] ≈ 4.75
        @test cp[2] ≈ 4.75
        @test cp[3] ≈ 4.75
    end
    
    @testset "Reference database coverage" begin
        # Test that reference data exists for common elements
        @test !isnothing(CellBase.lookup_reference_state(:Cu))
        @test !isnothing(CellBase.lookup_reference_state(:Fe))
        @test !isnothing(CellBase.lookup_reference_state(:Al))
        @test !isnothing(CellBase.lookup_reference_state(:Si))
        @test !isnothing(CellBase.lookup_reference_state(:Ti))
        @test !isnothing(CellBase.lookup_reference_state(:W))
        @test !isnothing(CellBase.lookup_reference_state(:Au))
        
        # Test elements with no data
        @test isnothing(CellBase.lookup_reference_state(:Pm))  # Promethium
        @test isnothing(CellBase.lookup_reference_state(:Og))  # Oganesson
    end
    
    @testset "Error handling" begin
        # Missing lattice constant and no reference data
        @test_throws ArgumentError bulk("Pm")  # No reference data
        
        # Unknown structure
        @test_throws ArgumentError bulk("Cu", "unknown", a=3.0)
        
        # Wrong number of atoms for structure
        @test_throws ArgumentError bulk("Cu", "zincblende", a=3.0)
        
        # Conflicting c and c/a
        @test_throws ArgumentError bulk("Mg", "hcp", a=3.0, c=5.0, covera=1.6)
        
        # Cubic not supported for hcp
        @test_throws ArgumentError bulk("Mg", "hcp", a=3.0, cubic=true)
        
        # Molecular structures not supported
        @test_throws ArgumentError bulk("H")  # Hydrogen is diatomic
        @test_throws ArgumentError bulk("He")  # Helium is atomic
    end
    
    @testset "string_to_symbols parser" begin
        # Test the internal parser
        syms = CellBase.string_to_symbols("MgO")
        @test syms == [:Mg, :O]
        
        syms = CellBase.string_to_symbols("H2O")
        @test syms == [:H, :H, :O]
        
        syms = CellBase.string_to_symbols("CaTiO3")
        @test syms == [:Ca, :Ti, :O, :O, :O]
        
        # Test invalid element validation
        @test_throws ArgumentError CellBase.string_to_symbols("XxYz")  # Invalid symbols
        @test_throws ArgumentError CellBase.string_to_symbols("MgQ")   # Invalid: Q is not an element
    end
end