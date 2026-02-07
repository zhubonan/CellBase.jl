#=
Tests for STRU file I/O
=#

using CellBase
using Test

const BOHR_TO_ANGSTROM = 0.52917721092

@testset "STRU" begin
    this_dir = splitpath(@__FILE__)[1:end-1]

    @testset "Read Direct coordinates" begin
        fpath = joinpath(this_dir..., "Si.stru")
        cell = CellBase.read_stru(fpath)

        @test species(cell) == [:Si, :Si]
        @test length(species(cell)) == 2

        lat = cellmat(cell)
        expected_lat = 10.2 * BOHR_TO_ANGSTROM * [0.5 0.5 0.0; 0.5 0.0 0.5; 0.0 0.5 0.5]
        @test lat ≈ expected_lat atol = 1e-6

        pos = positions(cell)
        expected_pos = expected_lat * [0.0 0.25; 0.0 0.25; 0.0 0.25]
        @test pos ≈ expected_pos atol = 1e-6

        @test get(cell.metadata, :lattice_constant_bohr, nothing) == 10.2
    end

    @testset "Read Cartesian coordinates" begin
        fpath = joinpath(this_dir..., "Si_cart.stru")
        cell = CellBase.read_stru(fpath)

        @test species(cell) == [:Si, :Si]
        @test length(species(cell)) == 2

        lat = cellmat(cell)
        expected_lat = 10.2 * BOHR_TO_ANGSTROM * [0.5 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        @test lat ≈ expected_lat atol = 1e-6

        pos = positions(cell)
        expected_pos = [0.0 (2.70 * 10.2 * BOHR_TO_ANGSTROM); 0.0 0.0; 0.0 0.0]
        @test pos ≈ expected_pos atol = 1e-5
    end

    @testset "Read Cartesian_angstrom" begin
        fpath = joinpath(this_dir..., "Si_ang.stru")
        cell = CellBase.read_stru(fpath)

        @test species(cell) == [:Si, :Si]
        @test length(species(cell)) == 2

        lat = cellmat(cell)
        expected_lat = 10.2 * BOHR_TO_ANGSTROM * [0.5 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        @test lat ≈ expected_lat atol = 1e-6

        pos = positions(cell)
        expected_pos = [0.0 5.4; 0.0 0.0; 0.0 0.0]
        @test pos ≈ expected_pos atol = 1e-5
    end

    @testset "Write-Read round-trip" begin
        lat = Lattice(5.43, 5.43, 5.43)
        cell_species = [:Si, :O, :O]
        pos = [0.0 1.36 2.72; 0.0 1.36 1.36; 0.0 0.0 2.72]
        original_cell = Cell(lat, cell_species, pos)

        mktempdir() do tempd
            fpath = joinpath(tempd, "test.stru")
            write_stru(fpath, original_cell)

            read_cell = read_stru(fpath)

            @test length(species(read_cell)) == length(species(original_cell))
            @test Set(species(read_cell)) == Set(species(original_cell))
            @test cellmat(read_cell) ≈ cellmat(original_cell) atol = 1e-6

            pos_orig = positions(sort(original_cell))
            pos_read = positions(sort(read_cell))

            @test pos_read ≈ pos_orig atol = 1e-4
        end
    end

    @testset "Write structure" begin
        lat = Lattice(5.43, 5.43, 5.43)
        cell_species = [:Si, :O]
        pos = [0.0 1.36; 0.0 1.36; 0.0 0.0]
        cell = Cell(lat, cell_species, pos)

        mktemp() do path, io
            write_stru(io, cell)
            seekstart(io)
            lines = readlines(io)

            @test any(x -> contains(x, "ATOMIC_SPECIES"), lines)
            @test any(x -> contains(x, "LATTICE_CONSTANT"), lines)
            @test any(x -> contains(x, "LATTICE_VECTORS"), lines)
            @test any(x -> contains(x, "ATOMIC_POSITIONS"), lines)
            @test any(x -> contains(x, "Direct"), lines)

            close(io)
        end
    end
end
