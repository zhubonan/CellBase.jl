#=
Tests for IO related codes
=#

using CellBase
using Test

@testset "IO" begin
    this_dir = splitpath(@__FILE__)[1:end-1]
    @testset "SHELX" begin
        fpath = joinpath(this_dir..., "lno.res")
        cell = CellBase.read_res(fpath)
        @test species(cell) == [:Li, :Li, :O, :O, :O, :O, :Ni, :Ni]
        @test cellpar(cell)[1] == 2.73724
        @test cellpar(cell)[2] == 6.02753

        # Read packed res
        fpath = joinpath(this_dir..., "lno.pack.res")
        cells = CellBase.read_res_many(fpath)
        @test length(cells) == 2
        @test cellpar(cells[1])[1] == 2.73724
        @test cellpar(cells[1])[2] == 6.02753

        mktempdir() do tempd
            cell.metadata[:comments] = ["COMMENT1", "COMMENT2"]
            cell.metadata[:info] = Dict("A" => 2, "B" => 3)
            write_res(joinpath(tempd, "test.res"), cell)
            written = readlines(joinpath(tempd, "test.res"))
            @test "REM COMMENT1" in written
            @test "REM COMMENT2" in written
            @test "REM A: 2" in written
            @test "REM Composition: Li 2.0 Ni 2.0 O 4.0" in written
        end
    end

    @testset "CELL" begin
        cell = CellBase.read_cell(joinpath(this_dir..., "Si2.cell"))
        @test species(cell) == [:Si, :Si]
        @test cellmat(cell)[1] ≈ 2.69546455 atol = 1e-5

        # Test writing CEll file
        pos = rand(3, 3)
        latt = rand(3, 3)
        tmp, f = mktemp()
        CellBase.CellIO.write_cell(f, latt, pos, [:H, :H, :H])
        close(f)
        celltmp = CellBase.read_cell(tmp)
        @test all(cellmat(celltmp) .== latt)
        @test all(positions(celltmp) .== pos)
        rm(tmp)

        CellBase.write_cell(tmp, celltmp)
        celltmp2 = CellBase.read_cell(tmp)
        @test all(cellmat(celltmp2) .== cellmat(celltmp))
    end

    @testset "POSCAR" begin
        cell = Cell(Lattice(10, 10, 10), [:H, :O, :H], [[0, 0.0, 0], [1, 0, 0], [2, 0, 0]])

        mktemp() do path, io
            write_poscar(io, cell)
            seekstart(io)
            lines = readlines(io)
            @test any(x -> startswith(x, "   H  O"), lines)
            @test any(x -> startswith(x, "   2  1"), lines)
            @test any(x -> startswith(x, "Direct"), lines)
            @test any(x -> startswith(x, "    0.00"), lines)
            @test any(x -> startswith(x, "    0.20"), lines)
            @test any(x -> startswith(x, "    0.10"), lines)
            @test any(x -> startswith(x, "     10.00"), lines)
            close(io)
            write_poscar(path, cell)
            @test any(x -> startswith(x, "    0.20"), lines)
            @test any(x -> startswith(x, "    0.10"), lines)
            @test any(x -> startswith(x, "     10.00"), lines)
            newcell = read_poscar(path)
            @test cellmat(newcell) == cellmat(cell)
            @test species(newcell) == [:H, :H, :O]
            @test positions(newcell)[1, 1] == 0.0
            @test positions(newcell)[1, 2] == 2.0
            @test positions(newcell)[1, 3] == 1.0
        end
    end

    @testset "Hyperdimensional writers" begin
        hyper = Cell(
            Lattice(10.0, 10.0, 10.0),
            [:H, :He],
            [
                1.0 2.0
                3.0 4.0
                5.0 6.0
                0.25 -0.5
            ],
        )
        hyper.arrays[:forces] = [
            0.1 0.2
            0.3 0.4
            0.5 0.6
        ]
        hyper.metadata[:label] = "hyper"

        mktemp() do path, io
            @test_logs (:warn, r"write_poscar drops auxiliary dimensions") write_poscar(io, hyper)
            seekstart(io)
            lines = readlines(io)
            @test any(x -> occursin("0.1000000000000000", x), lines)
            @test any(x -> occursin("0.2000000000000000", x), lines)
        end

        mktempdir() do tempd
            @test_logs (:warn, r"write_res drops auxiliary dimensions") write_res(joinpath(tempd, "hyper.res"), hyper)
            @test_logs (:warn, r"write_cell drops auxiliary dimensions") CellBase.write_cell(joinpath(tempd, "hyper.cell"), hyper)
        end

        mktemp() do path, io
            CellBase.push_xyz!(io, hyper)
            seekstart(io)
            lines = readlines(io)
            @test occursin("extra_dim_1:R:1", lines[2])
            @test parse(Float64, split(lines[3])[end]) ≈ 0.25
            @test parse(Float64, split(lines[4])[end]) ≈ -0.5
        end

        mktemp() do path, io
            CellBase.write_xyz(path, [hyper])
            xyz_cells = CellBase.read_xyz(path)
            @test length(xyz_cells) == 1
            @test size(positions(xyz_cells[1])) == (4, 2)
            @test positions(xyz_cells[1])[1:3, :] == positions(hyper)[1:3, :]
            @test positions(xyz_cells[1])[4, :] == positions(hyper)[4, :]
            @test CellBase.array(xyz_cells[1], :forces) == CellBase.array(hyper, :forces)
            @test CellBase.metadata(xyz_cells[1])[:label] == "hyper"
        end
    end

    @testset "CASTEP" begin
        snapshots =
            CellBase.read_castep(joinpath(this_dir..., "Fe.castep"), only_first=false)
        @test length(snapshots) == 11
        snap = snapshots[1]
        @test snap.forces[1] ≈ -3.20036
    end

    include("test_stru.jl")

end
