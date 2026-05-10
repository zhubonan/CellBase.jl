#=
Code for input/output

Format:

- SHELX: read and write in AIRSS style
- CELL: read and write
- XYZ: read and write via ExtXYZ.jl
- .castep: read only, including energies and forces
- STRU: read and write ABACUS STRU format
=#

include("io_cell.jl")
include("io_res.jl")
include("io_xyz.jl")
include("io_sheap.jl")
include("io_dotcastep.jl")
include("io_poscar.jl")
include("io_stru.jl")

export read_res, write_res, read_cell, write_xyz, read_stru, write_stru
