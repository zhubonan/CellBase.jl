module CellBase
import AtomsBase
import AtomsBase: periodicity
const AB = AtomsBase

# Re-export from submodules
export bulk, versioninfo

include("elemental_data.jl")
include("mathutils.jl")
include("minkowski.jl")
include("site.jl")
include("lattice.jl")
include("periodic.jl")
include("cell.jl")
include("composition.jl")
include("neighbour.jl")
include("spg.jl")
include("reference_data.jl")
include("build.jl")
include("io/io.jl")
include("external/atomsbase.jl")

end # module
