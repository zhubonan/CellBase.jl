using LinearAlgebra

export bulk

"""
    string_to_symbols(name::String)

Parse a chemical formula string into a vector of element symbols.

# Examples
```julia
string_to_symbols("MgO")  # Returns [:Mg, :O]
string_to_symbols("NaCl") # Returns [:Na, :Cl]
```
"""
function string_to_symbols(name::String)
    # Simple parser for chemical formulas like "MgO", "H2O", "CaTiO3"
    symbols = Symbol[]
    i = 1
    while i <= length(name)
        # Get element symbol (first letter uppercase, optional second lowercase)
        if i < length(name) && islowercase(name[i+1])
            sym = Symbol(name[i:i+1])
            i += 2
        else
            sym = Symbol(name[i:i])
            i += 1
        end
        
        # Validate element symbol against PeriodicTable
        try
            _ = elements[sym]
        catch
            throw(ArgumentError("Invalid element symbol '$(sym)' in formula '$(name)'"))
        end
        
        # Parse count (if any)
        count = 0
        while i <= length(name) && isdigit(name[i])
            count = count * 10 + parse(Int, name[i])
            i += 1
        end
        count = count == 0 ? 1 : count
        
        # Add atoms
        for _ in 1:count
            push!(symbols, sym)
        end
    end
    return symbols
end

"""
    @build_unsupported_cell(want, have)

Create an error message for unsupported cell types.
"""
macro build_unsupported_cell(want, have)
    return :(throw(ArgumentError("Cannot create " * string($(esc(want))) * " cell for " * string($(esc(have))) * " structure")))
end

"""
    bulk(
        name::String,
        crystalstructure::Union{String,Symbol}="auto";
        a::Union{Real,Nothing}=nothing,
        b::Union{Real,Nothing}=nothing,
        c::Union{Real,Nothing}=nothing,
        alpha::Union{Real,Nothing}=nothing,
        covera::Union{Real,Nothing}=nothing,
        u::Union{Real,Nothing}=nothing,
        orthorhombic::Bool=false,
        cubic::Bool=false,
    ) -> Cell

Create bulk crystal structures.

Crystal structure and lattice constant(s) will be guessed if not provided.

# Arguments
- `name`: Chemical symbol or symbols as in `"MgO"` or `"NaCl"`. For single 
  elements, reference data is automatically used when crystalstructure="auto".
- `crystalstructure`: Structure type as a String or Symbol. Use "auto" to 
  automatically detect from reference database for single elements.
  Supported: sc, fcc, bcc, hcp, diamond, tetragonal, bct, rhombohedral, 
  orthorhombic, monoclinic, cubic, zincblende, rocksalt, cesiumchloride, 
  fluorite, wurtzite.
- `a`: Lattice constant. If not specified, uses value from reference database.
- `b`: Lattice constant.
- `c`: Lattice constant.
- `alpha`: Angle in degrees for rhombohedral lattice.
- `covera`: c/a ratio used for hcp/tetragonal. Default is ideal ratio: sqrt(8/3).
- `u`: Internal coordinate for Wurtzite structure.
- `orthorhombic`: Construct orthorhombic unit cell instead of primitive cell.
- `cubic`: Construct cubic unit cell if possible.

# Examples
```julia
# Auto-detect from reference database
bulk("Cu")                    # FCC copper, a=3.61Å
bulk("Fe")                    # BCC iron, a=2.87Å  
bulk("Mg")                    # HCP magnesium, a=3.21Å, c/a=1.624

# Explicit structure specification
bulk("Cu", "fcc", a=3.65)
bulk("Si", "diamond", a=5.43)
bulk("NaCl", "rocksalt", a=5.64)

# Cell type options
bulk("Cu", "fcc", cubic=true)  # 4-atom cubic cell
```
"""
function bulk(
    name::String,
    crystalstructure::Union{String,Symbol}="auto";
    a::Union{Real,Nothing}=nothing,
    b::Union{Real,Nothing}=nothing,
    c::Union{Real,Nothing}=nothing,
    alpha::Union{Real,Nothing}=nothing,
    covera::Union{Real,Nothing}=nothing,
    u::Union{Real,Nothing}=nothing,
    orthorhombic::Bool=false,
    cubic::Bool=false,
)
    # Parse chemical formula
    symbols = string_to_symbols(name)
    natoms = length(symbols)
    
    # Handle auto-detection of crystalstructure for single elements
    ref_data = nothing
    if crystalstructure == "auto" || crystalstructure == :auto
        if natoms == 1
            # Look up in reference database
            ref_data = lookup_reference_state(symbols[1])
            if isnothing(ref_data)
                throw(ArgumentError("No reference data for element '$(symbols[1])'. " *
                                   "Please specify crystalstructure explicitly."))
            end
            
            crystalstructure = ref_data[:symmetry]
            
            # Check for special cases that need basis
            if haskey(ref_data, :basis) && isnothing(ref_data[:basis])
                throw(ArgumentError("Element '$(symbols[1])' has crystal structure '$crystalstructure' " *
                                   "which requires a basis that is not yet implemented. " *
                                   "Please specify a different structure manually."))
            end
            
            # Check for molecular structures
            if crystalstructure in ["diatom", "atom"]
                throw(ArgumentError("Element '$(symbols[1])' is a $(crystalstructure) in its reference state. " *
                                   "Please specify a crystal structure for bulk creation."))
            end
        else
            throw(ArgumentError("Cannot auto-detect crystal structure for compound '$name' with $(natoms) atoms. " *
                               "Please specify crystalstructure parameter."))
        end
    end
    
    crystalstructure = Symbol(lowercase(String(crystalstructure)))
    
    # Check for conflicting specifications
    if !isnothing(covera) && !isnothing(c)
        throw(ArgumentError("Don't specify both c and c/a!"))
    end
    
    # Apply reference data for lattice constants if not specified
    if !isnothing(ref_data)
        # Get lattice constant a from reference if not specified
        if isnothing(a) && haskey(ref_data, :a)
            a = Float64(ref_data[:a])
        end
        
        # Get c/a ratio from reference if applicable
        if crystalstructure in [:hcp, :wurtzite] && isnothing(covera) && haskey(ref_data, :covera)
            covera = Float64(ref_data[:covera])
        end
        
        # Get b/a ratio for orthorhombic
        if crystalstructure == :orthorhombic && isnothing(b) && haskey(ref_data, :bovera) && !isnothing(a)
            b = Float64(ref_data[:bovera]) * a
        end
        
        # Get c/a ratio for orthorhombic/tetragonal
        if crystalstructure in [:orthorhombic, :tetragonal, :bct] && isnothing(c) && isnothing(covera) && haskey(ref_data, :covera) && !isnothing(a)
            c = Float64(ref_data[:covera]) * a
        end
        
        # Get alpha for rhombohedral
        if crystalstructure == :rhombohedral && isnothing(alpha) && haskey(ref_data, :alpha)
            alpha = Float64(ref_data[:alpha])
        end
    end
    
    # Verify lattice constant
    if isnothing(a)
        throw(ArgumentError("Lattice constant 'a' must be specified or available in reference data"))
    end
    a = Float64(a)
    
    # Number of atoms for each structure type (primitive cell)
    structure_natoms = Dict(
        :sc => 1, :fcc => 1, :bcc => 1,
        :tetragonal => 1, :bct => 1,
        :hcp => 2,
        :rhombohedral => 1,
        :orthorhombic => 1,
        :monoclinic => 1,
        :mcl => 1,
        :cubic => 1,
        :diamond => 2,
        :zincblende => 2, :rocksalt => 2, :cesiumchloride => 2,
        :fluorite => 3, :wurtzite => 4
    )
    
    if !haskey(structure_natoms, crystalstructure)
        throw(ArgumentError("Unknown crystal structure: $crystalstructure"))
    end
    
    expected_natoms = structure_natoms[crystalstructure]
    
    # Check if we need to replicate atoms for single-element structures
    if natoms == 1 && expected_natoms > 1
        # Replicate the single symbol for structures that need multiple atoms
        if crystalstructure in [:hcp, :diamond]
            symbols = repeat(symbols, expected_natoms)
            natoms = length(symbols)
        end
    end
    
    # Wurtzite needs to expand 2 atoms to 4 (2 of each type)
    if crystalstructure == :wurtzite && natoms == 2
        symbols = [symbols[1], symbols[1], symbols[2], symbols[2]]
        natoms = 4
    end
    
    if natoms != expected_natoms
        throw(ArgumentError("Structure $crystalstructure requires $expected_natoms atom(s), " *
                           "but '$name' has $natoms atom(s)"))
    end
    
    # Handle c/a ratio for hcp and wurtzite
    if crystalstructure in [:hcp, :wurtzite]
        if cubic
            @build_unsupported_cell "cubic" crystalstructure
        end
        
        if !isnothing(c)
            covera = c / a
        elseif isnothing(covera)
            covera = sqrt(8 / 3)  # Ideal c/a ratio
        end
    end
    
    # Build the structure
    if orthorhombic
        return _orthorhombic_bulk(symbols, crystalstructure, a, covera, u)
    elseif cubic && crystalstructure in [:bcc, :cesiumchloride]
        return _orthorhombic_bulk(symbols, crystalstructure, a, covera, u)
    elseif cubic && crystalstructure != :sc
        return _cubic_bulk(symbols, crystalstructure, a)
    else
        return _primitive_bulk(symbols, crystalstructure, a, b, c, alpha, covera, u)
    end
end

"""
    _primitive_bulk(symbols, structure, a, b, c, alpha, covera, u)

Build primitive cell structures.
"""
function _primitive_bulk(symbols, structure, a, b, c, alpha, covera, u)
    T = Float64
    
    if structure == :sc
        # Simple cubic
        lattice = Lattice(a, a, a)
        positions = zeros(T, 3, 1)
        return Cell(lattice, symbols, positions)
        
    elseif structure == :fcc
        # FCC - primitive cell
        b = a / 2
        lattice = Lattice([
            0.0  b    b;
            b    0.0  b;
            b    b    0.0
        ])
        positions = zeros(T, 3, 1)
        return Cell(lattice, symbols, positions)
        
    elseif structure == :bcc
        # BCC - primitive cell
        b = a / 2
        lattice = Lattice([
            -b    b    b;
            b    -b    b;
            b    b    -b
        ])
        positions = zeros(T, 3, 1)
        return Cell(lattice, symbols, positions)
        
    elseif structure == :hcp
        # HCP - primitive cell with 2 atoms
        sqrt3 = sqrt(T(3))
        c = covera * a
        # Lattice vectors as columns
        lattice = Lattice([
            a       -a/2         0.0;
            0.0     a*sqrt3/2    0.0;
            0.0     0.0          c
        ])
        # Two atoms at (0, 0, 0) and (1/3, 2/3, 1/2) in fractional coordinates
        positions = zeros(T, 3, 2)
        positions[:, 1] = [0.0, 0.0, 0.0]
        positions[:, 2] = lattice.matrix * [1/3, 2/3, 1/2]
        return Cell(lattice, symbols, positions)
        
    elseif structure == :diamond
        # Diamond = two interpenetrating FCC lattices
        return _build_diamond(symbols, a)
        
    elseif structure == :zincblende
        # Zincblende (ZnS structure)
        return _build_zincblende(symbols, a)
        
    elseif structure == :rocksalt
        # Rocksalt (NaCl structure)
        return _build_rocksalt(symbols, a)
        
    elseif structure == :cesiumchloride
        # Cesium chloride structure
        return _build_cesiumchloride(symbols, a)
        
    elseif structure == :fluorite
        # Fluorite (CaF2 structure)
        return _build_fluorite(symbols, a)
        
    elseif structure == :wurtzite
        # Wurtzite structure
        return _build_wurtzite(symbols, a, covera, u)
        
    elseif structure == :tetragonal
        # Tetragonal cell
        c_val = isnothing(c) ? a : Float64(c)
        lattice = Lattice(a, a, c_val)
        positions = zeros(T, 3, 1)
        return Cell(lattice, symbols, positions)
        
    elseif structure == :bct
        # Body-centered tetragonal
        c_val = isnothing(c) ? a : Float64(c)
        lattice = Lattice(a, a, c_val)
        positions = zeros(T, 3, 1)
        return Cell(lattice, symbols, positions)
        
    elseif structure == :orthorhombic
        # Orthorhombic cell
        b_val = isnothing(b) ? a : Float64(b)
        c_val = isnothing(c) ? a : Float64(c)
        lattice = Lattice(a, b_val, c_val)
        positions = zeros(T, 3, 1)
        return Cell(lattice, symbols, positions)
        
    elseif structure == :rhombohedral
        # Rhombohedral cell
        return _build_rhombohedral(symbols, a, alpha)
        
    elseif structure == :cubic
        # Generic cubic cell
        lattice = Lattice(a, a, a)
        positions = zeros(T, 3, 1)
        return Cell(lattice, symbols, positions)
        
    elseif structure == :monoclinic
        throw(ArgumentError("Monoclinic structure not yet implemented"))
        
    else
        throw(ArgumentError("Structure $structure not yet implemented"))
    end
end

"""
    _build_diamond(symbols, a)

Build diamond structure (primitive cell).
"""
function _build_diamond(symbols, a)
    T = Float64
    # Diamond structure requires 2 atoms (same element)
    # Note: For single elements, bulk() automatically replicates the symbol
    if length(symbols) != 2
        throw(ArgumentError("Diamond structure requires 2 atoms. " *
                           "For single elements, use bulk(\"Si\", \"diamond\") and the symbol will be auto-replicated."))
    end
    
    b = a / 2
    lattice = Lattice([
        0.0  b    b;
        b    0.0  b;
        b    b    0.0
    ])
    
    # Two atoms at (0, 0, 0) and (1/4, 1/4, 1/4) in conventional cell
    positions = zeros(T, 3, 2)
    positions[:, 1] = [0.0, 0.0, 0.0]
    positions[:, 2] = [a/4, a/4, a/4]
    
    return Cell(lattice, symbols, positions)
end

"""
    _build_zincblende(symbols, a)

Build zincblende structure (primitive cell).
"""
function _build_zincblende(symbols, a)
    T = Float64
    if length(symbols) != 2
        throw(ArgumentError("Zincblende structure requires 2 different atoms (e.g., 'ZnS')"))
    end
    
    b = a / 2
    lattice = Lattice([
        0.0  b    b;
        b    0.0  b;
        b    b    0.0
    ])
    
    # Two atoms at (0, 0, 0) and (1/4, 1/4, 1/4)
    positions = zeros(T, 3, 2)
    positions[:, 1] = [0.0, 0.0, 0.0]
    positions[:, 2] = [a/4, a/4, a/4]
    
    return Cell(lattice, symbols, positions)
end

"""
    _build_rocksalt(symbols, a)

Build rocksalt (NaCl) structure.
"""
function _build_rocksalt(symbols, a)
    T = Float64
    if length(symbols) != 2
        throw(ArgumentError("Rocksalt structure requires 2 different atoms (e.g., 'NaCl')"))
    end
    
    b = a / 2
    lattice = Lattice([
        0.0  b    b;
        b    0.0  b;
        b    b    0.0
    ])
    
    # Two atoms: one at origin, other at face center (1/2, 0, 0)
    positions = zeros(T, 3, 2)
    positions[:, 1] = [0.0, 0.0, 0.0]
    positions[:, 2] = [a/2, 0.0, 0.0]
    
    return Cell(lattice, symbols, positions)
end

"""
    _build_cesiumchloride(symbols, a)

Build cesium chloride structure.
"""
function _build_cesiumchloride(symbols, a)
    T = Float64
    if length(symbols) != 2
        throw(ArgumentError("Cesium chloride structure requires 2 different atoms (e.g., 'CsCl')"))
    end
    
    lattice = Lattice(a, a, a)
    
    # Two atoms: one at origin, other at body center (1/2, 1/2, 1/2)
    positions = zeros(T, 3, 2)
    positions[:, 1] = [0.0, 0.0, 0.0]
    positions[:, 2] = [a/2, a/2, a/2]
    
    return Cell(lattice, symbols, positions)
end

"""
    _build_fluorite(symbols, a)

Build fluorite (CaF2) structure.
"""
function _build_fluorite(symbols, a)
    T = Float64
    if length(symbols) != 3
        throw(ArgumentError("Fluorite structure requires 3 atoms (e.g., 'CaFF')"))
    end
    
    b = a / 2
    lattice = Lattice([
        0.0  b    b;
        b    0.0  b;
        b    b    0.0
    ])
    
    # Ca at (0,0,0), F at (1/4,1/4,1/4) and (3/4,3/4,3/4)
    positions = zeros(T, 3, 3)
    positions[:, 1] = [0.0, 0.0, 0.0]
    positions[:, 2] = [a/4, a/4, a/4]
    positions[:, 3] = [3*a/4, 3*a/4, 3*a/4]
    
    return Cell(lattice, symbols, positions)
end

"""
    _build_wurtzite(symbols, a, covera, u)

Build wurtzite structure.
"""
function _build_wurtzite(symbols, a, covera, u)
    T = Float64
    
    # Default u parameter
    if isnothing(u)
        u = T(0.375)
    end
    
    sqrt3 = sqrt(T(3))
    c = covera * a
    
    # Lattice vectors as columns
    lattice = Lattice([
        a       -a/2         0.0;
        0.0     a*sqrt3/2    0.0;
        0.0     0.0          c
    ])
    
    # 4 atoms: 2 of each type
    positions = zeros(T, 3, 4)
    # Type A
    positions[:, 1] = lattice.matrix * [0.0, 0.0, 0.0]
    positions[:, 2] = lattice.matrix * [1/3, 2/3, 0.5]
    # Type B
    positions[:, 3] = lattice.matrix * [0.0, 0.0, u]
    positions[:, 4] = lattice.matrix * [1/3, 2/3, 0.5 + u]
    
    # Symbols should already be expanded to 4 elements [A, A, B, B]
    if length(symbols) != 4
        throw(ArgumentError("Wurtzite structure requires 4 atoms (2 of each type, e.g., 'ZnO')"))
    end
    
    return Cell(lattice, symbols, positions)
end

"""
    _build_rhombohedral(symbols, a, alpha)

Build rhombohedral structure (primitive cell).
"""
function _build_rhombohedral(symbols, a, alpha)
    T = Float64
    
    if isnothing(alpha)
        # Default to 60 degrees if not specified
        alpha = T(60.0)
    end
    
    # Build rhombohedral lattice
    # All sides equal a, angles alpha, beta, gamma
    lattice = Lattice(a, a, a, alpha, alpha, alpha)
    
    # Single atom at origin for simple rhombohedral
    positions = zeros(T, 3, 1)
    
    return Cell(lattice, symbols, positions)
end

"""
    _orthorhombic_bulk(symbols, structure, a, covera, u)

Build orthorhombic conventional cell structures.
"""
function _orthorhombic_bulk(symbols, structure, a, covera, u)
    T = Float64
    
    if structure == :fcc
        # Orthorhombic FCC
        b = a / sqrt(T(2))
        lattice = Lattice(b, b, a)
        positions = zeros(T, 3, 2)
        positions[:, 1] = [0.0, 0.0, 0.0]
        positions[:, 2] = [b/2, b/2, a/2]
        # Both atoms are the same for single element FCC
        if length(symbols) == 1
            symbols = [symbols[1], symbols[1]]
        end
        return Cell(lattice, symbols, positions)
        
    elseif structure == :bcc
        # Orthorhombic BCC (conventional cubic)
        lattice = Lattice(a, a, a)
        positions = zeros(T, 3, 2)
        positions[:, 1] = [0.0, 0.0, 0.0]
        positions[:, 2] = [a/2, a/2, a/2]
        if length(symbols) == 1
            symbols = [symbols[1], symbols[1]]
        end
        return Cell(lattice, symbols, positions)
        
    elseif structure == :hcp
        # Orthorhombic HCP
        sqrt3 = sqrt(T(3))
        c = covera * a
        lattice = Lattice(a, a*sqrt3, c)
        positions = zeros(T, 3, 4)
        positions[:, 1] = [0.0, 0.0, 0.0]
        positions[:, 2] = [a/2, a*sqrt3/2, 0.0]
        positions[:, 3] = [a/2, a*sqrt3/6, c/2]
        positions[:, 4] = [0.0, 2*a*sqrt3/3, c/2]
        # Replicate symbols
        expanded_symbols = repeat(symbols, 2)
        return Cell(lattice, expanded_symbols, positions)
        
    elseif structure == :diamond
        # Orthorhombic diamond (via zincblende)
        if length(symbols) == 1
            return _orthorhombic_bulk([symbols[1], symbols[1]], :zincblende, a, covera, u)
        else
            return _orthorhombic_bulk(symbols, :zincblende, a, covera, u)
        end
        
    elseif structure == :zincblende
        # Orthorhombic zincblende
        b = a / sqrt(T(2))
        lattice = Lattice(b, b, a)
        positions = zeros(T, 3, 4)
        positions[:, 1] = [0.0, 0.0, 0.0]
        positions[:, 2] = [b/2, 0.0, a/4]
        positions[:, 3] = [b/2, b/2, a/2]
        positions[:, 4] = [0.0, b/2, 3*a/4]
        # Replicate symbols: 2 of each type
        expanded_symbols = [symbols[1], symbols[1], symbols[2], symbols[2]]
        return Cell(lattice, expanded_symbols, positions)
        
    elseif structure == :rocksalt
        # Orthorhombic rocksalt
        b = a / sqrt(T(2))
        lattice = Lattice(b, b, a)
        positions = zeros(T, 3, 4)
        positions[:, 1] = [0.0, 0.0, 0.0]
        positions[:, 2] = [b/2, b/2, 0.0]
        positions[:, 3] = [b/2, b/2, a/2]
        positions[:, 4] = [0.0, 0.0, a/2]
        # Replicate symbols: 2 of each type
        expanded_symbols = [symbols[1], symbols[1], symbols[2], symbols[2]]
        return Cell(lattice, expanded_symbols, positions)
        
    elseif structure == :cesiumchloride
        # Orthorhombic cesium chloride (same as primitive)
        lattice = Lattice(a, a, a)
        positions = zeros(T, 3, 2)
        positions[:, 1] = [0.0, 0.0, 0.0]
        positions[:, 2] = [a/2, a/2, a/2]
        return Cell(lattice, symbols, positions)
        
    elseif structure == :wurtzite
        # Orthorhombic wurtzite
        if isnothing(u)
            u = T(0.375)
        end
        sqrt3 = sqrt(T(3))
        c = covera * a
        lattice = Lattice(a, a*sqrt3, c)
        positions = zeros(T, 3, 8)
        # First layer A atoms
        positions[:, 1] = [0.0, 0.0, 0.0]
        positions[:, 2] = [a/2, a*sqrt3/2, 0.0]
        # B atoms
        positions[:, 3] = [0.0, a/sqrt3, c*(0.5-u)]
        positions[:, 4] = [0.0, 0.0, c*(1-u)]
        # Second layer A atoms
        positions[:, 5] = [a/2, a*sqrt3/6, c/2]
        positions[:, 6] = [0.0, 2*a/sqrt3, c/2]
        # B atoms
        positions[:, 7] = [a/2, a*sqrt3/2, c*(0.5-u)]
        positions[:, 8] = [a/2, a*sqrt3/6, c*(1-u)]
        
        # Replicate symbols: 4 of each type
        expanded_symbols = repeat(symbols, 4)
        return Cell(lattice, expanded_symbols, positions)
        
    else
        throw(ArgumentError("Cannot create orthorhombic cell for $structure structure"))
    end
end

"""
    _cubic_bulk(symbols, structure, a)

Build cubic conventional cell structures.
"""
function _cubic_bulk(symbols, structure, a)
    T = Float64
    
    if structure == :fcc
        # Cubic FCC conventional cell
        lattice = Lattice(a, a, a)
        positions = zeros(T, 3, 4)
        positions[:, 1] = [0.0, 0.0, 0.0]
        positions[:, 2] = [0.0, a/2, a/2]
        positions[:, 3] = [a/2, 0.0, a/2]
        positions[:, 4] = [a/2, a/2, 0.0]
        # Replicate symbols
        expanded_symbols = repeat(symbols, 4)
        return Cell(lattice, expanded_symbols, positions)
        
    elseif structure == :diamond
        # Cubic diamond = zincblende with same element
        if length(symbols) == 1
            return _cubic_bulk([symbols[1], symbols[1]], :zincblende, a)
        else
            return _cubic_bulk(symbols, :zincblende, a)
        end
        
    elseif structure == :zincblende
        # Cubic zincblende
        lattice = Lattice(a, a, a)
        positions = zeros(T, 3, 8)
        # FCC lattice points for type 1
        positions[:, 1] = [0.0, 0.0, 0.0]
        positions[:, 2] = [0.0, a/2, a/2]
        positions[:, 3] = [a/2, 0.0, a/2]
        positions[:, 4] = [a/2, a/2, 0.0]
        # FCC lattice points for type 2 (shifted by a/4, a/4, a/4)
        positions[:, 5] = [a/4, a/4, a/4]
        positions[:, 6] = [a/4, 3*a/4, 3*a/4]
        positions[:, 7] = [3*a/4, a/4, 3*a/4]
        positions[:, 8] = [3*a/4, 3*a/4, a/4]
        # Replicate symbols: 4 of each type
        expanded_symbols = repeat(symbols, 4)
        return Cell(lattice, expanded_symbols, positions)
        
    elseif structure == :rocksalt
        # Cubic rocksalt
        lattice = Lattice(a, a, a)
        positions = zeros(T, 3, 8)
        # FCC lattice points for type 1
        positions[:, 1] = [0.0, 0.0, 0.0]
        positions[:, 2] = [0.0, a/2, a/2]
        positions[:, 3] = [a/2, 0.0, a/2]
        positions[:, 4] = [a/2, a/2, 0.0]
        # FCC lattice points for type 2 (shifted by a/2, 0, 0)
        positions[:, 5] = [a/2, 0.0, 0.0]
        positions[:, 6] = [a/2, a/2, a/2]
        positions[:, 7] = [0.0, 0.0, a/2]
        positions[:, 8] = [0.0, a/2, 0.0]
        # Replicate symbols: 4 of each type
        expanded_symbols = repeat(symbols, 4)
        return Cell(lattice, expanded_symbols, positions)
        
    elseif structure == :fluorite
        # Cubic fluorite
        lattice = Lattice(a, a, a)
        positions = zeros(T, 3, 12)
        # Ca positions (FCC)
        positions[:, 1] = [0.0, 0.0, 0.0]
        positions[:, 2] = [0.0, a/2, a/2]
        positions[:, 3] = [a/2, 0.0, a/2]
        positions[:, 4] = [a/2, a/2, 0.0]
        # F positions (at 1/4 and 3/4)
        positions[:, 5] = [a/4, a/4, a/4]
        positions[:, 6] = [a/4, 3*a/4, 3*a/4]
        positions[:, 7] = [3*a/4, a/4, 3*a/4]
        positions[:, 8] = [3*a/4, 3*a/4, a/4]
        positions[:, 9] = [3*a/4, 3*a/4, 3*a/4]
        positions[:, 10] = [3*a/4, a/4, a/4]
        positions[:, 11] = [a/4, 3*a/4, a/4]
        positions[:, 12] = [a/4, a/4, 3*a/4]
        
        expanded_symbols = [repeat([symbols[1]], 4); repeat([symbols[2]], 4); repeat([symbols[3]], 4)]
        return Cell(lattice, expanded_symbols, positions)
        
    else
        throw(ArgumentError("Cannot create cubic cell for $structure structure"))
    end
end
