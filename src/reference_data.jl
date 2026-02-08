#=
Elemental Reference States Database

Data from ASE (Atomic Simulation Environment), originally from:
Ashcroft and Mermin, "Solid State Physics" (standard solid state physics reference)

Each entry is a Dict containing:
- symmetry: Crystal structure type
- a: Lattice constant in Angstroms
- covera: c/a ratio for hcp/tetragonal/bct (optional)
- bovera: b/a ratio for orthorhombic (optional)  
- alpha: Angle in degrees for rhombohedral (optional)
- magmom_per_atom: Magnetic moment per atom (optional)

Special structures:
- "diatom": Diatomic molecule (d = bond length)
- "atom": Single atom (no lattice constant)
=#

const REFERENCE_STATES = Dict{Symbol, Union{Nothing, Dict{Symbol, Any}}}(
    # X - dummy
    :X => nothing,
    
    # H - Hydrogen (diatomic molecule)
    :H => Dict(:symmetry => "diatom", :d => 0.74),
    
    # He - Helium (atom)
    :He => Dict(:symmetry => "atom"),
    
    # Li - Lithium (BCC)
    :Li => Dict(:symmetry => "bcc", :a => 3.49),
    
    # Be - Beryllium (HCP)
    :Be => Dict(:symmetry => "hcp", :a => 2.29, :covera => 1.567),
    
    # B - Boron (tetragonal, complex basis)
    :B => Dict(:symmetry => "tetragonal", :a => 8.73, :covera => 0.576, :basis => nothing),
    
    # C - Carbon (diamond)
    :C => Dict(:symmetry => "diamond", :a => 3.57),
    
    # N - Nitrogen (diatomic)
    :N => Dict(:symmetry => "diatom", :d => 1.10),
    
    # O - Oxygen (diatomic)
    :O => Dict(:symmetry => "diatom", :d => 1.21),
    
    # F - Fluorine (diatomic)
    :F => Dict(:symmetry => "diatom", :d => 1.42),
    
    # Ne - Neon (FCC)
    :Ne => Dict(:symmetry => "fcc", :a => 4.43),
    
    # Na - Sodium (BCC)
    :Na => Dict(:symmetry => "bcc", :a => 4.23),
    
    # Mg - Magnesium (HCP)
    :Mg => Dict(:symmetry => "hcp", :a => 3.21, :covera => 1.624),
    
    # Al - Aluminum (FCC)
    :Al => Dict(:symmetry => "fcc", :a => 4.05),
    
    # Si - Silicon (diamond)
    :Si => Dict(:symmetry => "diamond", :a => 5.43),
    
    # P - Phosphorus (cubic, complex)
    :P => Dict(:symmetry => "cubic", :a => 7.17, :basis => nothing),
    
    # S - Sulfur (orthorhombic)
    :S => Dict(:symmetry => "orthorhombic", :a => 10.47, :bovera => 1.229, :covera => 2.339, :basis => nothing),
    
    # Cl - Chlorine (orthorhombic)
    :Cl => Dict(:symmetry => "orthorhombic", :a => 6.24, :bovera => 0.718, :covera => 1.324, :basis => nothing),
    
    # Ar - Argon (FCC)
    :Ar => Dict(:symmetry => "fcc", :a => 5.26),
    
    # K - Potassium (BCC)
    :K => Dict(:symmetry => "bcc", :a => 5.23),
    
    # Ca - Calcium (FCC)
    :Ca => Dict(:symmetry => "fcc", :a => 5.58),
    
    # Sc - Scandium (HCP)
    :Sc => Dict(:symmetry => "hcp", :a => 3.31, :covera => 1.594),
    
    # Ti - Titanium (HCP)
    :Ti => Dict(:symmetry => "hcp", :a => 2.95, :covera => 1.588),
    
    # V - Vanadium (BCC)
    :V => Dict(:symmetry => "bcc", :a => 3.02),
    
    # Cr - Chromium (BCC)
    :Cr => Dict(:symmetry => "bcc", :a => 2.88),
    
    # Mn - Manganese (cubic, complex)
    :Mn => Dict(:symmetry => "cubic", :a => 8.89, :basis => nothing),
    
    # Fe - Iron (BCC, magnetic)
    :Fe => Dict(:symmetry => "bcc", :a => 2.87, :magmom_per_atom => 2.3),
    
    # Co - Cobalt (HCP, magnetic)
    :Co => Dict(:symmetry => "hcp", :a => 2.51, :covera => 1.622, :magmom_per_atom => 1.2),
    
    # Ni - Nickel (FCC, magnetic)
    :Ni => Dict(:symmetry => "fcc", :a => 3.52, :magmom_per_atom => 0.6),
    
    # Cu - Copper (FCC)
    :Cu => Dict(:symmetry => "fcc", :a => 3.61),
    
    # Zn - Zinc (HCP)
    :Zn => Dict(:symmetry => "hcp", :a => 2.66, :covera => 1.856),
    
    # Ga - Gallium (orthorhombic)
    :Ga => Dict(:symmetry => "orthorhombic", :a => 4.51, :bovera => 1.001, :covera => 1.695, :basis => nothing),
    
    # Ge - Germanium (diamond)
    :Ge => Dict(:symmetry => "diamond", :a => 5.66),
    
    # As - Arsenic (rhombohedral)
    :As => Dict(:symmetry => "rhombohedral", :a => 4.13, :alpha => 54.10),
    
    # Se - Selenium (HCP, complex)
    :Se => Dict(:symmetry => "hcp", :a => 4.36, :covera => 1.136, :basis => nothing),
    
    # Br - Bromine (orthorhombic)
    :Br => Dict(:symmetry => "orthorhombic", :a => 6.67, :bovera => 0.672, :covera => 1.307, :basis => nothing),
    
    # Kr - Krypton (FCC)
    :Kr => Dict(:symmetry => "fcc", :a => 5.72),
    
    # Rb - Rubidium (BCC)
    :Rb => Dict(:symmetry => "bcc", :a => 5.59),
    
    # Sr - Strontium (FCC)
    :Sr => Dict(:symmetry => "fcc", :a => 6.08),
    
    # Y - Yttrium (HCP)
    :Y => Dict(:symmetry => "hcp", :a => 3.65, :covera => 1.571),
    
    # Zr - Zirconium (HCP)
    :Zr => Dict(:symmetry => "hcp", :a => 3.23, :covera => 1.593),
    
    # Nb - Niobium (BCC)
    :Nb => Dict(:symmetry => "bcc", :a => 3.30),
    
    # Mo - Molybdenum (BCC)
    :Mo => Dict(:symmetry => "bcc", :a => 3.15),
    
    # Tc - Technetium (HCP)
    :Tc => Dict(:symmetry => "hcp", :a => 2.74, :covera => 1.604),
    
    # Ru - Ruthenium (HCP)
    :Ru => Dict(:symmetry => "hcp", :a => 2.70, :covera => 1.584),
    
    # Rh - Rhodium (FCC)
    :Rh => Dict(:symmetry => "fcc", :a => 3.80),
    
    # Pd - Palladium (FCC)
    :Pd => Dict(:symmetry => "fcc", :a => 3.89),
    
    # Ag - Silver (FCC)
    :Ag => Dict(:symmetry => "fcc", :a => 4.09),
    
    # Cd - Cadmium (HCP)
    :Cd => Dict(:symmetry => "hcp", :a => 2.98, :covera => 1.886),
    
    # In - Indium (BCT - body-centered tetragonal)
    :In => Dict(:symmetry => "bct", :a => 4.59 / sqrt(2), :covera => 1.076 * sqrt(2)),
    
    # Sn - Tin (BCT, beta-Sn/white tin)
    :Sn => Dict(:symmetry => "bct", :a => 5.82, :covera => 0.546),
    
    # Sb - Antimony (rhombohedral)
    :Sb => Dict(:symmetry => "rhombohedral", :a => 4.51, :alpha => 57.60),
    
    # Te - Tellurium (HCP, complex)
    :Te => Dict(:symmetry => "hcp", :a => 4.45, :covera => 1.330, :basis => nothing),
    
    # I - Iodine (orthorhombic)
    :I => Dict(:symmetry => "orthorhombic", :a => 7.27, :bovera => 0.659, :covera => 1.347, :basis => nothing),
    
    # Xe - Xenon (FCC)
    :Xe => Dict(:symmetry => "fcc", :a => 6.20),
    
    # Cs - Cesium (BCC)
    :Cs => Dict(:symmetry => "bcc", :a => 6.05),
    
    # Ba - Barium (BCC)
    :Ba => Dict(:symmetry => "bcc", :a => 5.02),
    
    # La - Lanthanum (HCP)
    :La => Dict(:symmetry => "hcp", :a => 3.75, :covera => 1.619),
    
    # Ce - Cerium (FCC)
    :Ce => Dict(:symmetry => "fcc", :a => 5.16),
    
    # Pr - Praseodymium (HCP)
    :Pr => Dict(:symmetry => "hcp", :a => 3.67, :covera => 1.614),
    
    # Nd - Neodymium (HCP)
    :Nd => Dict(:symmetry => "hcp", :a => 3.66, :covera => 1.614),
    
    # Pm - Promethium (no data)
    :Pm => nothing,
    
    # Sm - Samarium (rhombohedral, complex)
    :Sm => Dict(:symmetry => "rhombohedral", :a => 9.00, :alpha => 23.13),
    
    # Eu - Europium (BCC)
    :Eu => Dict(:symmetry => "bcc", :a => 4.61),
    
    # Gd - Gadolinium (HCP)
    :Gd => Dict(:symmetry => "hcp", :a => 3.64, :covera => 1.588),
    
    # Tb - Terbium (HCP)
    :Tb => Dict(:symmetry => "hcp", :a => 3.60, :covera => 1.581),
    
    # Dy - Dysprosium (HCP)
    :Dy => Dict(:symmetry => "hcp", :a => 3.59, :covera => 1.573),
    
    # Ho - Holmium (HCP)
    :Ho => Dict(:symmetry => "hcp", :a => 3.58, :covera => 1.570),
    
    # Er - Erbium (HCP)
    :Er => Dict(:symmetry => "hcp", :a => 3.56, :covera => 1.570),
    
    # Tm - Thulium (HCP)
    :Tm => Dict(:symmetry => "hcp", :a => 3.54, :covera => 1.570),
    
    # Yb - Ytterbium (FCC)
    :Yb => Dict(:symmetry => "fcc", :a => 5.49),
    
    # Lu - Lutetium (HCP)
    :Lu => Dict(:symmetry => "hcp", :a => 3.51, :covera => 1.585),
    
    # Hf - Hafnium (HCP)
    :Hf => Dict(:symmetry => "hcp", :a => 3.20, :covera => 1.582),
    
    # Ta - Tantalum (BCC)
    :Ta => Dict(:symmetry => "bcc", :a => 3.31),
    
    # W - Tungsten (BCC)
    :W => Dict(:symmetry => "bcc", :a => 3.16),
    
    # Re - Rhenium (HCP)
    :Re => Dict(:symmetry => "hcp", :a => 2.76, :covera => 1.615),
    
    # Os - Osmium (HCP)
    :Os => Dict(:symmetry => "hcp", :a => 2.74, :covera => 1.579),
    
    # Ir - Iridium (FCC)
    :Ir => Dict(:symmetry => "fcc", :a => 3.84),
    
    # Pt - Platinum (FCC)
    :Pt => Dict(:symmetry => "fcc", :a => 3.92),
    
    # Au - Gold (FCC)
    :Au => Dict(:symmetry => "fcc", :a => 4.08),
    
    # Hg - Mercury (rhombohedral)
    :Hg => Dict(:symmetry => "rhombohedral", :a => 2.99, :alpha => 70.45),
    
    # Tl - Thallium (HCP)
    :Tl => Dict(:symmetry => "hcp", :a => 3.46, :covera => 1.599),
    
    # Pb - Lead (FCC)
    :Pb => Dict(:symmetry => "fcc", :a => 4.95),
    
    # Bi - Bismuth (rhombohedral)
    :Bi => Dict(:symmetry => "rhombohedral", :a => 4.75, :alpha => 57.14),
    
    # Po - Polonium (SC)
    :Po => Dict(:symmetry => "sc", :a => 3.35),
    
    # At - Astatine (no data)
    :At => nothing,
    
    # Rn - Radon (no data)
    :Rn => nothing,
    
    # Fr - Francium (no data)
    :Fr => nothing,
    
    # Ra - Radium (no data)
    :Ra => nothing,
    
    # Ac - Actinium (FCC)
    :Ac => Dict(:symmetry => "fcc", :a => 5.31),
    
    # Th - Thorium (FCC)
    :Th => Dict(:symmetry => "fcc", :a => 5.08),
    
    # Pa - Protactinium (tetragonal)
    :Pa => Dict(:symmetry => "tetragonal", :a => 3.92, :covera => 0.825),
    
    # U - Uranium (orthorhombic)
    :U => Dict(:symmetry => "orthorhombic", :a => 2.85, :bovera => 1.736, :covera => 2.056),
    
    # Np - Neptunium (orthorhombic)
    :Np => Dict(:symmetry => "orthorhombic", :a => 4.72, :bovera => 1.035, :covera => 1.411),
    
    # Pu - Plutonium (monoclinic)
    :Pu => Dict(:symmetry => "monoclinic"),
    
    # Am - Americium (no data)
    :Am => nothing,
    
    # Cm - Curium (no data)
    :Cm => nothing,
    
    # Bk - Berkelium (no data)
    :Bk => nothing,
    
    # Cf - Californium (no data)
    :Cf => nothing,
    
    # Es - Einsteinium (no data)
    :Es => nothing,
    
    # Fm - Fermium (no data)
    :Fm => nothing,
    
    # Md - Mendelevium (no data)
    :Md => nothing,
    
    # No - Nobelium (no data)
    :No => nothing,
    
    # Lr - Lawrencium (no data)
    :Lr => nothing,
    
    # Rf - Rutherfordium (no data)
    :Rf => nothing,
    
    # Db - Dubnium (no data)
    :Db => nothing,
    
    # Sg - Seaborgium (no data)
    :Sg => nothing,
    
    # Bh - Bohrium (no data)
    :Bh => nothing,
    
    # Hs - Hassium (no data)
    :Hs => nothing,
    
    # Mt - Meitnerium (no data)
    :Mt => nothing,
    
    # Ds - Darmstadtium (no data)
    :Ds => nothing,
    
    # Rg - Roentgenium (no data)
    :Rg => nothing,
    
    # Cn - Copernicium (no data)
    :Cn => nothing,
    
    # Nh - Nihonium (no data)
    :Nh => nothing,
    
    # Fl - Flerovium (no data)
    :Fl => nothing,
    
    # Mc - Moscovium (no data)
    :Mc => nothing,
    
    # Lv - Livermorium (no data)
    :Lv => nothing,
    
    # Ts - Tennessine (no data)
    :Ts => nothing,
    
    # Og - Oganesson (no data)
    :Og => nothing,
)

"""
    lookup_reference_state(symbol::Symbol)

Look up the reference crystal structure data for an element.
Returns nothing if element not found or has no data.
"""
function lookup_reference_state(symbol::Symbol)
    return get(REFERENCE_STATES, symbol, nothing)
end

"""
    get_supported_structures()

Return a list of structure types supported by the bulk() function.
"""
function get_supported_structures()
    return ["sc", "fcc", "bcc", "hcp", "diamond", 
            "tetragonal", "bct", "rhombohedral", 
            "orthorhombic", "monoclinic", "cubic",
            "zincblende", "rocksalt", "cesiumchloride", 
            "fluorite", "wurtzite"]
end
