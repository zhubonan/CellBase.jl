#=
For handling compositions
=#
import Base
export Composition, formula, atomic_weight

"""
Type representing a composition

The Composition type supports index access:

```julia-repl
julia> comp = Composition(:CO2)
Composition(:CO2)
julia> comp[:C] 
1
```

Internally, the speices and quantities are stored as two vectors, yet the interface
mimics that of a `Dict{Sybmol,Float64}`.

"""
struct Composition
    species::Vector{Symbol}
    counts::Vector{Float64}
    function Composition(species, counts)
        length(species) != length(counts) && error("Size mismatch")
        sidx = sortperm(species)
        new(species[sidx], counts[sidx])
    end
end


"""
    Composition(pairs::Pair{Symbol, T}...) where T

Return `Composition` from pairs.
"""
function Composition(pairs::Pair{Symbol,T}...) where {T}
    x = Array{Symbol}(undef, length(pairs))
    y = Array{Float64}(undef, length(pairs))
    i = 1
    for (a, b) in pairs
        x[i] = a
        y[i] = b
        i += 1
    end
    Composition(x, y)
end

Base.hash(t::Composition) = Base.hash((t.species, t.counts))

function formula(t::Composition)
    args = []
    for (s, c) in zip(t.species, t.counts)
        push!(args, s)
        if round(c) == c
            # Omitted 1, e.g. C1O2 -> CO2
            Int(c) == 1 && continue
            push!(args, Symbol(Int(c)))
        else
            push!(args, Symbol(c))
        end
    end
    Symbol(args...)
end

Base.string(t::Composition) = string(formula(t))

formula(t::Cell) = formula(Composition(t))

function Base.show(io::IO, ::MIME"text/plain", o::Composition)
    print(io, "Composition(:$(formula(o)))")
end

Base.show(io::IO, o::Composition) = Base.show(io, MIME("text/plain"), o)

"""
    Composition(comp::Dict)

Return `Composition` from dictionary
"""
function Composition(comp::Dict)
    Composition(collect(comp)...)
end


function Base.getindex(t::Composition, key::Symbol)
    idx = findfirst(x -> x == key, t.species)
    isnothing(idx) && return 0.0
    t.counts[idx]
end

function Base.setindex!(t::Composition, value, key::Symbol)
    idx = findfirst(x -> x == key, t.species)
    if isnothing(idx)
        push!(t.species, key)
        push!(t.counts, value)
        #Sort the speices and indices
        perm = sortperm(t.species)
        news = t.species[perm]
        newc = t.counts[perm]
        t.species .= news
        t.counts .= newc
    else
        t.counts[idx] = value
    end
end


Base.haskey(t::Composition, key) = key in t.species
Base.keys(t::Composition) = t.species
Base.pairs(t::Composition) = [a => b for (a, b) in zip(t.species, t.counts)]
Base.iterate(t::Composition) = Base.iterate(Base.pairs(t))
Base.iterate(t::Composition, i) = Base.iterate(Base.pairs(t), i)
Base.length(t::Composition) = Base.length(t.species)
function Base.:(==)(t::Composition, t2::Composition)
    (t.species == t2.species) && (t.counts == t2.counts)
end
function Base.isequal(t::Composition, t2::Composition)
    (t.species == t2.species) && (t.counts == t2.counts)
end
Base.reduce(t::Composition) = t / gcd(Int.(t.counts)...)

"""
    reduced_and_factor(t::Composition)

Return reduced formula and reduction factor
"""
function reduced_and_factor(t::Composition)
    n = nform(t)
    t / n, n
end

nform(t::Composition) = gcd(Int.(t.counts))
nform(t::Cell) = nform(Composition(t))

"""
    Composition(cell::Cell)

Return the composition of a Cell
"""
function Composition(cell::Cell)
    sp_array = species(cell)
    unique_sp = sort(unique(sp_array))
    num_atoms = Array{Int}(undef, size(unique_sp))
    for i in 1:length(unique_sp)
        num_atoms[i] = count(x -> x == unique_sp[i], sp_array)
    end
    Composition(unique_sp, num_atoms)
end

"""
    Composition(string::AbstractString)

Construct composition from a string.
"""
function Composition(string::AbstractString)
    Composition(parse_formula_with_bracket(string))
end

Composition(symbol::Symbol) = Composition(string(symbol))

"""
    parse_formula(formula)

Parse a formula into a dictionary.
"""
function parse_formula(formula, factor=1.0)
    reg = r"([A-Z][a-z]*)\s*([\.\d]*)"
    output = Dict{Symbol,Float64}()
    for match in eachmatch(reg, formula)
        num = 1.0
        match[2] != "" && (num = parse(Float64, match[2]))
        x = get(output, Symbol(match[1]), 0.0)
        x += num * factor
        output[Symbol(match[1])] = x
    end
    return output
end

"""
    expand_bracket(formula)

Return a formula with bracket expanded.
"""
function expand_bracket(formula)
    reg_bracket = r"\(([^\(\)]+)\)\s*([\.eE\d]*)"
    m = match(reg_bracket, formula)
    if !isnothing(m)
        # Expanding is needed
        if m[2] == ""
            factor = 1.0
        else
            factor = parse(Float64, m[2])
        end
        # Parse the contant inside the bracket
        inner = parse_formula(m[1], factor)
        expanded_sym = join(["$(sym)$(num)" for (sym, num) in inner])
        formula = replace(formula, m.match => expanded_sym)
        return expand_bracket(formula)
    end
    return formula
end

"""
    parse_formula_with_bracket(formula)

Parse formula and expand the bracket if necessary.
"""
function parse_formula_with_bracket(formula)
    expanded = expand_bracket(formula)
    return parse_formula(expanded)
end


function Base.:+(a::Composition, b::Composition)
    all_keys = unique(vcat(keys(a), keys(b)))
    nums = zeros(Float64, length(all_keys))
    for (i, key) in enumerate(all_keys)
        nums[i] = a[key] + b[key]
    end
    Composition(all_keys, nums)
end

function Base.:-(a::Composition, b::Composition)
    nums = copy(a.counts)
    out = deepcopy(a)
    for (i, key) in enumerate(b.species)
        if key in a.species
            out[key] = a[key] - b[key]
        else
            throw(ErrorException("$(a) does not contain  $(key)"))
        end
    end
    out
end

Base.:*(a::Composition, b::Real) = Composition(a.species, a.counts .* b)
Base.:/(a::Composition, b::Real) = Composition(a.species, a.counts ./ b)

atomic_weight(comp::Composition) = sum(elements[key].atomic_mass * value for (key, value) in pairs(comp))

function Base.contains(a::Composition, b::Composition)
    akeys = keys(a)
    for (s, v) in pairs(b)
        if !(s in akeys)
            return false
        end
        if a[s] < v
            return false
        end
    end
    return true
end