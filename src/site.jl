#=
Representation of single atomic Sites

The `Site` type is deprecated.
=#

using Printf
using StaticArrays
export Site, distance_between, distance_squared_between

"""
Site{T}

Represent a site with absolute coordinates.

Coordinates can be accessed with :x, :y, :z.

# Example

```julia-repl
julia> site = Site(Float64[1,2,3], 1, :H)
Site{Vector{Float64}}:
    H    1    1.00000    2.00000    3.00000

julia> site.x
1.0

julia> cell=Cell(Lattice(10, 10, 10), [:H], [[1.,2., 3.]])
Cell with 1 ions
Lattice: 
  10.000      0.000      0.000
   0.000     10.000      0.000
   0.000      0.000     10.000
Sites: 
   H     1.00000     2.00000     3.00000
```

Updating the position of the `Site` object will update that of the `Cell` object.
```julia-repl
julia> site = cell[1]
Site{SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}:
    H    1    1.00000    2.00000    3.00000

julia> site.position[1] 
1.0

julia> site.position[1] = 2.;

3×1 Matrix{Float64}:
 2.0
 2.0
 3.0
```

"""
struct Site{T}
    position::T
    index::Int
    symbol::Symbol
end

"Initialise a site from only position and symbol"
function Site(pos, symbol)
    Site(pos, 0, symbol)
end

function Base.show(io::IO, s::Site{T}) where {T}
    print(
        io,
        @sprintf "Site{%s}: %5s %4d at %10.5f %10.5f %10.5f" T s.symbol s.index s.position[1] s.position[2] s.position[3]
    )
end

function Base.show(io::IO, ::MIME"text/plain", s::Site{T}) where {T}
    print(io, @sprintf "Site{%s}:\n" T)
    print(
        io,
        @sprintf "%5s %4d %10.5f %10.5f %10.5f" s.symbol s.index s.position[1] s.position[2] s.position[3]
    )
end

position(s::Site) = s.position
index(s::Site) = s.index
symbol(s::Site) = s.symbol
function Base.getproperty(s::Site, sym::Symbol)
    if sym == :x
        return s.position[1]
    elseif sym == :y
        return s.position[2]
    elseif sym == :z
        return s.position[3]
    else
        return getfield(s, sym)
    end
end

Base.propertynames(s::Site) = (Base.fieldnames(Site)..., :x, :y, :z)

distance_squared_between(s1::Site, s2::Site) = sum((s1.position .- s2.position) .^ 2)
distance_squared_between(s1::AbstractVector, s2::AbstractVector) = sum((s1 - s2) .^ 2)
distance_squared_between(s1::AbstractVector, s2::AbstractVector, shift) =
    sum((s1 - s2 .- shift) .^ 2)
distance_squared_between(s1::Site, s2::Site, shift) =
    sum((s1.position .- s2.position .- shift) .^ 2)

distance_between(s1, s2, shift) = sqrt(distance_squared_between(s1, s2, shift))
distance_between(s1, s2) = sqrt(distance_squared_between(s1, s2))


"Unit vector from site 1 to shifted site 2"
function unit_vector_between(s1::Site, s2::Site, shift)
    v = s2.position .- s1.position .+ shift
    v ./ norm(v)
end
