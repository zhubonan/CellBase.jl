#=
Utility modules - contains common routines
=#
using LinearAlgebra

export rotation_matrix, rotation_matrix_align, axis_from_string

const dgrd = π / 180.0
const rddg = 1 / dgrd

function angle(x::AbstractVector, y::AbstractVector, radian=false)
    lx = norm(x)
    ly = norm(y)
    dotp = dot(x, y)
    cos_ang = dotp / lx / ly
    ang = acos(cos_ang)
    if !radian
        ang *= 180.0 / pi
    end
    return ang
end


"""
Convert cell vectors to cell parameters

Returns an static array of the cell parameters
"""
function vec2cellpar(va::AbstractVector, vb::AbstractVector, vc::AbstractVector)
    a = norm(va)
    b = norm(vb)
    c = norm(vc)

    α = angle(vb, vc)
    β = angle(va, vc)
    γ = angle(va, vb)
    return [a, b, c, α, β, γ]
end

function vec2cellpar(cell::AbstractMatrix)
    va = cell[:, 1]
    vb = cell[:, 2]
    vc = cell[:, 3]
    return vec2cellpar(va, vb, vc)
end


function _cellpar_trig(α, β, γ)
    # For the orthorhombic cell case
    eps = 1.4e-14   # Error tolorance
    # Alpha
    if abs(abs(α) - 0.5π) < eps
        cos_alpha = 0.0
    else
        cos_alpha = cos(α)
    end
    # β
    if abs(abs(β) - 0.5π) < eps
        cos_beta = 0.0
    else
        cos_beta = cos(β)
    end

    if abs(γ - 0.5π) < eps
        cos_gamma = 0.0
        sin_gamma = 1.0
    elseif abs(γ + 0.5π) < eps
        cos_gamma = 0.0
        sin_gamma = -1.0
    else
        cos_gamma = cos(γ)
        sin_gamma = sin(γ)
    end
    return cos_alpha, cos_beta, cos_gamma, sin_gamma
end

radian(a, b, c) = dgrd * a, dgrd * b, dgrd * c

"Check if cell parameters are valid"
function isvalidcellpar(a, b, c, α, β, γ; degree=true)
    if degree
        α, β, γ = radian(α, β, γ)
    end
    cos_alpha, cos_beta, cos_gamma, sin_gamma = _cellpar_trig(α, β, γ)
    cx = cos_beta
    cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    tmp = 1 - cx * cx - cy * cy
    tmp >= 0.0
end


"convert cell parameters to column vectors"
function cellpar2mat(a, b, c, α, β, γ; degree=true)

    if degree
        α, β, γ = radian(α, β, γ)
    end
    a_direction = [1.0, 0.0, 0.0]
    ab_normal = [0.0, 0.0, 1.0]

    cos_alpha, cos_beta, cos_gamma, sin_gamma = _cellpar_trig(α, β, γ)

    # Buildce cell vectors
    va = [a, 0.0, 0.0]
    vb = [cos_gamma * b, sin_gamma * b, 0.0]
    cx = cos_beta
    cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    tmp = 1.0 - cx * cx - cy * cy
    @assert tmp > 0 "Cell parameters are not valid"
    cz = sqrt(tmp)
    vc = [cx * c, cy * c, cz * c]

    return hcat(va, vb, vc)
end

"Compute volume from cell parameters"
function volume(a, b, c, α, β, γ; degree=true)
    if degree
        α, β, γ = radian(α, β, γ)
    end
    a * b * c * sqrt(1 + 2 * cos(α) * cos(β) * cos(γ) - cos(α)^2 - cos(β)^2 - cos(γ)^2)
end


"""
    rotation_matrix(angle::Real, axis::AbstractVector) -> Matrix{Float64}

Compute the 3×3 rotation matrix for rotating by `angle` (in degrees) around `axis`.
Uses the Rodrigues rotation formula.

The rotation matrix R satisfies: v_rotated = R * v

# Arguments
- `angle`: Rotation angle in degrees
- `axis`: 3D rotation axis vector (will be normalized)

# Returns
3×3 rotation matrix

# Examples
```julia
R = rotation_matrix(90, [0, 0, 1])  # 90 degree rotation around z-axis
R = rotation_matrix(45, [1, 1, 1])  # 45 degree rotation around [111] direction
```
"""
function rotation_matrix(angle::Real, axis::AbstractVector)
    # Normalize axis
    ax = normalize(axis)
    
    # Convert angle to radians
    θ = deg2rad(angle)
    c = cos(θ)
    s = sin(θ)
    
    # Rodrigues rotation formula: R = I*cos(θ) + (1-cos(θ))*aa^T + sin(θ)*[a×]
    # where [a×] is the cross-product matrix
    
    # Cross-product matrix K
    K = [0    -ax[3]  ax[2];
         ax[3]  0    -ax[1];
        -ax[2] ax[1]  0]
    
    # Rotation matrix: R = I + sin(θ)*K + (1-cos(θ))*K^2
    I_mat = Matrix{Float64}(I, 3, 3)
    R = I_mat + s * K + (1 - c) * (K * K)
    
    return R
end


"""
    rotation_matrix_align(v1::AbstractVector, v2::AbstractVector) -> Matrix{Float64}

Compute the rotation matrix that aligns vector `v1` with vector `v2`.

# Arguments
- `v1`: Source vector (will be rotated to align with v2)
- `v2`: Target vector direction

# Returns
3×3 rotation matrix
"""
function rotation_matrix_align(v1::AbstractVector, v2::AbstractVector)
    # Normalize vectors
    u1 = normalize(v1)
    u2 = normalize(v2)
    
    # Check if vectors are parallel
    dot_prod = dot(u1, u2)
    
    # Parallel (same direction)
    if dot_prod > 1 - 1e-10
        return Matrix{Float64}(I, 3, 3)
    end
    
    # Anti-parallel (opposite direction)
    if dot_prod < -1 + 1e-10
        # Find perpendicular axis
        if abs(u1[1]) < abs(u1[2])
            axis = normalize(cross(u1, [1.0, 0.0, 0.0]))
        else
            axis = normalize(cross(u1, [0.0, 1.0, 0.0]))
        end
        return rotation_matrix(180, axis)
    end
    
    # General case: angle = acos(dot(u1, u2))
    # Rotation axis is perpendicular to both
    axis = cross(u1, u2)
    angle_rad = acos(clamp(dot_prod, -1.0, 1.0))
    angle_deg = rad2deg(angle_rad)
    
    return rotation_matrix(angle_deg, axis)
end


"""
    axis_from_string(s::Union{String,Symbol}) -> Vector{Float64}

Convert a string/symbol axis specification to a vector.

Supported formats:
- "x", :x → [1, 0, 0]
- "y", :y → [0, 1, 0]
- "z", :z → [0, 0, 1]
- "-x", :-x → [-1, 0, 0]
- "-y", :-y → [0, -1, 0]
- "-z", :-z → [0, 0, -1]

# Examples
```julia
axis_from_string("x")    # [1, 0, 0]
axis_from_string(:z)     # [0, 0, 1]
axis_from_string("-y")   # [0, -1, 0]
```
"""
function axis_from_string(s::Union{String,Symbol})
    str = lowercase(string(s))
    
    # Handle negative prefix
    negative = startswith(str, "-")
    axis_char = negative ? str[2:end] : str
    
    vec = if axis_char == "x"
        [1.0, 0.0, 0.0]
    elseif axis_char == "y"
        [0.0, 1.0, 0.0]
    elseif axis_char == "z"
        [0.0, 0.0, 1.0]
    else
        throw(ArgumentError("Invalid axis specification: '$s'. Use 'x', 'y', 'z', '-x', '-y', or '-z'."))
    end
    
    return negative ? -vec : vec
end
