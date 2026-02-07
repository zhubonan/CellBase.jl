#= 
functions for handling periodic conditions
=#
using LinearAlgebra
using StaticArrays

"""
    shift_vectors(lattice::AbstractMatrix, shift1, shift2, shift3)

Compute the shift vectors needed to include all parts of the lattice within a
cut off radius.
"""
function shift_vectors(
    lattice::AbstractMatrix,
    s1max,
    s2max,
    s3max,
    s1min=-s1max,
    s2min=-s2max,
    s3min=-s3max,
)
    a1, a2, a3, =
        SVector{3}(lattice[:, 1]), SVector{3}(lattice[:, 2]), SVector{3}(lattice[:, 3])
    nshifts = (s1max - s1min + 1) * (s2max - s2min + 1) * (s3max - s3min + 1)
    shift_vectors = SVector{3,Float64}[]

    for s3 = s3min:s3max
        for s2 = s2min:s2max
            for s1 = s1min:s1max
                push!(shift_vectors, s1 .* a1 .+ s2 .* a2 .+ s3 .* a3)
            end  # s1
        end  # s2
    end  # s3
    shift_vectors
end

"Compute shift vectors given lattice matrix and cut off radius"
function shift_vectors(lattice::AbstractMatrix, rc::Real; safe=false)
    shifts = max_shifts(lattice, rc; safe=safe)
    shift_vectors(lattice, shifts...)
end

"""
    shift_indices(lattice::AbstractArray, rc::Real)

Compute a vector containing shift indices
"""
function shift_indices(lattice::AbstractArray, rc::Real)
    shift1max, shift2max, shift3max = max_shifts(lattice, rc)
    shift_indices(shift1max, shift2max, shift3max)
end

function shift_indices(shift1::Int, shift2::Int, shift3::Int)
    nshifts = (shift1 * 2 + 1) * (shift2 * 2 + 1) * (shift3 * 2 + 1)
    idx = Array{Int}(undef, (3, nshifts))
    idx = SVector{3,Float64}[]
    for s1 = -shift1:shift1, s2 = -shift2:shift2, s3 = -shift3:shift3
        push!(idx, SA[s1, s2, s3])
    end
    idx
end


"""
Compute the require periodicities required to fill a cut off sphere of a certain radius
"""
function max_shifts(lattice::AbstractArray, rc::Real; safe=false)

    # Not suer why this, probably a overkill
    if safe
        diag_vec = sum(lattice, dims=2)
        diag = sqrt(dot(diag_vec, diag_vec))
        rc += diag / 2
    end

    invl = inv(lattice)
    b1, b2, b3 = invl[:, 1], invl[:, 2], invl[:, 3]

    shift1max = ceil(Int, rc * norm(b1))
    shift2max = ceil(Int, rc * norm(b2))
    shift3max = ceil(Int, rc * norm(b3))
    return shift1max, shift2max, shift3max
end
