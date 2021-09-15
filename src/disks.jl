using Parameters: @with_kw

abstract type AccretionDisk end
abstract type OpticallyThinDisk <: AccretionDisk end

"""
Geometric disk with no opacity.
Does not consider relativistic, doppler or any other effects.
"""
@with_kw struct GeometricDisk <: AccretionDisk

    "Inclination of the disk into the plane."
    α::Float64 = π / 50

    "Inclination of the disk around the axis of view."
    β::Float64 = 0.0

    "Inner disk radius (arbitrary units)."
    r_inner::Float64 = 12.0

    "Outer disk radius (arbitrary units)."
    r_outer::Float64 = 44.0

end


"""
Disk with Gaussian profile for opacity.
"""
@with_kw struct GaussianThinDisk <: OpticallyThinDisk

    "Inclination of the disk into the plane."
    α::Float64 = π / 50

    "Inclination of the disk around the axis of view."
    β::Float64 = 0.0

    "Inner disk radius (arbitrary units)."
    r_inner::Float64 = 12.0

    "Outer disk radius (arbitrary units)."
    r_outer::Float64 = 44.0

    "Singularity the disk is bound to."
    s::Singularity

    "The standard deviation on the Gaussian disk profile"
    σ::Float64 = 4

    "Mean of the Gaussian disk profile."
    μ::Float64 = (r_outer - r_inner) / 2 + r_inner
end


function profile(disk::GaussianThinDisk, x)
    begin
        if !(disk.r_inner < x < disk.r_outer)
            return 0
        else
            y = (1 / disk.σ * sqrt(2 * π)) * exp(-0.5 * ((x - disk.μ) / disk.σ)^2)
            return max(y, 0)
        end
    end
end

@doc raw"""

    orthonormal(s::Singularity, r, dt, dr, dθ, dϕ)

Transform given coordinates into the orthonormal reference frame of a given singularity.

Returns `(dt, dr, dθ, dϕ)` in orthonormal frame.
"""
function orthonormal(s::Singularity, r, dt, dr, dθ, dϕ)
    
    (
        sqrt(r / 2 * s.m + r) * dt,
        2 * s.m / (sqrt(r * (2 * s.m + r))) * dt + sqrt((2 * s.m + r) / r) * dr,
        dθ * r,
        dϕ * r,
    )
end


function intersection(
    disk::GaussianThinDisk,
    coord1::NamedTuple,
    coord2::NamedTuple,
    r,
)::Float64
    if disk.r_inner < r < disk.r_outer
        if coord1.z != coord2.z # prevent divide by zero
            if 0 <= coord1.z / (coord1.z - coord2.z) < 1
                # calculate relevant components of the four velocity
                dt, dx, dy, dz = (coord2.t - coord1.t, coord2.x - coord1.x, coord2.y - coord1.y, coord2.z - coord1.z)
                
                # transform to spherical
                dr, dθ, dϕ = dxyz2spher(coord2.x, coord2.y, coord2.z, dx, dy, dz)
                
                # shift to orthonormal reference frame
                dt, dr, dθ, dϕ = orthonormal(disk.s, r, dt, dr, dθ, dϕ)

                # calculate width of θ
                θ_width = atan(dθ, sqrt(dr^2 + dϕ^2))

                # thickness with ~magic~ constant
                return  30 * profile(disk, r) / sin(abs(θ_width))
            end
        end
    end
    return 0
end

function intersection(
    disk::GeometricDisk,
    coord1::NamedTuple,
    coord2::NamedTuple,
    r,
)::Float64
    if disk.r_inner < r < disk.r_outer
        if coord1.z != coord2.z # prevent divide by zero
            if 0 <= coord1.z / (coord1.z - coord2.z) < 1
                return 255
            end
        end
    end
    return 0
end

@views function intersection(disk::AccretionDisk, geodesic::AbstractArray, β)
    (t, x, y, z) = geodesic[:, begin]

    # transform into plane of disk
    x, y, z = inclination_transform(x, y, z, disk.α, disk.β + β)

    sum(
        i -> begin
            (t2, x2, y2, z2) = geodesic[:, i]

            r = sqrt(x2^2 + y2^2 + z2^2)
            x2, y2, z2 = inclination_transform(x2, y2, z2, disk.α, disk.β + β)

            # use named tuples
            value = intersection(
                disk,
                (t = t, x = x, y = y, z = z),
                (t = t2, x = x2, y = y2, z = z2),
                r,
            )

            # update state
            t, x, y, z = t2, x2, y2, z2

            value
        end,
        2:size(geodesic)[2],
    )
end

# module exports
export GeometricDisk
