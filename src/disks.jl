using Parameters: @with_kw

abstract type AccretionDisk end

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


function intersection(disk::GeometricDisk, geodesic::AbstractArray, β)
    (t, x, y, z) = geodesic[:, begin]

    # transform into plane of disk
    x, y, z = inclination_transform(x, y, z, disk.α, disk.β + β)

    sum(
        i -> begin
            value = 0
            (t2, x2, y2, z2) = geodesic[:, i]

            r = sqrt(x2^2 + y2^2 + z2^2)
            x2, y2, z2 = inclination_transform(x2, y2, z2, disk.α, disk.β + β)

            if disk.r_inner < r < disk.r_outer
                if z != z2 # prevent divide by zero
                    if 0 <= z / (z - z2) < 1
                        value = 255
                    end
                end
            end

            # update state
            t, x, y, z = t2, x2, y2, z2

            value
        end,
        2:size(geodesic)[2],
    )
end

# module exports
export GeometricDisk
