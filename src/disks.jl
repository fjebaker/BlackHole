using Parameters: @with_kw

abstract type AccretionDisk end

"""
Geometric disk with no opacity.
Does not consider relativistic, doppler or any other effects.
"""
@with_kw struct GeometricDisk <: AccretionDisk

    "Inclination of the disk into the plane."
    α::Float64 = π/50

    "Inclination of the disk around the axis of view."
    β::Float64 = 0.0

    "Inner disk radius (arbitrary units)."
    r_inner::Float64 = 12.0

    "Outer disk radius (arbitrary units)."
    r_outer::Float64 = 44.0

end
