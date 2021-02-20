module BlackHole

abstract type Singularity end
abstract type AccretionDisk end

"""
Includes are order sensitive
"""

include("geodesics.jl")

include("eddingtonfinkelstein.jl")
include("accretiondisk.jl")
include("render.jl")

#export EddingtonFinkelstein, calcgeodesics, plot!, GeometricDisk, intersection, renderdisk
end