module BlackHole

abstract type Singularity end
abstract type AccretionDisk end

"""

s = EddingtonFinkelstein(2.0)

geodesics = calcgeodesics(s, num=1000, Δϕ=0.005)

disk = GeometricDisk(
    α=π/3, 
    β=0.0, 
    rinner=12,
    router=44
)

image = renderdisk(
    disk,
    geodesics, 
    height=720, 
    width=1080,
    fov_index=200    
)

"""

include("geodesics.jl")

include("eddingtonfinkelstein.jl")
include("accretiondisk.jl")
include("render.jl")

export EddingtonFinkelstein, calcgeodesics, plot!, GeometricDisk, intersection, renderdisk
end