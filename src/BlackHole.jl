module BlackHole

abstract type Singularity end
abstract type AccretionDisk end

"""
Includes are order sensitive
"""

include("geodesics.jl")

include("eddingtonfinkelstein.jl")

include("accretiondisk.jl")
include("disks/geometric.jl")
# include("disks/opticallythin.jl")

include("render.jl")

end