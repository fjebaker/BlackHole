"""
BlackHole rendering library.
"""

module BlackHole

include("metrics.jl")

include("disks.jl")
include("utils.jl")


include("geodesics.jl")

include("gpu.jl")
include("rendering.jl")

end
