"""
BlackHole rendering library.
"""

module BlackHole

include("disks.jl")
include("utils.jl")

include("metrics.jl")

include("geodesics.jl")

include("gpu.jl")
include("rendering.jl")

end
