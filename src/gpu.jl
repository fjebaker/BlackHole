using CUDA

@doc raw"""
    function kernel_index_calculator(
        oimg::AbstractArray{Int},
        β_store::AbstractArray{<:Number},
        width::Int,
        height::Int,
        fov::Int,
        num_geodesics::Int
    )

CUDA kernel for calculating which index of the geodesic array to use in each pixel. The index
calculated for a given `y, x` position is stored in `oimg[y, x]`. `β_store[y, x]` also stores
the corresponding pixel's angle around the `x` axis, used for calculating inclination into the
accretion disk in the second kernel (`kernel_geodesic_render`).

`fov` is a "field of view" index, used to adjust which pixel maps to the end of the geodesic
array (`num_geodesics`). As a rough guideline, `fov=1000` maps the pixel at distance `r=1000`
to the last index of the geodesic array.
"""
function kernel_index_calculator(
    oimg::AbstractArray{Int},
    β_store::AbstractArray{<:Number},
    width::Int,
    height::Int,
    fov::Int,
    num_geodesics::Int,
)
    # which pixel does this thread work on
    i = threadIdx().x + blockDim().x * (blockIdx().x - 1)

    # guard
    if i > width * height
        return nothing
    end

    # calculate position in image
    frame = (i - 1)
    rowᵢ = (frame) % height + 1
    colᵢ = (frame) ÷ height + 1

    # x, y from center of image
    x = colᵢ - 1 - (width ÷ 2)
    y = (rowᵢ - 1 - (height ÷ 2)) * -1

    # angle around the unit circle, positive x axis is 0
    β = atan(y, x) # atan2

    # distance from center of image
    r = sqrt(x^2 + y^2)

    # use this distance to calculate which index of the geodesic to assign to the pixel
    geoindex = 1 + trunc(Int, min(r / fov * num_geodesics, num_geodesics - 1))

    # store
    oimg[i] = geoindex
    β_store[i] = β

    return nothing
end

@doc raw"""
    kernel_geodesic_render(
        indeximage::AbstractArray{T},
        β_store::AbstractArray{<:Number},
        geo_matrix,::AbstractArray{<:Number},
        width::Int, height::Int,
        step_num::Int,
        disk::GeometricDisk,
        max_index::Int
    ) where T <: Int

CUDA kernel for computing the integrated paths of geodesics through an accretion disk.
`indeximage` is the result from calling `kernel_index_calculator`, and is a matrix containing
the geodesic index for each pixel.

`max_index` is the maximal index of the geodesics matrix, i.e. `size(geo_matrix, 3)`.

To be launch with `threads = 512 blocks = length(indeximage)`, until better topology is
implemented.

# Extended help

This kernel works by using one compute block per pixel, and calculating the difference
between two points in spacetime along the geodesic per thread. The thread then calculates
whether the line length between these points intersects the accretion disk, and sets
a value in a shared memory buffer: in essense, each thread calculates
```jl
    geodesics[:, thread + 1, geoindex] - geodesics[:, thread, geoindex]
```
transformed by the angles of the accretion disk. This allows the intersection to be calculated
from the ``z`` component:
``
0 \leq \frac{z_i}{z_i - z_{i+1}} < 1
``

Internally, the offsets are mapped with the column major formula
```
n x m x p:
    arr[i, j, k] == arr[(n*m)*(k-1) + n*(j-1) + i]
```
to avoid `getindex` function calls.

At the end of the kernel, a simple reduce algorithm is implemented to sum over the shared
array. The first thread of each block then stores the value back into `indeximage` as
the output of the kernel.

Better topology for this kernel would have the y dimension continue to calculate the same
pixel, incase there are more geodesic integration steps (`step_num`) than threads can
be spawned. Due to the size of this kernel, only about `512` threads can be spawned on my
hardware before the device runs out of resources.
"""
function kernel_geodesic_render(
    indeximage::AbstractArray{T},
    β_store::AbstractArray{<:Number},
    geo_matrix::AbstractArray{<:Number},
    width::Int,
    height::Int,
    step_num::Int,
    disk::GeometricDisk,
    max_index::Int,
) where {T<:Int}

    # work out what this thread needs to do
    elements = blockDim().x
    thread = threadIdx().x

    px_index = blockIdx().x

    shared = @cuStaticSharedMem(T1, (2048,))
    # initialise data
    @inbounds shared[thread] = 0

    # guard
    if px_index > width * height || (thread + 1) > step_num
        return nothing
    end

    # retrieve geodesic index from array
    @inbounds geoindex = indeximage[px_index]

    # guard
    if geoindex ≤ 0 || geoindex > max_index
        return nothing
    end

    # retrieve beta calculation from array
    @inbounds β = β_store[px_index]

    α = disk.α

    # offset into the array from the geoindex
    offset = 4 * step_num * (geoindex - 1)

    ind1::Int = (4 * (thread - 1)) + 1 + offset # thread
    ind2::Int = (4 * thread) + 1 + offset # thread + 1

    # retrieve distant points first, so we can check if we're actually in the disk
    x2 = geo_matrix[ind2+1]
    y2 = geo_matrix[ind2+2]
    z2 = geo_matrix[ind2+3]

    radius = sqrt(x2^2 + y2^2 + z2^2)

    if disk.r_inner < radius < disk.r_outer
        # potential intersection

        # precalc
        cosβ = cos(β)
        cosα = cos(α)
        sinβ = sin(β)
        sinα = sin(α)

        x1 = geo_matrix[ind1+1]
        y1 = geo_matrix[ind1+2]
        z1 = geo_matrix[ind1+3]

        # transform into the plane of the disk
        # only calculate components we need
        z1 = z1 * cosβ + y1 * sinβ
        z1 = z1 * cosα + x1 * sinα

        z2 = z2 * cosβ + y2 * sinβ
        z2 = z2 * cosα + x2 * sinα

        # check if intersection
        if z1 != z2
            if 0 ≤ z1 / (z1 - z2) < 1
                @inbounds shared[thread] = 255
            end
        end
    end

    # reduce algorithm implementation
    d = 1
    while d < elements
        # sync threads in reduce
        sync_threads()

        ind = (2 * d * (thread - 1) + 1)

        if ind ≤ elements && ind + d ≤ elements
            # sum over shared values
            @inbounds shared[ind] = shared[ind] + shared[ind+d]
        end

        d *= 2
    end

    # if leading thread, save value to output
    if thread == 1
        @inbounds indeximage[px_index] = shared[1]
    end

    return nothing
end



@doc raw"""
    renderdisk(
        ::Val{:gpu},
        disk::AccretionDisk,
        geodesics::AbstractArray{<:Number},
        width::Int,
        height::Int,
        fov::Int,
    )

Dispatch method for rendering the disk using the GPU.

# Extended help
    
The kernel launch configuration is hard-coded
into this function, which launches two kernels: `kernel_index_calculator` and `kernel_geodesic_render`.
See the documentation of these individual functions for more details.
"""
function renderdisk(
    ::Val{:gpu},
    disk::AccretionDisk,
    geodesics::AbstractArray{<:Number},
    width::Int,
    height::Int,
    fov::Int,
)
    β_store = CUDA.zeros(Float64, (height, width))
    outimage = similar(β_store, Int64)

    CUDA.@sync @cuda threads = 1024 blocks = cld(length(outimage), 1024) kernel_index_calculator(
        outimage,
        β_store,
        width,
        height,
        fov,
        size(geodesics, 3),
    )

    β_store .+= disk.β

    CUDA.@sync @cuda threads = 512 blocks = length(outimage) kernel_geodesic_render(
        outimage,
        β_store,
        geodesics,
        width,
        height,
        size(geodesics, 2),
        disk,
        size(geodesics, 3),
    )

    cat(outimage, outimage, outimage; dims = 3)
end
