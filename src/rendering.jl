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


function render_column!(
    column_view::AbstractArray{<:Number},
    x::Int,
    geodesics::AbstractArray{<:Number},
    disk::AccretionDisk,
    fov::Int,
    height::Int,
)
    column_view .= 0 # zero out the column

    mid_height = height ÷ 2
    n = size(geodesics)[3]

    for y = 1:(height-1)
        β = atan(y - mid_height, x) # atan2
        r = sqrt(x^2 + (y - mid_height)^2) # distance from middle of the image
        index = convert(Int, round(1 + r / fov * n))

        if index < n
            value = intersection(disk, geodesics[:, :, index], β)
            @. column_view[height-y, :] = truncator(value)
        end
    end
end


@doc raw"""
    renderdisk(
        ::Val{:cpu},
        geodesics::AbstractArray{<:Number},
        width::Int,
        height::Int,
        fov::Int,
        disk::AccretionDisk
    )

Dispatch method for rendering the disk using the CPU.
"""
function renderdisk(
    ::Val{:cpu},
    geodesics::AbstractArray{<:Number},
    width::Int,
    height::Int,
    fov::Int,
    disk::AccretionDisk,
)
    image_out = zeros(UInt8, (height, width, 3))

    mid_width = width ÷ 2

    for x = (-mid_width):(mid_width-1)
        render_column!(
            @view(image_out[:, x+mid_width+1, :]),
            x,
            geodesics,
            disk,
            fov,
            height,
        )
    end

    image_out
end

@doc raw"""
    renderdisk(
        geodesics::AbstractArray{<:Number},
        width::Int,
        height::Int,
        fov::Int,
        disk::AccretionDisk
    )

Render an accretion disk `disk` given an array of pre-calculated geodesics `geodesics` for the disk.

Returns an array of dimensions `width` by `height`, representing the rendered image. 

This method checks whether CUDA is installed, and dispatches a GPU kernel accordingly. Otherwise, the
rendering functions are executed in parallel on the CPU. 
"""
function renderdisk(
    geodesics::AbstractArray{<:Number},
    width::Int,
    height::Int,
    fov::Int,
    disk::AccretionDisk,
)
    if CUDA.has_cuda_gpu()
        return renderdisk(:gpu, geodesics, width, height, fov, disk)
    else
        return renderdisk(:cpu, geodesics, width, height, fov, disk)
    end
end

renderdisk(s::Symbol, args...; kwargs...) = renderdisk(Val(s), args...; kwargs...)
