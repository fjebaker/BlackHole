function render_column!(
    column_view,
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
            value = intersection(disk, @view(geodesics[:, :, index]), β)
            @. column_view[height-y, :] = truncator(value)
        end
    end
end


@doc raw"""
    renderdisk(
        ::Val{:cpu},
        disk::AccretionDisk,
        geodesics::AbstractArray{<:Number},
        width::Int,
        height::Int,
        fov::Int,
    )

Dispatch method for rendering the disk using the CPU.
"""
function renderdisk(
    ::Val{:cpu},
    disk::AccretionDisk,
    geodesics::AbstractArray{<:Number},
    width::Int,
    height::Int,
    fov::Int,
)
    image_out = zeros(UInt8, (height, width, 3))

    mid_width = width ÷ 2

    Threads.@threads for x = (-mid_width):(mid_width-1)
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
        disk::AccretionDisk,
        geodesics::AbstractArray{<:Number};
        width::Int = 480,
        height::Int = 720,
        fov::Int = 200
    )

Render an accretion disk `disk` given an array of pre-calculated geodesics `geodesics` for the disk.

Returns an array of dimensions `width` by `height`, representing the rendered image. 

This method checks whether CUDA is installed, and dispatches a GPU kernel accordingly. Otherwise, the
rendering functions are executed in parallel on the CPU. 
"""
function renderdisk(
    disk::AccretionDisk,
    geodesics::AbstractArray{<:Number};
    width::Int = 720,
    height::Int = 480,
    fov::Int = 200,
)
    if CUDA.has_cuda_gpu()
        return renderdisk(:gpu, disk, geodesics, width, height, fov)
    else
        return renderdisk(:cpu, disk, geodesics, width, height, fov)
    end
end

renderdisk(s::Symbol, args...; kwargs...) = renderdisk(Val(s), args...; kwargs...)

export renderdisk
