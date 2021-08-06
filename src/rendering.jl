
function rrenderdisk(
    ::Val{:cpu},
    geodesics::AbstractArray{<:Number},
    width::Int,
    height::Int,
    fov::Int,
    disk::AccretionDisk
)

end

function renderdisk(
    geodesics::AbstractArray{<:Number},
    width::Int,
    height::Int,
    fov::Int,
    disk::AccretionDisk
)
    if CUDA.has_cuda_gpu()
        return renderdisk(:gpu, geodesics, width, height, fov, disk)
    else
        return renderdisk(:cpu, geodesics, width, height, fov, disk)
    end
end

renderdisk(s::Symbol, args...; kwargs...) = renderdisk(Val(s), args... ; kwargs...)
