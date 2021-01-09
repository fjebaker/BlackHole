

function truncator(px)
    convert(
        UInt8, 
        trunc(min(max(0, px), 255))
    )
end

function render_column(
    d::AccretionDisk, 
    geodesics::AbstractArray{<:Geodesic},
    height::Int,
    fov_index::Int,
    x::Int
    )

    res = zeros(Float64, (height, 3))
    mid = height ÷ 2 # integer division
    n = length(geodesics)

    for y in 1:height-1
        β = atan(y-mid, x) # atan2
        r = sqrt(x^2 + (y-mid)^2) # distance from middle of image
        i = convert(
            Int,
            trunc(1 + r / fov_index * n)
        )
        if i < n
            res[height-y, :] .= intersection(geodesics[i], d, β)
        end
    end
    res
end

function renderdisk(
    d::AccretionDisk,
    geodesics::AbstractArray{<:Geodesic}
    ;
    height::Int=720,
    width::Int=1080,
    fov_index::Int=200
    )

    data = zeros(UInt8, (height, width, 3))

    mid = width ÷ 2 # integer division

    for x in (-mid):(mid - 1)
        col = render_column(
            d, geodesics, height, fov_index, x
        )

        data[:, x + mid + 1, :] .= truncator.(col)
    end

    data
end