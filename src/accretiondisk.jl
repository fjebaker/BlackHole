using Parameters: @with_kw

""" Helper Functions """

function transform(x, y, z, α, β)
    y, z = (
        y * cos(β) - z * sin(β),
        z * cos(β) + y * sin(β)
    )

    (
        x * cos(α) - z * sin(α), 
        y, 
        z * cos(α) + x * sin(α)
    )
end

function calcintersect(g, d::AccretionDisk, β, intensity::Function)
    t, x, y, z = g[:, begin]
    z = 0 # rounding errors in ODE solve
    x, y, z = transform(x, y, z, d.α, d.β + β)
    
    mapfoldl(
        (index) -> begin
            t2, x2, y2, z2 = g[:, index]
            z2 = 0
            radius = sqrt(x2^2 + y2^2 + z2^2)
            x2, y2, z2 = transform(x2, y2, z2, d.α, d.β + β)
            
            ret = intensity(
                t, x, y, z, t2, x2, y2, z2, radius
            )

            # update state
            x, y, z = x2, y2, z2
            
            ret
        end,
        +,
        2:size(g)[2],
        init=zeros(Float64, 3)
    )
end


""" Exported Functions """

intersection(g, d::AccretionDisk) = error("Not implemented for disk type.")

export intersection