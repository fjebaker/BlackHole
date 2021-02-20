using Parameters: @with_kw

""" Structure Definitions """

@with_kw struct GeometricDisk <: AccretionDisk
    α::Float64
    β::Float64  
    rinner::Float64 
    router::Float64
end

@with_kw struct OpticallyThinDisk <: AccretionDisk
    α::Float64
    β::Float64  
    rinner::Float64 
    router::Float64

    profile::Function
end


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
    x, y, z = g[:, begin]
    z = 0 # rounding errors in ODE solve
    x, y, z = transform(x, y, z, d.α, d.β + β)
    
    mapfoldl(
        (index) -> begin
            x2, y2, z2 = g[:, index]
            z2 = 0
            radius = sqrt(x2^2 + y2^2 + z2^2)
            x2, y2, z2 = transform(x2, y2, z2, d.α, d.β + β)
            
            ret = intensity(
                x, y, z, x2, y2, z2, radius
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

function intersection(g, d::GeometricDisk, β)

    intensity(x, y, z, x2, y2, z2, r) = begin 
        ret = zeros(Float64, 3)
        if d.rinner < r < d.router
            if z != z2 # no divide by zero thanks
                if 0 <= z / (z-z2) < 1
                    ret[1] = 255
                    ret[2] = 255
                    ret[3] = 255
                end
            end
        end
        ret
    end

    calcintersect(g, d, β, intensity)
end


export GeometricDisk, intersection