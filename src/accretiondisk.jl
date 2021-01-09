using Parameters: @with_kw

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


function transform(x, y, z, d::AccretionDisk)
    y, z = (
        y * cos(d.β) - z * sin(d.β),
        z * cos(d.β) + y * sin(d.β)
    )

    (
        x * cos(d.α) - z * sin(d.α), 
        y, 
        z * cos(d.α) + x * sin(d.α)
    )
end

function calcintersect(g, d::AccretionDisk, intensity::Function)
    x, y, z = g.curve[:, begin]
    z = 0 # rounding errors in ODE solve
    x, y, z = transform(x, y, z, d)
    
    mapfoldl(
        (index) -> begin
            x2, y2, z2 = g.curve[:, index]
            z2 = 0
            radius = sqrt(x2^2 + y2^2 + z2^2)
            x2, y2, z2 = transform(x2, y2, z2, d)
            
            ret = intensity(
                x, y, z, x2, y2, z2, radius
            )

            # update state
            x, y, z = x2, y2, z2
            
            ret
        end,
        +,
        2:size(g.curve)[2],
        init=zeros(Float64, 3)
    )
end

intersection(g, d::AccretionDisk) = error("Not implemented for disk type.")

function intersection(g::Geodesic, d::GeometricDisk)

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

    calcintersect(g, d, intensity)
end