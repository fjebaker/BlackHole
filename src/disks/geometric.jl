using Parameters: @with_kw

""" Definitions """

@with_kw struct GeometricDisk <: AccretionDisk
    α::Float64
    β::Float64  
    rinner::Float64 
    router::Float64
end

function intersection(g, d::GeometricDisk, β)

    intensity(t, x, y, z, t2, x2, y2, z2, r) = begin 
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