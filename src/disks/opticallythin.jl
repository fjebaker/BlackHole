using Parameters: @with_kw

""" Structure Definitions """

abstract type OpticallyThinDisk <: AccretionDisk end

@with_kw struct ThinDisk <: OpticallyThinDisk
    α::Float64
    β::Float64  
    rinner::Float64 
    router::Float64

    s::Singularity
    profile::Function
end

@with_kw struct GaussianThinDisk <: OpticallyThinDisk
    α::Float64
    β::Float64  
    rinner::Float64 
    router::Float64

    s::Singularity
    σ::Float64
end

""" Helper functions """

function orthonormal(r, s::Singularity, dt, dr, dθ, dϕ)
    """ transform given coordinates into the orthonormal reference frame """
    (
        sqrt(r / 2 * s.m + r) * dt,
        2 * s.m / (sqrt(r * (2 * s.m + r))) * dt + sqrt((2 * s.m + r) / r) * dr,
        dθ * r,
        dϕ * r
    )
end


""" Exported Methods """

profile(x, d::OpticallyThinDisk) = d.profile(x, d)

# Gaussian profile specifics
profile(x, d::GaussianThinDisk) = begin
    μ = (d.router - d.rinner) / 2 + d.rinner
    if !(d.rinner < x < d.router)
        return 0
    else
        y = ( 1 / d.σ * sqrt(2 * π)) * exp(-0.5 * ((x - μ) / d.σ)^2)
        return max(y, 0)
    end
end

function intersection(g, d::OpticallyThinDisk, β)

    intensity(t, x, y, z, t2, x2, y2, z2, r) = begin
        ret = zeros(Float64, 3)
        if d.rinner < r < d.router
            if z != z2 # no divide by zero thanks
                if 0 <= z / (z-z2) < 1
                    
                    # relevant components of four velocity
                    dt, dx, dy, dz = (t2 - t, x2 - x, y2 - y, z2 - z)
                    
                    dr, dθ, dϕ = xyz2spher(x2, y2, z2, dx, dy, dz)
                    
                    dt, dr, dθ, dϕ = orthonormal(
                        r, s, dt, dr, dθ, dϕ
                    )

                    # angle calculation
                    θ_width = atan(
                        dθ, sqrt(dr^2 + dϕ^2)
                    )

                    # thickness, with ~magic~ constant
                    thickness = 3 * profile(r, d) / sin(abs(θ_width))

                    ret[1] = thickness
                    ret[2] = thickness
                    ret[3] = thickness
                end
            end
        end
        ret
    end

    calcintersect(g, d, β, intensity)
end


export ThinDisk, GaussianThinDisk, profile, intersection