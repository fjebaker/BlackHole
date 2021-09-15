

"""
Super type of all space time metrics.
"""
abstract type Singularity end

geodesic(s::Singularity, pos, vel) = error("Not implemented for this type of Singularity.")
constrain(s::Singularity, r, dϕ) = error("Not implemented for this type of Singularity.")

"""
Eddington-Finkelstein metric describing a Black Hole in spherically symmetric space.
"""
struct EddingtonFinkelstein <: Singularity
    "Mass (arbitrary units)."
    m::Float64
end

"""
Christoffel components for Eddington-Finkelstein, and components of the geodesic equation.

`pos` is a position tuple/array, with t,r,θ,ϕ indexed by 1,2,3,4.
"""
module cc_edfink
ttt(m, pos) = 2.0 * m^2 / pos[2]^3
ttr(m, pos) = (2.0 * m^2 + m * pos[2]) / pos[2]^3
trr(m, pos) = 2.0 * (m^2 + m * pos[2]) / pos[2]^3
tθθ(m, pos) = -2.0 * m
tϕϕ(m, pos) = -2.0 * m * sin(pos[3])^2

rtt(m, pos) = -(2.0 * m^2 - m * pos[2]) / pos[2]^3
rtr(m, pos) = -2.0 * m^2 / pos[2]^3
rrr(m, pos) = -(2.0 * m^2 + m * pos[2]) / pos[2]^3
rθθ(m, pos) = 2.0 * m - pos[2]
rϕϕ(m, pos) = (2.0 * m - pos[2]) * sin(pos[3])^2

θrθ(m, pos) = 1.0 / pos[2]
θϕϕ(m, pos) = -cos(pos[3]) * sin(pos[3])

ϕrϕ(m, pos) = 1.0 / pos[2]
ϕθϕ(m, pos) = cos(pos[3]) / sin(pos[3])

"""
Euler Lagrange μ = t component
"""
function el_t(m::Real, pos, vel)::Float64
    -1.0 * (
        ttt(m, pos) * vel[1]^2 +
        trr(m, pos) * vel[2]^2 +
        tθθ(m, pos) * vel[3]^2 +
        2.0 * ttr(m, pos) * vel[1] * vel[2] +
        tϕϕ(m, pos) * vel[4]^2
    )
end

"""
Euler Lagrange μ = r component
"""
function el_r(m::Real, pos, vel)::Float64
    -1.0 * (
        rrr(m, pos) * vel[2]^2 +
        2.0 * rtr(m, pos) * vel[1] * vel[2] +
        rtt(m, pos) * vel[1]^2 +
        rϕϕ(m, pos) * vel[4]^2
    )
end

"""
Euler Lagrange μ = θ component
"""
function el_θ(m::Real, pos, vel)::Float64
    -1.0 * (2.0 * θrθ(m, pos) * vel[2] * vel[3] + θϕϕ(m, pos) * vel[4]^2)
end

"""
Euler Lagrange μ = ϕ component
"""
function el_ϕ(m::Real, pos, vel)::Float64
    -2.0 * (ϕrϕ(m, pos) * vel[2] * vel[4] + ϕθϕ(m, pos) * vel[3] * vel[4])
end

end

@doc raw"""

    geodesic(s::EddingtonFinkelstein, pos, vel)

Components of the Euler-Lagrange equation for the Eddington-Finkelstein metric.

Returns `(t, r, θ, ϕ)` components.
"""
function geodesic(s::EddingtonFinkelstein, pos, vel)
    (
        cc_edfink.el_t(s.m, pos, vel),
        cc_edfink.el_r(s.m, pos, vel),
        cc_edfink.el_θ(s.m, pos, vel),
        cc_edfink.el_ϕ(s.m, pos, vel),
    )
end

@doc raw"""

    constrain(s::EddingtonFinkelstein, r, dϕ)::Float64

Calculates the time component for a light-like 4-vector defined by `(t, r, θ=0, ϕ)`.
"""
function constrain(s::EddingtonFinkelstein, r, dϕ)::Float64
    (sqrt(-(2.0 * s.m * r - r^2) * dϕ^2 + 1) * r - 2.0 * s.m) / (2.0 * s.m - r)
end


# module exports
export Singularity, EddingtonFinkelstein
