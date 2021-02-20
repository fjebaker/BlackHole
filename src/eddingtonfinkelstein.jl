module __ef
    """ 
    Christoffel Components and Geodesic Equation for Eddington Finkelstein metric
    """
    ttt(m, pos) = 2 * m^2 / pos[2]^3
    ttr(m, pos) = (2 * m^2 + m * pos[2]) / pos[2]^3
    trr(m, pos) = 2 * (m^2 + m * pos[2]) / pos[2]^3
    tθθ(m, pos) = - 2 * m
    tϕϕ(m, pos) = - 2 * m * sin(pos[3])^2

    rtt(m, pos) = - (2 * m^2 - m * pos[2]) / pos[2]^3
    rtr(m, pos) = - 2 * m^2 / pos[2]^3
    rrr(m, pos) = - (2 * m^2 + m * pos[2]) / pos[2]^3
    rθθ(m, pos) = 2 * m - pos[2]
    rϕϕ(m, pos) = (2 * m - pos[2]) * sin(pos[3])^2

    θrθ(m, pos) = 1/pos[2]
    θϕϕ(m, pos) = - cos(pos[3]) * sin(pos[3])

    ϕrϕ(m, pos) = 1/pos[2]
    ϕθϕ(m, pos) = cos(pos[3]) / sin(pos[3])

    function geodesic(pos, vel, m)
        """ Returns components of the Euler-Lagrange (Geodesic) equation """
        (
            # μ=t
            - 1 * (
                ttt(m, pos) * vel[1]^2 
                + trr(m, pos) * vel[2]^2 
                + tθθ(m, pos) * vel[3]^2 
                + 2 * ttr(m, pos) * vel[1] * vel[2]
                + tϕϕ(m, pos) * vel[4]^2
            ),
            # μ=r
            - 1 * (
                rrr(m, pos) * vel[2]^2
                + 2 * rtr(m, pos) * vel[1] * vel[2]
                + rtt(m, pos) * vel[1]^2
                + rϕϕ(m, pos) * vel[4]^2
            ),
            # μ=θ
            - 1 * (
                2 * θrθ(m, pos) * vel[2] * vel[3]
                + θϕϕ(m, pos) * vel[4]^2
            ),
            # μ=ϕ
            - 2 * (
                ϕrϕ(m, pos) * vel[2] * vel[4]
                + ϕθϕ(m, pos) * vel[3] * vel[4]
            )
        )
    end

    function constrain(m, r, dϕ)
        """ Calculates time component such that vector (t, r, θ=0, ϕ) is light-like """
        (sqrt(
                -(2 * m * r - r^2) * dϕ^2 + 1
            ) * r - 2 * m
        ) / (2 * m - r)
    end

end

struct EddingtonFinkelstein <: Singularity
    m::Float64
end

function getgeodesic(s::EddingtonFinkelstein)
    Geodesic{EddingtonFinkelstein}
end

function getintfunc(s::EddingtonFinkelstein)
    getintfunc(__ef.geodesic)
end

function getconstraint(s::EddingtonFinkelstein)
    __ef.constrain
end

