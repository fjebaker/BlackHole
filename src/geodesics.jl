using DifferentialEquations

import Base: axes, getindex, size

struct Geodesic{T}
    curve::Array{Float64,2}

    function Geodesic{T}(sol, λ_span) where {T<:Singularity}
        begin
            """ convert to x,y,z and instantiate """

            new(
                hcat(
                    map(
                        (λ) -> begin
                            t, r, θ, ϕ = sol(λ)
                            Float64[
                                t, # store time
                                r * sin(θ) * cos(ϕ),
                                r * sin(θ) * sin(ϕ),
                                r * cos(θ),
                            ]
                        end,
                        λ_span[1]:λ_span[2],
                    )...,
                ),
            )
        end
    end
end

axes(g::Geodesic, i::Int) = axes(g.curve, i)
getindex(g::Geodesic, i...) = getindex(g.curve, i...)
size(g::Geodesic) = size(g.curve)

getgeodesic(s::Singularity) = error("Not defined for $s !")
getintfunc(s::Singularity) = error("Not defined for $s !")
getconstraint(s::Singularity) = error("Not defined for $s !")

function getintfunc(geodesic::Function)::Function
    (du, u, m, λ) -> begin
        x, v = u[:, 1], u[:, 2]
        du[:, 1] .= v # velocity
        du[:, 2] .= geodesic(x, v, m)
    end
end

function observer(u, t, integrator)::Bool
    r = u[:, 1][2]
    r <= 5.0 || abs(integrator.uprev[:, 1][2] - r) > 50
end
function stopintegration!(integrator)
    terminate!(integrator)
end

function calcgeodesics(s::Singularity; num::Int=1000, Δϕ::Float64=0.005)
    i_x = Vector{Float64}([0, 100, π / 2, 0])
    i_v = Vector{Float64}([0, -1, 0, 0])
    λ_span = (0.0, 300.0)

    factory = getgeodesic(s)
    stepper! = getintfunc(s)
    constrain = getconstraint(s)

    mapper = (i) -> begin
        dϕ = i * Δϕ / 1000
        i_v[1] = constrain(s.m, i_x[2], dϕ)
        i_v[4] = convert(Float64, dϕ)

        charter = DiscreteCallback(observer, stopintegration!)
        prob = ODEProblem(stepper!, hcat(i_x, i_v), λ_span, s.m)

        sol = solve(prob; callback=charter)
        factory(sol, λ_span)
    end

    map(mapper, 0:num)
end

export calcgeodesics, plot!
