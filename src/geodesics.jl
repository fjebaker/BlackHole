using StaticArrays
using DifferentialEquations

@doc raw"""

    observer(u, t, integrator)::Bool

Discrete callback function which ensures that the geodesics are smooth, by returning `false`
if the integrator gets too close to the event horizon, or if the next proposed step is
arbitrarily distant from the previous.

To prevent bad results, this observer ensures the geodesic does not leave a defined chart
of "good" spacetime values, by terminating the integration proceedure early.
"""
@views function observer(u, t, integrator)::Bool
    r = u[:, 1][2]
    r <= 10.0 || abs(integrator.uprev[:, 1][2] - r) > 50.0 # make sure the radius change is gradual
end


@doc raw"""

    diffprob!(du, u, s::Singularity, λ)

ODEProblem integration function. Here `λ` is the proper time step used by the integrator.

`u` unpacks to position 4 vector `x = u[:, 1]` and velocity 4 vector `v = u[:, 2]`.
"""
function diffprob!(du, u, s::Singularity, λ)
    x, v = u[:, 1], u[:, 2]
    du[:, 1] .= v
    du[:, 2] .= geodesic(s, x, v)
end

@doc raw"""

    unpack(solutions::Vector{OrdinaryDiffEq.ODECompositeSolution})::Array{Float64,3}

Method for unpacking the integration solutions into a standardised 3-dimensional array.
The array is guarunteed to give the same number of indexes per solution, by interpolating
a fixed number of points (fixed to `640` for technical reasons related to GPU warp 
configurations) between the start and end time of the integrator.

The ouput array has dimensions `4 x 640 x length(solutions)`.
"""
function unpack(solutions::Vector{OrdinaryDiffEq.ODECompositeSolution})::Array{Float64,3}

    # maximum integration length so we can standardize size
    # in the current implementation, this is closely tied to the GPU maximum threads
    # that the kernel can run, given a thread per integration step. Only change if you
    # know what you're doing.
    maxlength::Int64 = 640

    # allocate output matrix
    output::Array{Float64,3} = zeros(Float64, (4, maxlength, length(solutions)))

    # use lazy foreach in J1.6
    Iterators.foreach(enumerate(solutions)) do (ind, sol)

        # range of time to generate from
        tmin, tmax = sol.t[1], sol.t[end]

        i = 1
        for tstep in range(tmin, tmax, length = maxlength)
            # transform coordinates
            output[:, i, ind] .= spher2xyz(sol(tstep))
            i += 1
        end

    end

    output
end

@doc raw"""

    function calcgeodesics(
        s::Singularity;
        num::Int=1500,
        Δϕ::Float64=0.004,
    )::Array{Float64,3}

Calculate integrated geodesic curves for a given singularity.

Performs `num` geodesic traces, each differing in initial velocity 4-vector by an angle `Δϕ`
in the `x, y` plane.

Note, the method for extracting the results from the integration solutions, and the actual
method of integration, is quite tentative at the moment, and requires some investigation.

If `num` is too low, it can lead to aliasing problems in the final render.

If `Δϕ` is too low, the integrator calculates light paths which swing around and terminate
in odd locations, leading to "jets" streaming out of the black hole into the plane of the disk
(where there really should not be any jets).

# Extended help

The root cause of the integration issue lies probably in how the integration chart was adapted
with the `observer` function. A more discriminating observer may improve the renderings
significantly.
"""
function calcgeodesics(
    s::Singularity;
    num::Int = 1300,
    Δϕ::Float64 = 0.004,
)::Array{Float64,3}
    x₀ = SA_F64[0.0, 100.0, π/2.0, 0.0]  # static array
    λs = (0.0, 500.0)

    solutions = Vector{OrdinaryDiffEq.ODECompositeSolution}(undef, num)

    for i = 1:num
        dϕ = i * Δϕ / 1000.0
        v₀ = SA_F64[constrain(s, x₀[2], dϕ), -1.0, 0.0, dϕ]

        # configure the chart through a discrete callback
        chartbounds = DiscreteCallback(observer, terminate!)

        # inital input
        Matrix{Float64}(hcat(x₀, v₀))

        prob = ODEProblem(diffprob!, Matrix{Float64}(hcat(x₀, v₀)), λs, s)

        sol = solve(prob, callback = chartbounds)
        solutions[i] = sol
    end

    unpack(solutions)
end


# module exports
export calcgeodesics
