module AMICI

import Pkg
#Pkg.add(["OrdinaryDiffEq", "Zygote", "SciMLSensitivity"])

import Base.@kwdef

export ExpData, Model, run_simulation

@kwdef struct ExpData
    id::String
    ts::Vector{Float64}
    k::Vector{Float64}
    sigmay::Matrix{Float64}
    my::Matrix{Float64}
    pscale::Vector{String}
end

@kwdef struct Model
    xdot::Function
    w::Function
    x0::Function
    x_solver::Function
    x_rdata::Function
    tcl::Function
    y::Function
    sigmay::Function
    Jy::Function
end

import Pkg
#Pkg.add(["OrdinaryDiffEq", "Zygote", "SciMLSensitivity", "ModelingToolkit", "Symbolics", "SciMLBase", "Sundials"])


using OrdinaryDiffEq: solve, KenCarp4, FBDF, QNDF, Rosenbrock23, TRBDF2, Rodas4, RadauIIA5
using Sundials: CVODE_BDF
using SciMLBase: ODEProblem
using SciMLSensitivity: remake

function run_simulation(p::Vector{Float64}, model::Model, prob::ODEProblem, edata::ExpData, solver::String)::Tuple{Float64, Vector{Vector{Float64}}, Vector{Vector{Float64}}}
    up = map(transform_scale, zip(p, edata.pscale))

    # initialization & prepare
    x0 = model.x0(up, edata.k)
    tcl = model.tcl(x0, up, edata.k)

    if solver == "KenCarp4"
        solv = KenCarp4()
    elseif solver == "FBDF"
        solv = FBDF()
    elseif solver == "QNDF"
        solv = QNDF()
    elseif solver == "Rosenbrock23"
        solv = Rosenbrock23()
    elseif solver == "TRBDF2"
        solv = TRBDF2()
    elseif solver == "Rodas4"
        solv = Rodas4()
    elseif solver == "RadauIIA5"
        solv = RadauIIA5()
    elseif solver == "CVODE_BDF"
        solv = CVODE_BDF(linear_solver=:KLU)
    else
        solv = QNDF()
    end

    # simulate
    _prob = remake(prob,
                   p=[up; edata.k; tcl],
                   u0=model.x_solver(x0),
                   tspan=[min(edata.ts[begin], 0.0), edata.ts[end]])
    time = @elapsed sol = solve(_prob, solv, saveat=edata.ts,
                                abstol=1e-16, reltol=1e-8)
    println("julia ($solver) simulation time $time [s]")

    # observables & loss
    nt = length(edata.ts)
    x = sol.u
    y = map(it -> model.y(edata.ts[it], x[it], up, edata.k, tcl),
            range(; stop = nt))
    sigmay = map(it -> model.sigmay(y[it], up, edata.k),
                 range(; stop = nt))
    llhs = map(it -> model.Jy(y[it], edata.my[it, :], sigmay[it]),
               range(; stop = length(edata.ts)))
    return - sum(sum(llhs)), x, y
end


function transform_scale(t::Tuple{Float64, String})
    p, pscale = t
    if pscale == "none"
        return p
    elseif pscale == "ln"
        return exp(p)
    elseif pscale == "log10"
        return 10 ^ p
    end
end

end