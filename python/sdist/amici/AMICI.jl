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
#Pkg.add(["OrdinaryDiffEq", "Zygote", "SciMLSensitivity", "ModelingToolkit", "Symbolics", "SciMLBase"])


using OrdinaryDiffEq: solve, KenCarp4
using SciMLBase: ODEProblem
using SciMLSensitivity: remake

function run_simulation(p::Vector{Float64}, model::Model, prob::ODEProblem, edata::ExpData)::Tuple{Float64, Vector{Vector{Float64}}, Vector{Vector{Float64}}}
    up = map(transform_scale, zip(p, edata.pscale))

    # initialization & prepare
    x0 = model.x0(up, edata.k)
    tcl = model.tcl(x0, up, edata.k)

    # simulate
    _prob = remake(prob,
                   p=[up; edata.k; tcl],
                   u0=model.x_solver(x0),
                   tspan=[min(edata.ts[begin], 0.0), edata.ts[end]])
    time = @elapsed sol = solve(_prob, KenCarp4(), saveat=edata.ts)
    println("julia simulation time $time [s]")

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