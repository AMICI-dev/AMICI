module AMICI

import Pkg
Pkg.add(["OrdinaryDiffEq", "Zygote", "SciMLSensitivity"])

export ExpData, Model

struct ExpData
    ts::Vector{Float64}
    k::Vector{Float64}
    sigmay::Matrix{Float64}
    my::Matrix{Float64}
    pscale::Vector{String}
end

struct Model
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
Pkg.add(["OrdinaryDiffEq", "Zygote", "SciMLSensitivity"])

using DifferentialEquations: ODEProblem, solve, KenCarp4
using SciMLSensitivity: remake

function get_solver()

function run_simulation(p::Vector{Float64}, model::Model, prob::ODEProblem, edata::ExpData)
    up = map(transform_scale, zip(p, edata.pscale))

    # initialization & prepare
    x0 = model.x0(up, edata.k)
    tcl = model.tcl(x0, up, edata.k)

    # simulate
    _prob = remake(prob, p=hstack(up, edata.k, tcl), u0=model.x_solver(x0))
    x = solve(_prob, KenCarp4(), saveat=edata.ts)

    # observables & loss
    y = map(it -> model.y(edata.ts[it], x[it, :] , up, edata.k, tcl),
            range(; stop = length(edata.ts)))
    sigmay = mapslices(yi -> model.sigmay(yi, up, k), y, dims=1)
    llhs = map(it -> model.Jy(y[it,:], edata.my[it,:], sigmay[it, :]),
               range(; stop = length(edata.ts)))
    return - sum(filter(!isnan, llhs)
end


function transform_scale(p::Float64, pscale::String)
    if pscale == 'lin'
        return p
    elseif pscale == 'ln'
        return exp(p)
    elseif pscale == 'log10'
        return 10 ^ p
    end
end