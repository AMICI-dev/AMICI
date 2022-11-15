module TPL_MODULE_NAME

Pkg.add(["Symbolics", "OrdinaryDiffEq", "ModelingToolkit"])
export model, prob

using Symbolics

function xdot(u::Vector{Symbolics.Num}, p::Vector{Symbolics.Num}, t::Symbolics.Num)::Vector{Symbolics.Num}
    _p = p[begin:TPL_NP]
    k = p[(TPL_NP+1):(TPL_NP + TPL_NK)]
    tcl = p[(TPL_NP+TPL_NK+1):end]

    TPL_X_SYMS = u
    TPL_P_SYMS = _p
    TPL_K_SYMS = k
    TPL_TCL_SYMS = tcl

TPL_W_EQ

TPL_XDOT_EQ

    return TPL_XDOT_RET
end

function w(t::Float64, x::Vector{Float64}, p::Vector{Float64}, k::Vector{Float64}, tcl::Vector{Float64})::Vector{Float64}
    TPL_X_SYMS = x
    TPL_P_SYMS = p
    TPL_K_SYMS = k
    TPL_TCL_SYMS = tcl

TPL_W_EQ

    return TPL_W_RET
end

function x0(p::Vector{Float64}, k::Vector{Float64})::Vector{Float64}
    TPL_P_SYMS = p
    TPL_K_SYMS = k

TPL_X0_EQ

    return TPL_X0_RET
end

function x_solver(x::Vector{Float64})
    TPL_X_RDATA_SYMS = x

TPL_X_SOLVER_EQ

    return TPL_X_SOLVER_RET
end

function x_rdata(x::Vector{Float64}, tcl::Vector{Float64})::Vector{Float64}
    TPL_X_SYMS = x
    TPL_TCL_SYMS = tcl

TPL_X_RDATA_EQ

    return TPL_X_RDATA_RET
end

function tcl(x::Vector{Float64}, p::Vector{Float64}, k::Vector{Float64})::Vector{Float64}
    TPL_X_RDATA_SYMS = x
    TPL_P_SYMS = p
    TPL_K_SYMS = k

TPL_TOTAL_CL_EQ

    return TPL_TOTAL_CL_RET
end

function y(t::Float64, x::Vector{Float64}, p::Vector{Float64}, k::Vector{Float64}, tcl::Vector{Float64})::Vector{Float64}

    TPL_X_SYMS = x
    TPL_P_SYMS = p
    TPL_K_SYMS = k
    TPL_W_SYMS = w(t, x, p, k, tcl)

TPL_Y_EQ

    return TPL_Y_RET
end

function sigmay(y::Vector{Float64}, p::Vector{Float64}, k::Vector{Float64})::Vector{Float64}
    TPL_Y_SYMS = y
    TPL_P_SYMS = p
    TPL_K_SYMS = k

TPL_SIGMAY_EQ

    return TPL_SIGMAY_RET
end

function Jy(y::Vector{Float64}, my::Vector{Float64}, sigmay::Vector{Float64})::Vector{Float64}
    TPL_Y_SYMS = y
    TPL_MY_SYMS = my
    TPL_SIGMAY_SYMS = sigmay

TPL_JY_EQ

    return TPL_JY_RET
end

using AMICI: Model

model = Model(
    xdot=xdot,
    w=w,
    x0=x0,
    x_solver=x_solver,
    x_rdata=x_rdata,
    tcl=tcl,
    y=y,
    sigmay=sigmay,
    Jy=Jy,
)

using OrdinaryDiffEq: ODEProblem
using ModelingToolkit: modelingtoolkitize

tspan_ref = (0, 1)
x0_ref = zeros(Float64, TPL_NX)
p_ref = ones(Float64, TPL_NP + TPL_NK + TPL_NTCL)

pre_prob = ODEProblem(
    xdot,
    x0_ref,
    tspan_ref,
    p_ref,
)
time = @elapsed sys = modelingtoolkitize(pre_prob)
println("toolkitize finished after $time [s]")

prob = ODEProblem(
    sys,
    [],
    tspan_ref,
    jac=true,
    sparse=true
)

end