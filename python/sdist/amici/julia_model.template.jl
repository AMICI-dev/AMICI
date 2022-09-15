module TPL_MODEL_NAME_model

using AMICI

export model

function xdot(u::Array{Float64}, p::Array{Float64}, t::Float64)
    _p = p[begin:TPL_NP]
    k = p[(TPL_NP+1):(TPL_NP + TPL_NK)]
    tcl = p[(TPL_NP+TPL_NK+1):end]

    TPL_X_SYMS = u
    TPL_P_SYMS = _p
    TPL_K_SYMS = k
    TPL_TCL_SYMS = tcl
    TPL_W_SYMS = w(t, u, p, k, tcl)

TPL_XDOT_EQ

    return TPL_XDOT_RET
end

function w(t::Float64, x::Array{Float64}, p::Array{Float64}, k::Array{Float64}, tcl::Array{Float64})
    TPL_X_SYMS = x
    TPL_P_SYMS = p
    TPL_K_SYMS = k
    TPL_TCL_SYMS = tcl

TPL_W_EQ

    return TPL_W_RET
end

function x0(p::Array{Float64}, k::Array{Float64})
    TPL_P_SYMS = p
    TPL_K_SYMS = k

TPL_X0_EQ

    return TPL_X0_RET
end

function x_solver(x::Array{Float64})
    TPL_X_RDATA_SYMS = x

TPL_X_SOLVER_EQ

    return TPL_X_SOLVER_RET
end

function x_rdata(x::Array{Float64}, tcl::Array{Float64})
    TPL_X_SYMS = x
    TPL_TCL_SYMS = tcl

TPL_X_RDATA_EQ

    return TPL_X_RDATA_RET
end

function tcl(x::Array{Float64}, p::Array{Float64}, k::Array{Float64})
    TPL_X_RDATA_SYMS = x
    TPL_P_SYMS = p
    TPL_K_SYMS = k

TPL_TOTAL_CL_EQ

    return TPL_TOTAL_CL_RET
end

function y(t::Float64, x::Array{Float64}, p::Array{Float64}, k::Array{Float64}, tcl::Array{Float64})

    TPL_X_SYMS = x
    TPL_P_SYMS = p
    TPL_K_SYMS = k
    TPL_W_SYMS = w(t, x, p, k, tcl)

TPL_Y_EQ

    return TPL_Y_RET
end

function sigmay(y::Array{Float64}, p::Array{Float64}, k::Array{Float64})
    TPL_Y_SYMS = y
    TPL_P_SYMS = p
    TPL_K_SYMS = k

TPL_SIGMAY_EQ

    return TPL_SIGMAY_RET
end

function Jy(y::Array{Float64}, my::Array{Float64}, sigmay::Array{Float64})
    TPL_Y_SYMS = y
    TPL_MY_SYMS = my
    TPL_SIGMAY_SYMS = sigmay

TPL_JY_EQ

    TPL_JY_RET
end

model = Model(
    xdot::Function = xdot,
    w::Function = w,
    x0::Function = x0,
    x_solver::Function = x_solver,
    x_rdata::Function = x_rdata,
    tcl::Function = tcl,
    y::Function = y,
    sigmay::Function = sigmay,
    Jy::Function = Jy,
)

end
