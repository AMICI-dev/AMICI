#include "amici/rdata.h"

#include "amici/misc.h"
#include "amici/model.h"
#include "amici/edata.h"
#include "amici/symbolic_functions.h"
#include "amici/solver.h"
#include "amici/exception.h"
#include "amici/steadystateproblem.h"
#include "amici/backwardproblem.h"
#include "amici/forwardproblem.h"

#include <cstring>

namespace amici {

ReturnData::ReturnData(Solver const &solver, const Model &model)
    : ReturnData(
          model.getTimepoints(), model.np(), model.nk(), model.nx_rdata,
          model.nx_solver, model.nxtrue_rdata, model.ny, model.nytrue, model.nz,
          model.nztrue, model.ne, model.nJ, model.nplist(), model.nMaxEvent(),
          model.nt(), solver.getNewtonMaxSteps(), model.nw,
          model.getParameterScale(), model.o2mode, solver.getSensitivityOrder(),
          static_cast<SensitivityMethod>(solver.getSensitivityMethod())) {}

ReturnData::ReturnData(std::vector<realtype> ts, int np, int nk, int nx,
                       int nx_solver, int nxtrue, int ny, int nytrue, int nz,
                       int nztrue, int ne, int nJ, int nplist, int nmaxevent,
                       int nt, int newton_maxsteps, int nw,
                       std::vector<ParameterScaling> pscale,
                       SecondOrderMode o2mode, SensitivityOrder sensi,
                       SensitivityMethod sensi_meth)
    : ts(std::move(ts)), np(np), nk(nk), nx(nx), nx_solver(nx_solver),
      nxtrue(nxtrue), ny(ny), nytrue(nytrue), nz(nz), nztrue(nztrue), ne(ne),
      nJ(nJ), nplist(nplist), nmaxevent(nmaxevent), nt(nt), nw(nw),
      newton_maxsteps(newton_maxsteps), pscale(std::move(pscale)),
      o2mode(o2mode), sensi(sensi), sensi_meth(sensi_meth), x_solver(nx_solver),
      sx_solver(nx_solver, nplist), x_rdata(nx), sx_rdata(nx, nplist),
      nroots(ne) {
    xdot.resize(nx_solver, getNaN());

    J.resize(nx_solver * nx_solver, getNaN());

    // initialize with 0.0, so we only need to write non-zero values
    z.resize(nmaxevent * nz, 0.0);
    sigmaz.resize(nmaxevent * nz, 0.0);

    rz.resize(nmaxevent * nz, 0.0);
    x.resize(nt * nx, 0.0);
    y.resize(nt * ny, 0.0);
    sigmay.resize(nt * ny, 0.0);
    w.resize(nt * nw, 0.0);

    preeq_numsteps.resize(3, 0);
    preeq_numlinsteps.resize(newton_maxsteps * 2, 0);
    posteq_numsteps.resize(3, 0);
    posteq_numlinsteps.resize(newton_maxsteps * 2, 0);

    if (nt > 0) {
        numsteps.resize(nt, 0);
        numrhsevals.resize(nt, 0);
        numerrtestfails.resize(nt, 0);
        numnonlinsolvconvfails.resize(nt, 0);
        order.resize(nt, 0);

        if (sensi_meth == SensitivityMethod::adjoint &&
            sensi >= SensitivityOrder::first) {
            numstepsB.resize(nt, 0);
            numrhsevalsB.resize(nt, 0);
            numerrtestfailsB.resize(nt, 0);
            numnonlinsolvconvfailsB.resize(nt, 0);
        }
    }

    x0.resize(nx, getNaN());

    x_ss.resize(nx, getNaN());

    res.resize(nt * nytrue, 0.0);

    llh = getNaN();
    chi2 = getNaN();
    if (sensi >= SensitivityOrder::first) {
        sllh.resize(nplist, getNaN());
        sx0.resize(nx * nplist, getNaN());
        sx_ss.resize(nx * nplist, getNaN());

        if (sensi_meth == SensitivityMethod::forward ||
            sensi >= SensitivityOrder::second) {
            // for second order we can fill in from the augmented states
            sx.resize(nt * nx * nplist, 0.0);
            sy.resize(nt * ny * nplist, 0.0);
            sz.resize(nmaxevent * nz * nplist, 0.0);
            srz.resize(nmaxevent * nz * nplist, 0.0);

            FIM.resize(nplist * nplist, 0.0);
            sres.resize(nt * nytrue * nplist, 0.0);
        }

        ssigmay.resize(nt * ny * nplist, 0.0);
        ssigmaz.resize(nmaxevent * nz * nplist, 0.0);
        if (sensi >= SensitivityOrder::second) {
            s2llh.resize(nplist * (nJ - 1), getNaN());
            if (sensi_meth == SensitivityMethod::forward)
                s2rz.resize(nmaxevent * nztrue * nplist * nplist, 0.0);
        }
    }
}

void ReturnData::processPreEquilibration(SteadystateProblem const &preeq,
                                         Model &model) {
    ModelContext mc(&model);
    readSimulationState(preeq.getSimulationState(), model);
    model.fx_rdata(x_rdata, x_solver);
    x_ss = x_rdata.getVector();
    if (sensi >= SensitivityOrder::first) {
        model.fsx_rdata(sx_rdata, sx_solver);
        for (int ip = 0; ip < nplist; ip++)
            std::copy_n(sx_rdata.data(ip), nx, &sx_ss.at(ip * nx));
    }
    /* Get cpu time for Newton solve in seconds */
    preeq_cpu_time = preeq.getCPUTime() / 1000;
    preeq_status = static_cast<int>(preeq.getNewtonStatus());
    preeq_wrms = preeq.getResidualNorm();
    if (preeq.getNewtonStatus() == NewtonStatus::newt_sim)
        preeq_t = preeq.getSteadyStateTime();
    writeSlice(preeq.getNumSteps(), preeq_numsteps);
    writeSlice(preeq.getNumLinSteps(), preeq_numlinsteps);
}

void ReturnData::processPostEquilibration(SteadystateProblem const &posteq,
                                          Model &model, ExpData const *edata) {
    ModelContext mc(&model);
    for (int it = 0; it < nt; it++) {
        auto t = model.getTimepoint(it);
        if (std::isinf(t)) {
            readSimulationState(posteq.getSimulationState(), model);
            getDataOutput(it, model, edata);
        }
    }
    /* Get cpu time for Newton solve in seconds */
    posteq_cpu_time = posteq.getCPUTime() / 1000;
    posteq_status = static_cast<int>(posteq.getNewtonStatus());
    posteq_wrms = posteq.getResidualNorm();
    if (posteq.getNewtonStatus() == NewtonStatus::newt_sim)
        preeq_t = posteq.getSteadyStateTime();
    writeSlice(posteq.getNumSteps(), preeq_numsteps);
    writeSlice(posteq.getNumLinSteps(), preeq_numlinsteps);
}

void ReturnData::processForwardProblem(ForwardProblem const &fwd, Model &model,
                                       ExpData const *edata) {
    ModelContext mc(&model);
    if (edata)
        initializeObjectiveFunction();

    model.fx_rdata(x_rdata, fwd.getInitialState());
    x0 = x_rdata.getVector();

    if (computingFSA()) {
        model.fsx_rdata(sx_rdata, fwd.getInitialStateSensitivity());
        for (int ip = 0; ip < nplist; ip++)
            std::copy_n(sx_rdata.data(ip), nx, &sx0.at(ip * nx));
    }

    // process timpoint data
    for (int it = 0; it <= fwd.getTimePointCounter(); it++) {
        readSimulationState(fwd.getSimulationStateTimepoint(it), model);
        getDataOutput(it, model, edata);
    }
    // check for integration failure but consider postequilibration
    for (int it = fwd.getTimePointCounter() + 1; it < nt; it++)
        if (!std::isinf(model.getTimepoint(it)))
            invalidate(it);

    // process event data
    if (nz > 0) {
        auto rootidx = fwd.getRootIndexes();
        for (int iroot = 0; iroot <= fwd.getEventCounter(); iroot++) {
            readSimulationState(fwd.getSimulationStateEvent(iroot), model);
            getEventOutput(iroot, t, rootidx.at(iroot), model, edata);
        }
    }
}

void ReturnData::getDataOutput(int it, Model &model, ExpData const *edata) {
    model.fx_rdata(x_rdata, x_solver);
    std::copy_n(x_rdata.data(), nx, &x.at(it * nx));
    model.getExpression(slice(w, it, nw), ts[it], x_solver);

    model.getObservable(slice(y, it, ny), ts[it], x_solver);
    model.getObservableSigma(slice(sigmay, it, ny), it, edata);
    if (edata) {
        model.addObservableObjective(llh, it, x_solver, *edata);
        fres(it, *edata);
        fchi2(it);
    }

    if (sensi >= SensitivityOrder::first && nplist > 0) {

        model.getObservableSigmaSensitivity(slice(ssigmay, it, nplist * ny), it,
                                            edata);

        if (sensi_meth == SensitivityMethod::forward) {
            getDataSensisFSA(it, model, edata);
        } else {
            if (edata) {
                model.addPartialObservableObjectiveSensitivity(
                    sllh, s2llh, it, x_solver, *edata);
            }
        }
    }
}

void ReturnData::getDataSensisFSA(int it, Model &model, ExpData const *edata) {
    model.fsx_rdata(sx_rdata, sx_solver);
    for (int ip = 0; ip < nplist; ip++) {
        std::copy_n(sx_rdata.data(ip), nx, &sx.at((it * nplist + ip) * nx));
    }

    model.getObservableSensitivity(slice(sy, it, nplist * ny), ts[it], x_solver,
                                   sx_solver);

    if (edata) {
        model.addObservableObjectiveSensitivity(sllh, s2llh, it, x_solver,
                                                sx_solver, *edata);
        fsres(it, *edata);
        fFIM(it);
    }
}

void ReturnData::getEventOutput(int iroot, realtype t, std::vector<int> rootidx,
                                Model &model, ExpData const *edata) {

    for (int ie = 0; ie < ne; ie++) {
        if (rootidx.at(ie) != 1 || nroots.at(ie) >= nmaxevent)
            continue;

        /* get event output */
        model.getEvent(slice(z, nroots.at(ie), nz), ie, t, x_solver);
        /* if called from fillEvent at last timepoint,
         then also get the root function value */
        if (t == model.getTimepoint(nt - 1))
            model.getEventRegularization(slice(rz, nroots.at(ie), nz), ie, t,
                                         x_solver);

        if (edata) {
            model.getEventSigma(slice(sigmaz, nroots.at(ie), nz), ie,
                                nroots.at(ie), t, edata);
            model.addEventObjective(llh, ie, nroots.at(ie), t, x_solver,
                                    *edata);

            /* if called from fillEvent at last timepoint,
               add regularization based on rz */
            if (t == model.getTimepoint(nt - 1))
                model.addEventObjectiveRegularization(llh, ie, nroots.at(ie), t,
                                                      x_solver, *edata);
        }

        if (sensi >= SensitivityOrder::first) {
            if (sensi_meth == SensitivityMethod::forward) {
                getEventSensisFSA(iroot, ie, t, model, edata);
            } else {
                if (edata)
                    model.addPartialEventObjectiveSensitivity(
                        sllh, s2llh, ie, nroots.at(ie), t, x_solver, *edata);
            }
        }
        nroots.at(ie)++;
    }
}

void ReturnData::getEventSensisFSA(int iroot, int ie, realtype t, Model &model,
                                   ExpData const *edata) {
    if (t == model.getTimepoint(nt - 1)) {
        // call from fillEvent at last timepoint
        model.getUnobservedEventSensitivity(
            slice(sz, nroots.at(ie), nz * nplist), ie);
        model.getEventRegularizationSensitivity(
            slice(srz, nroots.at(ie), nz * nplist), ie, t, x_solver, sx_solver);
    } else {
        model.getEventSensitivity(slice(sz, nroots.at(ie), nz * nplist), ie, t,
                                  x_solver, sx_solver);
    }

    if (edata) {
        model.addEventObjectiveSensitivity(sllh, s2llh, ie, nroots.at(ie), t,
                                           x_solver, sx_solver, *edata);
    }
}

void ReturnData::processBackwardProblem(ForwardProblem const &fwd,
                                        BackwardProblem const &bwd,
                                        Model &model) {
    std::vector<realtype> llhS0(model.nJ * model.nplist(), 0.0);
    auto xB = bwd.getAdjointState();
    auto xQB = bwd.getAdjointQuadrature();
    auto sx0 = fwd.getInitialStateSensitivity();

    for (int iJ = 0; iJ < model.nJ; iJ++) {
        if (iJ == 0) {
            for (int ip = 0; ip < model.nplist(); ++ip) {
                llhS0[ip] = 0.0;
                for (int ix = 0; ix < model.nxtrue_solver; ++ix) {
                    llhS0[ip] += xB[ix] * sx0.at(ix, ip);
                }
            }
        } else {
            for (int ip = 0; ip < model.nplist(); ++ip) {
                llhS0[ip + iJ * model.nplist()] = 0.0;
                for (int ix = 0; ix < model.nxtrue_solver; ++ix) {
                    llhS0[ip + iJ * model.nplist()] +=
                        xB[ix + iJ * model.nxtrue_solver] * sx0.at(ix, ip) +
                        xB[ix] * sx0.at(ix + iJ * model.nxtrue_solver, ip);
                }
            }
        }
    }

    for (int iJ = 0; iJ < model.nJ; iJ++) {
        for (int ip = 0; ip < model.nplist(); ip++) {
            if (iJ == 0) {
                sllh.at(ip) -= llhS0[ip] + xQB[ip * model.nJ];
            } else {
                s2llh.at(iJ - 1 + ip * (model.nJ - 1)) -=
                    llhS0[ip + iJ * model.nplist()] + xQB[iJ + ip * model.nJ];
            }
        }
    }
}

void ReturnData::processSolver(Solver const &solver) {

    cpu_time = solver.getCpuTime();

    writeSlice(solver.getNumSteps(), numsteps);
    writeSlice(solver.getNumRhsEvals(), numrhsevals);
    writeSlice(solver.getNumErrTestFails(), numerrtestfails);
    writeSlice(solver.getNumNonlinSolvConvFails(), numnonlinsolvconvfails);
    writeSlice(solver.getLastOrder(), order);

    cpu_timeB = solver.getCpuTimeB();
    writeSlice(solver.getNumStepsB(), numstepsB);
    writeSlice(solver.getNumRhsEvalsB(), numrhsevalsB);
    writeSlice(solver.getNumErrTestFailsB(), numerrtestfailsB);
    writeSlice(solver.getNumNonlinSolvConvFailsB(), numnonlinsolvconvfailsB);
}

void ReturnData::storeJacobianAndDerivativeInReturnData(
    ForwardProblem const &fwd, Model &model) {
    // dont use SimulationState here since we explicitely want to evaluate at
    // final timepoint
    auto t = fwd.getTime();
    auto x = fwd.getState();
    auto dx = fwd.getStateDerivative();

    storeJacobianAndDerivativeInReturnData(t, x, dx, model);
}

void ReturnData::storeJacobianAndDerivativeInReturnData(
    SteadystateProblem const &posteq, Model &model) {
    auto t = std::numeric_limits<realtype>::infinity();
    auto x = posteq.getState();
    auto dx = AmiVector(model.nx_solver);

    storeJacobianAndDerivativeInReturnData(t, x, dx, model);
}

void ReturnData::readSimulationState(SimulationState const &state,
                                     Model &model) {
    x_solver = state.x;
    if (computingFSA())
        sx_solver = state.sx;
    t = state.t;
    model.setModelState(state.state);
}

void ReturnData::storeJacobianAndDerivativeInReturnData(realtype t,
                                                        const AmiVector &x,
                                                        const AmiVector &dx,
                                                        Model &model) {
    AmiVector xdot(nx_solver);
    model.fxdot(t, x, dx, xdot);
    this->xdot = xdot.getVector();

    SUNMatrixWrapper J(SUNMatrixWrapper(nx_solver, nx_solver));
    model.fJ(t, 0.0, x, dx, xdot, J.get());
    // CVODES uses colmajor, so we need to transform to rowmajor
    for (int ix = 0; ix < model.nx_solver; ix++) {
        for (int jx = 0; jx < model.nx_solver; jx++) {
            this->J[ix * model.nx_solver + jx] =
                J.data()[ix + model.nx_solver * jx];
        }
    }
}

void ReturnData::invalidate(const int it_start) {
    if (it_start >= nt)
        return

    invalidateLLH();
    invalidateSLLH();

    for (int it = it_start; it < nt; it++){
        for (int ix = 0; ix < nx; ix++)
            x.at(ix + nx * it) = getNaN();
        for (int iy = 0; iy < ny; iy++)
            y.at(iy + ny * it) = getNaN();
        for (int iw = 0; iw < nw; iw++)
            w.at(iw + nw * it) = getNaN();
    }

    if (!sx.empty()) {
        for (int it = it_start; it < nt; it++){
            for (int ip = 0; ip < nplist; ip++) {
                for (int ix = 0; ix < nx; ix++)
                    sx.at(ix + nx*(ip + it*nplist)) = getNaN();
            }
        }
    }
    if(!sy.empty()) {
        for (int it = it_start; it < nt; it++){
            for (int ip = 0; ip < nplist; ip++) {
                for (int iy = 0; iy < ny; iy++)
                    sy.at(iy + ny*(ip + it*nplist)) = getNaN();
            }
        }
    }
}

void ReturnData::invalidateLLH()
{
    llh = getNaN();
    chi2 = getNaN();
}

void ReturnData::invalidateSLLH() {
    std::fill(sllh.begin(), sllh.end(), getNaN());
    std::fill(s2llh.begin(), s2llh.end(), getNaN());
}

void ReturnData::applyChainRuleFactorToSimulationResults(const Model &model) {
    // chain-rule factor: multiplier for am_p
    std::vector<realtype> coefficient(nplist, 1.0);
    std::vector<realtype> pcoefficient(nplist, 1.0);

    std::vector<realtype> unscaledParameters = model.getParameters();
    unscaleParameters(unscaledParameters, model.getParameterScale(), unscaledParameters);

    std::vector<realtype> augcoefficient(np, 1.0);

    if (sensi == SensitivityOrder::second && o2mode == SecondOrderMode::full) {
        for (int ip = 0; ip < np; ++ip) {
            switch (pscale[ip]) {
            case ParameterScaling::log10:
                augcoefficient.at(ip) = unscaledParameters.at(ip) * log(10);
                break;
            case ParameterScaling::ln:
                augcoefficient.at(ip) = unscaledParameters.at(ip);
                break;
            case ParameterScaling::none:
                break;
            }
        }
    }

    for (int ip = 0; ip < nplist; ++ip) {
        switch (pscale[model.plist(ip)]) {
        case ParameterScaling::log10:
            coefficient.at(ip) = log(10.0);
            pcoefficient.at(ip) = unscaledParameters.at(model.plist(ip)) * log(10);
            break;
        case ParameterScaling::ln:
            pcoefficient.at(ip) = unscaledParameters.at(model.plist(ip));
            break;
        case ParameterScaling::none:
            break;
        }
    }

    if (sensi >= SensitivityOrder::first) {
        // recover first order sensitivies from states for adjoint sensitivity
        // analysis
        if (sensi == SensitivityOrder::second && o2mode == SecondOrderMode::full) {
            if (sensi_meth == SensitivityMethod::adjoint) {
                for (int ip = 0; ip < nplist; ++ip)
                    for (int ix = 0; ix < nxtrue; ++ix)
                        for (int it = 0; it < nt; ++it)
                            sx.at(ix + nxtrue*(ip + it*nplist)) =
                                x.at(it * nx + nxtrue + ip * nxtrue + ix);

                for (int ip = 0; ip < nplist; ++ip)
                    for (int iy = 0; iy < nytrue; ++iy)
                        for (int it = 0; it < nt; ++it)
                            sy.at(iy + nytrue*(ip + it*nplist)) =
                                y.at(it * ny + nytrue + ip * nytrue + iy);

                for (int ip = 0; ip < nplist; ++ip)
                    for (int iz = 0; iz < nztrue; ++iz)
                        for (int it = 0; it < nt; ++it)
                            sz.at(iz + nztrue*(ip + it*nplist)) =
                                z.at(it * nz + nztrue + ip * nztrue + iz);
            }
        }

        for (int ip = 0; ip < nplist; ++ip)
            sllh.at(ip) *= pcoefficient.at(ip);

        if(!sres.empty())
            for (int iyt = 0; iyt < nytrue*nt; ++iyt)
                for (int ip = 0; ip < nplist; ++ip)
                    sres.at((iyt * nplist + ip)) *= pcoefficient.at(ip);

        if(!FIM.empty())
            for (int ip = 0; ip < nplist; ++ip)
                for (int jp = 0; jp < nplist; ++jp)
                    FIM.at(jp + ip * nplist) *= pcoefficient.at(ip)*pcoefficient.at(jp);

#define chainRule(QUANT, IND1, N1T, N1, IND2, N2)                             \
    if (!s##QUANT.empty())                                                    \
        for (int IND1 = 0; (IND1) < (N1T); ++(IND1))                          \
            for (int ip = 0; ip < nplist; ++ip)                               \
                for (int IND2 = 0; (IND2) < (N2); ++(IND2)) {                 \
                    s##QUANT.at(((IND2)*nplist + ip)*(N1) + (IND1)) *=        \
                        pcoefficient.at(ip);                                  \
                }

        chainRule(x, ix, nxtrue, nx, it, nt);
        chainRule(y, iy, nytrue, ny, it, nt);
        chainRule(sigmay, iy, nytrue, ny, it, nt);
        chainRule(z, iz, nztrue, nz, ie, nmaxevent);
        chainRule(sigmaz, iz, nztrue, nz, ie, nmaxevent);
        chainRule(rz, iz, nztrue, nz, ie, nmaxevent);
        chainRule(x0, ix, nxtrue, nx, it, 1);
    }

    if (o2mode == SecondOrderMode::full) { // full
        if (!s2llh.empty() && !sllh.empty()) {
            for (int ip = 0; ip < nplist; ++ip) {
                for (int iJ = 1; iJ < nJ; ++iJ) {
                    s2llh[ip * nplist + (iJ - 1)] *=
                            pcoefficient.at(ip) * augcoefficient[iJ - 1];
                    if (model.plist(ip) == iJ - 1)
                        s2llh[ip * nplist + (iJ - 1)] +=
                                sllh.at(ip) * coefficient.at(ip);
                }
            }
        }

#define s2ChainRule(QUANT, IND1, N1T, N1, IND2, N2)                           \
    if (!s##QUANT.empty())                                                    \
        for (int ip = 0; ip < nplist; ++ip)                                   \
            for (int iJ = 1; iJ < nJ; ++iJ)                                   \
                for (int IND1 = 0; IND1 < N1T; ++IND1)                        \
                    for (int IND2 = 0; IND2 < N2; ++IND2) {                   \
                    s##QUANT.at((IND2*nplist + ip)*N1 + IND1 + iJ*N1T) *=     \
                        pcoefficient.at(ip) * augcoefficient[iJ - 1];         \
                        if (model.plist(ip) == iJ - 1)                        \
                            s##QUANT.at((IND2*nplist + ip)*N1 + IND1 + iJ*N1T) +=   \
                                s##QUANT.at((IND2*nplist + ip)*N1 + IND1) *   \
                                    coefficient[ip];                          \
}


        s2ChainRule(x, ix, nxtrue, nx, it, nt);
        s2ChainRule(y, iy, nytrue, ny, it, nt);
        s2ChainRule(sigmay, iy, nytrue, ny, it, nt);
        s2ChainRule(z, iz, nztrue, nz, ie, nmaxevent);
        s2ChainRule(sigmaz, iz, nztrue, nz, ie, nmaxevent);
        s2ChainRule(rz, iz, nztrue, nz, ie, nmaxevent);
    }

    if (o2mode == SecondOrderMode::directional) { // directional
            for (int ip = 0; ip < nplist; ++ip) {
                s2llh.at(ip) *= pcoefficient.at(ip);
                s2llh.at(ip) += model.k()[nk - nplist + ip] * sllh.at(ip) /
                             unscaledParameters[model.plist(ip)];
        }

#define s2vecChainRule(QUANT, IND1, N1T, N1, IND2, N2)                        \
    if (!s##QUANT.empty())                                                    \
        for (int ip = 0; ip < nplist; ++ip)                                   \
            for (int IND1 = 0; IND1 < N1T; ++IND1)                            \
                for (int IND2 = 0; IND2 < N2; ++IND2) {                       \
                    s##QUANT.at((IND2*nplist + ip)*N1 + IND1 + N1T) *=        \
                        pcoefficient.at(ip);                                  \
                    s##QUANT.at((IND2*nplist + ip)*N1 + IND1 + N1T) +=        \
                        model.k()[nk - nplist + ip] *                         \
                        s##QUANT.at((IND2*nplist + ip)*N1 + IND1) /           \
                        unscaledParameters[model.plist(ip)];                  \
                }

        s2vecChainRule(x, ix, nxtrue, nx, it, nt);
        s2vecChainRule(y, iy, nytrue, ny, it, nt);
        s2vecChainRule(sigmay, iy, nytrue, ny, it, nt);
        s2vecChainRule(z, iz, nztrue, nz, ie, nmaxevent);
        s2vecChainRule(sigmaz, iz, nztrue, nz, ie, nmaxevent);
        s2vecChainRule(rz, iz, nztrue, nz, ie, nmaxevent);
    }
}

void ReturnData::initializeObjectiveFunction()
{
    llh = 0.0;
    chi2 = 0.0;
    std::fill(sllh.begin(),sllh.end(), 0.0);
    std::fill(s2llh.begin(),s2llh.end(), 0.0);
}

void ReturnData::fres(const int it, const ExpData &edata) {
    if ( res.empty())
        return;

    auto observedData = edata.getObservedDataPtr(it);
    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * edata.nytrue();
        int iyt = iy + it * ny;
        if (!edata.isSetObservedData(it, iy))
            continue;
        res.at(iyt_true) =
        (y.at(iyt) - observedData[iy]) / sigmay.at(iyt);
    }
}

void ReturnData::fchi2(const int it) {
    if (res.empty())
        return;

    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * nytrue;
        chi2 += pow(res.at(iyt_true), 2);
    }
}

void ReturnData::fsres(const int it, const ExpData &edata) {
    if (sres.empty())
        return;

    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * edata.nytrue();
        int iyt = iy + it * ny;
        if (!edata.isSetObservedData(it, iy))
            continue;
        for (int ip = 0; ip < nplist; ++ip) {
            sres.at(iyt_true * nplist + ip) =
            sy.at(iy + ny * (ip + it * nplist)) /
            sigmay.at(iyt);
        }
    }
}

void ReturnData::fFIM(const int it) {
    if (sres.empty())
        return;

    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * nytrue;
        for (int ip = 0; ip < nplist; ++ip) {
            for (int jp = 0; jp < nplist; ++jp) {
                FIM.at(ip + nplist * jp) += sres.at(iyt_true * nplist + ip) *
                    sres.at(iyt_true * nplist + jp);
            }
        }
    }
}

ModelContext::ModelContext(Model *model)
    : model(model), original_state(model->getModelState()) {}

ModelContext::~ModelContext()
{
    restore();
}

void ModelContext::restore()
{
    model->setModelState(original_state);
}


} // namespace amici
