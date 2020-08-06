#include "amici/rdata.h"

#include "amici/backwardproblem.h"
#include "amici/edata.h"
#include "amici/exception.h"
#include "amici/forwardproblem.h"
#include "amici/misc.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/steadystateproblem.h"
#include "amici/symbolic_functions.h"

#include <cstring>

namespace amici {

ReturnData::ReturnData(Solver const &solver, const Model &model)
    : ReturnData(model.getTimepoints(), model.np(), model.nk(), model.nx_rdata,
                 model.nx_solver, model.nxtrue_rdata, model.nx_solver_reinit,
                 model.ny, model.nytrue, model.nz, model.nztrue, model.ne,
                 model.nJ, model.nplist(), model.nMaxEvent(), model.nt(),
                 solver.getNewtonMaxSteps(), model.nw,
                 model.getParameterScale(), model.o2mode,
                 solver.getSensitivityOrder(), solver.getSensitivityMethod(),
                 solver.getReturnDataReportingMode()) {}

ReturnData::ReturnData(std::vector<realtype> ts, int np, int nk, int nx,
                       int nx_solver, int nxtrue, int nx_solver_reinit, int ny, int nytrue, int nz,
                       int nztrue, int ne, int nJ, int nplist, int nmaxevent,
                       int nt, int newton_maxsteps, int nw,
                       std::vector<ParameterScaling> pscale,
                       SecondOrderMode o2mode, SensitivityOrder sensi,
                       SensitivityMethod sensi_meth, RDataReporting rdrm)
    : ts(std::move(ts)), np(np), nk(nk), nx(nx), nx_solver(nx_solver),
      nxtrue(nxtrue), nx_solver_reinit(nx_solver_reinit), ny(ny), nytrue(nytrue), nz(nz), nztrue(nztrue), ne(ne),
      nJ(nJ), nplist(nplist), nmaxevent(nmaxevent), nt(nt), nw(nw),
      newton_maxsteps(newton_maxsteps), pscale(std::move(pscale)),
      o2mode(o2mode), sensi(sensi), sensi_meth(sensi_meth),
      rdata_reporting(rdrm), x_solver_(nx_solver), sx_solver_(nx_solver, nplist),
      x_rdata_(nx), sx_rdata_(nx, nplist), nroots_(ne) {

    switch (rdata_reporting) {
    case RDataReporting::full:
        initializeFullReporting();
        break;

    case RDataReporting::residuals:
        initializeResidualReporting();
        break;

    case RDataReporting::likelihood:
        initializeLikelihoodReporting();
        break;
    }
}

void ReturnData::initializeLikelihoodReporting() {
    llh = getNaN();
    chi2 = getNaN();
    if (sensi >= SensitivityOrder::first) {
        sllh.resize(nplist, getNaN());
        if (sensi >= SensitivityOrder::second)
            s2llh.resize(nplist * (nJ - 1), getNaN());

        if (sensi_meth == SensitivityMethod::forward ||
            sensi >= SensitivityOrder::second)
            FIM.resize(nplist * nplist, 0.0);
    }
}

void ReturnData::initializeResidualReporting() {
    y.resize(nt * ny, 0.0);
    sigmay.resize(nt * ny, 0.0);
    res.resize(nt * nytrue, 0.0);
    if ((sensi_meth == SensitivityMethod::forward &&
         sensi >= SensitivityOrder::first)
        || sensi >= SensitivityOrder::second) {

        sy.resize(nt * ny * nplist, 0.0);
        ssigmay.resize(nt * ny * nplist, 0.0);
        sres.resize(nt * nytrue * nplist, 0.0);
    }
}

void ReturnData::initializeFullReporting() {

    initializeLikelihoodReporting();
    initializeResidualReporting();

    xdot.resize(nx_solver, getNaN());

    J.resize(nx_solver * nx_solver, getNaN());

    // initialize with 0.0, so we only need to write non-zero values
    z.resize(nmaxevent * nz, 0.0);
    sigmaz.resize(nmaxevent * nz, 0.0);

    rz.resize(nmaxevent * nz, 0.0);
    x.resize(nt * nx, 0.0);
    w.resize(nt * nw, 0.0);

    preeq_numsteps.resize(3, 0);
    preeq_status.resize(3, SteadyStateStatus::not_run);
    posteq_numsteps.resize(3, 0);
    posteq_status.resize(3, SteadyStateStatus::not_run);

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

    if (sensi >= SensitivityOrder::first) {
        sx0.resize(nx * nplist, getNaN());
        sx_ss.resize(nx * nplist, getNaN());

        if (sensi_meth == SensitivityMethod::forward ||
            sensi >= SensitivityOrder::second) {
            // for second order we can fill in from the augmented states
            sx.resize(nt * nx * nplist, 0.0);
            sz.resize(nmaxevent * nz * nplist, 0.0);
            srz.resize(nmaxevent * nz * nplist, 0.0);
        }

        ssigmay.resize(nt * ny * nplist, 0.0);
        ssigmaz.resize(nmaxevent * nz * nplist, 0.0);
        if (sensi >= SensitivityOrder::second &&
            sensi_meth == SensitivityMethod::forward)
            s2rz.resize(nmaxevent * nztrue * nplist * nplist, 0.0);
    }
}

void ReturnData::processSimulationObjects(SteadystateProblem const *preeq,
                                          ForwardProblem const *fwd,
                                          BackwardProblem const *bwd,
                                          SteadystateProblem const *posteq,
                                          Model &model, Solver const &solver,
                                          ExpData const *edata) {
    ModelContext mc(&model);

    processSolver(solver);

    if (preeq)
        processPreEquilibration(*preeq, model);

    if (fwd)
        processForwardProblem(*fwd, model, edata);
    else
        invalidate(0);

    if (posteq)
        processPostEquilibration(*posteq, model, edata);

    if (fwd && !posteq)
        storeJacobianAndDerivativeInReturnData(*fwd, model);
    else if (posteq)
        storeJacobianAndDerivativeInReturnData(*posteq, model);

    if (fwd && bwd)
        processBackwardProblem(*fwd, *bwd, preeq, model);
    else if (solver.computingASA())
        invalidateSLLH();

    applyChainRuleFactorToSimulationResults(model);
}

void ReturnData::processPreEquilibration(SteadystateProblem const &preeq,
                                         Model &model) {
    readSimulationState(preeq.getFinalSimulationState(), model);

    if (!x_ss.empty()) {
        model.fx_rdata(x_rdata_, x_solver_);
        writeSlice(x_rdata_, x_ss);
    }
    if (!sx_ss.empty() && sensi >= SensitivityOrder::first) {
        model.fsx_rdata(sx_rdata_, sx_solver_);
        for (int ip = 0; ip < nplist; ip++)
            writeSlice(sx_rdata_[ip], slice(sx_ss, ip, nx));
    }
    /* Get cpu time for Newton solve in milliseconds */
    preeq_cpu_time = preeq.getCPUTime();
    preeq_cpu_timeB = preeq.getCPUTimeB();
    preeq_numstepsB = preeq.getNumStepsB();
    preeq_wrms = preeq.getResidualNorm();
    preeq_status = preeq.getSteadyStateStatus();
    if (preeq_status[1] == SteadyStateStatus::success)
        preeq_t = preeq.getSteadyStateTime();
    if (!preeq_numsteps.empty())
        writeSlice(preeq.getNumSteps(), preeq_numsteps);
    if (!preeq.getNumLinSteps().empty() && !preeq_numlinsteps.empty()) {
        preeq_numlinsteps.resize(newton_maxsteps * 2, 0);
        writeSlice(preeq.getNumLinSteps(), preeq_numlinsteps);
    }
}

void ReturnData::processPostEquilibration(SteadystateProblem const &posteq,
                                          Model &model, ExpData const *edata) {
    for (int it = 0; it < nt; it++) {
        auto t = model.getTimepoint(it);
        if (std::isinf(t)) {
            readSimulationState(posteq.getFinalSimulationState(), model);
            getDataOutput(it, model, edata);
        }
    }
    /* Get cpu time for Newton solve in milliseconds */
    posteq_cpu_time = posteq.getCPUTime();
    posteq_cpu_timeB = posteq.getCPUTimeB();
    posteq_numstepsB = posteq.getNumStepsB();
    posteq_wrms = posteq.getResidualNorm();
    posteq_status = posteq.getSteadyStateStatus();
    if (posteq_status[1] == SteadyStateStatus::success)
        posteq_t = posteq.getSteadyStateTime();
    if (!posteq_numsteps.empty())
        writeSlice(posteq.getNumSteps(), posteq_numsteps);
    if (!posteq.getNumLinSteps().empty() && !posteq_numlinsteps.empty()) {
        posteq_numlinsteps.resize(newton_maxsteps * 2, 0);
        writeSlice(posteq.getNumLinSteps(), posteq_numlinsteps);
    }
}

void ReturnData::processForwardProblem(ForwardProblem const &fwd, Model &model,
                                       ExpData const *edata) {
    if (edata)
        initializeObjectiveFunction();

    auto initialState = fwd.getInitialSimulationState();
    if (initialState.x.getLength() == 0)
        return; // if x wasn't set forward problem failed during initialization
    readSimulationState(initialState, model);

    if (!x0.empty()) {
        model.fx_rdata(x_rdata_, x_solver_);
        writeSlice(x_rdata_, x0);
    }

    if (!sx0.empty()) {
        model.fsx_rdata(sx_rdata_, sx_solver_);
        for (int ip = 0; ip < nplist; ip++)
            writeSlice(sx_rdata_[ip], slice(sx0, ip, nx));
    }

    // process timpoint data
    realtype tf = fwd.getFinalTime();
    for (int it = 0; it < model.nt(); it++) {
        if (model.getTimepoint(it) <= tf) {
            readSimulationState(fwd.getSimulationStateTimepoint(it), model);
            getDataOutput(it, model, edata);
        } else {
            // check for integration failure but consider postequilibration
            if (!std::isinf(model.getTimepoint(it)))
                invalidate(it);
        }
    }

    // process event data
    if (nz > 0) {
        auto rootidx = fwd.getRootIndexes();
        for (int iroot = 0; iroot <= fwd.getEventCounter(); iroot++) {
            readSimulationState(fwd.getSimulationStateEvent(iroot), model);
            getEventOutput(iroot, t_, rootidx.at(iroot), model, edata);
        }
    }
}

void ReturnData::getDataOutput(int it, Model &model, ExpData const *edata) {
    if (!x.empty()) {
        model.fx_rdata(x_rdata_, x_solver_);
        writeSlice(x_rdata_, slice(x, it, nx));
    }
    if (!w.empty())
        model.getExpression(slice(w, it, nw), ts[it], x_solver_);
    if (!y.empty())
        model.getObservable(slice(y, it, ny), ts[it], x_solver_);
    if (!sigmay.empty())
        model.getObservableSigma(slice(sigmay, it, ny), it, edata);

    if (edata) {
        if (!isNaN(llh))
            model.addObservableObjective(llh, it, x_solver_, *edata);
        fres(it, model, *edata);
        fchi2(it);
    }

    if (sensi >= SensitivityOrder::first && nplist > 0) {

        if (!ssigmay.empty())
            model.getObservableSigmaSensitivity(slice(ssigmay, it, nplist * ny),
                                                it, edata);

        if (sensi_meth == SensitivityMethod::forward) {
            getDataSensisFSA(it, model, edata);
        } else if (edata && !sllh.empty()) {
            model.addPartialObservableObjectiveSensitivity(
                sllh, s2llh, it, x_solver_, *edata);
        }
    }
}

void ReturnData::getDataSensisFSA(int it, Model &model, ExpData const *edata) {
    if (!sx.empty()) {
        model.fsx_rdata(sx_rdata_, sx_solver_);
        for (int ip = 0; ip < nplist; ip++) {
            writeSlice(sx_rdata_[ip],
                       slice(sx, it * nplist + ip, nx));
        }
    }

    if (!sy.empty()) {
        model.getObservableSensitivity(slice(sy, it, nplist * ny), ts[it],
                                       x_solver_, sx_solver_);
    }

    if (edata) {
        if (!sllh.empty())
            model.addObservableObjectiveSensitivity(sllh, s2llh, it, x_solver_,
                                                    sx_solver_, *edata);
        fsres(it, model, *edata);
        fFIM(it, model, *edata);
    }
}

void ReturnData::getEventOutput(int iroot, realtype t, std::vector<int> rootidx,
                                Model &model, ExpData const *edata) {

    for (int ie = 0; ie < ne; ie++) {
        if (rootidx.at(ie) != 1 || nroots_.at(ie) >= nmaxevent)
            continue;

        /* get event output */
        if (!z.empty())
            model.getEvent(slice(z, nroots_.at(ie), nz), ie, t, x_solver_);
        /* if called from fillEvent at last timepoint,
         then also get the root function value */
        if (t == model.getTimepoint(nt - 1))
            if (!rz.empty())
                model.getEventRegularization(slice(rz, nroots_.at(ie), nz), ie,
                                             t, x_solver_);

        if (edata) {
            if (!sigmaz.empty())
                model.getEventSigma(slice(sigmaz, nroots_.at(ie), nz), ie,
                                    nroots_.at(ie), t, edata);
            if (!isNaN(llh))
                model.addEventObjective(llh, ie, nroots_.at(ie), t, x_solver_,
                                        *edata);

            /* if called from fillEvent at last timepoint,
               add regularization based on rz */
            if (t == model.getTimepoint(nt - 1) && !isNaN(llh)) {
                model.addEventObjectiveRegularization(
                    llh, ie, nroots_.at(ie), t, x_solver_, *edata);
            }
        }

        if (sensi >= SensitivityOrder::first) {
            if (sensi_meth == SensitivityMethod::forward) {
                getEventSensisFSA(iroot, ie, t, model, edata);
            } else if (edata && !sllh.empty()) {
                model.addPartialEventObjectiveSensitivity(
                    sllh, s2llh, ie, nroots_.at(ie), t, x_solver_, *edata);
            }
        }
        nroots_.at(ie)++;
    }
}

void ReturnData::getEventSensisFSA(int iroot, int ie, realtype t, Model &model,
                                   ExpData const *edata) {
    if (t == model.getTimepoint(nt - 1)) {
        // call from fillEvent at last timepoint
        if (!sz.empty())
            model.getUnobservedEventSensitivity(
                slice(sz, nroots_.at(ie), nz * nplist), ie);
        if (!srz.empty())
            model.getEventRegularizationSensitivity(
                slice(srz, nroots_.at(ie), nz * nplist), ie, t, x_solver_,
                sx_solver_);
    } else if (!sz.empty()) {
        model.getEventSensitivity(slice(sz, nroots_.at(ie), nz * nplist), ie,
                                  t, x_solver_, sx_solver_);
    }

    if (edata && !sllh.empty()) {
        model.addEventObjectiveSensitivity(sllh, s2llh, ie, nroots_.at(ie),
                                           t, x_solver_, sx_solver_, *edata);
    }
}

void ReturnData::processBackwardProblem(ForwardProblem const &fwd,
                                        BackwardProblem const &bwd,
                                        SteadystateProblem const *preeq,
                                        Model &model) {
    if (sllh.empty())
        return;
    readSimulationState(fwd.getInitialSimulationState(), model);

    std::vector<realtype> llhS0(model.nJ * model.nplist(), 0.0);
    auto xB = bwd.getAdjointState();
    auto xQB = bwd.getAdjointQuadrature();

    if (preeq && preeq->hasQuadrature()) {
        handleSx0Backward(model, *preeq, xQB);
    } else {
        handleSx0Forward(model, llhS0, xB);
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

void ReturnData::handleSx0Backward(const Model &model,
                                   SteadystateProblem const &preeq,
                                   AmiVector &xQB) const {
    /* If preequilibration is run in adjoint mode, the scalar product of sx0
       with its adjoint counterpart (see handleSx0Forward()) is not necessary:
       the actual simulation is "extended" by the preequilibration time.
       At initialization (at t=-inf), the adjoint state is in steady state (= 0)
       and so is the scalar product. Instead of the scalar product, the
       quadratures xQB from preequilibration contribute to the gradient
       (see example notebook on equilibration for further documentation). */
    auto xQBpreeq = preeq.getAdjointQuadrature();
    for (int ip = 0; ip < model.nplist(); ++ip)
        xQB[ip] += xQBpreeq[ip];
}

void ReturnData::handleSx0Forward(const Model &model,
                                  std::vector<realtype> &llhS0,
                                  AmiVector &xB) const {
    /* If preequilibration is run in forward mode or is not needed, then adjoint
       sensitivity analysis still needs the state sensitivities at t=0 (sx0),
       to compute the gradient. For each parameter, the scalar product of sx0
       with its adjoint counterpart contributes to the gradient. */
    for (int iJ = 0; iJ < model.nJ; iJ++) {
        if (iJ == 0) {
            for (int ip = 0; ip < model.nplist(); ++ip) {
                llhS0[ip] = 0.0;
                for (int ix = 0; ix < model.nxtrue_solver; ++ix) {
                    llhS0[ip] += xB[ix] * sx_solver_.at(ix, ip);
                }
            }
        } else {
            for (int ip = 0; ip < model.nplist(); ++ip) {
                llhS0[ip + iJ * model.nplist()] = 0.0;
                for (int ix = 0; ix < model.nxtrue_solver; ++ix) {
                    llhS0[ip + iJ * model.nplist()] +=
                        xB[ix + iJ * model.nxtrue_solver] * sx_solver_.at(ix, ip) +
                        xB[ix] * sx_solver_.at(ix + iJ * model.nxtrue_solver, ip);
                }
            }
        }
    }
}

void ReturnData::processSolver(Solver const &solver) {

    cpu_time = solver.getCpuTime();
    if (!numsteps.empty())
        numsteps = solver.getNumSteps();
    if (!numrhsevals.empty())
        numrhsevals = solver.getNumRhsEvals();
    if (!numerrtestfails.empty())
        numerrtestfails = solver.getNumErrTestFails();
    if (!numnonlinsolvconvfails.empty())
        numnonlinsolvconvfails = solver.getNumNonlinSolvConvFails();
    if (!order.empty())
        order = solver.getLastOrder();

    cpu_timeB = solver.getCpuTimeB();
    if (!numstepsB.empty())
        numstepsB = solver.getNumStepsB();
    if (!numrhsevalsB.empty())
        numrhsevalsB = solver.getNumRhsEvalsB();
    if (!numerrtestfailsB.empty())
        numerrtestfailsB = solver.getNumErrTestFailsB();
    if (!numnonlinsolvconvfailsB.empty())
        numnonlinsolvconvfailsB = solver.getNumNonlinSolvConvFailsB();
}

void ReturnData::readSimulationState(SimulationState const &state,
                                     Model &model) {
    x_solver_ = state.x;
    dx_solver_ = state.dx;
    if (computingFSA() || state.t == model.t0())
        sx_solver_ = state.sx;
    t_ = state.t;
    model.setModelState(state.state);
}

void ReturnData::invalidate(const int it_start) {
    if (it_start >= nt)
        return;

    invalidateLLH();
    invalidateSLLH();

    if (!x.empty())
        std::fill(x.begin() + nx * it_start, x.end(), getNaN());
    if (!y.empty())
        std::fill(y.begin() + ny * it_start, y.end(), getNaN());
    if (!w.empty())
        std::fill(w.begin() + nw * it_start, w.end(), getNaN());

    if (!sx.empty())
        std::fill(sx.begin() + nx * nplist * it_start, sx.end(), getNaN());
    if (!sy.empty())
        std::fill(sy.begin() + ny * nplist * it_start, sy.end(), getNaN());
}

void ReturnData::invalidateLLH() {
    llh = getNaN();
    chi2 = getNaN();
}

void ReturnData::invalidateSLLH() {
    if (!sllh.empty()) {
        std::fill(sllh.begin(), sllh.end(), getNaN());
        std::fill(s2llh.begin(), s2llh.end(), getNaN());
    }
}

void ReturnData::applyChainRuleFactorToSimulationResults(const Model &model) {
    // chain-rule factor: multiplier for am_p
    std::vector<realtype> coefficient(nplist, 1.0);
    std::vector<realtype> pcoefficient(nplist, 1.0);

    std::vector<realtype> unscaledParameters = model.getParameters();
    unscaleParameters(unscaledParameters, model.getParameterScale(),
                      unscaledParameters);

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
            pcoefficient.at(ip) =
                unscaledParameters.at(model.plist(ip)) * log(10);
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
        if (sensi == SensitivityOrder::second
            && o2mode == SecondOrderMode::full
            && sensi_meth == SensitivityMethod::adjoint) {
            if (!sx.empty() && !x.empty())
                for (int ip = 0; ip < nplist; ++ip)
                    for (int ix = 0; ix < nxtrue; ++ix)
                        for (int it = 0; it < nt; ++it)
                            sx.at(ix + nxtrue * (ip + it * nplist)) =
                                x.at(it * nx + nxtrue + ip * nxtrue + ix);

            if (!sy.empty() && !y.empty())
                for (int ip = 0; ip < nplist; ++ip)
                    for (int iy = 0; iy < nytrue; ++iy)
                        for (int it = 0; it < nt; ++it)
                            sy.at(iy + nytrue * (ip + it * nplist)) =
                                y.at(it * ny + nytrue + ip * nytrue + iy);

            if (!sz.empty() && !z.empty())
                for (int ip = 0; ip < nplist; ++ip)
                    for (int iz = 0; iz < nztrue; ++iz)
                        for (int it = 0; it < nt; ++it)
                            sz.at(iz + nztrue * (ip + it * nplist)) =
                                z.at(it * nz + nztrue + ip * nztrue + iz);

        }

        if (!sllh.empty())
            for (int ip = 0; ip < nplist; ++ip)
                sllh.at(ip) *= pcoefficient.at(ip);


        if (!sres.empty())
            for (int iyt = 0; iyt < nytrue * nt; ++iyt)
                for (int ip = 0; ip < nplist; ++ip)
                    sres.at((iyt * nplist + ip)) *= pcoefficient.at(ip);

        if(!FIM.empty())
            for (int ip = 0; ip < nplist; ++ip)
                for (int jp = 0; jp < nplist; ++jp)
                    FIM.at(jp + ip * nplist) *=
                        pcoefficient.at(ip)*pcoefficient.at(jp);

#define chainRule(QUANT, IND1, N1T, N1, IND2, N2)                              \
    if (!s##QUANT.empty())                                                     \
        for (int IND1 = 0; (IND1) < (N1T); ++(IND1))                           \
            for (int ip = 0; ip < nplist; ++ip)                                \
                for (int IND2 = 0; (IND2) < (N2); ++(IND2)) {                  \
                    s##QUANT.at(((IND2)*nplist + ip) * (N1) + (IND1)) *=       \
                        pcoefficient.at(ip);                                   \
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

#define s2ChainRule(QUANT, IND1, N1T, N1, IND2, N2)                            \
    if (!s##QUANT.empty())                                                     \
        for (int ip = 0; ip < nplist; ++ip)                                    \
            for (int iJ = 1; iJ < nJ; ++iJ)                                    \
                for (int IND1 = 0; IND1 < N1T; ++IND1)                         \
                    for (int IND2 = 0; IND2 < N2; ++IND2) {                    \
                        s##QUANT.at((IND2 * nplist + ip) * N1 + IND1 +         \
                                    iJ * N1T) *=                               \
                            pcoefficient.at(ip) * augcoefficient[iJ - 1];      \
                        if (model.plist(ip) == iJ - 1)                         \
                            s##QUANT.at((IND2 * nplist + ip) * N1 + IND1 +     \
                                        iJ * N1T) +=                           \
                                s##QUANT.at((IND2 * nplist + ip) * N1 +        \
                                            IND1) *                            \
                                coefficient[ip];                               \
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

#define s2vecChainRule(QUANT, IND1, N1T, N1, IND2, N2)                         \
    if (!s##QUANT.empty())                                                     \
        for (int ip = 0; ip < nplist; ++ip)                                    \
            for (int IND1 = 0; IND1 < N1T; ++IND1)                             \
                for (int IND2 = 0; IND2 < N2; ++IND2) {                        \
                    s##QUANT.at((IND2 * nplist + ip) * N1 + IND1 + N1T) *=     \
                        pcoefficient.at(ip);                                   \
                    s##QUANT.at((IND2 * nplist + ip) * N1 + IND1 + N1T) +=     \
                        model.k()[nk - nplist + ip] *                          \
                        s##QUANT.at((IND2 * nplist + ip) * N1 + IND1) /        \
                        unscaledParameters[model.plist(ip)];                   \
                }

        s2vecChainRule(x, ix, nxtrue, nx, it, nt);
        s2vecChainRule(y, iy, nytrue, ny, it, nt);
        s2vecChainRule(sigmay, iy, nytrue, ny, it, nt);
        s2vecChainRule(z, iz, nztrue, nz, ie, nmaxevent);
        s2vecChainRule(sigmaz, iz, nztrue, nz, ie, nmaxevent);
        s2vecChainRule(rz, iz, nztrue, nz, ie, nmaxevent);
    }
}

void ReturnData::initializeObjectiveFunction() {
    if (rdata_reporting == RDataReporting::likelihood ||
        rdata_reporting == RDataReporting::full) {
        llh = 0.0;
        std::fill(sllh.begin(), sllh.end(), 0.0);
        std::fill(s2llh.begin(), s2llh.end(), 0.0);
    }
    if (rdata_reporting == RDataReporting::residuals ||
        rdata_reporting == RDataReporting::full)
        chi2 = 0.0;
}

static realtype fres(realtype y, realtype my, realtype sigma_y) {
    return (y - my) / sigma_y;
}

static realtype fsres(realtype y, realtype sy, realtype my,
                      realtype sigma_y, realtype ssigma_y) {
    double r = fres(sy, 0.0, sigma_y);
    if (ssigma_y > 0)
        r += fres(y, my, sigma_y * sigma_y / ssigma_y);
    return r;
}

void ReturnData::fres(const int it, Model &model, const ExpData &edata) {
    if (res.empty())
        return;

    std::vector<realtype> y_it(ny, 0.0);
    model.getObservable(y_it, ts[it], x_solver_);

    std::vector<realtype> sigmay_it(ny, 0.0);
    model.getObservableSigma(sigmay_it, it, &edata);

    auto observedData = edata.getObservedDataPtr(it);
    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt = iy + it * edata.nytrue();
        if (!edata.isSetObservedData(it, iy))
            continue;
        res.at(iyt) = amici::fres(y_it.at(iy), observedData[iy],
                                  sigmay_it.at(iy));
    }
}

void ReturnData::fchi2(const int it) {
    if (res.empty() || isNaN(chi2))
        return;

    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * nytrue;
        chi2 += pow(res.at(iyt_true), 2);
    }
}

void ReturnData::fsres(const int it, Model &model, const ExpData &edata) {
    if (sres.empty())
        return;

    std::vector<realtype> y_it(ny, 0.0);
    model.getObservable(y_it, ts[it], x_solver_);
    std::vector<realtype> sy_it(ny * nplist, 0.0);
    model.getObservableSensitivity(sy_it, ts[it], x_solver_, sx_solver_);

    std::vector<realtype> sigmay_it(ny, 0.0);
    model.getObservableSigma(sigmay_it, it, &edata);
    std::vector<realtype> ssigmay_it(ny * nplist, 0.0);
    model.getObservableSigmaSensitivity(ssigmay_it, it, &edata);

    auto observedData = edata.getObservedDataPtr(it);
    for (int iy = 0; iy < nytrue; ++iy) {
        if (!edata.isSetObservedData(it, iy))
            continue;
        for (int ip = 0; ip < nplist; ++ip) {
            int idx = (iy + it * edata.nytrue()) * nplist + ip;
            sres.at(idx) = amici::fsres(y_it.at(iy), sy_it.at(iy + ny * ip),
                                        observedData[iy], sigmay_it.at(iy),
                                        ssigmay_it.at(iy + ny * ip));
        }
    }
}

void ReturnData::fFIM(int it, Model &model, const ExpData &edata) {
    if (FIM.empty())
        return;

    std::vector<realtype> y_it(ny, 0.0);
    model.getObservable(y_it, ts[it], x_solver_);
    std::vector<realtype> sy_it(ny * nplist, 0.0);
    model.getObservableSensitivity(sy_it, ts[it], x_solver_, sx_solver_);

    std::vector<realtype> sigmay_it(ny, 0.0);
    model.getObservableSigma(sigmay_it, it, &edata);
    std::vector<realtype> ssigmay_it(ny * nplist, 0.0);
    model.getObservableSigmaSensitivity(ssigmay_it, it, &edata);

    auto observedData = edata.getObservedDataPtr(it);
    for (int iy = 0; iy < nytrue; ++iy) {
        if (!edata.isSetObservedData(it, iy))
            continue;
        for (int ip = 0; ip < nplist; ++ip) {
            for (int jp = 0; jp < nplist; ++jp) {
                FIM.at(ip + nplist * jp) +=
                    amici::fsres(y_it.at(iy), sy_it.at(iy + ny * ip),
                                 observedData[iy], sigmay_it.at(iy),
                                 ssigmay_it.at(iy + ny * ip))
                    *
                    amici::fsres(y_it.at(iy), sy_it.at(iy + ny * jp),
                                 observedData[iy], sigmay_it.at(iy),
                                 ssigmay_it.at(iy + ny * jp));
            }
        }
    }
}

ModelContext::ModelContext(Model *model)
    : model_(model), original_state_(model->getModelState()) {}

ModelContext::~ModelContext() { restore(); }

void ModelContext::restore() { model_->setModelState(original_state_); }

} // namespace amici
