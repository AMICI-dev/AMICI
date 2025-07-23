#include "amici/rdata.h"

#include "amici/backwardproblem.h"
#include "amici/edata.h"
#include "amici/forwardproblem.h"
#include "amici/misc.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/steadystateproblem.h"
#include "amici/symbolic_functions.h"
#include "amici/vector.h"

#include <algorithm>
#include <cmath>

namespace amici {

ReturnData::ReturnData(Solver const& solver, Model const& model)
    : ReturnData(
          model.getTimepoints(),
          ModelDimensions(static_cast<ModelDimensions const&>(model)),
          model.nMaxEvent(), solver.getNewtonMaxSteps(),
          model.getParameterList(), model.getParameterScale(), model.o2mode,
          solver.getSensitivityOrder(), solver.getSensitivityMethod(),
          solver.getReturnDataReportingMode(), model.hasQuadraticLLH(),
          model.getAddSigmaResiduals(), model.getMinimumSigmaResiduals()
      ) {}

ReturnData::ReturnData(
    std::vector<realtype> ts_, ModelDimensions const& model_dimensions_,
    int nmaxevent_, int newton_maxsteps_, std::vector<int> plist_,
    std::vector<ParameterScaling> pscale_, SecondOrderMode o2mode_,
    SensitivityOrder sensi_, SensitivityMethod sensi_meth_,
    RDataReporting rdrm_, bool quadratic_llh_, bool sigma_res_,
    realtype sigma_offset_
)
    : ModelDimensions(model_dimensions_)
    , ts(std::move(ts_))
    , nx(nx_rdata)
    , nxtrue(nxtrue_rdata)
    , nplist(plist_.size())
    , nmaxevent(nmaxevent_)
    , nt(ts.size())
    , newton_maxsteps(newton_maxsteps_)
    , pscale(std::move(pscale_))
    , o2mode(o2mode_)
    , sensi(sensi_)
    , sensi_meth(sensi_meth_)
    , rdata_reporting(rdrm_)
    , sigma_res(sigma_res_)
    , plist(plist_)
    , sigma_offset(sigma_offset_)
    , nroots_(ne) {
    switch (rdata_reporting) {
    case RDataReporting::full:
        initializeFullReporting(quadratic_llh_);
        break;

    case RDataReporting::residuals:
        initializeResidualReporting(quadratic_llh_);
        break;

    case RDataReporting::likelihood:
        initializeLikelihoodReporting(quadratic_llh_);
        break;

    case RDataReporting::observables_likelihood:
        initializeObservablesLikelihoodReporting(quadratic_llh_);
        break;
    }
}

void ReturnData::initializeLikelihoodReporting(bool enable_fim) {
    llh = getNaN();
    chi2 = getNaN();
    if (sensi >= SensitivityOrder::first
        && sensi_meth != SensitivityMethod::none) {
        sllh.resize(nplist, getNaN());
        if (sensi >= SensitivityOrder::second)
            s2llh.resize(nplist * (nJ - 1), getNaN());

        if ((sensi_meth == SensitivityMethod::forward
             || sensi >= SensitivityOrder::second)
            && enable_fim)
            FIM.resize(nplist * nplist, 0.0);
    }
}

void ReturnData::initializeObservablesLikelihoodReporting(bool enable_fim) {
    initializeLikelihoodReporting(enable_fim);

    y.resize(nt * ny, 0.0);
    sigmay.resize(nt * ny, 0.0);

    if ((sensi_meth == SensitivityMethod::forward
         && sensi >= SensitivityOrder::first)
        || sensi >= SensitivityOrder::second) {

        sy.resize(nt * ny * nplist, 0.0);
        ssigmay.resize(nt * ny * nplist, 0.0);
    }
}

void ReturnData::initializeResidualReporting(bool enable_res) {
    y.resize(nt * ny, 0.0);
    sigmay.resize(nt * ny, 0.0);
    if (enable_res)
        res.resize((sigma_res ? 2 : 1) * nt * nytrue, 0.0);

    if ((sensi_meth == SensitivityMethod::forward
         && sensi >= SensitivityOrder::first)
        || sensi >= SensitivityOrder::second) {

        sy.resize(nt * ny * nplist, 0.0);
        ssigmay.resize(nt * ny * nplist, 0.0);
        if (enable_res)
            sres.resize((sigma_res ? 2 : 1) * nt * nytrue * nplist, 0.0);
    }
}

void ReturnData::initializeFullReporting(bool quadratic_llh) {

    initializeLikelihoodReporting(quadratic_llh);
    initializeResidualReporting(quadratic_llh);

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

        if (sensi_meth == SensitivityMethod::adjoint
            && sensi >= SensitivityOrder::first) {
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

        if (sensi_meth == SensitivityMethod::forward
            || sensi >= SensitivityOrder::second) {
            // for second order we can fill in from the augmented states
            sx.resize(nt * nx * nplist, 0.0);
            sz.resize(nmaxevent * nz * nplist, 0.0);
            srz.resize(nmaxevent * nz * nplist, 0.0);
        }

        ssigmaz.resize(nmaxevent * nz * nplist, 0.0);
        if (sensi >= SensitivityOrder::second
            && sensi_meth == SensitivityMethod::forward)
            s2rz.resize(nmaxevent * nztrue * nplist * nplist, 0.0);
    }
}

void ReturnData::processSimulationObjects(
    ForwardProblem const* fwd, BackwardProblem const* bwd, Model& model,
    Solver const& solver, ExpData const* edata
) {
    ModelContext mc(&model);

    processSolver(solver);

    SteadystateProblem const* preeq = nullptr;
    SteadystateProblem const* posteq = nullptr;
    SteadyStateBackwardProblem const* preeq_bwd = nullptr;
    SteadyStateBackwardProblem const* posteq_bwd = nullptr;
    if (fwd) {
        preeq = fwd->getPreequilibrationProblem();
        posteq = fwd->getPostequilibrationProblem();
        if (bwd) {
            preeq_bwd = bwd->getPreequilibrationBwdProblem();
            posteq_bwd = bwd->getPostequilibrationBwdProblem();
        }
    }
    if (preeq)
        processPreEquilibration(*preeq, preeq_bwd, model);

    if (fwd)
        processForwardProblem(*fwd, model, edata);
    else
        invalidate(0);

    if (posteq)
        processPostEquilibration(*posteq, posteq_bwd, model, edata);

    if (edata && !edata->fixedParametersPreequilibration.empty() && fwd
        && !fwd->was_preequilibrated() && preeq) {
        // failure during preequilibration
        storeJacobianAndDerivativeInReturnData(*preeq, model);
    } else if (fwd && !posteq)
        storeJacobianAndDerivativeInReturnData(*fwd, model);
    else if (posteq)
        storeJacobianAndDerivativeInReturnData(*posteq, model);

    if (fwd && bwd)
        processBackwardProblem(*fwd, *bwd, preeq, preeq_bwd, model);
    else if (solver.computingASA())
        invalidateSLLH();

    applyChainRuleFactorToSimulationResults(model);
}

void ReturnData::processPreEquilibration(
    SteadystateProblem const& preeq,
    SteadyStateBackwardProblem const* preeq_bwd, Model& model
) {
    auto const simulation_state = preeq.getFinalSimulationState();
    model.setModelState(simulation_state.mod);

    if (!x_ss.empty()) {
        model.fx_rdata(x_ss, simulation_state.sol.x);
    }
    if (!sx_ss.empty() && sensi >= SensitivityOrder::first) {
        model.fsx_rdata(sx_ss, simulation_state.sol.sx, simulation_state.sol.x);
    }
    /* Get cpu time for pre-equilibration in milliseconds */
    preeq_cpu_time = preeq.getCPUTime();
    if (preeq_bwd) {
        preeq_cpu_timeB = preeq_bwd->getCPUTimeB();
        preeq_numstepsB = preeq_bwd->getNumStepsB();
    }
    preeq_wrms = preeq.getResidualNorm();
    preeq_status = preeq.getSteadyStateStatus();
    if (preeq_status[1] == SteadyStateStatus::success)
        preeq_t = preeq.getSteadyStateTime();
    if (!preeq_numsteps.empty())
        writeSlice(preeq.getNumSteps(), preeq_numsteps);
}

void ReturnData::processPostEquilibration(
    SteadystateProblem const& posteq,
    SteadyStateBackwardProblem const* posteq_bwd, Model& model,
    ExpData const* edata
) {
    for (int it = 0; it < nt; it++) {
        auto const t = model.getTimepoint(it);
        if (std::isinf(t)) {
            auto const& simulation_state = posteq.getFinalSimulationState();
            model.setModelState(simulation_state.mod);
            getDataOutput(it, model, simulation_state.sol, edata);
        }
    }
    /* Get cpu time for Newton solve in milliseconds */
    posteq_cpu_time = posteq.getCPUTime();
    if (posteq_bwd) {
        posteq_cpu_timeB = posteq_bwd->getCPUTimeB();
        posteq_numstepsB = posteq_bwd->getNumStepsB();
    }
    posteq_wrms = posteq.getResidualNorm();
    posteq_status = posteq.getSteadyStateStatus();
    if (posteq_status[1] == SteadyStateStatus::success)
        posteq_t = posteq.getSteadyStateTime();
    if (!posteq_numsteps.empty())
        writeSlice(posteq.getNumSteps(), posteq_numsteps);
}

void ReturnData::processForwardProblem(
    ForwardProblem const& fwd, Model& model, ExpData const* edata
) {
    if (edata)
        initializeObjectiveFunction(model.hasQuadraticLLH());

    auto const& initialState = fwd.getInitialSimulationState();
    if (initialState.sol.x.getLength() == 0 && model.nx_solver > 0)
        return; // if x wasn't set forward problem failed during initialization

    model.setModelState(initialState.mod);

    if (!x0.empty()) {
        model.fx_rdata(x0, initialState.sol.x);
    }

    if (!sx0.empty()) {
        model.fsx_rdata(sx0, initialState.sol.sx, initialState.sol.x);
    }

    // process timepoint data
    realtype tf = fwd.getFinalTime();
    for (int it = 0; it < model.nt(); it++) {
        if (model.getTimepoint(it) <= tf) {
            auto const simulation_state = fwd.getSimulationStateTimepoint(it);
            model.setModelState(simulation_state.mod);
            getDataOutput(it, model, simulation_state.sol, edata);
        } else {
            // check for integration failure but consider postequilibration
            if (!std::isinf(model.getTimepoint(it)))
                invalidate(it);
        }
    }

    // process event data
    if (nz > 0) {
        auto const& discontinuities = fwd.getDiscontinuities();
        Expects(
            static_cast<int>(discontinuities.size())
            == fwd.getEventCounter() + 1
        );
        for (int iroot = 0; iroot <= fwd.getEventCounter(); iroot++) {
            auto const simulation_state = fwd.getSimulationStateEvent(iroot);
            model.setModelState(simulation_state.mod);
            getEventOutput(
                discontinuities.at(iroot).root_info, model,
                simulation_state.sol, edata
            );
        }
    }
}

void ReturnData::getDataOutput(
    int it, Model& model, SolutionState const& sol, ExpData const* edata
) {
    if (!x.empty()) {
        model.fx_rdata(slice(x, it, nx), sol.x);
    }
    if (!w.empty())
        model.getExpression(slice(w, it, nw), ts[it], sol.x);
    if (!y.empty())
        model.getObservable(slice(y, it, ny), ts[it], sol.x);
    if (!sigmay.empty())
        model.getObservableSigma(slice(sigmay, it, ny), it, edata);

    if (edata) {
        if (!isNaN(llh))
            model.addObservableObjective(llh, it, sol.x, *edata);
        fres(it, model, sol, *edata);
        fchi2(it, *edata);
    }

    if (sensi >= SensitivityOrder::first && nplist > 0) {

        if (sensi_meth == SensitivityMethod::forward) {
            getDataSensisFSA(it, model, sol, edata);
        } else if (edata && !sllh.empty()) {
            model.addPartialObservableObjectiveSensitivity(
                sllh, s2llh, it, sol.x, *edata
            );
        }

        if (!ssigmay.empty())
            model.getObservableSigmaSensitivity(
                slice(ssigmay, it, nplist * ny), slice(sy, it, nplist * ny), it,
                edata
            );
    }
}

void ReturnData::getDataSensisFSA(
    int it, Model& model, SolutionState const& sol, ExpData const* edata
) {
    if (!sx.empty()) {
        model.fsx_rdata(
            gsl::span<realtype>(&sx.at(it * nplist * nx), nplist * nx), sol.sx,
            sol.x
        );
    }

    if (!sy.empty()) {
        model.getObservableSensitivity(
            slice(sy, it, nplist * ny), ts[it], sol.x, sol.sx
        );
    }

    if (edata) {
        if (!sllh.empty())
            model.addObservableObjectiveSensitivity(
                sllh, s2llh, it, sol.x, sol.sx, *edata
            );
        fsres(it, model, sol, *edata);
        fFIM(it, model, sol, *edata);
    }
}

void ReturnData::getEventOutput(
    std::vector<int> const& rootidx, Model& model, SolutionState const& sol,
    ExpData const* edata
) {

    for (int ie = 0; ie < ne; ie++) {
        if (rootidx.at(ie) != 1 || nroots_.at(ie) >= nmaxevent)
            continue;

        /* get event output */
        if (!z.empty())
            model.getEvent(slice(z, nroots_.at(ie), nz), ie, sol.t, sol.x);
        /* if called from fillEvent at last timepoint,
         then also get the root function value */
        if (sol.t == model.getTimepoint(nt - 1))
            if (!rz.empty())
                model.getEventRegularization(
                    slice(rz, nroots_.at(ie), nz), ie, sol.t, sol.x
                );

        if (edata) {
            if (!sigmaz.empty())
                model.getEventSigma(
                    slice(sigmaz, nroots_.at(ie), nz), ie, nroots_.at(ie),
                    sol.t, edata
                );
            if (!isNaN(llh))
                model.addEventObjective(
                    llh, ie, nroots_.at(ie), sol.t, sol.x, *edata
                );

            /* if called from fillEvent at last timepoint,
               add regularization based on rz */
            if (sol.t == model.getTimepoint(nt - 1) && !isNaN(llh)) {
                model.addEventObjectiveRegularization(
                    llh, ie, nroots_.at(ie), sol.t, sol.x, *edata
                );
            }
        }

        if (sensi >= SensitivityOrder::first) {
            if (sensi_meth == SensitivityMethod::forward) {
                getEventSensisFSA(ie, model, sol, edata);
            } else if (edata && !sllh.empty()) {
                model.addPartialEventObjectiveSensitivity(
                    sllh, s2llh, ie, nroots_.at(ie), sol.t, sol.x, *edata
                );
            }
        }
        nroots_.at(ie)++;
    }
}

void ReturnData::getEventSensisFSA(
    int ie, Model& model, SolutionState const& sol, ExpData const* edata
) {
    if (sol.t == model.getTimepoint(nt - 1)) {
        // call from fillEvent at last timepoint
        if (!sz.empty())
            model.getUnobservedEventSensitivity(
                slice(sz, nroots_.at(ie), nz * nplist), ie
            );
        if (!srz.empty())
            model.getEventRegularizationSensitivity(
                slice(srz, nroots_.at(ie), nz * nplist), ie, sol.t, sol.x,
                sol.sx
            );
    } else if (!sz.empty()) {
        model.getEventSensitivity(
            slice(sz, nroots_.at(ie), nz * nplist), ie, sol.t, sol.x, sol.sx
        );
    }

    if (edata && !sllh.empty()) {
        model.addEventObjectiveSensitivity(
            sllh, s2llh, ie, nroots_.at(ie), sol.t, sol.x, sol.sx, *edata
        );
    }
}

void ReturnData::processBackwardProblem(
    ForwardProblem const& fwd, BackwardProblem const& bwd,
    SteadystateProblem const* preeq,
    SteadyStateBackwardProblem const* preeq_bwd, Model& model
) {
    if (sllh.empty())
        return;

    auto simulation_state = fwd.getInitialSimulationState();
    model.setModelState(simulation_state.mod);

    std::vector<realtype> llhS0(model.nJ * model.nplist(), 0.0);
    // Adjoint quadratures before pre-equilibration
    auto xQB = bwd.getAdjointQuadraturePrePreeq();

    if (preeq_bwd && preeq_bwd->hasQuadrature()) {
        // Pre-equilibration with ASA steady-state shortcuts

        // If pre-equilibration is run in adjoint mode, the scalar product of
        // sx0 with its adjoint counterpart (see handleSx0Forward()) is not
        // necessary:
        // the actual simulation is "extended" by the pre-equilibration time.
        // At initialization (at t=-inf), the adjoint state is in steady state
        // (= 0) and so is the scalar product.
        // Instead of the scalar product, the quadratures xQB from
        // pre-equilibration contribute to the gradient
        // (see example notebook on equilibration for further documentation).

        // Adjoint quadratures after pre-equilibration
        auto const& xQB_post_preeq = preeq_bwd->getAdjointQuadrature();
        for (int ip = 0; ip < model.nplist(); ++ip)
            xQB[ip] += xQB_post_preeq.at(ip);

        handleSx0Backward(
            model, preeq->getStateSensitivity(), bwd.getAdjointState(), llhS0
        );
    } else if (preeq
               && preeq->getSteadyStateStatus()[1]
                      != SteadyStateStatus::not_run) {
        // Pre-equilibration with ASA backward simulation

        // Adjoint quadratures after pre-equilibration
        xQB = bwd.getAdjointQuadrature();

        handleSx0Backward(
            model, preeq->getStateSensitivity(), bwd.getAdjointState(), llhS0
        );

    } else {
        // No pre-equilibration, or pre-equilibration with FSA
        handleSx0Forward(
            model, simulation_state.sol, llhS0, bwd.getAdjointStatePrePreeq()
        );
    }

    for (int iJ = 0; iJ < model.nJ; iJ++) {
        for (int ip = 0; ip < model.nplist(); ip++) {
            if (iJ == 0) {
                sllh.at(ip) -= llhS0[ip] + xQB[ip * model.nJ];
            } else {
                s2llh.at(iJ - 1 + ip * (model.nJ - 1))
                    -= llhS0[ip + iJ * model.nplist()]
                       + xQB[iJ + ip * model.nJ];
            }
        }
    }
}

void ReturnData::handleSx0Backward(
    Model const& model, AmiVectorArray const& sx0, AmiVector const& xB,
    std::vector<realtype>& llhS0
) const {

    // Add the contribution for sx0 from preequilibration. If backward
    // preequilibration was done by simulation due to a singular Jacobian,
    // xB is not necessarily 0 and we may get a non-zero contribution here.
    for (int ip = 0; ip < model.nplist(); ++ip) {
        llhS0[ip] = 0.0;
        for (int ix = 0; ix < model.nxtrue_solver; ++ix) {
            llhS0[ip] += xB.at(ix) * sx0.at(ix, ip);
        }
    }
}

void ReturnData::handleSx0Forward(
    Model const& model, SolutionState const& sol, std::vector<realtype>& llhS0,
    AmiVector const& xB
) const {
    /* If preequilibration is run in forward mode or is not needed, then adjoint
       sensitivity analysis still needs the state sensitivities at t=0 (sx0),
       to compute the gradient. For each parameter, the scalar product of sx0
       with its adjoint counterpart contributes to the gradient. */
    for (int iJ = 0; iJ < model.nJ; iJ++) {
        if (iJ == 0) {
            for (int ip = 0; ip < model.nplist(); ++ip) {
                llhS0[ip] = 0.0;
                for (int ix = 0; ix < model.nxtrue_solver; ++ix) {
                    llhS0[ip] += xB[ix] * sol.sx.at(ix, ip);
                }
            }
        } else {
            for (int ip = 0; ip < model.nplist(); ++ip) {
                llhS0[ip + iJ * model.nplist()] = 0.0;
                for (int ix = 0; ix < model.nxtrue_solver; ++ix) {
                    llhS0[ip + iJ * model.nplist()]
                        += xB[ix + iJ * model.nxtrue_solver] * sol.sx.at(ix, ip)
                           + xB[ix]
                                 * sol.sx.at(ix + iJ * model.nxtrue_solver, ip);
                }
            }
        }
    }
}

void ReturnData::processSolver(Solver const& solver) {
    using std::ranges::copy;

    cpu_time = solver.getCpuTime();

    if (!numsteps.empty()) {
        Expects(numsteps.size() >= solver.getNumSteps().size());
        // copy instead of assignment to ensure length `nt`
        // (vector from solver may be shorter in case of integration errors)
        copy(solver.getNumSteps(), numsteps.begin());
    }

    if (!numrhsevals.empty()) {
        copy(solver.getNumRhsEvals(), numrhsevals.begin());
    }

    if (!numerrtestfails.empty()) {
        copy(solver.getNumErrTestFails(), numerrtestfails.begin());
    }

    if (!numnonlinsolvconvfails.empty()) {
        copy(
            solver.getNumNonlinSolvConvFails(), numnonlinsolvconvfails.begin()
        );
    }

    if (!order.empty()) {
        copy(solver.getLastOrder(), order.begin());
    }

    cpu_timeB = solver.getCpuTimeB();

    if (!numstepsB.empty()) {
        copy(solver.getNumStepsB(), numstepsB.begin());
    }

    if (!numrhsevalsB.empty()) {
        copy(solver.getNumRhsEvalsB(), numrhsevalsB.begin());
    }

    if (!numerrtestfailsB.empty()) {
        copy(solver.getNumErrTestFailsB(), numerrtestfailsB.begin());
    }

    if (!numnonlinsolvconvfailsB.empty()) {
        copy(
            solver.getNumNonlinSolvConvFailsB(), numnonlinsolvconvfailsB.begin()
        );
    }
}

void ReturnData::invalidate(int const it_start) {
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
        std::ranges::fill(sllh, getNaN());
        std::ranges::fill(s2llh, getNaN());
    }
}

void ReturnData::applyChainRuleFactorToSimulationResults(Model const& model) {
    // chain-rule factor: multiplier for am_p
    std::vector<realtype> coefficient(nplist, 1.0);
    // Only apply `pcoefficient` to non-zero sensitivities. This allows
    // having unused parameters that are set to NAN, and still having finite
    // sensitivities. Otherwise, this would lead to NAN-sensitivities w.r.t.
    // to such parameters if they are log-scaled.
    std::vector<realtype> pcoefficient(nplist, 1.0);

    std::vector<realtype> unscaledParameters = model.getParameters();
    unscaleParameters(
        unscaledParameters, model.getParameterScale(), unscaledParameters
    );

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
            pcoefficient.at(ip)
                = unscaledParameters.at(model.plist(ip)) * log(10);
            break;
        case ParameterScaling::ln:
            pcoefficient.at(ip) = unscaledParameters.at(model.plist(ip));
            break;
        case ParameterScaling::none:
            break;
        }
    }

    if (sensi >= SensitivityOrder::first) {
        // recover first order sensitivities from states for adjoint sensitivity
        // analysis
        if (sensi == SensitivityOrder::second && o2mode == SecondOrderMode::full
            && sensi_meth == SensitivityMethod::adjoint) {
            if (!sx.empty() && !x.empty())
                for (int ip = 0; ip < nplist; ++ip)
                    for (int ix = 0; ix < nxtrue; ++ix)
                        for (int it = 0; it < nt; ++it)
                            sx.at(ix + nxtrue * (ip + it * nplist))
                                = x.at(it * nx + nxtrue + ip * nxtrue + ix);

            if (!sy.empty() && !y.empty())
                for (int ip = 0; ip < nplist; ++ip)
                    for (int iy = 0; iy < nytrue; ++iy)
                        for (int it = 0; it < nt; ++it)
                            sy.at(iy + nytrue * (ip + it * nplist))
                                = y.at(it * ny + nytrue + ip * nytrue + iy);

            if (!sz.empty() && !z.empty())
                for (int ip = 0; ip < nplist; ++ip)
                    for (int iz = 0; iz < nztrue; ++iz)
                        for (int it = 0; it < nt; ++it)
                            sz.at(iz + nztrue * (ip + it * nplist))
                                = z.at(it * nz + nztrue + ip * nztrue + iz);
        }

        if (!sllh.empty())
            for (int ip = 0; ip < nplist; ++ip)
                if (auto&& val = sllh.at(ip); val != 0.0)
                    val *= pcoefficient.at(ip);

        if (!sres.empty())
            for (int ires = 0; ires < gsl::narrow<int>(res.size()); ++ires)
                for (int ip = 0; ip < nplist; ++ip)
                    if (auto&& val = sres.at(ires * nplist + ip); val != 0.0)
                        val *= pcoefficient.at(ip);

        if (!FIM.empty())
            for (int ip = 0; ip < nplist; ++ip)
                for (int jp = 0; jp < nplist; ++jp)
                    if (auto&& val = FIM.at(jp + ip * nplist); val != 0.0)
                        val *= pcoefficient.at(ip) * pcoefficient.at(jp);

        // apply chain rule to sensitivities
        auto chain_rule = [&](auto& sens, int n1, int stride1, int n2) {
            if (sens.empty()) {
                return;
            }
            using index_type =
                typename std::remove_reference_t<decltype(sens)>::size_type;
            Expects(
                sens.size() == gsl::narrow<index_type>(n2 * nplist * stride1)
            );
            Expects(n1 <= stride1);
            Expects(pcoefficient.size() == gsl::narrow<index_type>(nplist));
            for (index_type i1 = 0; i1 < gsl::narrow<index_type>(n1); ++i1)
                for (index_type ip = 0; ip < gsl::narrow<index_type>(nplist);
                     ++ip)
                    for (index_type i2 = 0; i2 < gsl::narrow<index_type>(n2);
                         ++i2)
                        if (auto&& val
                            = sens.at((i2 * nplist + ip) * stride1 + i1);
                            val != 0.0)
                            val *= pcoefficient[ip];
        };

        chain_rule(sx, nxtrue, nx, nt);
        chain_rule(sy, nytrue, ny, nt);
        chain_rule(ssigmay, nytrue, ny, nt);
        chain_rule(sz, nztrue, nz, nmaxevent);
        chain_rule(ssigmaz, nztrue, nz, nmaxevent);
        chain_rule(srz, nztrue, nz, nmaxevent);
        chain_rule(sx0, nxtrue, nx, 1);
        chain_rule(sx_ss, nxtrue, nx, 1);
    }

    if (o2mode == SecondOrderMode::full) {
        if (!s2llh.empty() && !sllh.empty()) {
            for (int ip = 0; ip < nplist; ++ip) {
                for (int iJ = 1; iJ < nJ; ++iJ) {
                    auto&& val = s2llh[ip * nplist + (iJ - 1)];
                    if (val != 0.0)
                        val *= pcoefficient.at(ip) * augcoefficient[iJ - 1];
                    if (model.plist(ip) == iJ - 1)
                        val += sllh.at(ip) * coefficient.at(ip);
                }
            }
        }

        auto chain_rule = [&](auto& sens, int n1, int stride1, int n2) {
            if (sens.empty())
                return;

            using index_type =
                typename std::remove_reference_t<decltype(sens)>::size_type;
            Expects(
                sens.size() == gsl::narrow<index_type>(n2 * nplist * stride1)
            );
            Expects(n1 <= stride1);
            Expects(pcoefficient.size() == gsl::narrow<index_type>(nplist));
            Expects(coefficient.size() == gsl::narrow<index_type>(nplist));

            for (int ip = 0; ip < nplist; ++ip)
                for (int iJ = 1; iJ < nJ; ++iJ)
                    for (int i1 = 0; i1 < n1; ++i1)
                        for (int i2 = 0; i2 < n2; ++i2) {
                            auto idx
                                = (i2 * nplist + ip) * stride1 + i1 + iJ * n1;
                            auto&& val = sens.at(idx);
                            if (val != 0.0)
                                val *= pcoefficient.at(ip)
                                       * augcoefficient[iJ - 1];
                            if (model.plist(ip) == iJ - 1)
                                val += sens.at(
                                           (i2 * nplist + ip) * stride1 + i1
                                       )
                                       * coefficient[ip];
                        }
        };

        chain_rule(sx, nxtrue, nx, nt);
        chain_rule(sy, nytrue, ny, nt);
        chain_rule(ssigmay, nytrue, ny, nt);
        chain_rule(sz, nztrue, nz, nmaxevent);
        chain_rule(ssigmaz, nztrue, nz, nmaxevent);
        chain_rule(srz, nztrue, nz, nmaxevent);
    } else if (o2mode == SecondOrderMode::directional) {
        for (int ip = 0; ip < nplist; ++ip) {
            auto&& val = sllh.at(ip);
            if (val != 0.0)
                val *= pcoefficient.at(ip);
            val += model.k()[nk - nplist + ip] * sllh.at(ip)
                   / unscaledParameters[model.plist(ip)];
        }

        auto chain_rule = [&](auto& sens, int n1, int stride1, int n2) {
            if (sens.empty())
                return;

            using index_type =
                typename std::remove_reference_t<decltype(sens)>::size_type;
            Expects(
                sens.size() == gsl::narrow<index_type>(n2 * nplist * stride1)
            );
            Expects(n1 <= stride1);
            Expects(pcoefficient.size() == gsl::narrow<index_type>(nplist));

            for (int ip = 0; ip < nplist; ++ip)
                for (int i1 = 0; i1 < n1; ++i1)
                    for (int i2 = 0; i2 < n2; ++i2) {
                        auto idx = (i2 * nplist + ip) * stride1 + i1 + n1;
                        auto&& val = sens.at(idx);
                        if (val != 0.0)
                            val *= pcoefficient.at(ip);
                        val += model.k()[nk - nplist + ip]
                               * sens.at((i2 * nplist + ip) * stride1 + i1)
                               / unscaledParameters[model.plist(ip)];
                    }
        };

        chain_rule(sx, nxtrue, nx, nt);
        chain_rule(sy, nytrue, ny, nt);
        chain_rule(ssigmay, nytrue, ny, nt);
        chain_rule(sz, nztrue, nz, nmaxevent);
        chain_rule(ssigmaz, nztrue, nz, nmaxevent);
        chain_rule(srz, nztrue, nz, nmaxevent);
    }
}

void ReturnData::initializeObjectiveFunction(bool enable_chi2) {
    if (rdata_reporting == RDataReporting::likelihood
        || rdata_reporting == RDataReporting::observables_likelihood
        || rdata_reporting == RDataReporting::full) {
        llh = 0.0;
        std::ranges::fill(sllh, 0.0);
        std::ranges::fill(s2llh, 0.0);
    }
    if ((rdata_reporting == RDataReporting::residuals
         || rdata_reporting == RDataReporting::full)
        && enable_chi2)
        chi2 = 0.0;
}

static realtype
fres(realtype y, realtype my, realtype sigma_y, ObservableScaling scale) {
    switch (scale) {
    case ObservableScaling::lin:
        return (my - y) / sigma_y;
    case ObservableScaling::log:
        return (std::log(my) - std::log(y)) / sigma_y;
    case ObservableScaling::log10:
        return (std::log10(my) - std::log10(y)) / sigma_y;
    default:
        throw std::invalid_argument("only lin, log, log10 allowed.");
    }
}

static realtype fres_error(realtype sigma_y, realtype sigma_offset) {
    return sqrt(2 * log(sigma_y) + sigma_offset);
}

void ReturnData::fres(
    int const it, Model& model, SolutionState const& sol, ExpData const& edata
) {
    if (res.empty())
        return;

    std::vector<realtype> y_it(ny, 0.0);
    model.getObservable(y_it, ts[it], sol.x);

    std::vector<realtype> sigmay_it(ny, 0.0);
    model.getObservableSigma(sigmay_it, it, &edata);

    auto observedData = edata.getObservedDataPtr(it);
    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt = iy + it * edata.nytrue();
        if (!edata.isSetObservedData(it, iy))
            continue;

        res.at(iyt) = amici::fres(
            y_it.at(iy), observedData[iy], sigmay_it.at(iy),
            model.getObservableScaling(iy)
        );

        if (sigma_res)
            res.at(iyt + nt * nytrue)
                = fres_error(sigmay_it.at(iy), sigma_offset);
    }
}

void ReturnData::fchi2(int const it, ExpData const& edata) {
    if (res.empty() || isNaN(chi2))
        return;

    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * nytrue;
        chi2 += pow(res.at(iyt_true), 2);
        if (sigma_res && edata.isSetObservedData(it, iy))
            chi2 += pow(res.at(iyt_true + nt * nytrue), 2) - sigma_offset;
    }
}

static realtype fsres(
    realtype y, realtype sy, realtype my, realtype sigma_y, realtype ssigma_y,
    ObservableScaling scale
) {
    auto res = fres(y, my, sigma_y, scale);
    switch (scale) {
    case ObservableScaling::lin:
        return (-sy - ssigma_y * res) / sigma_y;
    case ObservableScaling::log:
        return (-sy / y - ssigma_y * res) / sigma_y;
    case ObservableScaling::log10:
        return (-sy / (y * std::log(10)) - ssigma_y * res) / sigma_y;
    default:
        throw std::invalid_argument("only lin, log, log10 allowed.");
    }
}
static realtype
fsres_error(realtype sigma_y, realtype ssigma_y, realtype sigma_offset) {
    return ssigma_y / (fres_error(sigma_y, sigma_offset) * sigma_y);
}

void ReturnData::fsres(
    int const it, Model& model, SolutionState const& sol, ExpData const& edata
) {
    if (sres.empty())
        return;

    std::vector<realtype> y_it(ny, 0.0);
    model.getObservable(y_it, ts[it], sol.x);
    std::vector<realtype> sy_it(ny * nplist, 0.0);
    model.getObservableSensitivity(sy_it, ts[it], sol.x, sol.sx);

    std::vector<realtype> sigmay_it(ny, 0.0);
    model.getObservableSigma(sigmay_it, it, &edata);
    std::vector<realtype> ssigmay_it(ny * nplist, 0.0);
    model.getObservableSigmaSensitivity(ssigmay_it, sy_it, it, &edata);

    auto observedData = edata.getObservedDataPtr(it);
    for (int iy = 0; iy < nytrue; ++iy) {
        if (!edata.isSetObservedData(it, iy))
            continue;
        for (int ip = 0; ip < nplist; ++ip) {
            int idx = (iy + it * edata.nytrue()) * nplist + ip;

            sres.at(idx) = amici::fsres(
                y_it.at(iy), sy_it.at(iy + ny * ip), observedData[iy],
                sigmay_it.at(iy), ssigmay_it.at(iy + ny * ip),
                model.getObservableScaling(iy)
            );

            if (sigma_res) {
                int idx_res
                    = (iy + it * edata.nytrue() + edata.nytrue() * edata.nt())
                          * nplist
                      + ip;
                sres.at(idx_res) = fsres_error(
                    sigmay_it.at(iy), ssigmay_it.at(iy + ny * ip), sigma_offset
                );
            }
        }
    }
}

void ReturnData::fFIM(
    int it, Model& model, SolutionState const& sol, ExpData const& edata
) {
    if (FIM.empty())
        return;

    std::vector<realtype> y_it(ny, 0.0);
    model.getObservable(y_it, ts[it], sol.x);
    std::vector<realtype> sy_it(ny * nplist, 0.0);
    model.getObservableSensitivity(sy_it, ts[it], sol.x, sol.sx);

    std::vector<realtype> sigmay_it(ny, 0.0);
    model.getObservableSigma(sigmay_it, it, &edata);
    std::vector<realtype> ssigmay_it(ny * nplist, 0.0);
    model.getObservableSigmaSensitivity(ssigmay_it, sy_it, it, &edata);

    /*
     * https://www.wolframalpha.com/input/?i=d%2Fdu+d%2Fdv+log%28s%28u%2Cv%29%29+%2B+0.5+*+%28r%28u%2Cv%29%2Fs%28u%2Cv%29%29%5E2
     * r = (y - m)
     * r_du = y_du
     * d/du(d/(dv)(log(s) + 0.5 (r/s)^2)) =
     * -(2*y_du*s_dv)*r/s^3 - (2*y_dv*s_du)*r/s^3 + y_du_dv*r/s^2
     * 22222222222222222222   2222222222222222222   #############
     *
     * + (y_dv*y_du)/s^2 + (3*s_dv*s_du)*r^2/s^4 - s_du_dv*r^2/s^3
     *   111111111111111   333333333333333333333   ***************
     *
     * - (s_dv*s_du)/s^2 + (s_du_dv)/s
     *                     +++++++++++
     *
     *
     * we compute this using fsres:
     * sres_u * sres_v = (y_du*s - s_du * r) * (y_dv*s - s_dv * r) / s^4
     * = y_du*y_dv/s^2 - (y_du*s_dv + y_dv*s_du)*r/s^3 + s_du*s_dv*r^2/s^4
     *   1111111111111   22222222222222222222222222222   33333333333333333
     *
     * r should be on the same order as s. We keep 1/s^2 and drop 1/s terms.
     * drop:
     * ********: r^2/s^3 term, typically zero anyways
     * ########: r/s^2 term, requires second order sensitivities
     * ++++++++: 1/s term, typically zero anyways
     *
     * keep:
     * 123123123: accounted for by sres*sres,
     * but
     * - (s_dv*s_du)/s^2 is unaccounted for and
     * - (s_dv*y_du + s_du*y_dv)*r/s^3 is missing from 2222 and
     * + 2*(s_du*s_dv)*r^2/s^4 is missing from 3333
     *
     * s_dv*y_du and s_du*y_dv are usually zero since we do not have parameters
     * that affect both observables and sigmas. Accordingly, it is hard to know
     * emprically whether these terms are important or not.
     *
     * This leaves us with
     * + (s_du*s_dv)(2*r^2-s^2)/s^4
     * which may be problematic, since this expression does not factorise and
     * may thus introduce directions of negative curvature.
     *
     * For the least squares trick, where error residuals
     * er = sqrt(log(s) + c), with sensitivity er_du = s_du/(2*s*er). This
     * would yield terms (s_du*s_dv)*(s^2/(4*er^2))/s^4.
     * These terms are guaranteed to yield positive curvature, but go to zero
     * in the limit c -> Infty.
     *
     * Empirically, simply taking this limit and dropping all missing terms,
     * works substantially better. This was evaluated using the fides optimizer
     * on the Boehm2014 Benchmark example.
     */

    auto observedData = edata.getObservedDataPtr(it);

    for (int iy = 0; iy < nytrue; ++iy) {
        if (!edata.isSetObservedData(it, iy))
            continue;
        auto y = y_it.at(iy);
        auto m = observedData[iy];
        auto s = sigmay_it.at(iy);
        auto os = model.getObservableScaling(iy);
        // auto r = amici::fres(y, m, s);
        for (int ip = 0; ip < nplist; ++ip) {
            auto dy_i = sy_it.at(iy + ny * ip);
            auto ds_i = ssigmay_it.at(iy + ny * ip);
            auto sr_i = amici::fsres(y, dy_i, m, s, ds_i, os);
            realtype sre_i = 0.0;
            if (sigma_res)
                sre_i = fsres_error(s, ds_i, sigma_offset);
            for (int jp = 0; jp < nplist; ++jp) {
                auto dy_j = sy_it.at(iy + ny * jp);
                auto ds_j = ssigmay_it.at(iy + ny * jp);
                auto sr_j = amici::fsres(y, dy_j, m, s, ds_j, os);
                FIM.at(ip + nplist * jp) += sr_i * sr_j;
                if (sigma_res) {
                    auto sre_j = fsres_error(s, ds_j, sigma_offset);
                    FIM.at(ip + nplist * jp) += sre_i * sre_j;
                }
                /*+ ds_i*ds_j*(2*pow(r/pow(s,2.0), 2.0) - 1/pow(s,2.0));*/
            }
        }
    }
}

ModelContext::ModelContext(Model* model)
    : model_(model)
    , original_state_(model->getModelState()) {}

ModelContext::~ModelContext() noexcept(false) { restore(); }

void ModelContext::restore() { model_->setModelState(original_state_); }

} // namespace amici
