#include "amici/backwardproblem.h"

#include "amici/edata.h"
#include "amici/exception.h"
#include "amici/forwardproblem.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/steadystateproblem.h"

namespace amici {

BackwardProblem::BackwardProblem(ForwardProblem& fwd)
    : model_(fwd.model)
    , solver_(fwd.solver)
    , edata_(fwd.edata)
    , t_(fwd.getTime())
    , sx0_(fwd.getStateSensitivity())
    , discs_(fwd.getDiscontinuities())
    , dJydx_(fwd.getAdjointUpdates(*model_, *edata_))
    , dJzdx_(fwd.getDJzdx())
    , preeq_problem_(fwd.getPreequilibrationProblem())
    , posteq_problem_(fwd.getPostequilibrationProblem())
    , ws_(model_, solver_)
    , simulator_(model_, solver_, &ws_) {}

void BackwardProblem::workBackwardProblem() {

    if (model_->nx_solver <= 0
        || solver_->getSensitivityOrder() < SensitivityOrder::first
        || solver_->getSensitivityMethod() != SensitivityMethod::adjoint
        || model_->nplist() == 0) {
        return;
    }

    handlePostequilibration();

    // If we have posteq, infinity timepoints were already treated
    int it = model_->nt() - 1;
    while (it >= 0 && std::isinf(model_->getTimepoint(it))) {
        --it;
    }

    // initialize state vectors, depending on postequilibration
    model_->initializeB(ws_.xB_, ws_.dxB_, ws_.xQB_, it < model_->nt() - 1);
    ws_.discs_ = discs_;
    ws_.nroots_ = compute_nroots(discs_, model_->ne, model_->nMaxEvent());
    simulator_.run(
        t_, model_->t0(), it, model_->getTimepoints(), &dJydx_, &dJzdx_
    );

    // handle presimulation
    if (edata_ && edata_->t_presim > 0) {
        ConditionContext cc(
            model_, edata_, FixedParameterContext::presimulation
        );
        solver_->runB(model_->t0() - edata_->t_presim);
        solver_->writeSolutionB(&t_, ws_.xB_, ws_.dxB_, ws_.xQB_, ws_.which);
    }

    // handle pre-equilibration
    if (preeq_problem_) {
        ConditionContext cc2(
            model_, edata_, FixedParameterContext::preequilibration
        );
        preeq_problem_->workSteadyStateBackwardProblem(
            *solver_, *model_, ws_.xB_, true
        );
    }
}

void BackwardProblem::handlePostequilibration() {
    if (!posteq_problem_) {
        return;
    }

    // initialize xB - only process the postequilibration timepoints
    for (int it = 0; it < model_->nt(); it++) {
        if (std::isinf(model_->getTimepoint(it))) {
            for (int ix = 0; ix < model_->nxtrue_solver; ix++)
                ws_.xB_[ix] += dJydx_[ix + it * model_->nx_solver];
        }
    }

    posteq_problem_->workSteadyStateBackwardProblem(
        *solver_, *model_, ws_.xB_, false
    );
    ws_.xQB_ = posteq_problem_->getEquilibrationQuadratures();
}

void EventHandlingBwdSimulator::handleEventB(
    Discontinuity const& disc, std::vector<realtype> const* dJzdx
) {

    // TODO: unused!
    auto x_in_event = AmiVector(disc.x_post);
    for (int iv = 0; iv < x_in_event.getLength(); iv++)
        x_in_event[iv] = 0.5 * (disc.x_post.at(iv) + disc.x_pre.at(iv));

    // TODO: unused!
    auto xdot_in_event = AmiVector(disc.xdot_post);
    for (int iv = 0; iv < xdot_in_event.getLength(); iv++)
        xdot_in_event[iv]
            = 0.5 * (disc.xdot_post.at(iv) + disc.xdot_pre.at(iv));

    auto delta_x = AmiVector(disc.x_post.getLength(), solver_->getSunContext());
    for (int iv = 0; iv < xdot_in_event.getLength(); iv++)
        delta_x[iv] = (disc.x_post.at(iv) - disc.x_pre.at(iv));

    model_->updateHeavisideB(disc.root_info.data());

    for (int ie = 0; ie < model_->ne; ie++) {

        if (disc.root_info[ie] != 1) {
            continue;
        }

        model_->addAdjointQuadratureEventUpdate(
            ws_->xQB_, ie, t_, disc.x_pre, ws_->xB_, disc.xdot_post,
            disc.xdot_pre, delta_x
        );
        model_->addAdjointStateEventUpdate(
            ws_->xB_, ie, t_, disc.x_pre, disc.xdot_post, disc.xdot_pre, delta_x
        );

        if (model_->nz > 0) {
            for (int ix = 0; ix < model_->nxtrue_solver; ++ix) {
                for (int iJ = 0; iJ < model_->nJ; ++iJ) {
                    ws_->xB_[ix + iJ * model_->nxtrue_solver] += (*dJzdx
                    )[iJ
                      + (ix + ws_->nroots_[ie] * model_->nx_solver)
                            * model_->nJ];
                }
            }
        }

        ws_->nroots_[ie]--;
    }
}

void EventHandlingBwdSimulator::handleDataPointB(
    int const it, std::vector<realtype> const* dJydx
) {
    // The solver wasn't reset yet, as xB_ is necessary for solver setup.
    // For initial time point (we are integrating backwards!), diagnosis needs
    // to be stored outside this function.
    if (it < model_->nt() - 1)
        solver_->storeDiagnosisB(ws_->which);

    Expects(dJydx != nullptr);
    for (int ix = 0; ix < model_->nxtrue_solver; ix++) {
        for (int iJ = 0; iJ < model_->nJ; iJ++)
            // we only need the 1:nxtrue_solver (not the nx_true) slice here!
            ws_->xB_[ix + iJ * model_->nxtrue_solver]
                += (*dJydx)[iJ + (ix + it * model_->nx_solver) * model_->nJ];
    }
}

realtype EventHandlingBwdSimulator::getTnext(int const it) {
    if (it < 0 && ws_->discs_.empty()) {
        throw AmiException(
            "No more timepoints (it=%d, ie=%d) available at %f. This should "
            "not happen, please report a bug including this stacktrace at "
            "https://github.com/AMICI-dev/AMICI/issues/new/choose",
            it, ws_->discs_.size(), this->t_
        );
    }

    if (!ws_->discs_.empty()
        && (it < 0 || ws_->discs_.back().time > model_->getTimepoint(it))) {
        double tdisc = ws_->discs_.back().time;
        return tdisc;
    }

    return model_->getTimepoint(it);
}

BwdSimWorkspace::BwdSimWorkspace(
    gsl::not_null<Model*> model, gsl::not_null<Solver const*> solver
)
    : model_(model)
    , xB_(model_->nx_solver, solver->getSunContext())
    , dxB_(model_->nx_solver, solver->getSunContext())
    , xQB_(model_->nJ * model_->nplist(), solver->getSunContext()) {}

void EventHandlingBwdSimulator::run(
    realtype t_start, realtype t_end, realtype it,
    std::vector<realtype> const& timepoints, std::vector<realtype> const* dJydx,
    std::vector<realtype> const* dJzdx
) {
    Expects(model_->nz == 0 || dJzdx != nullptr);
    Expects(t_start >= t_end);
    Expects(it < 0 || t_start >= timepoints[it]);
    Expects(it < 0 || t_end <= timepoints.front());

    t_ = t_start;

    if ((it >= 0 || !ws_->discs_.empty()) && timepoints[it] > t_end) {
        handleDataPointB(it, dJydx);
        solver_->setupB(
            &ws_->which, timepoints[it], model_, ws_->xB_, ws_->dxB_, ws_->xQB_
        );
        // for the initial datapoint, diagnosis needs to be stored after setup
        // as it is not called in handleDataPointB
        solver_->storeDiagnosisB(ws_->which);
        --it;

        while (it >= 0 || !ws_->discs_.empty()) {
            // check if next timepoint is a discontinuity or a data-point
            double tnext = getTnext(it);

            if (tnext < t_) {
                solver_->runB(tnext);
                solver_->writeSolutionB(
                    &t_, ws_->xB_, ws_->dxB_, ws_->xQB_, ws_->which
                );
            }

            // handle data-point
            if (it >= 0 && tnext == timepoints[it]) {
                handleDataPointB(it, dJydx);
                it--;
            }

            // handle discontinuity
            if (!ws_->discs_.empty() && tnext == ws_->discs_.back().time) {
                handleEventB(ws_->discs_.back(), dJzdx);
                ws_->discs_.pop_back();
            }

            // reinitialize state
            solver_->reInitB(ws_->which, t_, ws_->xB_, ws_->dxB_);
            solver_->quadReInitB(ws_->which, ws_->xQB_);
        }
    }

    // we still need to integrate from first datapoint to t_start
    if (t_ > t_end) {
        solver_->runB(t_end);
        solver_->writeSolutionB(
            &t_, ws_->xB_, ws_->dxB_, ws_->xQB_, ws_->which
        );
    }
}

} // namespace amici
