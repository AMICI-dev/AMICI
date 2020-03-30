#include "amici/backwardproblem.h"

#include "amici/model.h"
#include "amici/solver.h"
#include "amici/exception.h"
#include "amici/edata.h"
#include "amici/rdata.h"
#include "amici/forwardproblem.h"
#include "amici/steadystateproblem.h"
#include "amici/misc.h"

#include <cstring>

namespace amici {

BackwardProblem::BackwardProblem(const ForwardProblem &fwd,
                                 const SteadystateProblem *posteq):
    model(fwd.model),
    rdata(fwd.rdata),
    solver(fwd.solver),
    t(fwd.getTime()),
    llhS0(static_cast<decltype(llhS0)::size_type>(fwd.model->nJ *
                                                  fwd.model->nplist()), 0.0),
    xB(fwd.model->nx_solver),
    dxB(fwd.model->nx_solver),
    xQB(fwd.model->nJ*fwd.model->nplist()),
    x_disc(fwd.getStatesAtDiscontinuities()),
    xdot_disc(fwd.getRHSAtDiscontinuities()),
    xdot_old_disc(fwd.getRHSBeforeDiscontinuities()),
    sx0(fwd.getStateSensitivity()),
    nroots(fwd.getNumberOfRoots()),
    discs(fwd.getDiscontinuities()),
    rootidx(fwd.getRootIndexes()),
    dJydx(fwd.getDJydx()),
    dJzdx(fwd.getDJzdx()) {
        for (int it = 0; it < fwd.model->nt(); it++) {
            if (std::isinf(fwd.model->getTimepoint(it))) {
                if (!posteq)
                    throw AmiException("Model has non-finite timpoint but"
                    "postequilibration did not run");
                writeSlice(slice(posteq->getDJydx(), it,
                                 fwd.model->nx_solver * fwd.model->nJ),
                           slice(this->dJydx, it,
                                 fwd.model->nx_solver * fwd.model->nJ));
            }
        }
        
    }


void BackwardProblem::workBackwardProblem() {


    if (model->nx_solver <= 0 ||
        solver->getSensitivityOrder() < SensitivityOrder::first ||
        solver->getSensitivityMethod() != SensitivityMethod::adjoint ||
        model->nplist() == 0) {
        return;
    }
    
    int it = rdata->nt - 1;
    model->initializeB(xB, dxB, xQB);
    handleDataPointB(it);
    solver->setupB(&which, rdata->ts[it], model, xB, dxB, xQB);
    
    --it;

    while (it >= 0 || discs.size() >= 0) {

        /* check if next timepoint is a discontinuity or a data-point */
        double tnext = getTnext(it);

        if (tnext < t) {
            solver->runB(tnext);
            solver->writeSolutionB(&t, xB, dxB, xQB, this->which);
            solver->getDiagnosisB(it, rdata, this->which);
        }

        /* handle discontinuity */
        if (tnext > rdata->ts[it]) {
            handleEventB();
        }

        /* handle data-point */
        if (tnext == rdata->ts[it]) {
            handleDataPointB(it);
            it--;
        }

        /* reinit states */
        solver->reInitB(which, t, xB, dxB);
        solver->quadReInitB(which, xQB);
    }

    /* we still need to integrate from first datapoint to tstart */
    if (t > model->t0()) {
        /* solve for backward problems */
        solver->runB(model->t0());
        solver->writeSolutionB(&t, xB, dxB, xQB, this->which);
        solver->getDiagnosisB(0, rdata, this->which);
    }

    computeLikelihoodSensitivities();
    rdata->cpu_timeB = solver->getCpuTimeB();
}


void BackwardProblem::handleEventB() {
    auto rootidx = this->rootidx.back();
    this->rootidx.pop_back();
    
    auto x_disc = this->x_disc.back();
    this->x_disc.pop_back();
    
    auto xdot_disc = this->xdot_disc.back();
    this->xdot_disc.pop_back();
    
    auto xdot_old_disc = this->xdot_old_disc.back();
    this->xdot_old_disc.pop_back();
    
    for (int ie = 0; ie < model->ne; ie++) {

        if (rootidx[ie] == 0) {
            continue;
        }

        model->addAdjointQuadratureEventUpdate(xQB, ie, t, x_disc, xB,
                                               xdot_disc,
                                               xdot_old_disc);
        model->addAdjointStateEventUpdate(xB, ie, t, x_disc,
                                          xdot_disc,
                                          xdot_old_disc);

        for (int ix = 0; ix < model->nxtrue_solver; ++ix) {
            for (int iJ = 0; iJ < model->nJ; ++iJ) {
                if (model->nz > 0) {
                    xB[ix + iJ * model->nxtrue_solver] +=
                            dJzdx[iJ + ( ix + nroots[ie] * model->nx_solver )
                                  * model->nJ];
                }
            }
        }



        nroots[ie]--;
    }

    model->updateHeavisideB(rootidx.data());
}


void BackwardProblem::handleDataPointB(const int it) {
    for (int ix = 0; ix < model->nxtrue_solver; ix++) {
        for (int iJ = 0; iJ < model->nJ; iJ++)
            // we only need the 1:nxtrue_solver (not the nx_true) slice here!
            xB[ix + iJ * model->nxtrue_solver] +=
                dJydx[iJ + ( ix + it * model->nx_solver ) * model->nJ];
    }
}

realtype BackwardProblem::getTnext(const int it) {
    if (discs.size() > 0 && discs.back() > rdata->ts[it]) {
        double tdisc = discs.back();
        discs.pop_back();
        return tdisc;
    }

    return rdata->ts[it];
}


void BackwardProblem::computeLikelihoodSensitivities()
{
    for (int iJ = 0; iJ < model->nJ; iJ++) {
        if (iJ == 0) {
            for (int ip = 0; ip < model->nplist(); ++ip) {
                llhS0[ip] = 0.0;
                for (int ix = 0; ix < model->nxtrue_solver; ++ix) {
                    llhS0[ip] += xB[ix] * sx0.at(ix,ip);
                }
            }
        } else {
            for (int ip = 0; ip < model->nplist(); ++ip) {
                llhS0[ip + iJ * model->nplist()] = 0.0;
                for (int ix = 0; ix < model->nxtrue_solver; ++ix) {
                    llhS0[ip + iJ * model->nplist()] +=
                        xB[ix + iJ * model->nxtrue_solver] * sx0.at(ix,ip)+
                        xB[ix] * sx0.at(ix + iJ * model->nxtrue_solver,ip);
                }
            }
        }
    }

    for (int iJ = 0; iJ < model->nJ; iJ++) {
        for (int ip = 0; ip < model->nplist(); ip++) {
            if (iJ == 0) {
                rdata->sllh.at(ip) -= llhS0[ip] + xQB[ip*model->nJ];
            } else {
                rdata->s2llh.at(iJ - 1 + ip * (model->nJ - 1)) -=
                    llhS0[ip + iJ * model->nplist()] +
                    xQB[iJ + ip*model->nJ];
            }
        }
    }
}

} // namespace amici
