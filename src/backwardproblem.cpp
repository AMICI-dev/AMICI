#include "amici/backwardproblem.h"

#include "amici/model.h"
#include "amici/solver.h"
#include "amici/exception.h"
#include "amici/edata.h"
#include "amici/rdata.h"
#include "amici/forwardproblem.h"

#include <cstring>

namespace amici {
    
/** Construct backward problem from forward problem
  * @param fwd pointer to corresponding forward problem
  * @return new BackwardProblem instance
  */
BackwardProblem::BackwardProblem(const ForwardProblem *fwd) :
    llhS0(static_cast<decltype(llhS0)::size_type>(fwd->model->nJ * fwd->model->nplist()), 0.0),
    xB(fwd->model->nx_solver),
    dxB(fwd->model->nx_solver),
    xQB(fwd->model->nJ*fwd->model->nplist()),
    x_disc(fwd->getStatesAtDiscontinuities()),
    xdot_disc(fwd->getRHSAtDiscontinuities()),
    xdot_old_disc(fwd->getRHSBeforeDiscontinuities()),
    sx(fwd->getStateSensitivity()),
    nroots(fwd->getNumberOfRoots()),
    discs(fwd->getDiscontinuities()),
    irdiscs(fwd->getDiscontinuities()),
    rootidx(fwd->getRootIndexes()),
    dJydx(fwd->getDJydx()),
    dJzdx(fwd->getDJzdx())
    {
        t = fwd->getTime();
        model = fwd->model;
        solver = fwd->solver;
        rdata = fwd->rdata;
        iroot = fwd->getRootCounter();
    }


void BackwardProblem::workBackwardProblem() {
    /**
     * workBackwardProblem solves the backward problem. if adjoint
     * sensitivities are enabled this will also compute sensitivies
     * workForwardProblem should be called before this is function is called
     */

    if (model->nx_solver <= 0 || solver->getSensitivityOrder() < SensitivityOrder::first ||
        solver->getSensitivityMethod() != SensitivityMethod::adjoint || model->nplist() == 0) {
        return;
    }
    
    solver->setupAMIB(this, model);

    int it = rdata->nt - 2;
    --iroot;

    while (it >= 0 || iroot >= 0) {

        /* check if next timepoint is a discontinuity or a data-point */
        double tnext = getTnext(discs, iroot, it);

        if (tnext < t) {
            solver->solveB(tnext, AMICI_NORMAL);
            solver->getB(which, &t, &xB, &dxB);
            solver->getQuadB(which, &t, &xQB);
        }

        /* handle discontinuity */

        if (model->ne > 0 && rdata->nmaxevent > 0 && iroot >= 0) {
            if (tnext == discs.at(iroot)) {
                handleEventB(iroot);
                --iroot;
            }
        }

        /* handle data-point */
        if (tnext == rdata->ts[it]) {
            handleDataPointB(it);
            it--;
        }

        /* reinit states */
        solver->reInitB(which, t, &xB, &dxB);
        solver->quadReInitB(which, &xQB);
        solver->calcICB(which, t, &xB, &dxB);
    }

    /* we still need to integrate from first datapoint to tstart */
    if (t > model->t0()) {
        if (model->nx_solver > 0) {
            /* solve for backward problems */
            solver->solveB(model->t0(), AMICI_NORMAL);
            solver->getQuadB(which, &(t), &xQB);
            solver->getB(which, &(t), &xB, &dxB);
        }
    }
    

    for (int iJ = 0; iJ < model->nJ; iJ++) {
        if (iJ == 0) {
            for (int ip = 0; ip < model->nplist(); ++ip) {
                llhS0[ip] = 0.0;
                for (int ix = 0; ix < model->nxtrue_solver; ++ix) {
                    llhS0[ip] += xB[ix] * sx.at(ix,ip);
                }
            }
        } else {
            for (int ip = 0; ip < model->nplist(); ++ip) {
                llhS0[ip + iJ * model->nplist()] = 0.0;
                for (int ix = 0; ix < model->nxtrue_solver; ++ix) {
                    llhS0[ip + iJ * model->nplist()] +=
                        xB[ix + iJ * model->nxtrue_solver] * sx.at(ix,ip)+
                        xB[ix] * sx.at(ix + iJ * model->nxtrue_solver,ip);
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

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

void BackwardProblem::handleEventB(const int iroot) {
    /**
     * handleEventB executes everything necessary for the handling of events
     * for the backward problem
     *
     * @param iroot index of event @type int
     */

    for (int ie = 0; ie < model->ne; ie++) {

        if (rootidx[iroot * model->ne + ie] == 0) {
            continue;
        }

        model->fdeltaqB(ie, t, &x_disc[iroot],&xB,&xdot_disc[iroot], &xdot_old_disc[iroot]);
        model->fdeltaxB(ie, t, &x_disc[iroot],&xB,&xdot_disc[iroot], &xdot_old_disc[iroot]);

        for (int ix = 0; ix < model->nxtrue_solver; ++ix) {
            for (int iJ = 0; iJ < model->nJ; ++iJ) {
                xB[ix + iJ * model->nxtrue_solver] +=
                        model->deltaxB[ix + iJ * model->nxtrue_solver];
                if (model->nz > 0) {
                    xB[ix + iJ * model->nxtrue_solver] +=
                            dJzdx[iJ + ( ix + nroots[ie] * model->nx_solver ) * model->nJ];
                }
            }
        }

        for (int iJ = 0; iJ < model->nJ; ++iJ) {
            for (int ip = 0; ip < model->nplist(); ++ip) {
                xQB[ip + iJ * model->nplist()] +=
                        model->deltaqB[ip + iJ * model->nplist()];
            }
        }

        nroots[ie]--;
    }

    model->updateHeavisideB(&rootidx[iroot * model->ne]);
}

/**
 * handleDataPoint executes everything necessary for the handling of data
 * points for the backward problems
 *
 * @param it index of data point @type int
 */
void BackwardProblem::handleDataPointB(const int it) {
    for (int ix = 0; ix < model->nxtrue_solver; ix++) {
        for (int iJ = 0; iJ < model->nJ; iJ++)
            // we only need the 1:nxtrue_solver (not the nx_true) slice here!
            xB[ix + iJ * model->nxtrue_solver] +=
                dJydx[iJ + ( ix + it * model->nx_solver ) * model->nJ];
    }
    solver->getDiagnosisB(it, rdata, this->which);
}
    
/**
 * getTnext computes the next timepoint to integrate to. This is the maximum
 * of tdata and troot but also takes into account if it<0 or iroot<0 where
 * these expressions
 * do not necessarily make sense
 *
 * @param troot timepoint of next event @type realtype
 * @param iroot index of next event @type int
 * @param it index of next data point @type int
 * @param model pointer to model specification object @type Model
 * @return tnext next timepoint @type realtype
 */
realtype BackwardProblem::getTnext(std::vector<realtype> const& troot, const int iroot,
                                   const int it) {
    if (it < 0
            || (iroot >= 0 && model->ne > 0 && troot.at(iroot) > rdata->ts[it])) {
        return troot.at(iroot);
    }

    return rdata->ts[it];
}

} // namespace amici
