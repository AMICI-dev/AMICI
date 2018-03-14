#include "amici/backwardproblem.h"

#include "amici/model.h"
#include "amici/solver.h"
#include "amici/exception.h"
#include "amici/edata.h"
#include "amici/rdata.h"
#include "amici/forwardproblem.h"

#include <cstring>

namespace amici {
    
    /** default constructor
     * @param fwd pointer to corresponding forward problem
     * @return new BackwardProblem instance
     */
BackwardProblem::BackwardProblem(ForwardProblem *fwd) :
    llhS0(fwd->model->nJ*fwd->model->nplist(),0.0),
    xB(fwd->model->nx),
    xB_old(fwd->model->nx),
    dxB(fwd->model->nx),
    xQB(fwd->model->nJ*fwd->model->nplist()),
    xQB_old(fwd->model->nJ*fwd->model->nplist()),
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
     *
     */
    double tnext;

    if (model->nx <= 0 || rdata->sensi < AMICI_SENSI_ORDER_FIRST ||
        rdata->sensi_meth != AMICI_SENSI_ASA || model->nplist() == 0) {
        return;
    }
    
    solver->setupAMIB(this, model);

    int it = rdata->nt - 2;
    --iroot;
    while (it >= 0 || iroot >= 0) {

        /* check if next timepoint is a discontinuity or a data-point */
        tnext = getTnext(discs.data(), iroot, it);

        if (tnext < t) {
            solver->AMISolveB(tnext, AMICI_NORMAL);

            solver->AMIGetB(which, &t, &xB,
                                     &dxB);
            solver->AMIGetQuadB(which, &t, &xQB);
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
        solver->AMIReInitB(which, t, &xB, &dxB);

        solver->AMIQuadReInitB(which, &xQB);

        solver->AMICalcICB(which, t, &xB, &dxB);
    }

    /* we still need to integrate from first datapoint to tstart */
    if (t > model->t0()) {
        if (model->nx > 0) {
            /* solve for backward problems */
            solver->AMISolveB(model->t0(), AMICI_NORMAL);

            solver->AMIGetQuadB(which, &(t), &xQB);

            solver->AMIGetB(which, &(t), &xB, &dxB);
        }
    }
    

    for (int iJ = 0; iJ < model->nJ; iJ++) {
        if (iJ == 0) {
            for (int ip = 0; ip < model->nplist(); ++ip) {
                llhS0[ip] = 0.0;
                for (int ix = 0; ix < model->nxtrue; ++ix) {
                    llhS0[ip] += xB[ix] * sx.at(ix,ip);
                }
            }
        } else {
            for (int ip = 0; ip < model->nplist(); ++ip) {
                llhS0[ip + iJ * model->nplist()] = 0.0;
                for (int ix = 0; ix < model->nxtrue; ++ix) {
                    llhS0[ip + iJ * model->nplist()] +=
                        xB[ix + iJ * model->nxtrue] * sx.at(ix,ip)+
                        xB[ix] * sx.at(ix + iJ * model->nxtrue,ip);
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

void BackwardProblem::handleEventB(int iroot) {
    /**
     * handleEventB executes everything necessary for the handling of events
     * for the backward problem
     *
     * @param[out] iroot index of event @type int
     */

    /* store current values */
    xB_old = xB;
    xQB_old = xQB;

    for (int ie = 0; ie < model->ne; ie++) {

        if (rootidx[iroot * model->ne + ie] != 0) {

            model->fdeltaqB(ie, t, &x_disc[iroot],&xB,&xdot_disc[iroot], &xdot_old_disc[iroot]);

            model->fdeltaxB(ie, t, &x_disc[iroot],&xB,&xdot_disc[iroot], &xdot_old_disc[iroot]);

            for (int ix = 0; ix < model->nxtrue; ++ix) {
                for (int iJ = 0; iJ < model->nJ; ++iJ) {
                    xB[ix + iJ * model->nxtrue] +=
                        model->deltaxB[ix + iJ * model->nxtrue];
                    if (model->nz > 0) {
                        xB[ix + iJ * model->nxtrue] +=
                            dJzdx[iJ + ( ix + nroots[ie] * model->nx ) * model->nJ];
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
    }

    model->updateHeavisideB(&rootidx[iroot * model->ne]);
}

    /**
     * handleDataPoint executes everything necessary for the handling of data
     * points for the backward problems
     *
     * @param[in] it index of data point @type int
     */
void BackwardProblem::handleDataPointB(int it) {


    for (int ix = 0; ix < model->nxtrue; ix++) {
        for (int iJ = 0; iJ < model->nJ; iJ++)
            // we only need the 1:nxtrue slice here!
            xB[ix + iJ * model->nxtrue] +=
                dJydx[iJ + ( ix + it * model->nx ) * model->nJ];
    }
    solver->getDiagnosisB(it, rdata, this);
}
    
    /**
     * getTnext computes the next timepoint to integrate to. This is the maximum
     * of tdata and troot but also takes into account if it<0 or iroot<0 where
     * these expressions
     * do not necessarily make sense
     *
     * @param[in] troot timepoint of next event @type realtype
     * @param[in] iroot index of next event @type int
     * @param[in] it index of next data point @type int
     * @param[in] model pointer to model specification object @type Model
     * @return tnext next timepoint @type realtype
     */
realtype BackwardProblem::getTnext(const realtype *troot, const int iroot,
                                   const int it) {

    realtype tnext;

    if (it < 0) {
        tnext = troot[iroot];
    } else {
        if (iroot < 0) {
            tnext = rdata->ts[it];
        } else {
            if (model->ne > 0) {
                if (troot[iroot] > rdata->ts[it]) {
                    tnext = troot[iroot];
                } else {
                    tnext = rdata->ts[it];
                }
            } else {
                tnext = rdata->ts[it];
            }
        }
    }

    return (tnext);
}

} // namespace amici
