#include "amici/backwardproblem.h"

#include "amici/model.h"
#include "amici/solver.h"
#include "amici/exception.h"
#include "amici/edata.h"
#include "amici/rdata.h"
#include "amici/forwardproblem.h"
#include "amici/newton_solver.h"

#include <cstring>

namespace amici {

BackwardProblem::BackwardProblem(const ForwardProblem *fwd) :
    model(fwd->model),
    rdata(fwd->rdata),
    solver(fwd->solver),
    t(fwd->getTime()),
    llhS0(static_cast<decltype(llhS0)::size_type>(fwd->model->nJ * fwd->model->nplist()), 0.0),
    xB(fwd->model->nx_solver),
    dxB(fwd->model->nx_solver),
    xQB(fwd->model->nJ*fwd->model->nplist()),
    x_disc(fwd->getStatesAtDiscontinuities()),
    xdot_disc(fwd->getRHSAtDiscontinuities()),
    xdot_old_disc(fwd->getRHSBeforeDiscontinuities()),
    sx0(fwd->getStateSensitivity()),
    nroots(fwd->getNumberOfRoots()),
    discs(fwd->getDiscontinuities()),
    irdiscs(fwd->getDiscontinuities()),
    iroot(fwd->getRootCounter()),
    rootidx(fwd->getRootIndexes()),
    dJydx(fwd->getDJydx()),
    dJzdx(fwd->getDJzdx()) {}


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

    if (std::isinf(model->getTimepoint(it)))
    {
      computeIntegralForSteadyState();
      xB.reset();
    }
    else
    {
      solver->setupB(&which, rdata->ts[it], model, xB, dxB, xQB);

      --it;
      --iroot;

      while (it >= 0 || iroot >= 0) {

          /* check if next timepoint is a discontinuity or a data-point */
          double tnext = getTnext(discs, iroot, it);

          if (tnext < t) {
              solver->runB(tnext);
              solver->writeSolutionB(&t, xB, dxB, xQB, this->which);
              solver->getDiagnosisB(it, rdata, this->which);
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
    }

    computeLikelihoodSensitivities();
    rdata->cpu_timeB = solver->getCpuTimeB();
}


void BackwardProblem::handleEventB(const int iroot) {
    for (int ie = 0; ie < model->ne; ie++) {

        if (rootidx[iroot * model->ne + ie] == 0) {
            continue;
        }

        model->addAdjointQuadratureEventUpdate(xQB, ie, t, x_disc[iroot], xB,
                                               xdot_disc[iroot],
                                               xdot_old_disc[iroot]);
        model->addAdjointStateEventUpdate(xB, ie, t, x_disc[iroot],
                                          xdot_disc[iroot],
                                          xdot_old_disc[iroot]);

        for (int ix = 0; ix < model->nxtrue_solver; ++ix) {
            for (int iJ = 0; iJ < model->nJ; ++iJ) {
                if (model->nz > 0) {
                    xB[ix + iJ * model->nxtrue_solver] +=
                            dJzdx[iJ + ( ix + nroots[ie] * model->nx_solver ) * model->nJ];
                }
            }
        }



        nroots[ie]--;
    }

    model->updateHeavisideB(&rootidx[iroot * model->ne]);
}


void BackwardProblem::handleDataPointB(const int it) {
    for (int ix = 0; ix < model->nxtrue_solver; ix++) {
        for (int iJ = 0; iJ < model->nJ; iJ++)
            // we only need the 1:nxtrue_solver (not the nx_true) slice here!
            xB[ix + iJ * model->nxtrue_solver] +=
                dJydx[iJ + ( ix + it * model->nx_solver ) * model->nJ];
    }
}

realtype BackwardProblem::getTnext(std::vector<realtype> const& troot,
                                   const int iroot, const int it) {
    if (it < 0
            || (iroot >= 0 && model->ne > 0 && troot.at(iroot) > rdata->ts[it])) {
        return troot.at(iroot);
    }

    return rdata->ts[it];
}

void BackwardProblem::computeIntegralForSteadyState()
{
    amici::AmiVector x_steadystate{rdata->x};
    realtype T = std::numeric_limits<realtype>::infinity();

    // Compute the linear system JB*v = dJydy
    auto newtonSolver = NewtonSolver::getSolver(
        &T, &x_steadystate, solver->getLinearSolver(), model, rdata,
        solver->getNewtonMaxLinearSteps(), solver->getNewtonMaxSteps(),
        solver->getAbsoluteTolerance(), solver->getRelativeTolerance(),
        solver->getNewtonDampingFactorMode(),
        solver->getNewtonDampingFactorLowerBound());

    newtonSolver->prepareLinearSystemB(0, -1);
    newtonSolver->solveLinearSystem(xB);

    // Compute the inner product v*dxotdp
    if (model->pythonGenerated)
    {
        for (int ip=0; ip<model->nplist(); ++ip)
        {
            if (model->ndxdotdp_explicit > 0) {
                auto col = model->dxdotdp_explicit.indexptrs();
                auto row = model->dxdotdp_explicit.indexvals();
                auto data_ptr = model->dxdotdp_explicit.data();
                for (sunindextype iCol = col[model->plist(ip)];
                     iCol < col[model->plist(ip) + 1]; ++iCol)
                    xQB[ip] += xB[row[iCol]] * data_ptr[iCol];
            }

            if (model->ndxdotdp_implicit > 0) {
                auto col = model->dxdotdp_implicit.indexptrs();
                auto row = model->dxdotdp_implicit.indexvals();
                auto data_ptr = model->dxdotdp_implicit.data();
                for (sunindextype iCol = col[model->plist(ip)];
                     iCol < col[model->plist(ip) + 1]; ++iCol)
                    xQB[ip] += xB[row[iCol]] * data_ptr[iCol];
            }
        }
    }
    else
    {
      for (int ip=0; ip<model->nplist(); ++ip)
        xQB[ip] = N_VDotProd(xB.getNVector(), model->dxdotdp.getNVector(ip));
    }
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
