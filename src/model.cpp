#include "amici/model.h"

#include "amici/amici.h"

#include <cstring>
#include <cmath>
#include <typeinfo>

namespace amici {

/** Sensitivity of measurements y, total derivative
 * @param it timepoint index
 * @param rdata pointer to return data instance
 */
void Model::fsy(const int it, ReturnData *rdata) {
    // Compute sy = dydx * sx + dydp
    for (int ip = 0; ip < nplist(); ++ip) {
        for (int iy = 0; iy < ny; ++iy)
            // copy dydp to sy
            rdata->sy.at((it*nplist()+ip)*ny+iy) =
                    dydp.at(iy + ip * ny);
    }
    // compute sy = 1.0*dydx*sx + 1.0*sy
    // dydx A[ny,nx] * sx B[nx,nplist] = sy C[ny,nplist]
    //        M  K          K  N              M  N
    //        lda           ldb               ldc
    amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans, ny, nplist(), nx,
                1.0, dydx.data(), ny, getsx(it,rdata), nx, 1.0,
                &rdata->sy[it*nplist()*ny], ny);
}

/** Sensitivity of z at final timepoint (ignores sensitivity of timepoint),
 * total derivative
 * @param nroots event index
 * @param rdata pointer to return data instance
 */
void Model::fsz_tf(const int nroots, ReturnData *rdata) {
    // Compute sz = dzdx * sz + dzdp

    for (int ip = 0; ip < nplist(); ++ip) {
        for (int iz = 0; iz < nz; ++iz)
            // copy dydp to sy
            rdata->sz.at((nroots*nplist()+ip)*nz + iz) =
                    0.0;
    }
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy, total
 * derivative
 * @param it timepoint index
 * @param dJydx vector with values of state derivative of Jy
 * @param rdata pointer to return data instance
 */
void Model::fsJy(const int it, const std::vector<realtype> dJydx, ReturnData *rdata) {

    // Compute dJydx*sx for current 'it'
    // dJydx        rdata->nt x nJ x nx
    // sx           rdata->nt x nx x nplist()
    std::vector<realtype> multResult(nJ * nplist(), 0);

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans, nJ,
                nplist(), nx, 1.0, &dJydx.at(it*nJ*nx), nJ, getsx(it,rdata), nx, 0.0,
                multResult.data(), nJ);

    // multResult    nJ x nplist()
    // dJydp         nJ x nplist()
    // dJydxTmp      nJ x nx
    // sxTmp         nx x nplist()

    // sJy += multResult + dJydp
    for (int iJ = 0; iJ < nJ; ++iJ) {
        if (iJ == 0)
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->sllh.at(ip) -= multResult.at(ip * nJ) + dJydp.at(ip * nJ);
        else
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->s2llh.at((iJ - 1) + ip * (nJ - 1)) -=
                        multResult.at(iJ + ip * nJ) + dJydp.at(iJ + ip * nJ);
    }
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t.
 * parameters
 * @param it timepoint index
 * @param edata pointer to experimental data instance
 * @param rdata pointer to return data instance
 */
void Model::fdJydp(const int it, const ExpData *edata,
                   const ReturnData *rdata) {

    // dJydy         nytrue x nJ x ny
    // dydp          ny x nplist()
    // dJydp         nJ x nplist()
    getmy(it,edata);
    fill(dJydp.begin(),dJydp.end(),0.0);
    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (isNaN(my.at(iyt)))
            continue;

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, nplist(), ny, 1.0, &dJydy.at(iyt*nJ*ny), nJ, dydp.data(), ny,
                    1.0, dJydp.data(), nJ);

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, nplist(), ny, 1.0, &dJydsigma.at(iyt*nJ*ny), nJ,
                    dsigmaydp.data(), ny, 1.0, dJydp.data(), nJ);
    }
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t.
 * state variables
 * @param dJydx pointer to vector with values of state derivative of Jy
 * @param it timepoint index
 * @param edata pointer to experimental data instance
 * @param rdata pointer to return data instance
 */
void Model::fdJydx(std::vector<realtype> *dJydx, const int it, const ExpData *edata, const ReturnData *rdata) {

    // dJydy         nJ x ny x nytrue
    // dydx          ny x nx
    // dJydx         nJ x nx x nt
    getmy(it,edata);
    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (isNaN(my.at(iyt)))
            continue;
    // dJydy A[nyt,nJ,ny] * dydx B[ny,nx] = dJydx C[it,nJ,nx]
    //         slice                                slice
    //             M  K            K  N                M  N
    //             lda             ldb                 ldc
        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, nx, ny, 1.0, &dJydy.at(iyt*ny*nJ), nJ, dydx.data(), ny, 1.0,
                    &dJydx->at(it*nx*nJ), nJ);
    }
}

/** Sensitivity of event-resolved measurement negative log-likelihood Jz, total
 * derivative
 * @param nroots event index
 * @param dJzdx vector with values of state derivative of Jz
 * @param sx pointer to state sensitivities
 * @param rdata pointer to return data instance
 */
void Model::fsJz(const int nroots, const std::vector<realtype> dJzdx, AmiVectorArray *sx, ReturnData *rdata) {
    // sJz           nJ x nplist()
    // dJzdp         nJ x nplist()
    // dJzdx         nmaxevent x nJ x nx
    // sx            rdata->nt x nx x nplist()

    // Compute dJzdx*sx for current 'ie'
    // dJzdx        rdata->nt x nJ x nx
    // sx           rdata->nt x nx x nplist()

    std::vector<realtype> multResult(nJ * nplist(), 0);
    std::vector<realtype> sxTmp(rdata->nplist * nx, 0);
    for (int ip = 0; ip < rdata->nplist; ++ip) {
        for (int ix = 0; ix < nx; ++ix)
            sxTmp.at(ix + ip * nx) = sx->at(ix,ip);
    }

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans, nJ,
                nplist(), nx, 1.0, &dJzdx.at(nroots*nx*nJ), nJ, sxTmp.data(), nx, 1.0,
                multResult.data(), nJ);

    // sJy += multResult + dJydp
    for (int iJ = 0; iJ < nJ; ++iJ) {
        if (iJ == 0)
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->sllh.at(ip) -= multResult.at(ip * nJ) + dJzdp.at(ip * nJ);
        else
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->s2llh.at((iJ - 1) + ip * (nJ - 1)) -=
                        multResult.at(iJ + ip * nJ) + dJzdp.at(iJ + ip * nJ);
    }
}

/** Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t.
 * parameters
 * @param nroots event index
 * @param t current timepoint
 * @param edata pointer to experimental data instance
 * @param rdata pointer to return data instance
 */
void Model::fdJzdp(const int nroots, realtype t, const ExpData *edata,
                   const ReturnData *rdata) {
    // dJzdz         nJ x nz x nztrue
    // dJzdsigma     nJ x nz x nztrue
    // dzdp          nz x nplist()
    // dJzdp         nJ x nplist()

    getmz(nroots,edata);
    std::fill(dJzdp.begin(),dJzdp.end(),0.0);
    for (int izt = 0; izt < nztrue; ++izt) {
        if (isNaN(mz.at(izt)))
            continue;

        if (t < rdata->ts.at(rdata->ts.size() - 1)) {
            // with z
            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                        nJ, nplist(), nz, 1.0, &dJzdz.at(izt*nz*nJ), nJ, dzdp.data(), nz,
                        1.0, dJzdp.data(), nJ);
        } else {
            // with rz
            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                        AMICI_BLAS_NoTrans, nJ, nplist(), nz, 1.0,
                        &dJrzdsigma.at(izt*nz*nJ), nJ, dsigmazdp.data(), nz, 1.0,
                        dJzdp.data(), nJ);
            
            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                        nJ, nplist(), nz, 1.0, &dJrzdz.at(izt*nz*nJ), nJ, dzdp.data(), nz,
                        1.0, dJzdp.data(), nJ);
        }

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, nplist(), nz, 1.0, &dJzdsigma.at(izt*nz*nJ), nJ,
                    dsigmazdp.data(), nz, 1.0, dJzdp.data(), nJ);
    }
}

/** Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t.
 * state variables
 * @param dJzdx pointer to vector with values of state derivative of Jz
 * @param nroots event index
 * @param t current timepoint
 * @param edata pointer to experimental data instance
 * @param rdata pointer to return data instance
 */
void Model::fdJzdx(std::vector<realtype> *dJzdx, const int nroots, realtype t, const ExpData *edata, const ReturnData *rdata) {
    // dJzdz         nJ x nz x nztrue
    // dzdx          nz x nx
    // dJzdx         nJ x nx x nmaxevent
    getmz(nroots,edata);
    for (int izt = 0; izt < nztrue; ++izt) {
        if (isNaN(mz.at(izt)))
            continue;
        
        if (t < rdata->ts.at(rdata->ts.size() - 1)) {
            // z
            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                        AMICI_BLAS_NoTrans, nJ, nx, nz, 1.0, &dJzdz.at(izt*nz*nJ), nJ,
                        dzdx.data(), nz, 1.0, &dJzdx->at(nroots*nx*nJ), nJ);
        } else {
            // rz
            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                        AMICI_BLAS_NoTrans, nJ, nx, nz, 1.0, &dJrzdz.at(izt*nz*nJ), nJ,
                        drzdx.data(), nz, 1.0, &dJzdx->at(nroots*nx*nJ), nJ);
        }
    }
}

/** initialization of model properties
 * @param x pointer to state variables
 * @param dx pointer to time derivative of states (DAE only)
 */
void Model::initialize(AmiVector *x, AmiVector *dx) {

    initializeStates(x);
    
    fdx0(x, dx);
    
    initHeaviside(x,dx);
    
}

/** initialization of initial states
 * @param x pointer to state variables
 */
void Model::initializeStates(AmiVector *x) {

    if (x0data.empty()) {
        fx0(x);
    } else {
        for (int ix = 0; ix < nx; ix++) {
            (*x)[ix] = (realtype) x0data.at(ix);
        }
    }
}

/**
 * initHeaviside initialises the heaviside variables h at the intial time t0
 * heaviside variables activate/deactivate on event occurences
 * @param x pointer to state variables
 * @param dx pointer to time derivative of states (DAE only)
 */
void Model::initHeaviside(AmiVector *x, AmiVector *dx) {
    std::vector<realtype> rootvals(ne,0.0);
    froot(tstart, x, dx, rootvals.data());
    for (int ie = 0; ie < ne; ie++) {
        if (rootvals.at(ie) < 0) {
            h.at(ie) = 0.0;
        } else if (rootvals.at(ie) == 0) {
            throw AmiException("Simulation started in an event. This could lead to "
                               "unexpected results, aborting simulation! Please "
                               "specify an earlier simulation start via "
                               "options.t0");
        } else {
            h.at(ie) = 1.0;
        }
    }
}


void Model::initializeVectors()
{
    dsigmaydp.resize(ny * nplist(), 0.0);
    dsigmazdp.resize(nz * nplist(), 0.0);
    dJydp.resize(nJ * nplist(), 0.0);
    dJzdp.resize(nJ * nplist(), 0.0);
    deltasx.resize(nx * nplist(), 0.0);
    deltaqB.resize(nJ * nplist(), 0.0);
    dxdotdp.resize(nx * nplist(), 0.0);
    dzdp.resize(nz * nplist(), 0.0);
    drzdp.resize(nz * nplist(), 0.0);
    dydp.resize(ny * nplist(), 0.0);
    stau.resize(nplist(), 0.0);
}

/** Initial states
 * @param x pointer to state variables
 */
void Model::fx0(AmiVector *x) {
    x->reset();
    fx0(x->data(),tstart, unscaledParameters.data(),fixedParameters.data());
}

/** Initial value for initial state sensitivities
     * @param sx pointer to state sensitivity variables
     * @param x pointer to state variables
     **/
void Model::fsx0(AmiVectorArray *sx, const AmiVector *x) {
    sx->reset();
    for(int ip = 0; (unsigned)ip<plist_.size(); ip++)
        fsx0(sx->data(ip),tstart,x->data(), unscaledParameters.data(),fixedParameters.data(),plist_.at(ip));
}

/** Sensitivity of event timepoint, total derivative
     * @param t current timepoint
     * @param ie event index
     * @param x pointer to state variables
     * @param sx pointer to state sensitivity variables
     */
void Model::fstau(const realtype t, const int ie, const AmiVector *x, const AmiVectorArray *sx) {
    std::fill(stau.begin(),stau.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++){
        fstau(&stau.at(ip),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),sx->data(ip),plist_.at(ip),ie);
    }
}

/** Observables / measurements
     * @param it timepoint index
     * @param rdata pointer to return data instance
     */
void Model::fy(int it, ReturnData *rdata) {
    fy(&rdata->y.at(it*ny),rdata->ts.at(it),getx(it,rdata), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/** partial derivative of observables y w.r.t. model parameters p
     * @param it timepoint index
     * @param rdata pointer to return data instance
     */
void Model::fdydp(const int it, ReturnData *rdata) {
    std::fill(dydp.begin(),dydp.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++){
        fdydp(&dydp.at(ip*ny),rdata->ts.at(it),getx(it,rdata), unscaledParameters.data(),fixedParameters.data(),h.data(),plist_.at(ip));
    }
}

/** partial derivative of observables y w.r.t. state variables x
     * @param it timepoint index
     * @param rdata pointer to return data instance
     */
void Model::fdydx(const int it, ReturnData *rdata) {
    std::fill(dydx.begin(),dydx.end(),0.0);
    fdydx(dydx.data(),rdata->ts.at(it),getx(it,rdata), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/** Event-resolved output
     * @param nroots number of events for event index
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param rdata pointer to return data instance
     */
void Model::fz(const int nroots, const int ie, const realtype t, const AmiVector *x, ReturnData *rdata) {
    fz(&rdata->z.at(nroots*nz),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/** Sensitivity of z, total derivative
     * @param nroots number of events for event index
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param sx current state sensitivities
     * @param rdata pointer to return data instance
     */
void Model::fsz(const int nroots, const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata) {
    for(int ip = 0; (unsigned)ip < plist_.size();  ip++ ){
        fsz(&rdata->sz.at((nroots*nplist()+ip)*nz),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),sx->data(ip),plist_.at(ip));
    }
}

/** Event root function of events (equal to froot but does not include
     * non-output events)
     * @param nroots number of events for event index
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param rdata pointer to return data instance
     */
void Model::frz(const int nroots, const int ie, const realtype t, const AmiVector *x, ReturnData *rdata) {
    frz(&rdata->rz.at(nroots*nz),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/** Sensitivity of rz, total derivative
     * @param nroots number of events for event index
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param sx current state sensitivities
     * @param rdata pointer to return data instance
     */
void Model::fsrz(const int nroots, const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata) {
    for(int ip = 0; (unsigned)ip < plist_.size();  ip++ ){
        fsrz(&rdata->srz.at((nroots*nplist()+ip)*nz),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),sx->data(ip),plist_.at(ip));
    }
}

/** partial derivative of event-resolved output z w.r.t. to model parameters p
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
void Model::fdzdp(const realtype t, const int ie, const AmiVector *x) {
    std::fill(dzdp.begin(),dzdp.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++){
        fdzdp(dzdp.data(),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),plist_.at(ip));
    }
}

/** partial derivative of event-resolved output z w.r.t. to model states x
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
void Model::fdzdx(const realtype t, const int ie, const AmiVector *x) {
    std::fill(dzdx.begin(),dzdx.end(),0.0);
    fdzdx(dzdx.data(),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/** Sensitivity of event-resolved root output w.r.t. to model parameters p
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
void Model::fdrzdp(const realtype t, const int ie, const AmiVector *x) {
    std::fill(drzdp.begin(),drzdp.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++){
        fdrzdp(drzdp.data(),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),plist_.at(ip));
    }
}

/** Sensitivity of event-resolved measurements rz w.r.t. to model states x
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
void Model::fdrzdx(const realtype t, const int ie, const AmiVector *x) {
    std::fill(drzdx.begin(),drzdx.end(),0.0);
    fdrzdx(drzdx.data(),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/** State update functions for events
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
void Model::fdeltax(const int ie, const realtype t, const AmiVector *x,
                    const AmiVector *xdot, const AmiVector *xdot_old) {
    std::fill(deltax.begin(),deltax.end(),0.0);
    fdeltax(deltax.data(),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),ie,xdot->data(),xdot_old->data());
}

/** Sensitivity update functions for events, total derivative
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param sx current state sensitivity
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
void Model::fdeltasx(const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx,
                     const AmiVector *xdot, const AmiVector *xdot_old) {
    fw(t,x->getNVector());
    std::fill(deltasx.begin(),deltasx.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++)
        fdeltasx(&deltasx.at(nx*ip),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),w.data(),
                 plist_.at(ip),ie,xdot->data(),xdot_old->data(),sx->data(ip),&stau.at(ip));
}

/** Adjoint state update functions for events
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param xB current adjoint state
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
void Model::fdeltaxB(const int ie, const realtype t, const AmiVector *x, const AmiVector *xB,
                     const AmiVector *xdot, const AmiVector *xdot_old) {
    std::fill(deltaxB.begin(),deltaxB.end(),0.0);
    fdeltaxB(deltaxB.data(),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),ie,xdot->data(),xdot_old->data(),xB->data());
}

/** Quadrature state update functions for events
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param xB current adjoint state
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
void Model::fdeltaqB(const int ie, const realtype t, const AmiVector *x, const AmiVector *xB,
                     const AmiVector *xdot, const AmiVector *xdot_old) {
    std::fill(deltaqB.begin(),deltaqB.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++)
        fdeltaqB(deltaqB.data(),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                 plist_.at(ip),ie,xdot->data(),xdot_old->data(),xB->data());
}

/** Standard deviation of measurements
     * @param it timepoint index
     * @param edata pointer to experimental data instance
     * @param rdata pointer to return data instance
     */
void Model::fsigma_y(const int it, const ExpData *edata, ReturnData *rdata) {
    std::fill(sigmay.begin(),sigmay.end(),0.0);
    fsigma_y(sigmay.data(),rdata->ts.at(it), unscaledParameters.data(),fixedParameters.data());
    for (int iytrue = 0; iytrue < nytrue; iytrue++) {
        /* extract the value for the standard deviation, if the data value
             is NaN, use
             the parameter value. Store this value in the return struct */
        if(edata){
            if (!isNaN(edata->sigmay[it * edata->nytrue + iytrue])) {
                sigmay.at(iytrue) = edata->sigmay[it * edata->nytrue + iytrue];
            }
        }
        rdata->sigmay[it * rdata->ny + iytrue] = sigmay.at(iytrue);
    }
}

/** partial derivative of standard deviation of measurements w.r.t. model
     * @param it timepoint index
     * @param rdata pointer to return data instance
     */
void Model::fdsigma_ydp(const int it, const ReturnData *rdata) {
    std::fill(dsigmaydp.begin(),dsigmaydp.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++)
        fdsigma_ydp(dsigmaydp.data(),rdata->ts.at(it), unscaledParameters.data(),fixedParameters.data(),plist_.at(ip));
}

/** Standard deviation of events
     * @param t current timepoint
     * @param ie event index
     * @param nroots array with event numbers
     * @param edata pointer to experimental data instance
     * @param rdata pointer to return data instance
     */
void Model::fsigma_z(const realtype t, const int ie, const int *nroots,
                     const ExpData *edata, ReturnData *rdata) {
    std::fill(sigmaz.begin(),sigmaz.end(),0.0);
    fsigma_z(sigmaz.data(),t, unscaledParameters.data(),fixedParameters.data());
    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (z2event.at(iztrue) - 1 == ie) {
            if(edata) {
                if (!isNaN(edata->sigmaz[nroots[ie]*edata->nztrue + iztrue])) {
                    sigmaz.at(iztrue) = edata->sigmaz[nroots[ie]*edata->nztrue + iztrue];
                }
            }
            rdata->sigmaz[nroots[ie]*rdata->nz + iztrue] = sigmaz.at(iztrue);
        }
    }
}

/** Sensitivity of standard deviation of events measurements w.r.t. model parameters p
     * @param t current timepoint
     */
void Model::fdsigma_zdp(const realtype t) {
    std::fill(dsigmazdp.begin(),dsigmazdp.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++)
        fdsigma_zdp(dsigmazdp.data(),t, unscaledParameters.data(),fixedParameters.data(),plist_.at(ip));
}

/** negative log-likelihood of measurements y
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fJy(const int it, ReturnData *rdata, const ExpData *edata) {
    std::vector<realtype> nllh(nJ,0.0);
    getmy(it,edata);
    for(int iytrue = 0; iytrue < nytrue; iytrue++){
        if(!isNaN(my.at(iytrue))){
            std::fill(nllh.begin(),nllh.end(),0.0);
            fJy(nllh.data(),iytrue, unscaledParameters.data(),fixedParameters.data(),gety(it,rdata),sigmay.data(),my.data());
            rdata->llh -= nllh.at(0);
        }
    }
}

/** negative log-likelihood of event-resolved measurements z
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fJz(const int nroots, ReturnData *rdata, const ExpData *edata) {
    std::vector<realtype> nllh(nJ,0.0);
    getmz(nroots,edata);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            std::fill(nllh.begin(),nllh.end(),0.0);
            fJz(nllh.data(),iztrue, unscaledParameters.data(),fixedParameters.data(),getz(nroots,rdata),sigmaz.data(),mz.data());
            rdata->llh -= nllh.at(0);
        }
    }
}

/** regularization of negative log-likelihood with roots of event-resolved
     * measurements rz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fJrz(const int nroots, ReturnData *rdata, const ExpData *edata) {
    std::vector<realtype> nllh(nJ,0.0);
    getrz(nroots,rdata);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            std::fill(nllh.begin(),nllh.end(),0.0);
            fJrz(nllh.data(),iztrue, unscaledParameters.data(),fixedParameters.data(),getrz(nroots,rdata),sigmaz.data());
            rdata->llh -= nllh.at(0);
        }
    }
}

/** partial derivative of time-resolved measurement negative log-likelihood Jy
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fdJydy(const int it, const ReturnData *rdata,
                   const ExpData *edata) {
    getmy(it,edata);
    std::fill(dJydy.begin(),dJydy.end(),0.0);
    for(int iytrue = 0; iytrue < nytrue; iytrue++){
        if(!isNaN(my.at(iytrue))){
            fdJydy(&dJydy.at(iytrue*ny*nJ),iytrue, unscaledParameters.data(),fixedParameters.data(),gety(it,rdata),sigmay.data(),my.data());
        }
    }
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy
     * w.r.t. standard deviation sigma
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fdJydsigma(const int it, const ReturnData *rdata,
                       const ExpData *edata) {
    getmy(it,edata);
    std::fill(dJydsigma.begin(),dJydsigma.end(),0.0);
    for(int iytrue = 0; iytrue < nytrue; iytrue++){
        if(!isNaN(my.at(iytrue))){
            fdJydsigma(&dJydsigma.at(iytrue*ny*nJ),iytrue, unscaledParameters.data(),fixedParameters.data(),gety(it,rdata),sigmay.data(),my.data());
        }
    }
}

/** partial derivative of event measurement negative log-likelihood Jz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fdJzdz(const int nroots, const ReturnData *rdata,
                   const ExpData *edata) {
    getmz(nroots,edata);
    std::fill(dJzdz.begin(),dJzdz.end(),0.0);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            fdJzdz(&dJzdz.at(iztrue*nz*nJ),iztrue, unscaledParameters.data(),fixedParameters.data(),getz(nroots,rdata),sigmaz.data(),mz.data());
        }
    }
}

/** Sensitivity of event measurement negative log-likelihood Jz
     * w.r.t. standard deviation sigmaz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fdJzdsigma(const int nroots, const ReturnData *rdata,
                       const ExpData *edata) {
    getmz(nroots,edata);
    std::fill(dJzdsigma.begin(),dJzdsigma.end(),0.0);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            fdJzdsigma(&dJzdsigma.at(iztrue*nz*nJ),iztrue, unscaledParameters.data(),fixedParameters.data(),getz(nroots,rdata),sigmaz.data(),mz.data());
        }
    }
}

/** partial derivative of event measurement negative log-likelihood Jz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fdJrzdz(const int nroots, const ReturnData *rdata,
                    const ExpData *edata) {
    getmz(nroots,edata);
    std::fill(dJrzdz.begin(),dJrzdz.end(),0.0);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            fdJrzdz(&dJrzdz.at(iztrue*nz*nJ),iztrue, unscaledParameters.data(),fixedParameters.data(),getrz(nroots,rdata),sigmaz.data());
        }
    }
}

/** Sensitivity of event measurement negative log-likelihood Jz
     * w.r.t. standard deviation sigmaz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fdJrzdsigma(const int nroots,const ReturnData *rdata,
                        const ExpData *edata) {
    std::fill(dJrzdsigma.begin(),dJrzdsigma.end(),0.0);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            fdJrzdsigma(&dJrzdsigma.at(iztrue*nz*nJ),iztrue, unscaledParameters.data(),fixedParameters.data(),getrz(nroots,rdata),sigmaz.data());
        }
    }
}

/**
     * @brief Recurring terms in xdot
     * @param t timepoint
     * @param x Vector with the states
     */
void Model::fw(const realtype t, const N_Vector x) {
    std::fill(w.begin(),w.end(),0.0);
    fw(w.data(),t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/**
     * @brief Recurring terms in xdot, parameter derivative
     * @param t timepoint
     * @param x Vector with the states
     */
void Model::fdwdp(const realtype t, const N_Vector x) {
    fw(t,x);
    std::fill(dwdp.begin(),dwdp.end(),0.0);
    fdwdp(dwdp.data(),t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),w.data());
}

/**
     * @brief Recurring terms in xdot, state derivative
     * @param t timepoint
     * @param x Vector with the states
     */
void Model::fdwdx(const realtype t, const N_Vector x) {
    fw(t,x);
    std::fill(dwdx.begin(),dwdx.end(),0.0);
    fdwdx(dwdx.data(),t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),w.data());
}

/** create my slice at timepoint
     * @param it timepoint index
     * @param edata pointer to experimental data instance
     */
void Model::getmy(const int it, const ExpData *edata){
    if(edata) {
        for(int iytrue = 0; iytrue < nytrue; iytrue++){
            my.at(iytrue) = static_cast<realtype>(edata->my[it*edata->nytrue + iytrue]);
        }
    } else {
        for(int iytrue = 0; iytrue < nytrue; iytrue++){
            my.at(iytrue) = getNaN();
        }
    }
}

/** create x slice at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @return x x-slice from rdata instance
     */
const realtype *Model::getx(const int it, const ReturnData *rdata) const {
    return &rdata->x.at(it*nx);
}

/** create sx slice at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @return sx sx-slice from rdata instance
     */
const realtype *Model::getsx(const int it, const ReturnData *rdata) const {
    return &rdata->sx.at(it*nx*nplist());
}

/** create y slice at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @return y y-slice from rdata instance
     */
const realtype *Model::gety(const int it, const ReturnData *rdata) const {
    return &rdata->y.at(it*ny);
}


/** get current timepoint from index
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @return current timepoint
     */
realtype Model::gett(const int it, const ReturnData *rdata) const {
    return rdata->ts.at(it);
}

/** create mz slice at event
     * @param nroots event occurence
     * @param edata pointer to experimental data instance
     */
void Model::getmz(const int nroots, const ExpData *edata) {
    if(edata){
        for(int iztrue = 0; iztrue < nztrue; iztrue++){
            mz.at(iztrue) = static_cast<realtype>(edata->mz[nroots*edata->nztrue + iztrue]);
        }
    } else {
        for(int iztrue = 0; iztrue < nztrue; iztrue++){
            mz.at(iztrue) = getNaN();
        }
    }
}

/** create z slice at event
     * @param nroots event occurence
     * @param rdata pointer to return data instance
     * @return z slice
     */
const realtype *Model::getz(const int nroots, const ReturnData *rdata) const {
    return(&rdata->z.at(nroots*nz));
}

/** create rz slice at event
     * @param nroots event occurence
     * @param rdata pointer to return data instance
     * @return rz slice
     */
const realtype *Model::getrz(const int nroots, const ReturnData *rdata) const {
    return(&rdata->rz.at(nroots*nz));
}

/** create sz slice at event
     * @param nroots event occurence
     * @param ip sensitivity index
     * @param rdata pointer to return data instance
     * @return z slice
     */
const realtype *Model::getsz(const int nroots, const int ip, const ReturnData *rdata) const {
    return(&rdata->sz.at((nroots*nplist()+ip)*nz));
}

/** create srz slice at event
     * @param nroots event occurence
     * @param ip sensitivity index
     * @param rdata pointer to return data instance
     * @return rz slice
     */
const realtype *Model::getsrz(const int nroots, const int ip, const ReturnData *rdata) const {
    return(&rdata->srz.at((nroots*nplist()+ip)*nz));
}

int Model::checkFinite(const int N, const realtype *array, const char *fun) const
{
    auto result = amici::checkFinite(N, array, fun);

    if(result != AMICI_SUCCESS) {
        amici::checkFinite(ts.size(), ts.data(), "ts");
        amici::checkFinite(fixedParameters.size(), fixedParameters.data(), "k");
        amici::checkFinite(unscaledParameters.size(), unscaledParameters.data(), "p");
        amici::checkFinite(w.size(), w.data(), "w");
    }

    return result;
}

void Model::unscaleParameters(double *bufferUnscaled) const
{
    /**
         * unscaleParameters removes parameter scaling according to the parameter
         * scaling in pscale
         *
         * @param[out] bufferUnscaled unscaled parameters are written to the array
         * @type double
         *
         * @return status flag indicating success of execution @type int
         */
    for (int ip = 0; ip < np(); ++ip) {
        switch (pscale[ip]) {
        case AMICI_SCALING_LOG10:
            bufferUnscaled[ip] = pow(10, originalParameters[ip]);
            break;
        case AMICI_SCALING_LN:
            bufferUnscaled[ip] = exp(originalParameters[ip]);
            break;
        case AMICI_SCALING_NONE:
            bufferUnscaled[ip] = originalParameters[ip];
            break;
        }
    }
}

bool operator ==(const Model &a, const Model &b)
{
    if (typeid(a) != typeid(b))
            return false;

    return (a.nx == b.nx)
            && (a.nxtrue == b.nxtrue)
            && (a.ny == b.ny)
            && (a.nytrue == b.nytrue)
            && (a.nz == b.nz)
            && (a.nztrue == b.nztrue)
            && (a.ne == b.ne)
            && (a.nw == b.nw)
            && (a.ndwdx == b.ndwdx)
            && (a.ndwdp == b.ndwdp)
            && (a.nnz == b.nnz)
            && (a.nJ == b.nJ)
            && (a.ubw == b.ubw)
            && (a.lbw == b.lbw)
            && (a.o2mode == b.o2mode)
            && (a.z2event == b.z2event)
            && (a.idlist == b.idlist)
            && (a.h == b.h)
            && (a.unscaledParameters == b.unscaledParameters)
            && (a.originalParameters == b.originalParameters)
            && (a.fixedParameters == b.fixedParameters)
            && (a.plist_ == b.plist_)
            && (a.x0data == b.x0data)
            && (a.sx0data == b.sx0data)
            && (a.ts == b.ts)
            && (a.qpositivex == b.qpositivex)
            && (a.nmaxevent == b.nmaxevent)
            && (a.pscale == b.pscale)
            && (a.tstart == b.tstart);
}



} // namespace amici
