#include "include/amici.h"
#include "include/amici_model.h"
#include <cstring>
#include <include/edata.h>
#include <include/rdata.h>

namespace amici {

/** Sensitivity of measurements y, total derivative
 * @param it timepoint index
 * @param rdata pointer to return data instance
 */
void Model::fsy(const int it, ReturnData *rdata) {
    // Compute sy = dydx * sx + dydp
    getsx(it,rdata);
    for (int ip = 0; ip < rdata->nplist; ++ip) {
        for (int iy = 0; iy < ny; ++iy)
            // copy dydp to sy
            rdata->sy[ip * rdata->nt * ny + iy * rdata->nt + it] =
                dydp.at(iy + ip * ny);

        // compute sy = 1.0*dydx*sx + 1.0*sy
        amici_dgemv(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, ny, nx, 1.0,
                    dydx.data(), ny, sx.at(ip).data(), 1, 1.0,
                    &rdata->sy[it + ip * rdata->nt * ny], rdata->nt);
    }
}

/** Sensitivity of z at final timepoint (ignores sensitivity of timepoint),
 * total derivative
 * @param nroots event index
 * @param rdata pointer to return data instance
 */
void Model::fsz_tf(const int nroots, ReturnData *rdata) {
    // Compute sz = dzdx * sz + dzdp

    for (int ip = 0; ip < rdata->nplist; ++ip) {
        for (int iz = 0; iz < nz; ++iz)
            // copy dydp to sy
            rdata->sz[nroots + (iz + ip * nz) * rdata->nmaxevent] =
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
    // sx           rdata->nt x nx x rdata->nplist
    std::vector<realtype> multResult(nJ * rdata->nplist, 0);
    std::vector<realtype> sxTmp(rdata->nplist * nx, 0);
    
    for (int ix = 0; ix < nx; ++ix) {
        for (int ip = 0; ip < rdata->nplist; ++ip)
            sxTmp.at(ix + ip * nx) = static_cast<realtype>(rdata->sx[it + (ix + ip * nx) * rdata->nt]);
    }

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans, nJ,
                rdata->nplist, nx, 1.0, &dJydx.at(it*nJ*nx), nJ, sxTmp.data(), nx, 0.0,
                multResult.data(), nJ);

    // multResult    nJ x rdata->nplist
    // dJydp         nJ x rdata->nplist
    // dJydxTmp      nJ x nx
    // sxTmp         nx x rdata->nplist

    // sJy += multResult + dJydp
    for (int iJ = 0; iJ < nJ; ++iJ) {
        if (iJ == 0)
            for (int ip = 0; ip < rdata->nplist; ++ip)
                rdata->sllh[ip] -= static_cast<double>(multResult.at(ip * nJ) + dJydp.at(ip * nJ));
        else
            for (int ip = 0; ip < rdata->nplist; ++ip)
                rdata->s2llh[(iJ - 1) + ip * (nJ - 1)] -=
                    static_cast<realtype>(multResult.at(iJ + ip * nJ) + dJydp.at(iJ + ip * nJ));
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
    // dydp          ny x rdata->nplist
    // dJydp         nJ x rdata->nplist
    getmy(it,edata);
    fill(dJydp.begin(),dJydp.end(),0.0);
    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (isNaN(my.at(iyt)))
            continue;

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, ny, 1.0, &dJydy.at(iyt*nJ*ny), nJ, dydp.data(), ny,
                    1.0, dJydp.data(), nJ);

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, ny, 1.0, &dJydsigma.at(iyt*nJ*ny), nJ,
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
void Model::fsJz(const int nroots, const std::vector<realtype> dJzdx, AmiVectorArray *sx, const ReturnData *rdata) {
    // sJz           nJ x rdata->nplist
    // dJzdp         nJ x rdata->nplist
    // dJzdx         nmaxevent x nJ x nx
    // sx            rdata->nt x nx x rdata->nplist

    // Compute dJzdx*sx for current 'ie'
    // dJzdx        rdata->nt x nJ x nx
    // sx           rdata->nt x nx x rdata->nplist

    std::vector<realtype> multResult(nJ * rdata->nplist, 0);
    std::vector<realtype> sxTmp(rdata->nplist * nx, 0);
    for (int ip = 0; ip < rdata->nplist; ++ip) {
        for (int ix = 0; ix < nx; ++ix)
            sxTmp.at(ix + ip * nx) = sx->at(ix,ip);
    }

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans, nJ,
                rdata->nplist, nx, 1.0, &dJzdx.at(nroots*nx*nJ), nJ, sxTmp.data(), nx, 1.0,
                multResult.data(), nJ);

    // sJy += multResult + dJydp
    for (int iJ = 0; iJ < nJ; ++iJ) {
        if (iJ == 0)
            for (int ip = 0; ip < rdata->nplist; ++ip)
                rdata->sllh[ip] -= static_cast<double>(multResult.at(ip * nJ) + dJzdp.at(ip * nJ));
        else
            for (int ip = 0; ip < rdata->nplist; ++ip)
                rdata->s2llh[(iJ - 1) + ip * (nJ - 1)] -=
                    static_cast<double>(multResult.at(iJ + ip * nJ) + dJzdp.at(iJ + ip * nJ));
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
    // dzdp          nz x rdata->nplist
    // dJzdp         nJ x rdata->nplist

    getmz(nroots,edata);
    std::fill(dJzdp.begin(),dJzdp.end(),0.0);
    for (int izt = 0; izt < nztrue; ++izt) {
        if (isNaN(mz.at(izt)))
            continue;

        if (t < gett(rdata->nt - 1,rdata)) {
            // with z
            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                        nJ, rdata->nplist, nz, 1.0, &dJzdz.at(izt*nz*nJ), nJ, dzdp.data(), nz,
                        1.0, dJzdp.data(), nJ);
        } else {
            // with rz
            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                        AMICI_BLAS_NoTrans, nJ, rdata->nplist, nz, 1.0,
                        &dJrzdsigma.at(izt*nz*nJ), nJ, dsigmazdp.data(), nz, 1.0,
                        dJzdp.data(), nJ);
            
            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                        nJ, rdata->nplist, nz, 1.0, &dJrzdz.at(izt*nz*nJ), nJ, dzdp.data(), nz,
                        1.0, dJzdp.data(), nJ);
        }

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, nz, 1.0, &dJzdsigma.at(izt*nz*nJ), nJ,
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
        
        if (t < gett(rdata->nt - 1,rdata)) {
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
 * @param udata pointer to UserData instance
 */
void Model::initialize(AmiVector *x, AmiVector *dx, const UserData *udata) {

    initializeStates(x, udata);
    
    fdx0(x, dx);
    
    initHeaviside(x,dx, udata);
    
}

/** initialization of initial states
 * @param x pointer to state variables
 * @param udata pointer to UserData instance
 */
void Model::initializeStates(AmiVector *x, const UserData *udata) {

    if (udata->getInitialStates().empty()) {
        fx0(x, udata);
    } else {
        for (int ix = 0; ix < nx; ix++) {
            (*x)[ix] = (realtype) udata->getInitialStates().at(ix);
        }
    }
}

/**
 * initHeaviside initialises the heaviside variables h at the intial time t0
 * heaviside variables activate/deactivate on event occurences
 * @param x pointer to state variables
 * @param dx pointer to time derivative of states (DAE only)
 * @param udata pointer to UserData instance
 */
void Model::initHeaviside(AmiVector *x, AmiVector *dx, const UserData *udata) {
    std::vector<realtype> rootvals(ne,0.0);
    froot(udata->t0(), x, dx, rootvals.data());
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
    
    /** Initial states
     * @param x pointer to state variables
     * @param udata pointer to UserData instance
     */
    void Model::fx0(AmiVector *x, const UserData *udata) {
        x->reset();
        fx0(x->data(),udata->t0(),p.data(),k.data());
    };

    /** Initial value for initial state sensitivities
     * @param sx pointer to state sensitivity variables
     * @param x pointer to state variables
     * @param udata pointer to UserData instance
     **/
    void Model::fsx0(AmiVectorArray *sx, const AmiVector *x, const UserData *udata) {
        sx->reset();
        for(int ip = 0; ip<plist.size(); ip++)
            fsx0(sx->data(ip),udata->t0(),x->data(),p.data(),k.data(),plist.at(ip));
    }
    
    /** Sensitivity of event timepoint, total derivative
     * @param t current timepoint
     * @param ie event index
     * @param x pointer to state variables
     * @param sx pointer to state sensitivity variables
     */
    void Model::fstau(const realtype t, const int ie, const AmiVector *x, const AmiVectorArray *sx) {
        std::fill(stau.begin(),stau.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++){
            fstau(&stau.at(ip),t,x->data(),p.data(),k.data(),h.data(),sx->data(ip),plist.at(ip),ie);
        }
    }
    
    /** Observables / measurements
     * @param it timepoint index
     * @param rdata pointer to return data instance
     */
    void Model::fy(int it, ReturnData *rdata) {
        getx(it,rdata);
        std::fill(y.begin(),y.end(),0.0);
        fy(y.data(),gett(it,rdata),x.data(),p.data(),k.data(),h.data());
        for(int iy = 0; iy < ny; iy++)
            rdata->y[it + rdata->nt*iy] = static_cast<double>(y.at(iy));
    }
    
    /** partial derivative of observables y w.r.t. model parameters p
     * @param it timepoint index
     * @param rdata pointer to return data instance
     */
    void Model::fdydp(const int it, ReturnData *rdata) {
        getx(it,rdata);
        std::fill(dydp.begin(),dydp.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++){
            fdydp(&dydp.at(ip*ny),gett(it,rdata),x.data(),p.data(),k.data(),h.data(),plist.at(ip));
        }
    }
    
    /** partial derivative of observables y w.r.t. state variables x
     const UserData *udata
     * @param it timepoint index
     * @param rdata pointer to return data instance
     */
    void Model::fdydx(const int it, ReturnData *rdata) {
        const realtype t = gett(it,rdata);
        getx(it,rdata);
        std::fill(dydx.begin(),dydx.end(),0.0);
        fdydx(dydx.data(),t,x.data(),p.data(),k.data(),h.data());
    }
    
    /** Event-resolved output
     * @param nroots number of events for event index
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param rdata pointer to return data instance
     */
    void Model::fz(const int nroots, const int ie, const realtype t, const AmiVector *x, ReturnData *rdata) {
        std::vector<realtype> zreturn(nz,0.0);
        fz(zreturn.data(),ie,t,x->data(),p.data(),k.data(),h.data());
        for(int iz = 0; iz < nz; iz++) {
            if (z2event[iz] - 1 == ie)
                rdata->z[nroots+rdata->nmaxevent*iz] = static_cast<double>(zreturn.at(iz));
        }
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
        std::vector<realtype> szreturn(nz,0.0);
        for(int ip = 0; ip < plist.size();  ip++ ){
            std::fill(szreturn.begin(), szreturn.end(), 0.0);
            fsz(szreturn.data(),ie,t,x->data(),p.data(),k.data(),h.data(),sx->data(ip),plist.at(ip));
            for(int iz = 0; iz < nz; iz++) {
                if (z2event[iz] - 1 == ie)
                    rdata->sz[nroots+rdata->nmaxevent*(ip*nz + iz)] = static_cast<double>(szreturn.at(iz));
            }
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
        std::vector<realtype> rzreturn(nz,0.0);
        frz(rzreturn.data(),ie,t,x->data(),p.data(),k.data(),h.data());
        for(int iz = 0; iz < nz; iz++) {
            if (z2event[iz] - 1 == ie)
                rdata->rz[nroots+rdata->nmaxevent*iz] = static_cast<double>(rzreturn.at(iz));
        }
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
        std::vector<realtype> srzreturn(nz,0.0);
        for(int ip = 0; ip < plist.size();  ip++ ){
            std::fill(srzreturn.begin(), srzreturn.end(), 0.0);
            fsrz(srzreturn.data(),ie,t,x->data(),p.data(),k.data(),h.data(),sx->data(ip),plist.at(ip));
            for(int iz = 0; iz < nz; iz++) {
                if (z2event[iz] - 1 == ie)
                    rdata->srz[nroots+rdata->nmaxevent*(ip*nz + iz)] = static_cast<double>(srzreturn.at(iz));
            }
        }
    }
    
    /** partial derivative of event-resolved output z w.r.t. to model parameters p
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void Model::fdzdp(const realtype t, const int ie, const AmiVector *x) {
        std::fill(dzdp.begin(),dzdp.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++){
            fdzdp(dzdp.data(),ie,t,x->data(),p.data(),k.data(),h.data(),plist.at(ip));
        }
    }
    
    /** partial derivative of event-resolved output z w.r.t. to model states x
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void Model::fdzdx(const realtype t, const int ie, const AmiVector *x) {
        std::fill(dzdx.begin(),dzdx.end(),0.0);
        fdzdx(dzdx.data(),ie,t,x->data(),p.data(),k.data(),h.data());
    }
    
    /** Sensitivity of event-resolved root output w.r.t. to model parameters p
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void Model::fdrzdp(const realtype t, const int ie, const AmiVector *x) {
        std::fill(drzdp.begin(),drzdp.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++){
            fdrzdp(drzdp.data(),ie,t,x->data(),p.data(),k.data(),h.data(),plist.at(ip));
        }
    }
    
    /** Sensitivity of event-resolved measurements rz w.r.t. to model states x
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void Model::fdrzdx(const realtype t, const int ie, const AmiVector *x) {
        std::fill(drzdx.begin(),drzdx.end(),0.0);
        fdrzdx(drzdx.data(),ie,t,x->data(),p.data(),k.data(),h.data());
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
        fdeltax(deltax.data(),t,x->data(),p.data(),k.data(),h.data(),ie,xdot->data(),xdot_old->data());
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
        for(int ip = 0; ip < plist.size(); ip++)
            fdeltasx(&deltasx.at(nx*ip),t,x->data(),p.data(),k.data(),h.data(),w.data(),
                          plist.at(ip),ie,xdot->data(),xdot_old->data(),sx->data(ip),&stau.at(ip));
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
        fdeltaxB(deltaxB.data(),t,x->data(),p.data(),k.data(),h.data(),ie,xdot->data(),xdot_old->data(),xB->data());
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
        for(int ip = 0; ip < plist.size(); ip++)
            fdeltaqB(deltaqB.data(),t,x->data(),p.data(),k.data(),h.data(),
                          plist.at(ip),ie,xdot->data(),xdot_old->data(),xB->data());
    }
    
    /** Standard deviation of measurements
     * @param it timepoint index
     * @param edata pointer to experimental data instance
     * @param rdata pointer to return data instance
     */
    void Model::fsigma_y(const int it, const ExpData *edata, ReturnData *rdata) {
        std::fill(sigmay.begin(),sigmay.end(),0.0);
        fsigma_y(sigmay.data(),gett(it,rdata),p.data(),k.data());
        for (int iy = 0; iy < nytrue; iy++) {
            /* extract the value for the standard deviation, if the data value
             is NaN, use
             the parameter value. Store this value in the return struct */
            if(edata){
                if (!isNaN(edata->sigmay[iy * rdata->nt + it])) {
                    sigmay.at(iy) = edata->sigmay[iy * rdata->nt + it];
                }
            }
            rdata->sigmay[iy * rdata->nt + it] = static_cast<double>(sigmay.at(iy));
        }
    }
    
    /** partial derivative of standard deviation of measurements w.r.t. model
     * @param it timepoint index
     * @param rdata pointer to return data instance
     */
    void Model::fdsigma_ydp(const int it, const ReturnData *rdata) {
        std::fill(dsigmaydp.begin(),dsigmaydp.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++)
            fdsigma_ydp(dsigmaydp.data(),gett(it,rdata),p.data(),k.data(),plist.at(ip));
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
        fsigma_z(sigmaz.data(),t,p.data(),k.data());
        for (int iz = 0; iz < nztrue; iz++) {
            if (z2event.at(iz) - 1 == ie) {
                if(edata) {
                    if (!isNaN(edata->sigmaz[nroots[ie]+rdata->nmaxevent*iz])) {
                        sigmaz.at(iz) = edata->sigmaz[nroots[ie]+rdata->nmaxevent*iz];
                    }
                }
                rdata->sigmaz[nroots[ie]+rdata->nmaxevent * iz] = static_cast<double>(sigmaz.at(iz));
            }
        }
    }
    
    /** Sensitivity of standard deviation of events measurements w.r.t. model parameters p
     * @param t current timepoint
     */
    void Model::fdsigma_zdp(const realtype t) {
        std::fill(dsigmazdp.begin(),dsigmazdp.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++)
            fdsigma_zdp(dsigmazdp.data(),t,p.data(),k.data(),plist.at(ip));
    }
    
    /** negative log-likelihood of measurements y
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
    void Model::fJy(const int it, ReturnData *rdata, const ExpData *edata) {
        std::vector<realtype> nllh(nJ,0.0);
        gety(it,rdata);
        getmy(it,edata);
        for(int iytrue = 0; iytrue < nytrue; iytrue++){
            if(!isNaN(my.at(iytrue))){
                std::fill(nllh.begin(),nllh.end(),0.0);
                fJy(nllh.data(),iytrue,p.data(),k.data(),y.data(),sigmay.data(),my.data());
                rdata->llh[0] -= static_cast<double>(nllh.at(0));
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
        getz(nroots,rdata);
        getmz(nroots,edata);
        for(int iztrue = 0; iztrue < nztrue; iztrue++){
            if(!isNaN(mz.at(iztrue))){
                std::fill(nllh.begin(),nllh.end(),0.0);
                fJz(nllh.data(),iztrue,p.data(),k.data(),z.data(),sigmaz.data(),mz.data());
                rdata->llh[0] -= static_cast<double>(nllh.at(0));
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
                fJrz(nllh.data(),iztrue,p.data(),k.data(),rz.data(),sigmaz.data());
                rdata->llh[0] -= static_cast<double>(nllh.at(0));
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
        gety(it,rdata);
        getmy(it,edata);
        std::fill(dJydy.begin(),dJydy.end(),0.0);
        for(int iytrue = 0; iytrue < nytrue; iytrue++){
            if(!isNaN(my.at(iytrue))){
                fdJydy(&dJydy.at(iytrue*ny*nJ),iytrue,p.data(),k.data(),y.data(),sigmay.data(),my.data());
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
        gety(it,rdata);
        getmy(it,edata);
        std::fill(dJydsigma.begin(),dJydsigma.end(),0.0);
        for(int iytrue = 0; iytrue < nytrue; iytrue++){
            if(!isNaN(my.at(iytrue))){
                fdJydsigma(&dJydsigma.at(iytrue*ny*nJ),iytrue,p.data(),k.data(),y.data(),sigmay.data(),my.data());
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
        getz(nroots,rdata);
        getmz(nroots,edata);
        std::fill(dJzdz.begin(),dJzdz.end(),0.0);
        for(int iztrue = 0; iztrue < nztrue; iztrue++){
            if(!isNaN(mz.at(iztrue))){
                fdJzdz(&dJzdz.at(iztrue*nz*nJ),iztrue,p.data(),k.data(),z.data(),sigmaz.data(),mz.data());
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
        getz(nroots,rdata);
        getmz(nroots,edata);
        std::fill(dJzdsigma.begin(),dJzdsigma.end(),0.0);
        for(int iztrue = 0; iztrue < nztrue; iztrue++){
            if(!isNaN(mz.at(iztrue))){
                fdJzdsigma(&dJzdsigma.at(iztrue*nz*nJ),iztrue,p.data(),k.data(),z.data(),sigmaz.data(),mz.data());
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
        getrz(nroots,rdata);
        getmz(nroots,edata);
        std::fill(dJrzdz.begin(),dJrzdz.end(),0.0);
        for(int iztrue = 0; iztrue < nztrue; iztrue++){
            if(!isNaN(mz.at(iztrue))){
                fdJrzdz(&dJrzdz.at(iztrue*nz*nJ),iztrue,p.data(),k.data(),rz.data(),sigmaz.data());
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
        getrz(nroots,rdata);
        std::fill(dJrzdsigma.begin(),dJrzdsigma.end(),0.0);
        for(int iztrue = 0; iztrue < nztrue; iztrue++){
            if(!isNaN(mz.at(iztrue))){
                fdJrzdsigma(&dJrzdsigma.at(iztrue*nz*nJ),iztrue,p.data(),k.data(),rz.data(),sigmaz.data());
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
        fw(w.data(),t,N_VGetArrayPointer(x),p.data(),k.data(),h.data());
    }
    
    /**
     * @brief Recurring terms in xdot, parameter derivative
     * @param t timepoint
     * @param x Vector with the states
     */
    void Model::fdwdp(const realtype t, const N_Vector x) {
        fw(t,x);
        std::fill(dwdp.begin(),dwdp.end(),0.0);
        fdwdp(dwdp.data(),t,N_VGetArrayPointer(x),p.data(),k.data(),h.data(),w.data());
    }
    
    /**
     * @brief Recurring terms in xdot, state derivative
     * @param t timepoint
     * @param x Vector with the states
     */
    void Model::fdwdx(const realtype t, const N_Vector x) {
        fw(t,x);
        std::fill(dwdx.begin(),dwdx.end(),0.0);
        fdwdx(dwdx.data(),t,N_VGetArrayPointer(x),p.data(),k.data(),h.data(),w.data());
    }
    
    /** create my slice at timepoint
     * @param it timepoint index
     * @param edata pointer to experimental data instance
     */
    void Model::getmy(const int it, const ExpData *edata) {
        if(edata) {
            for(int iytrue = 0; iytrue < nytrue; iytrue++){
                my.at(iytrue) = static_cast<realtype>(edata->my[it + edata->nt*iytrue]);
            }
        } else {
            for(int iytrue = 0; iytrue < nytrue; iytrue++){
                my.at(iytrue) = getNaN();
            }
        }
    }
    
    /** create y slice at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data instance
     */
    void Model::gety(const int it, const ReturnData *rdata) {
        for(int iy = 0; iy < ny; iy++){
            y.at(iy) = static_cast<realtype>(rdata->y[it + rdata->nt*iy]);
        }
    }
    
    /** create x slice at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data instance
     */
    void Model::getx(const int it, const ReturnData *rdata) {
        for(int ix = 0; ix < nx; ix++){
            x.at(ix) = static_cast<realtype>(rdata->x[it + rdata->nt*ix]);
        }
    }
    
    /** create sx slice at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data instance
     */
    void Model::getsx(const int it, const ReturnData *rdata) {
        for(int ip = 0; ip < rdata->nplist; ip++) {
            for(int ix = 0; ix < nx; ix++){
                sx.at(ip).at(ix) = static_cast<realtype>(rdata->sx[(ip * nx + ix) * rdata->nt + it]);
            }
        }
    }
    
    /** get current timepoint from index
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @return current timepoint
     */
    const realtype Model::gett(const int it, const ReturnData *rdata) const {
        return static_cast<realtype>(rdata->ts[it]);
    }
    
    /** create mz slice at event
     * @param nroots event occurence
     * @param edata pointer to experimental data instance
     */
    void Model::getmz(const int nroots, const ExpData *edata) {
        if(edata){
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                mz.at(iztrue) = static_cast<realtype>(edata->mz[nroots + edata->nmaxevent*iztrue]);
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
     */
    void Model::getz(const int nroots, const ReturnData *rdata) {
        for(int iz = 0; iz < nz; iz++){
            z.at(iz) = static_cast<realtype>(rdata->z[nroots+rdata->nmaxevent*iz]);
        }
    }
    
    /** create rz slice at event
     * @param nroots event occurence
     * @param rdata pointer to return data instance
     */
    void Model::getrz(const int nroots, const ReturnData *rdata) {
        for(int iz = 0; iz < nz; iz++){
            rz.at(iz) = static_cast<realtype>(rdata->rz[nroots+rdata->nmaxevent*iz]);
        }
    }
    

} // namespace amici
