#include "include/amici.h"
#include "include/amici_model.h"
#include <cstring>
#include <include/edata.h>
#include <include/rdata.h>

namespace amici {

/** Sensitivity of measurements y, total derivative
 * @param it timepoint index @type int
 * @param[in,out] rdata pointer to return data object @type ReturnData
 */
void Model::fsy(const int it, ReturnData *rdata) {
    // Compute sy = dydx * sx + dydp

    for (int ip = 0; ip < rdata->nplist; ++ip) {
        for (int iy = 0; iy < ny; ++iy)
            // copy dydp to sy
            rdata->sy[ip * rdata->nt * ny + iy * rdata->nt + it] =
                dydp.at(iy + ip * ny);

        getsx(it,rdata);

        // compute sy = 1.0*dydx*sx + 1.0*sy
        amici_dgemv(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, ny, nx, 1.0,
                    dydx.data(), ny, sx.at(ip).data(), 1, 1.0,
                    &rdata->sy[it + ip * rdata->nt * ny], rdata->nt);
    }
}

/** Sensitivity of z at final timepoint (ignores sensitivity of timepoint),
 * total derivative
 * @param ie event index @type int
 * @param[in,out] rdata pointer to return data object @type ReturnData
 */
void Model::fsz_tf(const int nroots, ReturnData *rdata) {
    // Compute sz = dzdx * sz + dzdp

    for (int ip = 0; ip < rdata->nplist; ++ip) {
        for (int iz = 0; iz < nz; ++iz)
            // copy dydp to sy
            rdata->sz[nroots + (iz + ip * nz) * rdata->nmaxevent] =
                0;
    }
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy, total
 * derivative
 * @param it timepoint index @type int
 * @param[in,out] rdata pointer to return data object @type ReturnData
 */
    void Model::fsJy(const int it, const std::vector<double> dJydx, ReturnData *rdata) {

    // Compute dJydx*sx for current 'it'
    // dJydx        rdata->nt x nJ x nx
    // sx           rdata->nt x nx x rdata->nplist
    std::vector<double> multResult(nJ * rdata->nplist, 0);
    std::vector<double> sxTmp(rdata->nplist * nx, 0);
    
    for (int ix = 0; ix < nx; ++ix) {
        for (int ip = 0; ip < rdata->nplist; ++ip)
            sxTmp.at(ix + ip * nx) = rdata->sx[it + (ix + ip * nx) * rdata->nt];
        for (int iJ = 0; iJ < nJ; ++iJ)
            dJydxTmp.at(iJ + ix * nJ) =
                dJydx.at(it + (iJ + ix * nJ) * rdata->nt);
    }

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans, nJ,
                rdata->nplist, nx, 1.0, dJydxTmp.data(), nJ, sxTmp.data(), nx, 0.0,
                multResult.data(), nJ);

    // multResult    nJ x rdata->nplist
    // dJydp         nJ x rdata->nplist
    // dJydxTmp      nJ x nx
    // sxTmp         nx x rdata->nplist

    // sJy += multResult + dJydp
    for (int iJ = 0; iJ < nJ; ++iJ) {
        if (iJ == 0)
            for (int ip = 0; ip < rdata->nplist; ++ip)
                rdata->sllh[ip] -= multResult.at(ip * nJ) + dJydp.at(ip * nJ);
        else
            for (int ip = 0; ip < rdata->nplist; ++ip)
                rdata->s2llh[(iJ - 1) + ip * (nJ - 1)] -=
                    multResult.at(iJ + ip * nJ) + dJydp.at(iJ + ip * nJ);
    }
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t.
 * parameters
 * @param it timepoint index @type int
 * @param edata pointer to experimental data object @type ExpData
 * @param rdata pointer to return data object @type ReturnData
 */
void Model::fdJydp(const int it, const ExpData *edata,
                  const ReturnData *rdata) {

    // dJydy         nytrue x nJ x ny
    // dydp          ny x rdata->nplist
    // dJydp         nJ x rdata->nplist
    getmy(it,edata);
    fill(dJydp.begin(),dJydp.end(),0.0);
    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (amiIsNaN(my.at(iyt)))
            continue;

        // copy current (iyt) dJydy and dJydsigma slices
        // dJydyTmp     nJ x ny
        // dJydsigmaTmp nJ x ny
        for (int iJ = 0; iJ < nJ; ++iJ) {
            for (int iy = 0; iy < ny; ++iy) {
                dJydyTmp.at(iJ + iy * nJ) =
                    dJydy.at(iyt + (iJ + iy * nJ) * nytrue);
                dJydsigmaTmp.at(iJ + iy * nJ) =
                    dJydsigma.at(iyt + (iJ + iy * nJ) * nytrue);
            }
        }

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, ny, 1.0, dJydyTmp.data(), nJ, dydp.data(), ny,
                    1.0, dJydp.data(), nJ);

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, ny, 1.0, dJydsigmaTmp.data(), nJ,
                    dsigmaydp.data(), ny, 1.0, dJydp.data(), nJ);
    }
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t.
 * state variables
 * @param it timepoint index @type int
 * @param edata pointer to experimental data object @type ExpData
 */
void Model::fdJydx(std::vector<double> dJydx, const int it, const ExpData *edata, const ReturnData *rdata) {

    // dJydy         nytrue x nJ x ny
    // dydx          ny x nx
    // dJydx         rdata->nt x nJ x nx
    getmy(it,edata);
    std::vector<double> multResult(nJ * nx, 0);
    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (amiIsNaN(my.at(iyt)))
            continue;

        // copy current (iyt) dJydy slice
        // dJydyTmp     nJ x ny
        for (int iJ = 0; iJ < nJ; ++iJ)
            for (int iy = 0; iy < ny; ++iy)
                dJydyTmp.at(iJ + iy * nJ) =
                    dJydy.at(iyt + (iJ + iy * nJ) * nytrue);

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, nx, ny, 1.0, dJydyTmp.data(), nJ, dydx.data(), ny, 1.0,
                    multResult.data(), nJ);
    }
    for (int iJ = 0; iJ < nJ; ++iJ)
        for (int ix = 0; ix < nx; ++ix)
            dJydx.at(it + (iJ + ix * nJ) * rdata->nt)=
                multResult.at(iJ + ix * nJ);
}

/** Sensitivity of event-resolved measurement negative log-likelihood Jz, total
 * derivative
 * @param ie event index @type int
 * @param rdata pointer to return data object @type ReturnData
 */
void Model::fsJz(const int nroots, const std::vector<double> dJzdx, AmiVectorArray *sx, const ReturnData *rdata) {
    // sJz           nJ x rdata->nplist
    // dJzdp         nJ x rdata->nplist
    // dJzdx         nmaxevent x nJ x nx
    // sx            rdata->nt x nx x rdata->nplist

    // Compute dJzdx*sx for current 'ie'
    // dJzdx        rdata->nt x nJ x nx
    // sx           rdata->nt x nx x rdata->nplist

    std::vector<double> multResult(nJ * rdata->nplist, 0);
    std::vector<double> sxTmp(rdata->nplist * nx, 0);
    for (int ip = 0; ip < rdata->nplist; ++ip) {
        for (int ix = 0; ix < nx; ++ix)
            sxTmp.at(ix + ip * nx) = sx->at(ix,ip);
    }

    for (int ix = 0; ix < nx; ++ix)
        for (int iJ = 0; iJ < nJ; ++iJ)
            dJzdxTmp.at(iJ + ix * nJ) =
                dJzdx.at(nroots +
                             (iJ + ix * nJ) * rdata->nmaxevent);

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans, nJ,
                rdata->nplist, nx, 1.0, dJzdxTmp.data(), nJ, sxTmp.data(), nx, 1.0,
                multResult.data(), nJ);

    // sJy += multResult + dJydp
    for (int iJ = 0; iJ < nJ; ++iJ) {
        if (iJ == 0)
            for (int ip = 0; ip < rdata->nplist; ++ip)
                rdata->sllh[ip] -= multResult.at(ip * nJ) + dJzdp.at(ip * nJ);
        else
            for (int ip = 0; ip < rdata->nplist; ++ip)
                rdata->s2llh[(iJ - 1) + ip * (nJ - 1)] -=
                    multResult.at(iJ + ip * nJ) + dJzdp.at(iJ + ip * nJ);
    }
}

/** Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t.
 * parameters
 * @param ie event index @type int
 * @param edata pointer to experimental data object @type ExpData
 * @param rdata pointer to return data object @type ReturnData
 */
void Model::fdJzdp(const int nroots, realtype t, const ExpData *edata,
                  const ReturnData *rdata) {
    // dJzdz         nztrue x nJ x nz
    // dJzdsigma     nztrue x nJ x nz
    // dzdp          nz x rdata->nplist
    // dJzdp         nJ x rdata->nplist

    getmz(nroots,edata);
    std::fill(dJzdp.begin(),dJzdp.end(),0.0);
    for (int izt = 0; izt < nztrue; ++izt) {
        if (amiIsNaN(mz.at(izt)))
            continue;

        // copy current (izt) dJzdz and dJzdsigma slices
        // dJzdzTmp     nJ x nz
        // dJzdsigmaTmp nJ x nz

        if (t < rdata->ts[rdata->nt - 1]) {
            for (int iJ = 0; iJ < nJ; ++iJ) {
                for (int iz = 0; iz < nz; ++iz) {
                    dJzdzTmp.at(iJ + iz * nJ) =
                        dJzdz.at(izt + (iJ + iz * nJ) * nztrue);
                }
            }
        } else {
            for (int iJ = 0; iJ < nJ; ++iJ) {
                for (int iz = 0; iz < nz; ++iz) {
                    dJzdzTmp.at(iJ + iz * nJ) =
                        dJrzdz.at(izt + (iJ + iz * nJ) * nztrue);
                    dJrzdsigmaTmp.at(iJ + iz * nJ) =
                        dJrzdsigma.at(izt + (iJ + iz * nJ) * nztrue);
                }
            }
            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                        AMICI_BLAS_NoTrans, nJ, rdata->nplist, nz, 1.0,
                        dJrzdsigmaTmp.data(), nJ, dsigmazdp.data(), nz, 1.0,
                        dJzdp.data(), nJ);
        }

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, nz, 1.0, dJzdzTmp.data(), nJ, dzdp.data(), nz,
                    1.0, dJzdp.data(), nJ);

        for (int iJ = 0; iJ < nJ; ++iJ) {
            for (int iz = 0; iz < nz; ++iz) {
                dJzdsigmaTmp.at(iJ + iz * nJ) =
                    dJzdsigma.at(izt + (iJ + iz * nJ) * nztrue);
            }
        }

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, nz, 1.0, dJzdsigmaTmp.data(), nJ,
                    dsigmazdp.data(), nz, 1.0, dJzdp.data(), nJ);
    }
}

/** Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t.
 * state variables
 * @param ie event index @type int
 * @param edata pointer to experimental data object @type ExpData
 */
void Model::fdJzdx(std::vector<double> dJzdx, const int nroots, realtype t, const ExpData *edata, const ReturnData *rdata) {
    // dJzdz         nztrue x nJ x nz
    // dzdx          nz x nx
    // dJzdx         nmaxevent x nJ x nx
    getmz(nroots,edata);
    std::vector<double> multResult(nJ * nx, 0);
    for (int izt = 0; izt < nztrue; ++izt) {
        if (amiIsNaN(mz.at(izt)))
            continue;

        // copy current (izt) dJzdz slice
        // dJzdzTmp     nJ x nz
        if (t < rdata->ts[rdata->nt - 1]) {
            for (int iJ = 0; iJ < nJ; ++iJ)
                for (int iz = 0; iz < nz; ++iz)
                    dJzdzTmp.at(iJ + iz * nJ) =
                        dJzdz.at(izt + (iJ + iz * nJ) * nztrue);

            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                        AMICI_BLAS_NoTrans, nJ, nx, nz, 1.0, dJzdzTmp.data(), nJ,
                        dzdx.data(), nz, 1.0, multResult.data(), nJ);
        } else {
            for (int iJ = 0; iJ < nJ; ++iJ) {
                for (int iz = 0; iz < nz; ++iz) {
                    dJzdzTmp.at(iJ + iz * nJ) =
                        dJrzdz.at(izt + (iJ + iz * nJ) * nztrue);
                }
            }

            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                        AMICI_BLAS_NoTrans, nJ, nx, nz, 1.0, dJzdzTmp.data(), nJ,
                        drzdx.data(), nz, 1.0, multResult.data(), nJ);
        }
    }
    for (int iJ = 0; iJ < nJ; ++iJ)
        for (int ix = 0; ix < nx; ++ix)
            dJzdx.at(nroots +
                         (iJ + ix * nJ) * edata->nmaxevent) =
                multResult.at(iJ + ix * nJ);
}

/** initialization of model properties
 */
void Model::initialize(AmiVector *x, AmiVector *dx, const UserData *udata) {

    initializeStates(x, udata);
    
    fdx0(x, dx);
    
    initHeaviside(x,dx, udata);
    
}

/** initialization of initial states
 * @param x0data array with initial state values @type double
 */
void Model::initializeStates(AmiVector *x, const UserData *udata) {

    if (udata->getInitialStates().empty()) {
        fx0(x, udata);
    } else {
        for (int ix = 0; ix < nx; ix++) {
            x[ix] = (realtype) udata->getInitialStates().at(ix);
        }
    }
}

/**
 * initHeaviside initialises the heaviside variables h at the intial time t0
 * heaviside variables activate/deactivate on event occurences
 *
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
     **/
    void Model::fx0(AmiVector *x, const UserData *udata) {
        x->reset();
        model_x0(x->data(),udata->t0(),p.data(),k.data());
    };

    /** Initial value for initial state sensitivities
     **/
    void Model::fsx0(AmiVectorArray *sx, const AmiVector *x, const UserData *udata) {
        sx->reset();
        for(int ip = 0; ip<plist.size(); ip++)
            model_sx0(sx->data(ip),udata->t0(),x->data(),p.data(),k.data(),plist.at(ip));
    }
    
    /** Sensitivity of event timepoint, total derivative
     * @param ie event index
     */
    void Model::fstau(const int ie,const realtype t, const AmiVector *x, const AmiVectorArray *sx) {
        std::fill(stau.begin(),stau.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++){
            model_stau(stau.data(),t,x->data(),p.data(),k.data(),h.data(),sx->data(ip),plist.at(ip),ie);
        }
    }
    
    /** Observables / measurements
     * @param it timepoint index
     * @param rdata pointer to return data object
     */
    void Model::fy(int it, ReturnData *rdata) {
        getx(it,rdata);
        std::fill(y.begin(),y.end(),0.0);
        model_y(y.data(),gett(it,rdata),x.data(),p.data(),k.data(),h.data());
        for(int iy = 0; iy < ny; iy++)
            rdata->y[it + rdata->nt*iy] = y.at(iy);
    }
    
    /** partial derivative of observables y w.r.t. model parameters p
     */
    void Model::fdydp(const int it, ReturnData *rdata) {
        getx(it,rdata);
        std::fill(dydp.begin(),dydp.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++){
            model_dydp(&dydp.at(ip*ny),gett(it,rdata),x.data(),p.data(),k.data(),h.data(),plist.at(ip));
        }
    }
    
    /** partial derivative of observables y w.r.t. state variables x
     const UserData *udata
     */
    void Model::fdydx(const int it, ReturnData *rdata) {
        const realtype t = gett(it,rdata);
        getx(it,rdata);
        std::fill(dydx.begin(),dydx.end(),0.0);
        model_dydx(dydx.data(),t,x.data(),p.data(),k.data(),h.data());
    }
    
    /** Event-resolved output
     * @param nroots number of events for event index
     * @param rdata pointer to return data object
     */
    void Model::fz(const int nroots, const int ie, const realtype t, const AmiVector *x, ReturnData *rdata) {
        std::vector<double> zreturn(nz,0.0);
        model_z(zreturn.data(),ie,t,x->data(),p.data(),k.data(),h.data());
        for(int iz = 0; iz < nz; iz++) {
            rdata->z[nroots+rdata->nmaxevent*iz] = zreturn.at(iz);
        }
    }
    
    /** Sensitivity of z, total derivative
     * @param nroots number of events for event index
     * @param rdata pointer to return data object
     */
    void Model::fsz(const int nroots, const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata) {
        for(int ip = 0; ip < plist.size();  ip++ ){
            std::vector<double> szreturn(nz,0.0);
            model_sz(szreturn.data(),ie,t,x->data(),p.data(),k.data(),h.data(),sx->data(ip),plist.at(ip));
            for(int iz = 0; iz < nz; iz++) {
                rdata->sz[nroots+rdata->nmaxevent*(ip*nz + iz)] = szreturn.at(iz);
            }
        }
    }
    
    /** Event root function of events (equal to froot but does not include
     * non-output events)
     * @param nroots number of events for event index
     * @param rdata pointer to return data object
     */
    void Model::frz(const int nroots, const int ie, const realtype t, const AmiVector *x, ReturnData *rdata) {
        std::vector<double> rzreturn(nz,0.0);
        model_rz(rzreturn.data(),ie,t,x->data(),p.data(),k.data(),h.data());
        for(int iz = 0; iz < nz; iz++) {
            rdata->rz[nroots+rdata->nmaxevent*iz] = rzreturn.at(iz);
        }
    }
    
    /** Sensitivity of rz, total derivative
     * @param nroots number of events for event index
     * @param rdata pointer to return data object
     */
    void Model::fsrz(const int nroots, const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata) {
        for(int ip = 0; ip < plist.size();  ip++ ){
            std::vector<double> srzreturn(nz,0.0);
            model_srz(srzreturn.data(),ie,t,x->data(),p.data(),k.data(),h.data(),sx->data(ip),plist.at(ip));
            for(int iz = 0; iz < nz; iz++) {
                rdata->srz[nroots+rdata->nmaxevent*(ip*nz + iz)] = srzreturn.at(iz);
            }
        }
    }
    
    /** partial derivative of event-resolved output z w.r.t. to model parameters p
     */
    void Model::fdzdp(const realtype t, const int ie, const AmiVector *x) {
        std::fill(dzdp.begin(),dzdp.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++){
            model_dzdp(dzdp.data(),ie,t,x->data(),p.data(),k.data(),h.data(),plist.at(ip));
        }
    }
    
    /** partial derivative of event-resolved output z w.r.t. to model states x
     */
    void Model::fdzdx(const realtype t, const int ie, const AmiVector *x) {
        std::fill(dzdx.begin(),dzdx.end(),0.0);
        model_dzdx(dzdx.data(),ie,t,x->data(),p.data(),k.data(),h.data());
    }
    
    /** Sensitivity of event-resolved root output w.r.t. to model parameters p
     */
    void Model::fdrzdp(const realtype t, const int ie, const AmiVector *x) {
        std::fill(drzdp.begin(),drzdp.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++){
            model_drzdp(drzdp.data(),ie,t,x->data(),p.data(),k.data(),h.data(),plist.at(ip));
        }
    }
    
    /** Sensitivity of event-resolved measurements rz w.r.t. to model states x
     */
    void Model::fdrzdx(const realtype t, const int ie, const AmiVector *x) {
        std::fill(drzdx.begin(),drzdx.end(),0.0);
        model_drzdx(drzdx.data(),ie,t,x->data(),p.data(),k.data(),h.data());
    }
    
    /** State update functions for events
     * @param ie event index
     */
    void Model::fdeltax(const int ie, const realtype t, const AmiVector *x,
                         const AmiVector *xdot, const AmiVector *xdot_old) {
        std::fill(deltax.begin(),deltax.end(),0.0);
        model_deltax(deltax.data(),t,x->data(),p.data(),k.data(),h.data(),ie,xdot->data(),xdot_old->data());
    }
    
    /** Sensitivity update functions for events, total derivative
     * @param ie event index
     */
    void Model::fdeltasx(const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx,
                          const AmiVector *xdot, const AmiVector *xdot_old) {
        fw(t,x->getNVector());
        std::fill(deltasx.begin(),deltasx.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++)
            model_deltasx(&deltasx.at(nx*ip),t,x->data(),p.data(),k.data(),h.data(),w.data(),
                          plist.at(ip),ie,xdot->data(),xdot_old->data(),sx->data(ip),stau.data());
    }
    
    /** Adjoint state update functions for events
     * @param ie event index
     */
    void Model::fdeltaxB(const int ie, const realtype t, const AmiVector *x, const AmiVector *xB,
                          const AmiVector *xdot, const AmiVector *xdot_old) {
        std::fill(deltaxB.begin(),deltaxB.end(),0.0);
        model_deltaxB(deltaxB.data(),t,x->data(),p.data(),k.data(),h.data(),ie,xdot->data(),xdot_old->data(),xB->data());
    }
    
    /** Quadrature state update functions for events
     * @param ie event index
     */
    void Model::fdeltaqB(const int ie, const realtype t, const AmiVector *x, const AmiVector *xB,
                          const AmiVector *xdot, const AmiVector *xdot_old) {
        std::fill(deltaqB.begin(),deltaqB.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++)
            model_deltaqB(deltaqB.data(),t,x->data(),p.data(),k.data(),h.data(),
                          plist.at(ip),ie,xdot->data(),xdot_old->data(),xB->data());
    }
    
    /** Standard deviation of measurements
     */
    void Model::fsigma_y(const int it, const ExpData *edata, ReturnData *rdata) {
        std::fill(sigmay.begin(),sigmay.end(),0.0);
        model_sigma_y(sigmay.data(),gett(it,rdata),p.data(),k.data());
        for (int iy = 0; iy < nytrue; iy++) {
            /* extract the value for the standard deviation, if the data value
             is NaN, use
             the parameter value. Store this value in the return struct */
            if(edata){
                if (!amiIsNaN(edata->sigmay[iy * rdata->nt + it])) {
                    sigmay.at(iy) = edata->sigmay[iy * rdata->nt + it];
                }
            }
            rdata->sigmay[iy * rdata->nt + it] = sigmay.at(iy);
        }
    }
    
    /** partial derivative of standard deviation of measurements w.r.t. model
     */
    void Model::fdsigma_ydp(const int it, const ReturnData *rdata) {
        std::fill(dsigmaydp.begin(),dsigmaydp.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++)
            model_dsigma_ydp(dsigmaydp.data(),gett(it,rdata),p.data(),k.data(),plist.at(ip));
    }
    
    /** Standard deviation of events
     */
    void Model::fsigma_z(const realtype t, const int ie, const int *nroots,
                         const ExpData *edata, ReturnData *rdata) {
        std::fill(sigmaz.begin(),sigmaz.end(),0.0);
        model_sigma_z(sigmaz.data(),t,p.data(),k.data());
        for (int iz = 0; iz < nztrue; iz++) {
            if (z2event.at(iz) - 1 == ie) {
                if(edata) {
                    if (!amiIsNaN(edata->sigmaz[nroots[ie]+rdata->nmaxevent*iz])) {
                        sigmaz.at(iz) = edata->sigmaz[nroots[ie]+rdata->nmaxevent*iz];
                    }
                }
                rdata->sigmaz[nroots[ie]+rdata->nmaxevent * iz] = sigmaz.at(iz);
            }
        }
    }
    
    /** Sensitivity of standard deviation of events measurements w.r.t. model parameters p
     */
    void Model::fdsigma_zdp(const realtype t) {
        std::fill(dsigmazdp.begin(),dsigmazdp.end(),0.0);
        for(int ip = 0; ip < plist.size(); ip++)
            model_dsigma_zdp(dsigmazdp.data(),t,p.data(),k.data(),plist.at(ip));
    }
    
    /** negative log-likelihood of measurements y
     * @param it timepoint index
     * @param rdata pointer to return data object
     * @param edata pointer to experimental data object
     */
    void Model::fJy(const int it, ReturnData *rdata, const ExpData *edata) {
        std::vector<double> nllh(nJ,0.0);
        gety(it,rdata);
        getmy(it,edata);
        for(int iytrue = 0; iytrue < nytrue; iytrue++){
            if(!amiIsNaN(my.at(iytrue))){
                std::fill(nllh.begin(),nllh.end(),0.0);
                model_Jy(nllh.data(),iytrue,p.data(),k.data(),y.data(),sigmay.data(),my.data());
                for(int iJ = 0; iJ < nJ; iJ++){
                    rdata->llh[iJ] -= nllh.at(iJ);
                }
            }
        }
    }
    
    /** negative log-likelihood of event-resolved measurements z
     * @param nroots event index
     * @param rdata pointer to return data object
     * @param edata pointer to experimental data object
     */
    void Model::fJz(const int nroots, ReturnData *rdata, const ExpData *edata) {
        std::vector<double> nllh(nJ,0.0);
        getz(nroots,rdata);
        getmz(nroots,edata);
        for(int iztrue = 0; iztrue < nztrue; iztrue++){
            if(!amiIsNaN(mz.at(nroots))){
                std::fill(nllh.begin(),nllh.end(),0.0);
                model_Jz(nllh.data(),iztrue,p.data(),k.data(),z.data(),sigmaz.data(),mz.data());
                for(int iJ = 0; iJ < nJ; iJ++){
                    rdata->llh[iJ] -= nllh.at(iJ);
                }
            }
        }
    }
    
    /** regularization of negative log-likelihood with roots of event-resolved
     * measurements rz
     * @param nroots event index
     * @param rdata pointer to return data object
     * @param edata pointer to experimental data object
     */
    void Model::fJrz(const int nroots, ReturnData *rdata, const ExpData *edata) {
        std::vector<double> nllh(nJ,0.0);
        getrz(nroots,rdata);
        for(int iztrue = 0; iztrue < nztrue; iztrue++){
            if(!amiIsNaN(mz.at(nroots))){
                std::fill(nllh.begin(),nllh.end(),0.0);
                model_Jrz(nllh.data(),iztrue,p.data(),k.data(),rz.data(),sigmaz.data());
                for(int iJ = 0; iJ < nJ; iJ++){
                    rdata->llh[iJ] -= nllh.at(iJ);
                }
            }
        }
    }
    
    /** partial derivative of time-resolved measurement negative log-likelihood Jy
     * @param it timepoint index
     * @param rdata pointer to return data object
     * @param edata pointer to experimental data object
     */
    void Model::fdJydy(const int it, const ReturnData *rdata,
                        const ExpData *edata) {
        gety(it,rdata);
        getmy(it,edata);
        std::vector<double> dJydy_slice(nJ*nytrue, 0.0);
        std::fill(dJydy.begin(),dJydy.end(),0.0);
        for(int iytrue = 0; iytrue < nytrue; iytrue++){
            if(!amiIsNaN(my.at(iytrue))){
                std::fill(dJydy_slice.begin(),dJydy_slice.end(),0.0);
                model_dJydy(dJydy_slice.data(),iytrue,p.data(),k.data(),y.data(),sigmay.data(),my.data());
                // TODO: fix slicing here such that slicing is no longer necessary in sJy
                for(int iJ = 0; iJ < nJ; iJ++){
                    for(int iy = 0; iy < ny; iy++ )
                        dJydy.at(iytrue+(iJ+iy*nJ)) = dJydy_slice.at(iJ+iy*nJ);
                }
            }
        }
    }
    
    /** Sensitivity of time-resolved measurement negative log-likelihood Jy
     * w.r.t. standard deviation sigma
     * @param it timepoint index
     * @param rdata pointer to return data object
     * @param edata pointer to experimental data object
     */
    void Model::fdJydsigma(const int it, const ReturnData *rdata,
                            const ExpData *edata) {
        gety(it,rdata);
        getmy(it,edata);
        std::vector<double> dJydsigma_slice(nJ*nytrue, 0.0);
        std::fill(dJydsigma.begin(),dJydsigma.end(),0.0);
        for(int iytrue = 0; iytrue < nytrue; iytrue++){
            if(!amiIsNaN(my.at(iytrue))){
                std::fill(dJydsigma_slice.begin(),dJydsigma_slice.end(),0.0);
                model_dJydsigma(dJydsigma_slice.data(),iytrue,p.data(),k.data(),y.data(),sigmay.data(),my.data());
                // TODO: fix slicing here such that slicing is no longer necessary in sJy
                for(int iJ = 0; iJ < nJ; iJ++){
                    for(int iy = 0; iy < ny; iy++ )
                        dJydsigma.at(iytrue+(iJ+iy*nJ)) = dJydsigma_slice.at(iJ+iy*nJ);
                }
            }
        }
    }
    
    /** partial derivative of event measurement negative log-likelihood Jz
     * @param nroots event index
     * @param rdata pointer to return data object
     * @param edata pointer to experimental data object
     */
    void Model::fdJzdz(const int nroots, const ReturnData *rdata,
                        const ExpData *edata) {
        getz(nroots,rdata);
        getmz(nroots,edata);
        std::vector<double> dJzdz_slice(nJ*nztrue, 0.0);
        std::fill(dJzdz.begin(),dJzdz.end(),0.0);
        for(int iztrue = 0; iztrue < nztrue; iztrue++){
            if(!amiIsNaN(mz.at(iztrue))){
                std::fill(dJzdz_slice.begin(),dJzdz_slice.end(),0.0);
                model_dJzdz(dJzdz_slice.data(),iztrue,p.data(),k.data(),z.data(),sigmaz.data(),mz.data());
                // TODO: fix slicing here such that slicing is no longer necessary in sJz
                for(int iJ = 0; iJ < nJ; iJ++){
                    for(int iz = 0; iz < nz; iz++ )
                        dJzdz.at(iztrue+(iJ+iz*nJ)) = dJzdz_slice.at(iJ+iz*nJ);
                }
            }
        }
    }
    
    /** Sensitivity of event measurement negative log-likelihood Jz
     * w.r.t. standard deviation sigmaz
     * @param nroots event index
     * @param rdata pointer to return data object
     * @param edata pointer to experimental data object
     */
    void Model::fdJzdsigma(const int nroots, const ReturnData *rdata,
                            const ExpData *edata) {
        getz(nroots,rdata);
        getmz(nroots,edata);
        std::vector<double> dJzdsigma_slice(nJ*nztrue, 0.0);
        std::fill(dJzdsigma.begin(),dJzdsigma.end(),0.0);
        for(int iztrue = 0; iztrue < nztrue; iztrue++){
            if(!amiIsNaN(mz.at(iztrue))){
                std::fill(dJzdsigma_slice.begin(),dJzdsigma_slice.end(),0.0);
                model_dJzdsigma(dJzdsigma_slice.data(),iztrue,p.data(),k.data(),z.data(),sigmaz.data(),mz.data());
                // TODO: fix slicing here such that slicing is no longer necessary in sJz
                for(int iJ = 0; iJ < nJ; iJ++){
                    for(int iz = 0; iz < nz; iz++ )
                        dJzdsigma.at(iztrue+(iJ+iz*nJ)) = dJzdsigma_slice.at(iJ+iz*nJ);
                }
            }
        }
    }
    
    /** partial derivative of event measurement negative log-likelihood Jz
     * @param nroots event index
     * @param rdata pointer to return data object
     * @param edata pointer to experimental data object
     */
    void Model::fdJrzdz(const int nroots, const ReturnData *rdata,
                         const ExpData *edata) {
        getrz(nroots,rdata);
        getmz(nroots,edata);
        std::vector<double> dJrzdz_slice(nJ*nztrue, 0.0);
        std::fill(dJrzdz.begin(),dJrzdz.end(),0.0);
        for(int iztrue = 0; iztrue < nztrue; iztrue++){
            if(!amiIsNaN(mz.at(iztrue))){
                std::fill(dJrzdz_slice.begin(),dJrzdz_slice.end(),0.0);
                model_dJrzdz(dJrzdz_slice.data(),iztrue,p.data(),k.data(),rz.data(),sigmaz.data());
                // TODO: fix slicing here such that slicing is no longer necessary in sJz
                for(int iJ = 0; iJ < nJ; iJ++){
                    for(int iz = 0; iz < nz; iz++)
                        dJrzdz.at(iztrue+(iJ+iz*nJ)) = dJrzdz_slice.at(iJ+iz*nJ);
                }
            }
        }
    }
    
    /** Sensitivity of event measurement negative log-likelihood Jz
     * w.r.t. standard deviation sigmaz
     * @param nroots event index
     * @param rdata pointer to return data object
     * @param edata pointer to experimental data object
     */
    void Model::fdJrzdsigma(const int nroots,const ReturnData *rdata,
                             const ExpData *edata) {
        getrz(nroots,rdata);
        std::vector<double> dJrzdsigma_slice(nJ*nztrue, 0.0);
        std::fill(dJrzdsigma.begin(),dJrzdsigma.end(),0.0);
        for(int iztrue = 0; iztrue < nztrue; iztrue++){
            if(!amiIsNaN(mz.at(iztrue))){
                std::fill(dJrzdsigma_slice.begin(),dJrzdsigma_slice.end(),0.0);
                model_dJrzdsigma(dJrzdsigma_slice.data(),iztrue,p.data(),k.data(),rz.data(),sigmaz.data());
                // TODO: fix slicing here such that slicing is no longer necessary in sJz
                for(int iJ = 0; iJ < nJ; iJ++){
                    for(int iz = 0; iz < nz; iz++ )
                        dJrzdsigma.at(iztrue+(iJ+iz*nJ)) = dJrzdsigma_slice.at(iJ+iz*nJ);
                }
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
        model_w(w.data(),t,N_VGetArrayPointer(x),p.data(),k.data(),h.data());
    }
    
    /**
     * @brief Recurring terms in xdot, parameter derivative
     * @param t timepoint
     * @param x Vector with the states
     */
    void Model::fdwdp(const realtype t, const N_Vector x) {
        fw(t,x);
        std::fill(dwdp.begin(),dwdp.end(),0.0);
        model_dwdp(dwdp.data(),t,N_VGetArrayPointer(x),p.data(),k.data(),h.data(),w.data());
    }
    
    /**
     * @brief Recurring terms in xdot, state derivative
     * @param t timepoint
     * @param x Vector with the states @type N_Vector
     */
    void Model::fdwdx(const realtype t, const N_Vector x) {
        fw(t,x);
        std::fill(dwdx.begin(),dwdx.end(),0.0);
        model_dwdx(dwdx.data(),t,N_VGetArrayPointer(x),p.data(),k.data(),h.data(),w.data());
    }
    
    /** create my slice at timepoint
     * @param it timepoint index
     * @param edata pointer to experimental data object
     */
    void Model::getmy(const int it, const ExpData *edata) {
        if(edata) {
            for(int iytrue = 0; iytrue < nytrue; iytrue++){
                my.at(iytrue) = edata->my[it + edata->nt*iytrue];
            }
        } else {
            for(int iytrue = 0; iytrue < nytrue; iytrue++){
                my.at(iytrue) = amiGetNaN();
            }
        }
    }
    
    /** create y slice at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data object
     */
    void Model::gety(const int it, const ReturnData *rdata) {
        for(int iy = 0; iy < ny; iy++){
            y.at(iy) = rdata->y[it + rdata->nt*iy];
        }
    }
    
    /** create x slice at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data object
     */
    void Model::getx(const int it, const ReturnData *rdata) {
        for(int ix = 0; ix < nx; ix++){
            x.at(ix) = rdata->x[it + rdata->nt*ix];
        }
    }
    
    /** create sx slice at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data object
     */
    void Model::getsx(const int it, const ReturnData *rdata) {
        for(int ip = 0; ip < rdata->nplist; ip++) {
            for(int ix = 0; ix < nx; ix++){
                sx.at(ip).at(ix) = rdata->sx[it + rdata->nt*(ip+rdata->nplist*ix) ];
            }
        }
    }
    
    /** create t  at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data object
     */
    const realtype Model::gett(const int it, const ReturnData *rdata) const {
        return rdata->ts[it];
    }
    
    /** create mz slice at event
     * @param nroots event occurence
     * @param edata pointer to experimental data object
     */
    void Model::getmz(const int nroots, const ExpData *edata) {
        if(edata){
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                mz.at(iztrue) = edata->mz[nroots + edata->nmaxevent*iztrue];
            }
        } else {
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                mz.at(iztrue) = amiGetNaN();
            }
        }
    }
    
    /** create z slice at event
     * @param nroots event occurence
     * @param rdata pointer to return data object
     */
    void Model::getz(const int nroots, const ReturnData *rdata) {
        for(int iz = 0; iz < nz; iz++){
            z.at(iz) = rdata->z[nroots+rdata->nmaxevent*iz];
        }
    }
    
    /** create rz slice at event
     * @param nroots event occurence
     * @param rdata pointer to return data object
     */
    void Model::getrz(const int nroots, const ReturnData *rdata) {
        for(int iz = 0; iz < nz; iz++){
            rz.at(iz) = rdata->rz[nroots+rdata->nmaxevent*iz];
        }
    }
    

} // namespace amici
