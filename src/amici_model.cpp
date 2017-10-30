#include "include/amici_model.h"
#include <cstring>
#include <include/edata.h>
#include <include/rdata.h>
#include <include/tdata.h>
#include <include/udata.h>

namespace amici {


UserData Model::getUserData() const { return UserData(np, nk, nx); }

UserData *Model::getNewUserData() const { return new UserData(np, nk, nx); }

Model::~Model() {
    if (z2event)
        delete[] z2event;

    if (idlist)
        delete[] idlist;
}

/** Sensitivity of measurements y, total derivative
 * @param[in] it timepoint index @type int
 * @param[in] tdata pointer to temp data object @type TempData
 * @param[in,out] rdata pointer to return data object @type ReturnData
 */
void Model::fsy(const int it, const TempData *tdata, ReturnData *rdata) {
    // Compute sy = dydx * sx + dydp

    for (int ip = 0; ip < rdata->nplist; ++ip) {
        for (int iy = 0; iy < ny; ++iy)
            // copy dydp to sy
            rdata->sy[ip * rdata->nt * ny + iy * rdata->nt + it] =
                tdata->dydp[iy + ip * ny];

        realtype *sx_tmp = N_VGetArrayPointer(tdata->sx[ip]);

        // compute sy = 1.0*dydx*sx + 1.0*sy
        amici_dgemv(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, ny, nx, 1.0,
                    tdata->dydx, ny, sx_tmp, 1, 1.0,
                    &rdata->sy[it + ip * rdata->nt * ny], rdata->nt);
    }

    return;
}

/** Sensitivity of z at final timepoint (ignores sensitivity of timepoint),
 * total derivative
 * @param[in] ie event index @type int
 * @param[in] tdata pointer to temp data object @type TempData
 * @param[in,out] rdata pointer to return data object @type ReturnData
 */
void Model::fsz_tf(const int ie, const TempData *tdata, ReturnData *rdata) {
    // Compute sz = dzdx * sz + dzdp

    for (int ip = 0; ip < rdata->nplist; ++ip) {
        for (int iz = 0; iz < nz; ++iz)
            // copy dydp to sy
            rdata->sz[tdata->nroots[ie] + (iz + ip * nz) * rdata->nmaxevent] =
                0;
    }

    return;
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy, total
 * derivative
 * @param[in] it timepoint index @type int
 * @param[in] tdata pointer to temp data object @type TempData
 * @param[in,out] rdata pointer to return data object @type ReturnData
 */
void Model::fsJy(const int it, const TempData *tdata, ReturnData *rdata) {

    // Compute dJydx*sx for current 'it'
    // dJydx        rdata->nt x nJ x nx
    // sx           rdata->nt x nx x rdata->nplist
    std::vector<double> multResult;
    multResult.resize(nJ * rdata->nplist, 0);
    std::vector<double> sxTmp;
    sxTmp.resize(rdata->nplist * nx, 0);
    for (int ix = 0; ix < nx; ++ix) {
        for (int ip = 0; ip < rdata->nplist; ++ip)
            sxTmp[ix + ip * nx] = rdata->sx[it + (ix + ip * nx) * rdata->nt];
        for (int iJ = 0; iJ < nJ; ++iJ)
            dJydxTmp[iJ + ix * nJ] =
                tdata->dJydx[it + (iJ + ix * nJ) * rdata->nt];
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
                rdata->sllh[ip] -= multResult[ip * nJ] + tdata->dJydp[ip * nJ];
        else
            for (int ip = 0; ip < rdata->nplist; ++ip)
                rdata->s2llh[(iJ - 1) + ip * (nJ - 1)] -=
                    multResult[iJ + ip * nJ] + tdata->dJydp[iJ + ip * nJ];
    }
    return;
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t.
 * parameters
 * @param[in] it timepoint index @type int
 * @param[in,out] tdata pointer to temp data object @type TempData
 * @param[in] edata pointer to experimental data object @type ExpData
 * @param[in] rdata pointer to return data object @type ReturnData
 */
void Model::fdJydp(const int it, TempData *tdata, const ExpData *edata,
                  const ReturnData *rdata) {

    // dJydy         nytrue x nJ x ny
    // dydp          ny x rdata->nplist
    // dJydp         nJ x rdata->nplist
    memset(tdata->dJydp, 0, nJ * rdata->nplist * sizeof(double));
    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (amiIsNaN(edata->my[rdata->nt * iyt + it]))
            continue;

        // copy current (iyt) dJydy and dJydsigma slices
        // dJydyTmp     nJ x ny
        // dJydsigmaTmp nJ x ny
        for (int iJ = 0; iJ < nJ; ++iJ) {
            for (int iy = 0; iy < ny; ++iy) {
                dJydyTmp[iJ + iy * nJ] =
                    tdata->dJydy[iyt + (iJ + iy * nJ) * nytrue];
                dJydsigmaTmp[iJ + iy * nJ] =
                    tdata->dJydsigma[iyt + (iJ + iy * nJ) * nytrue];
            }
        }

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, ny, 1.0, dJydyTmp.data(), nJ, tdata->dydp, ny,
                    1.0, tdata->dJydp, nJ);

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, ny, 1.0, dJydsigmaTmp.data(), nJ,
                    tdata->dsigmaydp, ny, 1.0, tdata->dJydp, nJ);
    }
    return;
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t.
 * state variables
 * @param[in] it timepoint index @type int
 * @param[in,out] tdata pointer to temp data object @type TempData
 * @param[in] edata pointer to experimental data object @type ExpData
 */
void Model::fdJydx(const int it, TempData *tdata, const ExpData *edata) {

    // dJydy         nytrue x nJ x ny
    // dydx          ny x nx
    // dJydx         rdata->nt x nJ x nx
    
    std::vector<double> multResult;
    multResult.resize(nJ * nx, 0);
    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (amiIsNaN(edata->my[tdata->rdata->nt * iyt + it]))
            continue;

        // copy current (iyt) dJydy slice
        // dJydyTmp     nJ x ny
        for (int iJ = 0; iJ < nJ; ++iJ)
            for (int iy = 0; iy < ny; ++iy)
                dJydyTmp[iJ + iy * nJ] =
                    tdata->dJydy[iyt + (iJ + iy * nJ) * nytrue];

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, nx, ny, 1.0, dJydyTmp.data(), nJ, tdata->dydx, ny, 1.0,
                    multResult.data(), nJ);
    }
    for (int iJ = 0; iJ < nJ; ++iJ)
        for (int ix = 0; ix < nx; ++ix)
            tdata->dJydx[it + (iJ + ix * nJ) * tdata->rdata->nt] =
                multResult[iJ + ix * nJ];
    return;
}

/** Sensitivity of event-resolved measurement negative log-likelihood Jz, total
 * derivative
 * @param[in] ie event index @type int
 * @param[in,out] tdata pointer to temp data object @type TempData
 * @param[in] rdata pointer to return data object @type ReturnData
 */
void Model::fsJz(const int ie, TempData *tdata, const ReturnData *rdata) {
    // sJz           nJ x rdata->nplist
    // dJzdp         nJ x rdata->nplist
    // dJzdx         nmaxevent x nJ x nx
    // sx            rdata->nt x nx x rdata->nplist

    // Compute dJzdx*sx for current 'ie'
    // dJzdx        rdata->nt x nJ x nx
    // sx           rdata->nt x nx x rdata->nplist

    std::vector<double> multResult;
    multResult.resize(nJ * rdata->nplist, 0);
    std::vector<double> sxTmp;
    sxTmp.resize(rdata->nplist * nx, 0);
    realtype *sx_tmp;
    for (int ip = 0; ip < rdata->nplist; ++ip) {
        sx_tmp = NV_DATA_S(tdata->sx[ip]);
        if (!sx_tmp) {
            throw NullPointerException("sx_tmp");
        }
        for (int ix = 0; ix < nx; ++ix)
            sxTmp[ix + ip * nx] = sx_tmp[ix];
    }

    for (int ix = 0; ix < nx; ++ix)
        for (int iJ = 0; iJ < nJ; ++iJ)
            dJzdxTmp[iJ + ix * nJ] =
                tdata->dJzdx[tdata->nroots[ie] +
                             (iJ + ix * nJ) * rdata->nmaxevent];

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans, nJ,
                rdata->nplist, nx, 1.0, dJzdxTmp.data(), nJ, sxTmp.data(), nx, 1.0,
                multResult.data(), nJ);

    // sJy += multResult + dJydp
    for (int iJ = 0; iJ < nJ; ++iJ) {
        if (iJ == 0)
            for (int ip = 0; ip < rdata->nplist; ++ip)
                rdata->sllh[ip] -= multResult[ip * nJ] + tdata->dJzdp[ip * nJ];
        else
            for (int ip = 0; ip < rdata->nplist; ++ip)
                rdata->s2llh[(iJ - 1) + ip * (nJ - 1)] -=
                    multResult[iJ + ip * nJ] + tdata->dJzdp[iJ + ip * nJ];
    }
 

    return;
}

/** Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t.
 * parameters
 * @param[in] ie event index @type int
 * @param[in,out] tdata pointer to temp data object @type TempData
 * @param[in] edata pointer to experimental data object @type ExpData
 * @param[in] rdata pointer to return data object @type ReturnData
 */
void Model::fdJzdp(const int ie, TempData *tdata, const ExpData *edata,
                  const ReturnData *rdata) {
    // dJzdz         nztrue x nJ x nz
    // dJzdsigma     nztrue x nJ x nz
    // dzdp          nz x rdata->nplist
    // dJzdp         nJ x rdata->nplist

    memset(tdata->dJzdp, 0, nJ * rdata->nplist * sizeof(double));

    for (int izt = 0; izt < nztrue; ++izt) {
        if (amiIsNaN(edata->mz[tdata->nroots[ie] + izt * rdata->nmaxevent]))
            continue;

        // copy current (izt) dJzdz and dJzdsigma slices
        // dJzdzTmp     nJ x nz
        // dJzdsigmaTmp nJ x nz

        if (tdata->t < rdata->ts[rdata->nt - 1]) {
            for (int iJ = 0; iJ < nJ; ++iJ) {
                for (int iz = 0; iz < nz; ++iz) {
                    dJzdzTmp[iJ + iz * nJ] =
                        tdata->dJzdz[izt + (iJ + iz * nJ) * nztrue];
                }
            }
        } else {
            for (int iJ = 0; iJ < nJ; ++iJ) {
                for (int iz = 0; iz < nz; ++iz) {
                    dJzdzTmp[iJ + iz * nJ] =
                        tdata->dJrzdz[izt + (iJ + iz * nJ) * nztrue];
                    dJrzdsigmaTmp[iJ + iz * nJ] =
                        tdata->dJrzdsigma[izt + (iJ + iz * nJ) * nztrue];
                }
            }
            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                        AMICI_BLAS_NoTrans, nJ, rdata->nplist, nz, 1.0,
                        dJrzdsigmaTmp.data(), nJ, tdata->dsigmazdp, nz, 1.0,
                        tdata->dJzdp, nJ);
        }

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, nz, 1.0, dJzdzTmp.data(), nJ, tdata->dzdp, nz,
                    1.0, tdata->dJzdp, nJ);

        for (int iJ = 0; iJ < nJ; ++iJ) {
            for (int iz = 0; iz < nz; ++iz) {
                dJzdsigmaTmp[iJ + iz * nJ] =
                    tdata->dJzdsigma[izt + (iJ + iz * nJ) * nztrue];
            }
        }

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, nz, 1.0, dJzdsigmaTmp.data(), nJ,
                    tdata->dsigmazdp, nz, 1.0, tdata->dJzdp, nJ);
    }

    return;
}

/** Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t.
 * state variables
 * @param[in] ie event index @type int
 * @param[in,out] tdata pointer to temp data object @type TempData
 * @param[in] edata pointer to experimental data object @type ExpData
 */
void Model::fdJzdx(const int ie, TempData *tdata, const ExpData *edata) {
    // dJzdz         nztrue x nJ x nz
    // dzdx          nz x nx
    // dJzdx         nmaxevent x nJ x nx

    std::vector<double> multResult;
    multResult.resize(nJ * nx, 0);
    for (int izt = 0; izt < nztrue; ++izt) {
        if (amiIsNaN(
                edata->mz[tdata->nroots[ie] + izt * tdata->rdata->nmaxevent]))
            continue;

        // copy current (izt) dJzdz slice
        // dJzdzTmp     nJ x nz
        if (tdata->t < tdata->rdata->ts[tdata->rdata->nt - 1]) {
            for (int iJ = 0; iJ < nJ; ++iJ)
                for (int iz = 0; iz < nz; ++iz)
                    dJzdzTmp[iJ + iz * nJ] =
                        tdata->dJzdz[izt + (iJ + iz * nJ) * nztrue];

            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                        AMICI_BLAS_NoTrans, nJ, nx, nz, 1.0, dJzdzTmp.data(), nJ,
                        tdata->dzdx, nz, 1.0, multResult.data(), nJ);
        } else {
            for (int iJ = 0; iJ < nJ; ++iJ) {
                for (int iz = 0; iz < nz; ++iz) {
                    dJzdzTmp[iJ + iz * nJ] =
                        tdata->dJrzdz[izt + (iJ + iz * nJ) * nztrue];
                }
            }

            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                        AMICI_BLAS_NoTrans, nJ, nx, nz, 1.0, dJzdzTmp.data(), nJ,
                        tdata->drzdx, nz, 1.0, multResult.data(), nJ);
        }
    }
    for (int iJ = 0; iJ < nJ; ++iJ)
        for (int ix = 0; ix < nx; ++ix)
            tdata->dJzdx[tdata->nroots[ie] +
                         (iJ + ix * nJ) * tdata->rdata->nmaxevent] =
                multResult[iJ + ix * nJ];

    return;
}

/** initialization of model properties
 * @param[in] udata pointer to user data object @type UserData
 * @param[out] tdata pointer to temp data object @type TempData
 */
void Model::initialize(const UserData *udata, TempData *tdata) {
    if (nx < 1)
        return;
    
    initializeStates(udata->x0data, tdata);
    
    fdx0(tdata->x, tdata->dx, tdata);
    
    initHeaviside(tdata);
    
    return;
}

/** initialization of initial states
 * @param[in] x0data array with initial state values @type double
 * @param[out] tdata pointer to temp data object @type TempData
 */
void Model::initializeStates(const double *x0data, TempData *tdata) {
    if (nx < 1)
        return;

    if (!tdata->x)
        throw NullPointerException("tdata->x");

    if (!x0data) {
        fx0(tdata->x, tdata);
    } else {
        realtype *x_tmp = NV_DATA_S(tdata->x);
        if (!x_tmp)
            throw NullPointerException("x_tmp");

        for (int ix = 0; ix < nx; ix++) {
            x_tmp[ix] = (realtype)x0data[ix];
        }
    }

    return;
}

/**
 * initHeaviside initialises the heaviside variables h at the intial time t0
 * heaviside variables activate/deactivate on event occurences
 *
 * @param[out] tdata pointer to the temporary data struct @type TempData
 */
void Model::initHeaviside(TempData *tdata) {
    
    froot(tdata->t, tdata->x, tdata->dx, tdata->rootvals, tdata);

    for (int ie = 0; ie < ne; ie++) {
        if (tdata->rootvals[ie] < 0) {
            tdata->h[ie] = 0.0;
        } else if (tdata->rootvals[ie] == 0) {
            throw AmiException("Simulation started in an event. This could lead to "
                               "unexpected results, aborting simulation! Please "
                               "specify an earlier simulation start via "
                               "@amimodel.t0");
        } else {
            tdata->h[ie] = 1.0;
        }
    }
    return;
}


} // namespace amici
