#include "include/amici_model.h"
#include <cstring>
#include <include/edata.h>
#include <include/rdata.h>
#include <include/tdata.h>
#include <include/udata.h>

// int Model::fdx0(N_Vector x0, N_Vector dx0, void *user_data)
//{
//    UserData *udata = (UserData*) user_data;
//    realtype *x0_tmp = N_VGetArrayPointer(x0);

//    return fdx0(udata->k, x0_tmp);
//}

Model::~Model() {
    if (z2event)
        delete[] z2event;

    if (idlist)
        delete[] idlist;
}

int Model::fsy(int it, TempData *tdata, ReturnData *rdata) {
    // Compute sy = dydx * sx + dydp

    int status = AMICI_SUCCESS;

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

    return status;
}

int Model::fsz_tf(int ie, TempData *tdata, ReturnData *rdata) {
    // Compute sz = dzdx * sz + dzdp

    int status = AMICI_SUCCESS;

    for (int ip = 0; ip < rdata->nplist; ++ip) {
        for (int iz = 0; iz < nz; ++iz)
            // copy dydp to sy
            rdata->sz[tdata->nroots[ie] + (iz + ip * nz) * rdata->nmaxevent] =
                0;
    }

    return status;
}

int Model::fsJy(int it, TempData *tdata, ReturnData *rdata) {
    int status = AMICI_SUCCESS;

    // Compute dJydx*sx for current 'it'
    // dJydx        rdata->nt x nJ x nx
    // sx           rdata->nt x nx x rdata->nplist

    double *multResult = new double[nJ * rdata->nplist];
    double *dJydxTmp = new double[nJ * nx];
    double *sxTmp = new double[rdata->nplist * nx];
    for (int ix = 0; ix < nx; ++ix) {
        for (int ip = 0; ip < rdata->nplist; ++ip)
            sxTmp[ix + ip * nx] = rdata->sx[it + (ix + ip * nx) * rdata->nt];
        for (int iJ = 0; iJ < nJ; ++iJ)
            dJydxTmp[iJ + ix * nJ] =
                tdata->dJydx[it + (iJ + ix * nJ) * rdata->nt];
    }

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans, nJ,
                rdata->nplist, nx, 1.0, dJydxTmp, nJ, sxTmp, nx, 0.0,
                multResult, nJ);

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

    delete[] dJydxTmp;
    delete[] multResult;
    delete[] sxTmp;

    return (status);
}

int Model::fdJydp(int it, TempData *tdata, const ExpData *edata,
                  ReturnData *rdata) {

    int status = AMICI_SUCCESS;

    // dJydy         nytrue x nJ x ny
    // dydp          ny x rdata->nplist
    // dJydp         nJ x rdata->nplist

    memset(tdata->dJydp, 0, nJ * rdata->nplist * sizeof(double));

    realtype *dJydyTmp = new double[nJ * ny];
    realtype *dJydsigmaTmp = new double[nJ * ny];

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
                    nJ, rdata->nplist, ny, 1.0, dJydyTmp, nJ, tdata->dydp, ny,
                    1.0, tdata->dJydp, nJ);

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, ny, 1.0, dJydsigmaTmp, nJ,
                    tdata->dsigmaydp, ny, 1.0, tdata->dJydp, nJ);
    }
    delete[] dJydyTmp;
    delete[] dJydsigmaTmp;

    return (status);
}

int Model::fdJydx(int it, TempData *tdata, const ExpData *edata) {
    int status = AMICI_SUCCESS;

    // dJydy         nytrue x nJ x ny
    // dydx          ny x nx
    // dJydx         rdata->nt x nJ x nx

    realtype *dJydyTmp = new realtype[nJ * ny];
    realtype *multResult = new realtype[nJ * nx]();

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
                    nJ, nx, ny, 1.0, dJydyTmp, nJ, tdata->dydx, ny, 1.0,
                    multResult, nJ);
    }
    for (int iJ = 0; iJ < nJ; ++iJ)
        for (int ix = 0; ix < nx; ++ix)
            tdata->dJydx[it + (iJ + ix * nJ) * tdata->rdata->nt] =
                multResult[iJ + ix * nJ];

    delete[] dJydyTmp;
    delete[] multResult;

    return (status);
}

int Model::fsJz(int ie, TempData *tdata, const ReturnData *rdata) {
    int status = AMICI_SUCCESS;

    // sJz           nJ x rdata->nplist
    // dJzdp         nJ x rdata->nplist
    // dJzdx         nmaxevent x nJ x nx
    // sx            rdata->nt x nx x rdata->nplist

    // Compute dJzdx*sx for current 'ie'
    // dJzdx        rdata->nt x nJ x nx
    // sx           rdata->nt x nx x rdata->nplist

    realtype *multResult = new realtype[nJ * rdata->nplist]();
    realtype *dJzdxTmp = new realtype[nJ * nx];
    realtype *sxTmp = new realtype[rdata->nplist * nx];
    realtype *sx_tmp;
    for (int ip = 0; ip < rdata->nplist; ++ip) {
        sx_tmp = NV_DATA_S(tdata->sx[ip]);
        if (!sx_tmp)
            return AMICI_ERROR_FSA;
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
                rdata->nplist, nx, 1.0, dJzdxTmp, nJ, sxTmp, nx, 1.0,
                multResult, nJ);

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

    delete[] dJzdxTmp;
    delete[] multResult;
    delete[] sxTmp;

    return (status);
}

int Model::fdJzdp(int ie, TempData *tdata, const ExpData *edata,
                  ReturnData *rdata) {
    int status = AMICI_SUCCESS;

    // dJzdz         nztrue x nJ x nz
    // dJzdsigma     nztrue x nJ x nz
    // dzdp          nz x rdata->nplist
    // dJzdp         nJ x rdata->nplist

    memset(tdata->dJzdp, 0, nJ * rdata->nplist * sizeof(double));

    realtype *dJzdzTmp = new double[nJ * nz];
    realtype *dJzdsigmaTmp = new double[nJ * nz];
    realtype *dJrzdsigmaTmp = NULL;
    if (tdata->t == rdata->ts[rdata->nt - 1]) {
        dJrzdsigmaTmp = new double[nJ * nz];
    }

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
                        dJrzdsigmaTmp, nJ, tdata->dsigmazdp, nz, 1.0,
                        tdata->dJzdp, nJ);
        }

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, nz, 1.0, dJzdzTmp, nJ, tdata->dzdp, nz,
                    1.0, tdata->dJzdp, nJ);

        for (int iJ = 0; iJ < nJ; ++iJ) {
            for (int iz = 0; iz < nz; ++iz) {
                dJzdsigmaTmp[iJ + iz * nJ] =
                    tdata->dJzdsigma[izt + (iJ + iz * nJ) * nztrue];
            }
        }

        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nJ, rdata->nplist, nz, 1.0, dJzdsigmaTmp, nJ,
                    tdata->dsigmazdp, nz, 1.0, tdata->dJzdp, nJ);
    }
    delete[] dJzdzTmp;
    delete[] dJzdsigmaTmp;
    if (dJrzdsigmaTmp)
        delete[] dJrzdsigmaTmp;

    return (status);
}

int Model::fdJzdx(int ie, TempData *tdata, const ExpData *edata) {
    int status = AMICI_SUCCESS;

    // dJzdz         nztrue x nJ x nz
    // dzdx          nz x nx
    // dJzdx         nmaxevent x nJ x nx

    realtype *dJzdzTmp = new realtype[nJ * nz];
    realtype *multResult = new realtype[nJ * nx]();
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
                        AMICI_BLAS_NoTrans, nJ, nx, nz, 1.0, dJzdzTmp, nJ,
                        tdata->dzdx, nz, 1.0, multResult, nJ);
        } else {
            for (int iJ = 0; iJ < nJ; ++iJ) {
                for (int iz = 0; iz < nz; ++iz) {
                    dJzdzTmp[iJ + iz * nJ] =
                        tdata->dJrzdz[izt + (iJ + iz * nJ) * nztrue];
                }
            }

            amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                        AMICI_BLAS_NoTrans, nJ, nx, nz, 1.0, dJzdzTmp, nJ,
                        tdata->drzdx, nz, 1.0, multResult, nJ);
        }
    }
    for (int iJ = 0; iJ < nJ; ++iJ)
        for (int ix = 0; ix < nx; ++ix)
            tdata->dJzdx[tdata->nroots[ie] +
                         (iJ + ix * nJ) * tdata->rdata->nmaxevent] =
                multResult[iJ + ix * nJ];

    delete[] dJzdzTmp;
    delete[] multResult;

    return (status);
}

int Model::initialize(UserData *udata, TempData *tdata) {
    if (nx < 1)
        return AMICI_SUCCESS;

    int status;

    if ((status = initializeStates(udata->x0data, tdata)) != AMICI_SUCCESS)
        return status;

    if ((status = fdx0(tdata->x, tdata->dx, tdata)) != AMICI_SUCCESS)
        return status;

    if ((status = initHeaviside(tdata)) != AMICI_SUCCESS)
        return status;

    return AMICI_SUCCESS;
}

int Model::initializeStates(double *x0data, TempData *tdata) {
    if (nx < 1)
        return AMICI_SUCCESS;

    if (tdata->x == NULL)
        return AMICI_ERROR_TDATA;

    if (x0data == NULL) {
        if (fx0(tdata->x, tdata) != AMICI_SUCCESS)
            return AMICI_ERROR_MODEL;
    } else {
        realtype *x_tmp = NV_DATA_S(tdata->x);
        if (!x_tmp)
            return AMICI_ERROR_TDATA;

        for (int ix = 0; ix < nx; ix++) {
            x_tmp[ix] = (realtype)x0data[ix];
        }
    }

    return AMICI_SUCCESS;
}

int Model::initHeaviside(TempData *tdata) {
    /**
     * initHeaviside initialises the heaviside variables h at the intial time t0
     * heaviside variables activate/deactivate on event occurences
     *
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */

    int status = AMICI_SUCCESS;

    status = froot(tdata->t, tdata->x, tdata->dx, tdata->rootvals, tdata);
    if (status != AMICI_SUCCESS)
        return status;

    for (int ie = 0; ie < ne; ie++) {
        if (tdata->rootvals[ie] < 0) {
            tdata->h_udata[ie] = 0.0;
        } else if (tdata->rootvals[ie] == 0) {
            errMsgIdAndTxt("AMICI:mex:initHeaviside",
                           "Simulation started in an event. This could lead to "
                           "unexpected results, aborting simulation! Please "
                           "specify an earlier simulation start via "
                           "@amimodel.t0");
            return AMICI_ERROR_EVENT;
        } else {
            tdata->h_udata[ie] = 1.0;
        }
    }
    return status;
}
