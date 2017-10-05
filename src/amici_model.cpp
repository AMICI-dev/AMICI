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
 * @return status flag indicating successful execution @type int
 */
int Model::fsy(const int it, const TempData *tdata, ReturnData *rdata) {
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

/** Sensitivity of z at final timepoint (ignores sensitivity of timepoint),
 * total derivative
 * @param[in] ie event index @type int
 * @param[in] tdata pointer to temp data object @type TempData
 * @param[in,out] rdata pointer to return data object @type ReturnData
 * @return status flag indicating successful execution @type int
 */
int Model::fsz_tf(const int ie, const TempData *tdata, ReturnData *rdata) {
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

/** Sensitivity of time-resolved measurement negative log-likelihood Jy, total
 * derivative
 * @param[in] it timepoint index @type int
 * @param[in] tdata pointer to temp data object @type TempData
 * @param[in,out] rdata pointer to return data object @type ReturnData
 * @return status flag indicating successful execution @type int
 */
int Model::fsJy(const int it, const TempData *tdata, ReturnData *rdata) {
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

/** Sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t.
 * parameters
 * @param[in] it timepoint index @type int
 * @param[in,out] tdata pointer to temp data object @type TempData
 * @param[in] edata pointer to experimental data object @type ExpData
 * @param[in] rdata pointer to return data object @type ReturnData
 * @return status flag indicating successful execution @type int
 */
int Model::fdJydp(const int it, TempData *tdata, const ExpData *edata,
                  const ReturnData *rdata) {

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

/** Second order sensitivity of time-resolved measurement negative 
 * log-likelihood Jy w.r.t. parameters
 * @param[in] it timepoint index @type int
 * @param[in,out] tdata pointer to temp data object @type TempData
 * @param[in] edata pointer to experimental data object @type ExpData
 * @param[in] rdata pointer to return data object @type ReturnData
 * @return status flag indicating successful execution @type int
 */
int Model::fddJydpdp(const int it, TempData *tdata, const ExpData *edata,
                  const ReturnData *rdata) {
    
    int status = AMICI_SUCCESS;
    
    memset(tdata->ddJydpdp, 0, rdata->nplist * rdata->nplist * sizeof(double));
    
    realtype *sx_tmp;
    
    /* Temporary variables, think of how to do this more efficiently */
    realtype *sxTmp = new double[nx * rdata->nplist];
    realtype *syTmp = new double[ny * rdata->nplist];
    realtype *ddJy_tmp1 = new double[ny * ny];
    realtype *ddJy_tmp2 = new double[ny * rdata->nplist];
    realtype *ddJy_tmp3 = new double[rdata->nplist * rdata->nplist];
    realtype *ddJy_tmp4 = new double[nx * nx];
    realtype *ddJy_tmp5 = new double[nx * rdata->nplist];
    
    /*
     Short description:
     
     ddJydpdp =
         ddJy/dsds * ds/dp * ds/dp     -- Part 1a
       + dJy/ds * dds/dpdp             -- Part 1b
       + ddJy/dyds * ds/dp * sy +      -- Part 2
         (ddJy/dyds * ds/dp * sy)'
       + ddJy/dydy * sy * sy           -- Part 3
       + dJy/dy * ddy/dxdx * sx * sx   -- Part 4a
       + dJy/dy * ddy/dxdp * sx +      -- Part 4b
         (dJy/dy * ddy/dxdp * sx)'
       + dJy/dy * ddy/dpdp             -- Part 4c
     */
    
    
    for (int ip = 0; ip < np; ip++) {
        sx_tmp = N_VGetArrayPointer(tdata->sx[ip]);
        for (int ix = 0; ix < nx; ix++)
            sxTmp[ix + ip * nx] = sx_tmp[ix];
        for (int iy = 0; iy < ny; iy++)
            syTmp[iy + ip * ny] = rdata->sy[ip * rdata->nt * ny +
                                            iy * rdata->nt + it];
    }
    
    
    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (amiIsNaN(edata->my[rdata->nt * iyt + it]))
            continue;

        // Part 1a: ddJy_tmp1 = Slicing ddJy/dsds
        for (int iy = 0; iy < ny; ++iy)
            for (int jy = 0; jy < ny; ++jy)
                ddJy_tmp1[jy + iy * ny] =
                tdata->ddJydsigmadsigma[iyt + (jy + iy * ny) * nytrue];
        //          ddJy_tmp2 = ddJydsds * dsdp
        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    ny, rdata->nplist, ny, 1.0, ddJy_tmp1, ny, tdata->dsigmaydp, ny,
                    0.0, ddJy_tmp2, nJ);
        //          tdata->ddJydpdp += dsdp' * ddJy_tmp2
        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_Trans, AMICI_BLAS_NoTrans,
                    ny, rdata->nplist, ny, 1.0, ddJy_tmp2, rdata->nplist,
                    tdata->dsigmaydp, ny, 1.0, tdata->ddJydpdp, rdata->nplist);
        
        // Part 1b: tdata->ddJydpdp += ddJy_s2sigma
        for (int ip = 0; ip < rdata->nplist; ++ip)
            for (int jp = 0; jp < rdata->nplist; ++jp)
                tdata->ddJydpdp[jp + ip * rdata->nplist] +=
                tdata->ddJy_s2sigma[iyt + (jp + ip * rdata->nplist) * rdata->nplist];
        
        // Part 2:  ddJy_tmp2 = ddJydsdy * ds/dp
        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    ny, rdata->nplist, ny, 1.0, tdata->ddJydsigmady, nytrue*ny, tdata->dsigmaydp, ny,
                    0.0, ddJy_tmp2, ny);
        //          ddJy_tmp3 = sy' * ddJy_tmp2
        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_Trans, AMICI_BLAS_NoTrans,
                    ny, rdata->nplist, ny, 1.0, ddJy_tmp2, rdata->nplist,
                    syTmp, ny, 1.0, ddJy_tmp3, rdata->nplist);
        //          tdata->ddJydpdp += ddJy_tmp3 + ddJy_tmp3'
        for (int ip = 0; ip < np; ip++)
            for (int jp = 0; jp < np; jp++)
                tdata->ddJydpdp[jp + ip*rdata->nplist] +=
                ddJy_tmp3[jp + ip*rdata->nplist] +
                ddJy_tmp3[ip + jp*rdata->nplist];
        
        // Part 3:  ddJy_tmp2 = ddJy/dydy * sy
        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    ny, rdata->nplist, ny, 1.0, tdata->ddJydydy, nytrue*ny, syTmp, ny,
                    0.0, ddJy_tmp2, ny);
        //          fddJydpdp += sy' * ddJy_tmp2
        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_Trans, AMICI_BLAS_NoTrans,
                    ny, rdata->nplist, ny, 1.0, ddJy_tmp2, rdata->nplist, syTmp, ny,
                    1.0, tdata->ddJydpdp, rdata->nplist);
        
        // Part 4a: ddJy_tmp4 = dJy/dy * ddy/dxdx
        amici_dgemv(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                    ny, nx*nx, 1.0, tdata->ddydxdx, ny, tdata->dJydy, 1,
                    0.0, ddJy_tmp4, 1);
        //          ddJy_tmp5 = ddJy_tmp4 * sx
        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    nx, rdata->nplist, nx, 1.0, ddJy_tmp4, nx, sxTmp, rdata->nplist,
                    0.0, ddJy_tmp5, nx);
        //          tdata->ddJydpdp += sx' * ddJy_tmp5
        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_Trans, AMICI_BLAS_NoTrans,
                    nx, rdata->nplist, nx, 1.0, ddJy_tmp5, rdata->nplist, sxTmp, nx,
                    1.0, tdata->ddJydpdp, rdata->nplist);
        
        // Part 4b: ddJy_tmp5 = dJy/dy * ddy/dpdx
        amici_dgemv(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                   ny, nx*rdata->nplist, 1.0, tdata->ddydpdx, ny, tdata->dJydy, 1,
                   0.0, ddJy_tmp5, 1);
        //          ddJy_tmp3 = ddJy_tmp5 * sx
        amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans, AMICI_BLAS_NoTrans,
                    rdata->nplist, rdata->nplist, nx, 1.0, ddJy_tmp5, rdata->nplist,
                    sxTmp, nx, 0.0, ddJy_tmp3, rdata->nplist);
        //          tdata->ddJydpdp += ddJy_tmp3 + ddJy_tmp3'
        for (int ip = 0; ip < np; ip++)
            for (int jp = 0; jp < np; jp++)
                tdata->ddJydpdp[jp + ip*rdata->nplist] +=
                ddJy_tmp3[jp + ip*rdata->nplist] +
                ddJy_tmp3[ip + jp*rdata->nplist];
        
        // Part 4c: tdata->ddJydpdp += dJy/dy * ddy/dpdp
        amici_dgemv(AMICI_BLAS_ColMajor, AMICI_BLAS_NoTrans,
                    ny, rdata->nplist*rdata->nplist, 1.0, tdata->ddydpdp, ny, tdata->dJydy, 1,
                    1.0, tdata->ddJydpdp, 1);
    }
    
    delete[] sxTmp;
    delete[] syTmp;
    delete[] ddJy_tmp1;
    delete[] ddJy_tmp2;
    delete[] ddJy_tmp3;
    delete[] ddJy_tmp4;
    delete[] ddJy_tmp5;
    
    return (status);
}


/** Quadrature equations for second order adjoint sensitivities
 * for negative log-likelihood Jy w.r.t. parameters
 * @param[in] it timepoint index @type int
 * @param[in,out] tdata pointer to temp data object @type TempData
 * @param[in] edata pointer to experimental data object @type ExpData
 * @param[in] rdata pointer to return data object @type ReturnData
 * @return status flag indicating successful execution @type int
 */
int Model::fqBo2dot(realtype t, N_Vector x, N_Vector *sx, N_Vector xB,
                    N_Vector qBdot, void *user_data) {
    
    UserData *udata = (UserData*) user_data;
    int status = AMICI_SUCCESS;
    int np = udata->nplist;
    
    realtype *qBo2dot = N_VGetArrayPointer(qBdot);
    realtype *xB_tmp = N_VGetArrayPointer(xB);
    realtype *sx_tmp = N_VGetArrayPointer(sx[0]);
    
    /* Temporary variables, think of how to do this more efficiently */
    realtype *sxTmp = new double[nx * np];
    realtype *dJdx = new double[nx * nx * nx];
    realtype *dJdp = new double[nx * nx * np];
    realtype *ddfdpdp = new double[nx * np * np];
    realtype *qBo2_tmp1 = new double[nx * nx];
    realtype *qBo2_tmp2 = new double[nx * np];
    realtype *qBo2_tmp3 = new double[np * np];
    
    /* Copy sensitivities to make them BLAS usable */
    for (int ip = 0; ip < np; ip++) {
        sx_tmp = N_VGetArrayPointer(sx[ip]);
        for (int ix = 0; ix < nx; ix++)
            sxTmp[ix + ip * nx] = sx_tmp[ix];
    }
    
    /* Prepare quadrature fields to be read */
    status = fdJdx(t, x, x, dJdx, user_data);
    if (status != AMICI_SUCCESS)
        return status;
    status = fdJdp(t, x, x, dJdp, user_data);
    if (status != AMICI_SUCCESS)
        return status;
    status = fddxdotdpdp(t, x, x, ddfdpdp, user_data);
    if (status != AMICI_SUCCESS)
        return status;
    
    // Compute matrix xB' * dJdx
    for (int ix = 0; ix < nx; ix++) {
        
        // col ix of xBdJdxTmp' = xB' * dJdx
        amici_dgemv(AMICI_BLAS_ColMajor, AMICI_BLAS_Trans,
                    nx, nx, 1.0, dJdx,
                    nx, xB_tmp, 1,
                    0.0, qBo2_tmp1, nx);
        // increment pointer to next row
        dJdx += nx*nx;
    }
    
    // qBo2_part1 = qBo2_part1_1 * sx
    amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_Trans, AMICI_BLAS_NoTrans,
                nx, np, nx, 1.0, qBo2_tmp1, nx, sxTmp, nx,
                0.0, qBo2_tmp2, nx);
    
    // qBo2dot_tmp = sx' * qBo2_part1
    amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_Trans, AMICI_BLAS_NoTrans,
                np, np, nx, 1.0, sxTmp, np,
                qBo2_tmp2, nx, 0.0, qBo2dot, np);

    
    // Compute matrix xB' * dJdp
    for (int ip = 0; ip < udata->nplist; ip++) {
        /* col ip of qBo2_part2_1' = xB' * dJdp */
        amici_dgemv(AMICI_BLAS_ColMajor, AMICI_BLAS_Trans,
                    nx, nx, 1.0, dJdp,
                    nx, xB_tmp, 1,
                    0.0, qBo2_tmp2, nx);
        /* increment pointer to next row */
        dJdp += nx*nx;
    }
    
    // qBo2_part2 = qBo2_part2_1 * sx
    amici_dgemm(AMICI_BLAS_ColMajor, AMICI_BLAS_Trans, AMICI_BLAS_NoTrans,
                np, np, nx, 1.0, qBo2_tmp2, np, sxTmp, nx,
                0.0, qBo2_tmp3, np);
    
    // qBo2dot_tmp += qBo2_part2 + qBo2_part2'
    for (int ip = 0; ip < np; ip++)
        for (int jp = 0; jp < np; jp++)
            qBo2dot[jp + ip*udata->nplist] +=
                qBo2_tmp3[jp + ip*udata->nplist] +
                qBo2_tmp3[ip + jp*udata->nplist];
    
    // qBo2dot_tmp += xB' * ddfdpdp
    amici_dgemv(AMICI_BLAS_ColMajor, AMICI_BLAS_Trans,
                nx, np * np, 1.0, ddfdpdp,
                nx, xB_tmp, 1, 1.0, qBo2dot, 1);
    
    delete[] sxTmp;
    delete[] dJdx;
    delete[] dJdp;
    delete[] ddfdpdp;
    delete[] qBo2_tmp1;
    delete[] qBo2_tmp2;
    delete[] qBo2_tmp3;
    
    return status;
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t.
 * state variables
 * @param[in] it timepoint index @type int
 * @param[in,out] tdata pointer to temp data object @type TempData
 * @param[in] edata pointer to experimental data object @type ExpData
 * @return status flag indicating successful execution @type int
 */

int Model::fdJydx(const int it, TempData *tdata, const ExpData *edata) {
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

/** Sensitivity of event-resolved measurement negative log-likelihood Jz, total
 * derivative
 * @param[in] ie event index @type int
 * @param[in,out] tdata pointer to temp data object @type TempData
 * @param[in] rdata pointer to return data object @type ReturnData
 * @return status flag indicating successful execution @type int
 */
int Model::fsJz(const int ie, TempData *tdata, const ReturnData *rdata) {
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
        if (!sx_tmp) {
            status = AMICI_ERROR_FSA;
            goto freturn;
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

freturn:
    delete[] dJzdxTmp;
    delete[] multResult;
    delete[] sxTmp;

    return (status);
}

/** Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t.
 * parameters
 * @param[in] ie event index @type int
 * @param[in,out] tdata pointer to temp data object @type TempData
 * @param[in] edata pointer to experimental data object @type ExpData
 * @param[in] rdata pointer to return data object @type ReturnData
 * @return status flag indicating successful execution @type int
 */
int Model::fdJzdp(const int ie, TempData *tdata, const ExpData *edata,
                  const ReturnData *rdata) {
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

/** Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t.
 * state variables
 * @param[in] ie event index @type int
 * @param[in,out] tdata pointer to temp data object @type TempData
 * @param[in] edata pointer to experimental data object @type ExpData
 * @return status flag indicating successful execution @type int
 */
int Model::fdJzdx(const int ie, TempData *tdata, const ExpData *edata) {
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

/** initialization of model properties
 * @param[in] udata pointer to user data object @type UserData
 * @param[out] tdata pointer to temp data object @type TempData
 * @return status flag indicating success of execution @type int
 */
int Model::initialize(const UserData *udata, TempData *tdata) {
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

/** initialization of initial states
 * @param[in] x0data array with initial state values @type double
 * @param[out] tdata pointer to temp data object @type TempData
 * @return status flag indicating success of execution @type int
 */
int Model::initializeStates(const double *x0data, TempData *tdata) {
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

/**
 * initHeaviside initialises the heaviside variables h at the intial time t0
 * heaviside variables activate/deactivate on event occurences
 *
 * @param[out] tdata pointer to the temporary data struct @type TempData
 * @return status flag indicating success of execution @type int
 */
int Model::initHeaviside(TempData *tdata) {

    int status = AMICI_SUCCESS;

    status = froot(tdata->t, tdata->x, tdata->dx, tdata->rootvals, tdata);
    if (status != AMICI_SUCCESS)
        return status;

    for (int ie = 0; ie < ne; ie++) {
        if (tdata->rootvals[ie] < 0) {
            tdata->h[ie] = 0.0;
        } else if (tdata->rootvals[ie] == 0) {
            errMsgIdAndTxt("AMICI:mex:initHeaviside",
                           "Simulation started in an event. This could lead to "
                           "unexpected results, aborting simulation! Please "
                           "specify an earlier simulation start via "
                           "@amimodel.t0");
            return AMICI_ERROR_EVENT;
        } else {
            tdata->h[ie] = 1.0;
        }
    }
    return status;
}
