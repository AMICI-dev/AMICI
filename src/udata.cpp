#include "include/udata.h"
#include <stdio.h>
#include <cstring>


UserData::UserData(int np,
                   int nx, int nxtrue,
                   int nk,
                   int ny, int nytrue,
                   int nz, int nztrue,
                   int ne, int ng,
                   int nw, int ndwdx, int ndwdp, int nnz,
                   int ubw, int lbw,
                   AMI_parameter_scaling pscale,
                   AMI_o2mode o2mode) :
    np(np),
    nx(nx), nxtrue(nxtrue), nk(nk),
    ny(ny), nytrue(nytrue),
    nz(nz), nztrue(nztrue),
    ne(ne), ng(ng),
    nw(nw), ndwdx(ndwdx),
    ndwdp(ndwdp), nnz(nnz),
    ubw(ubw), lbw(lbw),
    pscale(pscale), o2mode(o2mode)
{

    qpositivex = NULL;
    plist = NULL;
    nplist = 0;
    nt = 0;
    p = NULL;
    k = NULL;
    ts = NULL;
    tstart = 0;
    pbar = NULL;
    xbar = NULL;
    idlist = NULL;
    sensi = AMI_SENSI_ORDER_NONE;
    atol = 1e-16;
    rtol = 1e-8;
    maxsteps = 0;
    ism = 1;
    nmaxevent = 10;

    sensi_meth = AMI_SENSI_FSA;
    linsol = 9;
    interpType = 1;
    lmm = 2;
    iter = 2;
    stldet = true;
    x0data = NULL;

    sx0data = NULL;
    ordering = 0;
    z2event = NULL;
    h = NULL;

    // J;
    dxdotdp = NULL;
    w = NULL;
    dwdx = NULL;
    dwdp = NULL;
    M = NULL;
    dfdx = NULL;
    stau = NULL;

    nan_dxdotdp = false;
    nan_J = false;
    nan_JSparse = false;
    nan_xdot = false;
    nan_xBdot = false;
    nan_qBdot = false;
}

UserData::~UserData()
{
#ifdef AMICI_WITHOUT_MATLAB
        if(qpositivex) delete[] qpositivex;
        if(p) delete[] p;
        if(k) delete[] k;
        if(ts) delete[] ts;
        if(pbar) delete[] pbar;
        if(xbar) delete[] xbar;
        if(idlist) delete[] idlist;
        if(x0data) delete[] x0data;
        if(sx0data) delete[] sx0data;
        if(z2event) delete[] z2event;
#endif
        if(plist) delete[] plist;
        if(h) delete[] h;
        if(dxdotdp) delete[] dxdotdp;
        if(w) delete[] w;
        if(dwdx) delete[] dwdx;
        if(dwdp) delete[] dwdp;
        if(M) delete[] M;
        if(dfdx) delete[] dfdx;
        if(stau) delete[] stau;
        if(J) SparseDestroyMat(J);
}

#ifdef AMICI_WITHOUT_MATLAB
void printUserData(UserData *udata) {
    printf("am_qpositivex: %p\n", udata->qpositivex);
    printf("am_plist: %p\n", udata->plist);
    printf("am_np: %d\n", udata->np);
    printf("am_ny: %d\n", udata->ny);
    printf("am_nytrue: %d\n", udata->nytrue);
    printf("am_nx: %d\n", udata->nx);
    printf("am_nxtrue: %d\n", udata->nxtrue);
    printf("am_nz: %d\n", udata->nz);
    printf("am_nztrue: %d\n", udata->nztrue);
    printf("am_ne: %d\n", udata->ne);
    printf("am_nt: %d\n", udata->nt);
    printf("am_ng: %d\n", udata->ng);
    printf("am_nw: %d\n", udata->nw);
    printf("am_ndwdx: %d\n", udata->ndwdx);
    printf("am_nnz: %d\n", udata->nnz);
    printf("am_nmaxevent: %d\n", udata->nmaxevent);
    printf("am_pscale: %d\n", (int)udata->pscale);
    printf("am_p: %p\n", udata->p);
    printf("am_k: %p\n", udata->k);
    printf("am_tstart: %e\n", udata->tstart);
    printf("am_ts: %p\n", udata->ts);
    printf("am_pbar: %p\n", udata->pbar);
    printf("am_xbar: %p\n", udata->xbar);
    printf("am_idlist: %p\n", udata->idlist);
    printf("am_sensi: %d\n", udata->sensi);
    printf("am_atol: %e\n", udata->atol);
    printf("am_rtol: %e\n", udata->rtol);
    printf("am_maxsteps: %d\n", udata->maxsteps);
    printf("am_ism: %d\n", udata->ism);
    printf("am_sensi_meth: %d\n", udata->sensi_meth);
    printf("am_linsol: %d\n", udata->linsol);
    printf("am_interpType: %d\n", udata->interpType);
    printf("am_lmm: %d\n", udata->lmm);
    printf("am_iter: %d\n", udata->iter);
    printf("am_stldet: %d\n", udata->stldet);
    printf("am_ubw: %d\n", udata->ubw);
    printf("am_lbw: %d\n", udata->lbw);
    printf("am_x0data: %p\n", udata->x0data);
    printf("am_sx0data: %p\n", udata->sx0data);
    printf("am_ordering: %d\n", udata->ordering);
    printf("am_z2event: %p\n", udata->z2event);
    printf("am_h: %p\n", udata->h);
}
#endif
