#include "include/udata.h"

#include <cstdio>
#include <cstring>

UserData::UserData() :
    np(0), nk(0),
    nx(0), nxtrue(0),
    ny(0), nytrue(0),
    nz(0), nztrue(0),
    ne(0), nw(0),
    ndwdx(0), ndwdp(0),
    nnz(0), nJ(0),
    ubw(0), lbw(0),
    o2mode(AMICI_O2MODE_NONE), pscale(AMICI_SCALING_NONE)
{
    init();
}

UserData::UserData(int np,
                   int nx, int nxtrue,
                   int nk,
                   int ny, int nytrue,
                   int nz, int nztrue,
                   int ne, int nJ,
                   int nw, int ndwdx, int ndwdp, int nnz,
                   int ubw, int lbw,
                   AMICI_parameter_scaling pscale,
                   AMICI_o2mode o2mode) :
    np(np), nk(nk),
    nx(nx), nxtrue(nxtrue),
    ny(ny), nytrue(nytrue),
    nz(nz), nztrue(nztrue),
    ne(ne), nw(nw),
    ndwdx(ndwdx), ndwdp(ndwdp),
    nnz(nnz),nJ(nJ),
    ubw(ubw), lbw(lbw),
    o2mode(o2mode), pscale(pscale)
{
    init();
}

/**
 * processUserData initializes fields of the udata struct
 *
 * @param[out] udata pointer to the user data struct @type UserData
 * @return void
 */
void UserData::initTemporaryFields()
{
    if (nx>0) {
        /* initialise temporary jacobian storage */
        J = SparseNewMat(nx,nx,nnz,CSC_MAT);
        M = new realtype[nx*nx]();
        dfdx = new realtype[nx*nx]();
    }
    if (sensi >= AMICI_SENSI_ORDER_FIRST) {
        /* initialise temporary dxdotdp storage */
        dxdotdp = new realtype[nx*nplist]();
    }
    if (ne>0) {
        /* initialise temporary stau storage */
        stau = new realtype[nplist]();
        h = new realtype[ne]();
    }

    w = new realtype[nw]();
    dwdx = new realtype[ndwdx]();
    dwdp = new realtype[ndwdp]();
    
}

void UserData::freeTemporaryFields()
{
    if(dxdotdp) delete[] dxdotdp;
    if(w) delete[] w;
    if(dwdx) delete[] dwdx;
    if(dwdp) delete[] dwdp;
    if(M) delete[] M;
    if(dfdx) delete[] dfdx;
    if(stau) delete[] stau;
    if(J) SparseDestroyMat(J);

    J = NULL;
    dxdotdp = NULL;
    w = NULL;
    dwdx = NULL;
    dwdp = NULL;
    M = NULL;
    dfdx = NULL;
    stau = NULL;
}

UserData::~UserData()
{
    freeTemporaryFields();

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
    if(plist) delete[] plist;
    if(h) delete[] h;
}

void UserData::init()
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
    sensi = AMICI_SENSI_ORDER_NONE;
    atol = 1e-16;
    rtol = 1e-8;
    maxsteps = 0;
    ism = 1;
    nmaxevent = 10;

    sensi_meth = AMICI_SENSI_FSA;
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

    nan_dxdotdp = false;
    nan_J = false;
    nan_JSparse = false;
    nan_xdot = false;
    nan_xBdot = false;
    nan_qBdot = false;

    J = NULL;
    dxdotdp = NULL;
    w = NULL;
    dwdx = NULL;
    dwdp = NULL;
    M = NULL;
    dfdx = NULL;
    stau = NULL;
}

void printUserData(UserData *udata) {
    printf("qpositivex: %p\n", udata->qpositivex);
    printf("plist: %p\n", udata->plist);
    printf("np: %d\n", udata->np);
    printf("ny: %d\n", udata->ny);
    printf("nytrue: %d\n", udata->nytrue);
    printf("nx: %d\n", udata->nx);
    printf("nxtrue: %d\n", udata->nxtrue);
    printf("nz: %d\n", udata->nz);
    printf("nztrue: %d\n", udata->nztrue);
    printf("ne: %d\n", udata->ne);
    printf("nt: %d\n", udata->nt);
    printf("nJ: %d\n", udata->nJ);
    printf("nw: %d\n", udata->nw);
    printf("ndwdx: %d\n", udata->ndwdx);
    printf("nnz: %d\n", udata->nnz);
    printf("nmaxevent: %d\n", udata->nmaxevent);
    printf("pscale: %d\n", (int)udata->pscale);
    printf("p: %p\n", udata->p);
    printf("k: %p\n", udata->k);
    printf("tstart: %e\n", udata->tstart);
    printf("ts: %p\n", udata->ts);
    printf("pbar: %p\n", udata->pbar);
    printf("xbar: %p\n", udata->xbar);
    printf("idlist: %p\n", udata->idlist);
    printf("sensi: %d\n", udata->sensi);
    printf("atol: %e\n", udata->atol);
    printf("rtol: %e\n", udata->rtol);
    printf("maxsteps: %d\n", udata->maxsteps);
    printf("ism: %d\n", udata->ism);
    printf("sensi_meth: %d\n", udata->sensi_meth);
    printf("linsol: %d\n", udata->linsol);
    printf("interpType: %d\n", udata->interpType);
    printf("lmm: %d\n", udata->lmm);
    printf("iter: %d\n", udata->iter);
    printf("stldet: %d\n", udata->stldet);
    printf("ubw: %d\n", udata->ubw);
    printf("lbw: %d\n", udata->lbw);
    printf("x0data: %p\n", udata->x0data);
    printf("sx0data: %p\n", udata->sx0data);
    printf("ordering: %d\n", udata->ordering);
    printf("z2event: %p\n", udata->z2event);
    printf("h: %p\n", udata->h);
}
