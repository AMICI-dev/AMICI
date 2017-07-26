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

int UserData::unscaleParameters()
{
    switch(pscale) {
        case AMICI_SCALING_LOG10:
            for(int ip = 0; ip < np; ++ip) {
                p[ip] = pow(10, p[ip]);
            }
            break;
        case AMICI_SCALING_LN:
            for(int ip = 0; ip < np; ++ip)
                p[ip] = exp(p[ip]);
            break;
        case AMICI_SCALING_NONE:
            //this should never be reached
            break;
    }
    return AMICI_SUCCESS;
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
    newton_maxsteps = 0;
    newton_maxlinsteps = 0;
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

void UserData::print()
{
    printf("qpositivex: %p\n", qpositivex);
    printf("plist: %p\n", plist);
    printf("nplist: %d\n", nplist);
    printf("np: %d\n", np);
    printf("ny: %d\n", ny);
    printf("nytrue: %d\n", nytrue);
    printf("nx: %d\n", nx);
    printf("nxtrue: %d\n", nxtrue);
    printf("nz: %d\n", nz);
    printf("nztrue: %d\n", nztrue);
    printf("ne: %d\n", ne);
    printf("nt: %d\n", nt);
    printf("nJ: %d\n", nJ);
    printf("nw: %d\n", nw);
    printf("ndwdx: %d\n", ndwdx);
    printf("nnz: %d\n", nnz);
    printf("nmaxevent: %d\n", nmaxevent);
    printf("pscale: %d\n", (int)pscale);
    printf("p: %p\n", p);
    printf("k: %p\n", k);
    printf("tstart: %e\n", tstart);
    printf("ts: %p\n", ts);
    printf("pbar: %p\n", pbar);
    printf("xbar: %p\n", xbar);
    printf("idlist: %p\n", idlist);
    printf("sensi: %d\n", sensi);
    printf("atol: %e\n", atol);
    printf("rtol: %e\n", rtol);
    printf("maxsteps: %d\n", maxsteps);
    printf("newton_maxsteps: %d\n", newton_maxsteps);
    printf("newton_maxlinsteps: %d\n", newton_maxlinsteps);
    printf("ism: %d\n", ism);
    printf("sensi_meth: %d\n", sensi_meth);
    printf("linsol: %d\n", linsol);
    printf("interpType: %d\n", interpType);
    printf("lmm: %d\n", lmm);
    printf("iter: %d\n", iter);
    printf("stldet: %d\n", stldet);
    printf("ubw: %d\n", ubw);
    printf("lbw: %d\n", lbw);
    printf("x0data: %p\n", x0data);
    printf("sx0data: %p\n", sx0data);
    printf("ordering: %d\n", ordering);
    printf("z2event: %p\n", z2event);
    printf("h: %p\n", h);
}
