#include "include/udata.h"
#include "include/amici_model.h"
#include <cstdio>
#include <cstring>

UserData::UserData()
{
    init();
}

int UserData::unscaleParameters(const Model *model, double *bufferUnscaled) const
{
    switch(pscale) {
        case AMICI_SCALING_LOG10:
            for(int ip = 0; ip < model->np; ++ip) {
                bufferUnscaled[ip] = pow(10, p[ip]);
            }
            break;
        case AMICI_SCALING_LN:
            for(int ip = 0; ip < model->np; ++ip)
                bufferUnscaled[ip] = exp(p[ip]);
            break;
        case AMICI_SCALING_NONE:
            for(int ip = 0; ip < model->np; ++ip)
                bufferUnscaled[ip] = p[ip];
            break;
    }

    return AMICI_SUCCESS;
}

UserData::~UserData()
{
    if(qpositivex) delete[] qpositivex;
    if(p) delete[] p;
    if(k) delete[] k;
    if(ts) delete[] ts;
    if(pbar) delete[] pbar;
    if(xbar) delete[] xbar;
    if(x0data) delete[] x0data;
    if(sx0data) delete[] sx0data;
    if(plist) delete[] plist;
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
    sensi = AMICI_SENSI_ORDER_NONE;
    atol = 1e-16;
    rtol = 1e-8;
    maxsteps = 0;
    newton_maxsteps = 0;
    newton_maxlinsteps = 0;
    newton_linsol = 7;
    newton_precon = 1;
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

}

void UserData::print()
{
    printf("qpositivex: %p\n", qpositivex);
    printf("plist: %p\n", plist);
    printf("nplist: %d\n", nplist);
    printf("nt: %d\n", nt);
    printf("nmaxevent: %d\n", nmaxevent);
    printf("p: %p\n", p);
    printf("k: %p\n", k);
    printf("tstart: %e\n", tstart);
    printf("ts: %p\n", ts);
    printf("pbar: %p\n", pbar);
    printf("xbar: %p\n", xbar);
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
    printf("x0data: %p\n", x0data);
    printf("sx0data: %p\n", sx0data);
    printf("ordering: %d\n", ordering);
}
