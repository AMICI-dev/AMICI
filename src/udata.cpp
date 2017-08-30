#include "include/udata.h"
#include <cstdio>
#include <cstring>


UserData::UserData(int np, int nk, int nx) : np(np), nk(nk), nx(nx)
{

}

UserData::UserData() : np(0), nk(0), nx(0)
{

}

int UserData::unscaleParameters(double *bufferUnscaled) const {
    /**
     * unscaleParameters removes parameter scaling according to the parameter scaling in pscale
     *
     * @param[out] bufferUnscaled unscaled parameters are written to the array @type double
     * 
     * @return status flag indicating success of execution @type int
     */
    switch (pscale) {
    case AMICI_SCALING_LOG10:
        for (int ip = 0; ip < np; ++ip) {
            bufferUnscaled[ip] = pow(10, p[ip]);
        }
        break;
    case AMICI_SCALING_LN:
        for (int ip = 0; ip < np; ++ip)
            bufferUnscaled[ip] = exp(p[ip]);
        break;
    case AMICI_SCALING_NONE:
        for (int ip = 0; ip < np; ++ip)
            bufferUnscaled[ip] = p[ip];
        break;
    }

    return AMICI_SUCCESS;
}

UserData::~UserData() {
    if (qpositivex)
        delete[] qpositivex;
    if (p)
        delete[] p;
    if (k)
        delete[] k;
    if (ts)
        delete[] ts;
    if (pbar)
        delete[] pbar;
    if (xbar)
        delete[] xbar;
    if (x0data)
        delete[] x0data;
    if (sx0data)
        delete[] sx0data;
    if (plist)
        delete[] plist;
}

void UserData::print() const {
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
