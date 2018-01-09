#include "include/udata.h"
#include <cstdio>
#include <cstring>
#include <algorithm>

namespace amici {

UserData::UserData(int np, int nk, int nx) : np(np), nk(nk), nx(nx) {
    // these fields must always be initialized, others are optional or can be set later
    k = new double[nk]();
    p = new double[np]();
    qpositivex = new double[nx];
    std::fill(qpositivex, qpositivex + nx, 1);
}

UserData::UserData() : np(0), nk(0), nx(0) {}

UserData::UserData(const UserData &other) : UserData(other.np, other.nk, other.nx)
{
    nmaxevent = other.nmaxevent;

    if(other.qpositivex) {
        std::copy(other.qpositivex, other.qpositivex + other.nx, qpositivex);
    }

    if(other.plist) {
        plist = new int[other.nplist];
        std::copy(other.plist, other.plist + other.nplist, plist);
    }

    nplist = other.nplist;
    nt = other.nt;

    if(other.p) {
        std::copy(other.p, other.p + other.np, p);
    }
    if(other.k) {
        std::copy(other.k, other.k + other.nk, k);
    }

    pscale = other.pscale;
    tstart = other.tstart;

    if(other.ts) {
        ts = new double[other.nt];
        std::copy(other.ts, other.ts + other.nt, ts);
    }
    if(other.pbar) {
        pbar = new double[other.nplist];
        std::copy(other.pbar, other.pbar + other.nplist, pbar);
    }
    if(other.xbar) {
        xbar = new double[other.nx];
        std::copy(other.xbar, other.xbar + other.nx, xbar);
    }

    sensi = other.sensi;
    atol = other.atol;
    rtol = other.rtol;
    maxsteps = other.maxsteps;
    quad_atol = other.quad_atol;
    quad_rtol = other.quad_rtol;
    maxstepsB = other.maxstepsB;
    newton_maxsteps = other.newton_maxsteps;
    newton_maxlinsteps = other.newton_maxlinsteps;
    newton_preeq = other.newton_preeq;
    newton_precon = other.newton_precon;
    ism = other.ism;
    sensi_meth = other.sensi_meth;
    linsol = other.linsol;
    interpType = other.interpType;
    lmm = other.lmm;
    iter = other.iter;
    stldet = other.stldet;

    if(other.x0data) {
        x0data = new double[other.nx];
        std::copy(other.x0data, other.x0data + other.nx, x0data);
    }

    if(other.sx0data) {
        sx0data = new double[other.nx * other.nplist];
        std::copy(other.sx0data, other.sx0data + other.nx * other.nplist, sx0data);
    }

    ordering = other.ordering;
    newton_precon = other.newton_precon;
    ism = other.ism;
}


void UserData::unscaleParameters(double *bufferUnscaled) const {
    /**
     * unscaleParameters removes parameter scaling according to the parameter
     * scaling in pscale
     *
     * @param[out] bufferUnscaled unscaled parameters are written to the array
     * @type double
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
}

void UserData::setTimepoints(const double *timepoints, int numTimepoints) {
    if (ts) {
        delete[] ts;
    }

    nt = numTimepoints;
    ts = new double[nt];
    memcpy(ts, timepoints, sizeof(double) * nt);
}

void UserData::setParameters(const double *parameters) {
    if(parameters)
        memcpy(p, parameters, sizeof(double) * np);
}

void UserData::setConstants(const double *constants) {
    if(constants)
        memcpy(k, constants, sizeof(double) * nk);
}

void UserData::setPlist(const double *plist, int nplist) {
    if (this->plist) {
        delete[] this->plist;
    }

    if(plist){
        this->nplist = nplist;
        this->plist = new int[nplist];
        
        for (int ip = 0; ip < nplist; ip++) {
            this->plist[ip] = (int)plist[ip];
        }
    } else {
        this->plist = nullptr;
        nplist = 0;
    }
}

void UserData::setPlist(const int *plist, int nplist) {
    if (this->plist) {
        delete[] this->plist;
    }

    if(plist){
        this->nplist = nplist;
        this->plist = new int[nplist];
        memcpy(this->plist, plist, sizeof(int) * nplist);
    } else {
        this->plist = nullptr;
        nplist = 0;
    }
}

void UserData::requireSensitivitiesForAllParameters()
{
    nplist = np;
    plist = new int[nplist];
    for (int i = 0; i < nplist; ++i)
        plist[i] = i;

}

void UserData::setPbar(const double *parameterScaling) {
    if (pbar) {
        delete[] pbar;
    }

    if(parameterScaling) {
        pbar = new double[nplist];
        memcpy(pbar, parameterScaling, sizeof(double) * nplist);
    } else {
        pbar = nullptr;
    }
}
    
void UserData::setXbar(const double *stateScaling) {
    if (xbar) {
        delete[] xbar;
    }
    
    if(stateScaling) {
        xbar = new double[nx];
        memcpy(xbar, stateScaling, sizeof(double) * nx);
    } else {
        xbar = nullptr;
    }
}

void UserData::setStateInitialization(const double *stateInitialization) {
    if (x0data) {
        delete[] x0data;
    }

    if(stateInitialization) {
        x0data = new double[nx];
        memcpy(x0data, stateInitialization, sizeof(double) * nx);
    } else {
        x0data = nullptr;
    }
}

void UserData::setSensitivityInitialization(
    const double *sensitivityInitialization) {
    if (sx0data) {
        delete[] sx0data;
    }

    if(sensitivityInitialization) {
        sx0data = new double[nx * nplist];
        memcpy(sx0data, sensitivityInitialization, sizeof(double) * nx * nplist);
    } else {
        sx0data = nullptr;
    }
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
    printf("quad_atol: %e\n", quad_atol);
    printf("quad_rtol: %e\n", quad_rtol);
    printf("maxstepsB: %d\n", maxstepsB);
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

} // namespace amici
