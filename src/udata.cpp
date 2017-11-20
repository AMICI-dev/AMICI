#include "include/udata.h"
#include "include/amici_exception.h"
#include <cstdio>
#include <cstring>
#include <algorithm>

namespace amici {

UserData::UserData(const int np, const int nk, const int nx) : sizex(nx) {
    // these fields must always be initialized, others are optional or can be set later
    konst.resize(nk);
    par.resize(np);
    unpar.resize(np);
    qpositivex.assign(nx,1);
}

UserData::UserData() : sizex(0) {
    konst.resize(0);
    par.resize(0);
    qpositivex.assign(0,1);
}

UserData::UserData(const UserData &other) : UserData(other.np(), other.nk(), other.nx())
{
    nmaxevent = other.nmaxevent;
    qpositivex = other.qpositivex;
    p_index = other.p_index;
    par = other.par;
    unpar = other.unpar;
    konst = other.konst;
    pscale = other.pscale;
    tstart = other.tstart;
    ts = other.ts;
    pbar = other.pbar;
    xbar = other.xbar;
    sensi = other.sensi;
    atol = other.atol;
    rtol = other.rtol;
    maxsteps = other.maxsteps;
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
    x0data = other.x0data;
    sx0data = other.sx0data;
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
        for (int ip = 0; ip < np(); ++ip) {
            bufferUnscaled[ip] = pow(10, par[ip]);
        }
        break;
    case AMICI_SCALING_LN:
        for (int ip = 0; ip < np(); ++ip)
            bufferUnscaled[ip] = exp(par[ip]);
        break;
    case AMICI_SCALING_NONE:
        for (int ip = 0; ip < np(); ++ip)
            bufferUnscaled[ip] = par[ip];
        break;
    }
}

void UserData::setTimepoints(const double * timepoints, int numTimepoints) {
    ts.resize(numTimepoints);
    memcpy(ts.data(), timepoints, sizeof(double) * numTimepoints);
}

void UserData::setParameters(const double * parameters) {
    memcpy(par.data(), parameters, sizeof(double) * np());
    unscaleParameters(unpar.data());
}

void UserData::setConstants(const double *constants) {
    memcpy(konst.data(), constants, sizeof(double) * nk());
}

void UserData::setPlist(const double *plist, int length) {
    if(!plist)
        throw AmiException("Provided plist was a nullptr, please provide a valid pointer.");
    p_index.resize(length);
    pbar.resize(nplist());
    std::fill(pbar.begin(),pbar.end(),1.0);
    for (int ip = 0; ip < length; ip++) {
        p_index[ip] = (int)plist[ip];
    }
}

void UserData::setPlist(const int *plist, int length) {
    if(!plist)
        throw AmiException("Provided plist was a nullptr, please provide a valid pointer.");
    p_index.resize(length);
    pbar.resize(nplist());
    std::fill(pbar.begin(),pbar.end(),1.0);
    memcpy(p_index.data(), plist, sizeof(int) * length);
}
    
void UserData::setPScale(const AMICI_parameter_scaling pscale) {
    this->pscale = pscale;
}

void UserData::requireSensitivitiesForAllParameters()
{
    p_index.resize(nplist());
    for (int i = 0; i < nplist(); ++i)
        p_index[i] = i;

}

void UserData::setPbar(const double *parameterScaling) {
    if(parameterScaling) {
        pbar.resize(nplist());
        memcpy(pbar.data(), parameterScaling, sizeof(double) * nplist());
    } else {
        pbar.clear();
    }
}
    
void UserData::setXbar(const double *stateScaling) {
    if(stateScaling){
        xbar.resize(nx());
        memcpy(xbar.data(), stateScaling, sizeof(double) * nx());
    } else {
        xbar.clear();
    }
}

void UserData::setStateInitialization(const double *stateInitialization) {
    if(stateInitialization) {
        x0data.resize(nx());
        memcpy(x0data.data(), stateInitialization, sizeof(double) * nx());
    } else {
        x0data.clear();
    }
}

void UserData::setSensitivityInitialization(
    const double *sensitivityInitialization) {
    if(sensitivityInitialization) {
        sx0data.resize(nx() * nplist());
        memcpy(sx0data.data(), sensitivityInitialization, sizeof(double) * nx() * nplist());
    } else {
        sx0data.clear();
    }
}

UserData::~UserData() {
}

void UserData::print() const {
    printf("qpositivex: %p\n", qpositivex.data());
    printf("plist: %p\n", p_index.data());
    printf("nplist: %d\n", nplist());
    printf("nt: %d\n", nt());
    printf("nmaxevent: %d\n", nmaxevent);
    printf("p: %p\n", par.data());
    printf("k: %p\n", konst.data());
    printf("tstart: %e\n", tstart);
    printf("ts: %p\n", ts.data());
    printf("pbar: %p\n", pbar.data());
    printf("xbar: %p\n", xbar.data());
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
    printf("x0data: %p\n", x0data.data());
    printf("sx0data: %p\n", sx0data.data());
    printf("ordering: %d\n", ordering);
}

} // namespace amici
