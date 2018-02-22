#include "include/rdata.h"
#include "include/amici_misc.h"
#include "include/amici_model.h"
#include "include/symbolic_functions.h"
#include "include/amici_solver.h"
#include <cstring>

namespace amici {

ReturnData::ReturnData()
    /**
     * @brief default constructor
     */
    : np(0), nk(0), nx(0), nxtrue(0), ny(0), nytrue(0), nz(0), nztrue(0), ne(0),
      nJ(0), nplist(0), nmaxevent(0), nt(0), newton_maxsteps(0),
      pscale(AMICI_SCALING_NONE), o2mode(AMICI_O2MODE_NONE),
      sensi(AMICI_SENSI_ORDER_NONE), sensi_meth(AMICI_SENSI_NONE) {}

ReturnData::ReturnData(Solver const& solver, const Model *model)
    : np(model->np()), nk(model->nk()), nx(model->nx), nxtrue(model->nxtrue),
      ny(model->ny), nytrue(model->nytrue), nz(model->nz),
      nztrue(model->nztrue), ne(model->ne), nJ(model->nJ),
      nplist(model->nplist()), nmaxevent(model->nMaxEvent()), nt(model->nt()),
      newton_maxsteps(solver.getNewtonMaxSteps()), pscale(model->getParameterScale()),
      o2mode(model->o2mode), sensi(solver.getSensitivityOrder()),
      sensi_meth(static_cast<AMICI_sensi_meth>(solver.getSensitivityMethod())) {
    /**
     * @brief constructor that uses information from model and solver to
     * appropriately initialize fields
     * @param solver solver
     * @param model pointer to model specification object @type Model
     * @param initializeFields flag to initialize arrays (needs to happen elsewhere if false) @type
     * bool
     */

    ts = model->getTimepoints();

    xdot.resize(nx, getNaN());

    J.resize(nx * nx, getNaN());

    z.resize(nmaxevent * nz, getNaN());
    sigmaz.resize(nmaxevent * nz, getNaN());
    sz.resize(nmaxevent * nz * nplist, getNaN());
    ssigmaz.resize(nmaxevent * nz * nplist, getNaN());
    rz.resize(nmaxevent * nz, getNaN());
    srz.resize(nmaxevent * nz * nplist, getNaN());
    s2rz.resize(nmaxevent * nz * nplist * nplist, getNaN());

    x.resize(nt * nx, getNaN());
    sx.resize(nt * nx * nplist, getNaN());
    
    y.resize(nt * model->ny, getNaN());
    sigmay.resize(nt * model->ny, getNaN());
    res.clear();
    sy.resize(nt * model->ny * nplist, getNaN());
    ssigmay.resize(nt * model->ny * nplist, getNaN());
    sres.clear();

    numsteps.resize(nt, getNaN());
    numstepsB.resize(nt, getNaN());
    numrhsevals.resize(nt, getNaN());
    numrhsevalsB.resize(nt, getNaN());
    numerrtestfails.resize(nt, getNaN());
    numerrtestfailsB.resize(nt, getNaN());
    numnonlinsolvconvfails.resize(nt, getNaN());
    numnonlinsolvconvfailsB.resize(nt, getNaN());
    order.resize(nt, getNaN());
    newton_numsteps.resize(2, getNaN());
    newton_numlinsteps.resize(newton_maxsteps*2, getNaN());

    x0.resize(nx, getNaN());
    sx0.resize(nx * nplist, getNaN());
    
    sllh.resize(nplist, getNaN());
    s2llh.resize(nplist * (model->nJ - 1), getNaN());
}

void ReturnData::invalidate(const realtype t) {
    /**
     * @brief routine to set likelihood, state variables, outputs and respective sensitivities to NaN
     * (typically after integration failure)
     * @param t time of integration failure
     */
    invalidateLLH();
    
    // find it corresponding to datapoint after integration failure
    int it_start;
    for (it_start = 0; it_start < nt; it_start++)
        if(ts.at(it_start)>t)
            break;
    
    for (int it = it_start; it < nt; it++){
        for (int ix = 0; ix < nx; ix++)
            x.at(ix * nt + it) = getNaN();
        for (int iy = 0; iy < ny; iy++)
            y.at(iy * nt + it) = getNaN();
        for (int ip = 0; ip < nplist; ip++) {
            for (int ix = 0; ix < nx; ix++)
                sx.at((ip*nx + ix) * nt + it) = getNaN();
            for (int iy = 0; iy < ny; iy++)
                sy.at((ip*ny + iy) * nt + it) = getNaN();
        }
    }
}
    
void ReturnData::invalidateLLH() {
    /**
     * @brief routine to set likelihood and respective sensitivities to NaN
     * (typically after integration failure)
     */

    llh = getNaN();
    chi2 = getNaN();
    std::fill(sllh.begin(),sllh.end(),getNaN());
    std::fill(s2llh.begin(),s2llh.end(),getNaN());
}

void ReturnData::applyChainRuleFactorToSimulationResults(const Model *model) {
    /**
     * @brief applies the chain rule to account for parameter transformation
     * in the sensitivities of simulation results
     * @param model Model from which the ReturnData was obtained
     */
    if (pscale == AMICI_SCALING_NONE)
        return;

    // chain-rule factor: multiplier for am_p
    realtype coefficient;

    std::vector<realtype> pcoefficient(nplist);
    std::vector<realtype> unscaledParameters(np);
    model->unscaleParameters(unscaledParameters.data());
    std::vector<realtype> augcoefficient(np);

    switch (pscale) {
        case AMICI_SCALING_LOG10:
            coefficient = log(10.0);
            for (int ip = 0; ip < nplist; ++ip)
                pcoefficient.at(ip) = unscaledParameters[model->plist(ip)] * log(10);
            if (sensi == AMICI_SENSI_ORDER_SECOND)
                if (o2mode == AMICI_O2MODE_FULL)
                    for (int ip = 0; ip < np; ++ip)
                        augcoefficient.at(ip) = unscaledParameters.at(ip) * log(10);
            break;
        case AMICI_SCALING_LN:
            coefficient = 1.0;
            for (int ip = 0; ip < nplist; ++ip)
                pcoefficient.at(ip) = unscaledParameters[model->plist(ip)];
            if (sensi == AMICI_SENSI_ORDER_SECOND)
                if (o2mode == AMICI_O2MODE_FULL)
                    for (int ip = 0; ip < np; ++ip)
                        augcoefficient.at(ip) = unscaledParameters.at(ip);
            break;
        case AMICI_SCALING_NONE:
            // this should never be reached
            break;
    }
    

    if (sensi >= AMICI_SENSI_ORDER_FIRST) {
        // recover first order sensitivies from states for adjoint sensitivity
        // analysis
        if (sensi == AMICI_SENSI_ORDER_SECOND) {
            if (sensi_meth == AMICI_SENSI_ASA) {
                for (int ip = 0; ip < nplist; ++ip)
                    for (int ix = 0; ix < nxtrue; ++ix)
                        for (int it = 0; it < nt; ++it)
                            sx.at((ip * nxtrue + ix) * nt + it) =
                                x.at((nxtrue + ip * nxtrue + ix) * nt +
                                  it);

                for (int ip = 0; ip < nplist; ++ip)
                    for (int iy = 0; iy < nytrue; ++iy)
                        for (int it = 0; it < nt; ++it)
                            sy.at((ip * nytrue + iy) * nt + it) =
                                y.at((nytrue + ip * nytrue + iy) * nt +
                                  it);

                for (int ip = 0; ip < nplist; ++ip)
                    for (int iz = 0; iz < nztrue; ++iz)
                        for (int it = 0; it < nt; ++it)
                            sz.at((ip * nztrue + iz) * nt + it) =
                                z.at((nztrue + ip * nztrue + iz) * nt +
                                  it);
            }
        }

        for (int ip = 0; ip < nplist; ++ip)
            sllh.at(ip) *= pcoefficient.at(ip);

#define chainRule(QUANT, IND1, N1T, N1, IND2, N2)                               \
    for (int ip = 0; ip < nplist; ++ip)                                         \
        for (int IND1 = 0; IND1 < N1T; ++IND1)                                  \
            for (int IND2 = 0; IND2 < N2; ++IND2) {                             \
                s##QUANT.at((ip * N1 + IND1) * N2 + IND2) *=                    \
                    pcoefficient.at(ip);                                        \
            }

        chainRule(x, ix, nxtrue, nx, it, nt);
        chainRule(y, iy, nytrue, ny, it, nt);
        chainRule(sigmay, iy, nytrue, ny, it, nt);
        chainRule(z, iz, nztrue, nz, ie, nmaxevent);
        chainRule(sigmaz, iz, nztrue, nz, ie, nmaxevent);
        chainRule(rz, iz, nztrue, nz, ie, nmaxevent);
        chainRule(x0, ix, nxtrue, nx, it, 1);
    }

    if (o2mode == AMICI_O2MODE_FULL) { // full
        for (int ip = 0; ip < nplist; ++ip) {
            for (int iJ = 1; iJ < nJ; ++iJ) {
                s2llh[ip * nplist + (iJ - 1)] *=
                    pcoefficient.at(ip) * augcoefficient[iJ - 1];
                if (model->plist(ip) == iJ - 1)
                    s2llh[ip * nplist + (iJ - 1)] +=
                        sllh.at(ip) * coefficient;
            }
        }

#define s2ChainRule(QUANT, IND1, N1T, N1, IND2, N2)                             \
    for (int ip = 0; ip < nplist; ++ip)                                         \
        for (int iJ = 1; iJ < nJ; ++iJ)                                         \
            for (int IND1 = 0; IND1 < N1T; ++IND1)                              \
                for (int IND2 = 0; IND2 < N2; ++IND2) {                         \
                    s##QUANT.at((ip * N1 + iJ * N1T + IND1) * N2 + IND2) *=     \
                        pcoefficient.at(ip) * augcoefficient[iJ - 1];           \
                    if (model->plist(ip) == iJ - 1)                             \
                        s##QUANT[(ip * N1 + iJ * N1T + IND1) * N2 +             \
                                 IND2] +=                                       \
                            s##QUANT.at((ip * N1 + IND1) * N2 + IND2) *         \
                            coefficient;                                        \
                }

        s2ChainRule(x, ix, nxtrue, nx, it, nt);
        s2ChainRule(y, iy, nytrue, ny, it, nt);
        s2ChainRule(sigmay, iy, nytrue, ny, it, nt);
        s2ChainRule(z, iz, nztrue, nz, ie, nmaxevent);
        s2ChainRule(sigmaz, iz, nztrue, nz, ie, nmaxevent);
        s2ChainRule(rz, iz, nztrue, nz, ie, nmaxevent);
    }

    if (o2mode == AMICI_O2MODE_DIR) { // directional
            for (int ip = 0; ip < nplist; ++ip) {
                s2llh.at(ip) *= pcoefficient.at(ip);
                s2llh.at(ip) += model->k()[nk - nplist + ip] * sllh.at(ip) /
                             unscaledParameters[model->plist(ip)];
        }

#define s2vecChainRule(QUANT, IND1, N1T, N1, IND2, N2)                          \
    for (int ip = 0; ip < nplist; ++ip)                                         \
        for (int IND1 = 0; IND1 < N1T; ++IND1)                                  \
            for (int IND2 = 0; IND2 < N2; ++IND2) {                             \
                s##QUANT.at((ip * N1 + N1T + IND1) * N2 + IND2) *=              \
                    pcoefficient.at(ip);                                        \
                s##QUANT.at((ip * N1 + N1T + IND1) * N2 + IND2) +=              \
                    model->k()[nk - nplist + ip] *                              \
                    s##QUANT.at((ip * N1 + IND1) * N2 + IND2) /                 \
                    unscaledParameters[model->plist(ip)];                       \
            }

        s2vecChainRule(x, ix, nxtrue, nx, it, nt);
        s2vecChainRule(y, iy, nytrue, ny, it, nt);
        s2vecChainRule(sigmay, iy, nytrue, ny, it, nt);
        s2vecChainRule(z, iz, nztrue, nz, ie, nmaxevent);
        s2vecChainRule(sigmaz, iz, nztrue, nz, ie, nmaxevent);
        s2vecChainRule(rz, iz, nztrue, nz, ie, nmaxevent);
    }
    return;
}

} // namespace amici
