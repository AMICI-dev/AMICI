#include "amici/rdata.h"

#include "amici/misc.h"
#include "amici/model.h"
#include "amici/symbolic_functions.h"
#include "amici/solver.h"
#include "amici/exception.h"

#include <cstring>

namespace amici {

ReturnData::ReturnData()
    /**
     * @brief default constructor
     */
    : np(0), nk(0), nx(0), nxtrue(0), ny(0), nytrue(0), nz(0), nztrue(0), ne(0),
      nJ(0), nplist(0), nmaxevent(0), nt(0), newton_maxsteps(0),
    pscale(std::vector<ParameterScaling>(0, ParameterScaling::none)), o2mode(SecondOrderMode::none),
      sensi(SensitivityOrder::none), sensi_meth(SensitivityMethod::none) {}

ReturnData::ReturnData(Solver const& solver, const Model *model)
    : ts(model->getTimepoints()),
      np(model->np()), nk(model->nk()), nx(model->nx), nxtrue(model->nxtrue),
      ny(model->ny), nytrue(model->nytrue), nz(model->nz),
      nztrue(model->nztrue), ne(model->ne), nJ(model->nJ),
      nplist(model->nplist()), nmaxevent(model->nMaxEvent()), nt(model->nt()),
      newton_maxsteps(solver.getNewtonMaxSteps()), pscale(model->getParameterScale()),
      o2mode(model->o2mode), sensi(solver.getSensitivityOrder()),
      sensi_meth(static_cast<SensitivityMethod>(solver.getSensitivityMethod()))
    {
    /**
     * @brief constructor that uses information from model and solver to
     * appropriately initialize fields
     * @param solver solver
     * @param model pointer to model specification object
     * bool
     */
          
    xdot.resize(nx, getNaN());

    J.resize(nx * nx, getNaN());

    // initialize with 0.0, so we only need to write non-zero values
    z.resize(nmaxevent * nz, 0.0);
    sigmaz.resize(nmaxevent * nz, 0.0);

    rz.resize(nmaxevent * nz, 0.0);
    x.resize(nt * nx, 0.0);
    y.resize(nt * model->ny, 0.0);
    sigmay.resize(nt * model->ny, 0.0);

    newton_numsteps.resize(2, 0);
    newton_numlinsteps.resize(newton_maxsteps*2, 0);

    if(nt>0) {
        numsteps.resize(nt, 0);
        numrhsevals.resize(nt, 0);
        numerrtestfails.resize(nt, 0);
        numnonlinsolvconvfails.resize(nt, 0);
        order.resize(nt, 0);
        
        if (sensi_meth == SensitivityMethod::adjoint && sensi >= SensitivityOrder::first) {
            numstepsB.resize(nt, 0);
            numrhsevalsB.resize(nt, 0);
            numerrtestfailsB.resize(nt, 0);
            numnonlinsolvconvfailsB.resize(nt, 0);
        }
    }

    x0.resize(nx, getNaN());
    
    res.resize(nt * model->nytrue, 0.0);

    llh = getNaN();
    chi2 = getNaN();
    if (sensi >= SensitivityOrder::first){
        sllh.resize(nplist, getNaN());
        sx0.resize(nx * nplist, getNaN());
        
        if (sensi_meth == SensitivityMethod::forward || sensi >= SensitivityOrder::second){
            // for second order we can fill in from the augmented states
            sx.resize(nt * nx * nplist, 0.0);
            sy.resize(nt * ny * nplist, 0.0);
            sz.resize(nmaxevent * nz * nplist, 0.0);
            srz.resize(nmaxevent * nz * nplist, 0.0);
            
            FIM.resize(nplist * nplist, 0.0);
            sres.resize(nt * model->nytrue * nplist, 0.0);
        }
        
        ssigmay.resize(nt * model->ny * nplist, 0.0);
        ssigmaz.resize(nmaxevent * nz * nplist, 0.0);
        if (sensi >= SensitivityOrder::second) {
            s2llh.resize(nplist * (model->nJ - 1), getNaN());
            if (sensi_meth == SensitivityMethod::forward)
                s2rz.resize(nmaxevent * nztrue * nplist * nplist, 0.0);
        }
    }
    
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
            x.at(ix + nx * it) = getNaN();
        for (int iy = 0; iy < ny; iy++)
            y.at(iy + ny * it) = getNaN();
    }

    if (!sx.empty()) {
        for (int it = it_start; it < nt; it++){
            for (int ip = 0; ip < nplist; ip++) {
                for (int ix = 0; ix < nx; ix++)
                    sx.at(ix + nx*(ip + it*nplist)) = getNaN();
            }
        }
    }
    if(!sy.empty()) {
        for (int it = it_start; it < nt; it++){
            for (int ip = 0; ip < nplist; ip++) {
                for (int iy = 0; iy < ny; iy++)
                    sy.at(iy + ny*(ip + it*nplist)) = getNaN();
            }
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

    // chain-rule factor: multiplier for am_p
    std::vector<realtype> coefficient(nplist, 1.0);

    std::vector<realtype> pcoefficient(nplist, 1.0);
    std::vector<realtype> unscaledParameters(np);
    model->unscaleParameters(unscaledParameters.data());
    std::vector<realtype> augcoefficient(np, 1.0);
    
    if (sensi == SensitivityOrder::second && o2mode == SecondOrderMode::full) {
        for (int ip = 0; ip < np; ++ip) {
            switch (pscale[ip]) {
            case ParameterScaling::log10:
                augcoefficient.at(ip) = unscaledParameters.at(ip) * log(10);
                break;
            case ParameterScaling::ln:
                augcoefficient.at(ip) = unscaledParameters.at(ip);
                break;
            case ParameterScaling::none:
                break;
            }
        }
    }

    for (int ip = 0; ip < nplist; ++ip) {
        switch (pscale[model->plist(ip)]) {
        case ParameterScaling::log10:
            coefficient.at(ip) = log(10.0);
            pcoefficient.at(ip) = unscaledParameters.at(model->plist(ip)) * log(10);
            break;
        case ParameterScaling::ln:
            pcoefficient.at(ip) = unscaledParameters.at(model->plist(ip));
            break;
        case ParameterScaling::none:
            break;
        }
    }

    if (sensi >= SensitivityOrder::first) {
        // recover first order sensitivies from states for adjoint sensitivity
        // analysis
        if (sensi == SensitivityOrder::second && o2mode == SecondOrderMode::full) {
            if (sensi_meth == SensitivityMethod::adjoint) {
                for (int ip = 0; ip < nplist; ++ip)
                    for (int ix = 0; ix < nxtrue; ++ix)
                        for (int it = 0; it < nt; ++it)
                            sx.at(ix + nxtrue*(ip + it*nplist)) =
                                x.at(it * nx + nxtrue + ip * nxtrue + ix);

                for (int ip = 0; ip < nplist; ++ip)
                    for (int iy = 0; iy < nytrue; ++iy)
                        for (int it = 0; it < nt; ++it)
                            sy.at(iy + nytrue*(ip + it*nplist)) =
                                y.at(it * ny + nytrue + ip * nytrue + iy);

                for (int ip = 0; ip < nplist; ++ip)
                    for (int iz = 0; iz < nztrue; ++iz)
                        for (int it = 0; it < nt; ++it)
                            sz.at(iz + nztrue*(ip + it*nplist)) =
                                z.at(it * nz + nztrue + ip * nztrue + iz);
            }
        }

        for (int ip = 0; ip < nplist; ++ip)
            sllh.at(ip) *= pcoefficient.at(ip);
        
        if(!sres.empty())
            for (int iyt = 0; iyt < nytrue*nt; ++iyt)
                for (int ip = 0; ip < nplist; ++ip)
                    sres.at((iyt * nplist + ip)) *= pcoefficient.at(ip);
        
        if(!FIM.empty())
            for (int ip = 0; ip < nplist; ++ip)
                for (int jp = 0; jp < nplist; ++jp)
                    FIM.at(jp + ip * nplist) *= pcoefficient.at(ip)*pcoefficient.at(jp);

#define chainRule(QUANT, IND1, N1T, N1, IND2, N2)                               \
    if (!s##QUANT.empty())                                                        \
        for (int IND1 = 0; IND1 < N1T; ++IND1)                                  \
            for (int ip = 0; ip < nplist; ++ip)                                 \
                for (int IND2 = 0; IND2 < N2; ++IND2) {                         \
                    s##QUANT.at((IND2*nplist + ip)*N1 + IND1) *=                \
                        pcoefficient.at(ip);                                    \
                }

        chainRule(x, ix, nxtrue, nx, it, nt);
        chainRule(y, iy, nytrue, ny, it, nt);
        chainRule(sigmay, iy, nytrue, ny, it, nt);
        chainRule(z, iz, nztrue, nz, ie, nmaxevent);
        chainRule(sigmaz, iz, nztrue, nz, ie, nmaxevent);
        chainRule(rz, iz, nztrue, nz, ie, nmaxevent);
        chainRule(x0, ix, nxtrue, nx, it, 1);
    }

    if (o2mode == SecondOrderMode::full) { // full
        if (!s2llh.empty() && !sllh.empty()) {
            for (int ip = 0; ip < nplist; ++ip) {
                for (int iJ = 1; iJ < nJ; ++iJ) {
                    s2llh[ip * nplist + (iJ - 1)] *=
                            pcoefficient.at(ip) * augcoefficient[iJ - 1];
                    if (model->plist(ip) == iJ - 1)
                        s2llh[ip * nplist + (iJ - 1)] +=
                                sllh.at(ip) * coefficient.at(ip);
                }
            }
        }

#define s2ChainRule(QUANT, IND1, N1T, N1, IND2, N2)                            \
    if (!s##QUANT.empty())                                                       \
        for (int ip = 0; ip < nplist; ++ip)                                    \
            for (int iJ = 1; iJ < nJ; ++iJ)                                    \
                for (int IND1 = 0; IND1 < N1T; ++IND1)                         \
                    for (int IND2 = 0; IND2 < N2; ++IND2) {                    \
                    s##QUANT.at((IND2*nplist + ip)*N1 + IND1 + iJ*N1T) *=      \
                        pcoefficient.at(ip) * augcoefficient[iJ - 1];          \
                        if (model->plist(ip) == iJ - 1)                        \
                            s##QUANT.at((IND2*nplist + ip)*N1 + IND1 + iJ*N1T) +=   \
                                s##QUANT.at((IND2*nplist + ip)*N1 + IND1) *    \
                                    coefficient[ip];                           \
}


        s2ChainRule(x, ix, nxtrue, nx, it, nt);
        s2ChainRule(y, iy, nytrue, ny, it, nt);
        s2ChainRule(sigmay, iy, nytrue, ny, it, nt);
        s2ChainRule(z, iz, nztrue, nz, ie, nmaxevent);
        s2ChainRule(sigmaz, iz, nztrue, nz, ie, nmaxevent);
        s2ChainRule(rz, iz, nztrue, nz, ie, nmaxevent);
    }

    if (o2mode == SecondOrderMode::directional) { // directional
            for (int ip = 0; ip < nplist; ++ip) {
                s2llh.at(ip) *= pcoefficient.at(ip);
                s2llh.at(ip) += model->k()[nk - nplist + ip] * sllh.at(ip) /
                             unscaledParameters[model->plist(ip)];
        }

#define s2vecChainRule(QUANT, IND1, N1T, N1, IND2, N2)                          \
    if (!s##QUANT.empty())                                                        \
        for (int ip = 0; ip < nplist; ++ip)                                     \
            for (int IND1 = 0; IND1 < N1T; ++IND1)                              \
                for (int IND2 = 0; IND2 < N2; ++IND2) {                         \
                    s##QUANT.at((IND2*nplist + ip)*N1 + IND1 + N1T) *=          \
                        pcoefficient.at(ip);                                    \
                    s##QUANT.at((IND2*nplist + ip)*N1 + IND1 + N1T) +=          \
                        model->k()[nk - nplist + ip] *                          \
                        s##QUANT.at((IND2*nplist + ip)*N1 + IND1) /             \
                        unscaledParameters[model->plist(ip)];                   \
                }

        s2vecChainRule(x, ix, nxtrue, nx, it, nt);
        s2vecChainRule(y, iy, nytrue, ny, it, nt);
        s2vecChainRule(sigmay, iy, nytrue, ny, it, nt);
        s2vecChainRule(z, iz, nztrue, nz, ie, nmaxevent);
        s2vecChainRule(sigmaz, iz, nztrue, nz, ie, nmaxevent);
        s2vecChainRule(rz, iz, nztrue, nz, ie, nmaxevent);
    }
}

void ReturnData::initializeObjectiveFunction()
{
    llh = 0.0;
    chi2 = 0.0;
    std::fill(sllh.begin(),sllh.end(), 0.0);
    std::fill(s2llh.begin(),s2llh.end(), 0.0);
}

} // namespace amici
