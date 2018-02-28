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

ReturnData::ReturnData(const Solver &solver, const Model *model)
    : ReturnData(solver, model, true) {
    /**
     * @brief constructor that uses information from model and solver to
     * appropriately initialize fields
     * @param[in] solver
     * @param[in] model pointer to model specification object @type Model
     */
}

ReturnData::ReturnData(const Solver &solver, const Model *model,
                       bool initializeFields)
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
     * @param[in] solver solver
     * @param[in] model pointer to model specification object @type Model
     * @param[in] initializeFields flag to initialize arrays (needs to happen elsewhere if false) @type
     * bool
     */

    if (initializeFields) {
        initFields();
        auto const t = model->getTimepoints();
        std::copy(t.begin(), t.end(), ts);
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
        if(ts[it_start]>t)
            break;
    
    for (int it = it_start; it < nt; it++){
        for (int ix = 0; ix < nx; ix++)
            x[ix * nt + it] = getNaN();
        for (int iy = 0; iy < ny; iy++)
            y[iy * nt + it] = getNaN();
        for (int ip = 0; ip < np; ip++) {
            if(sx)
                for (int ix = 0; ix < nx; ix++)
                    sx[(ip*nx + ix) * nt + it] = getNaN();
            if(sy)
                for (int iy = 0; iy < ny; iy++)
                    sy[(ip*ny + iy) * nt + it] = getNaN();
        }
    }
}
    
void ReturnData::invalidateLLH() {
    /**
     * @brief routine to set likelihood and respective sensitivities to NaN
     * (typically after integration failure)
     */
    if (llh)
        *llh = getNaN();
    
    if (sllh)
        setLikelihoodSensitivityFirstOrderNaN();
    
    if (s2llh)
        setLikelihoodSensitivitySecondOrderNaN();
}

void ReturnData::setLikelihoodSensitivityFirstOrderNaN() {
    /**
     * @brief routine to set first order sensitivities to NaN (typically after
     * integration failure)
     */
    fillArray(sllh, nplist, getNaN());
}

void ReturnData::setLikelihoodSensitivitySecondOrderNaN() {
    /**
     * @brief routine to set second order sensitivities to NaN (typically after
     * integration failure)
     */
    fillArray(s2llh, nplist * (nJ - 1), getNaN());
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
                pcoefficient[ip] = unscaledParameters[model->plist(ip)] * log(10);
            if (sensi == AMICI_SENSI_ORDER_SECOND)
                if (o2mode == AMICI_O2MODE_FULL)
                    for (int ip = 0; ip < np; ++ip)
                        augcoefficient[ip] = unscaledParameters[ip] * log(10);
            break;
        case AMICI_SCALING_LN:
            coefficient = 1.0;
            for (int ip = 0; ip < nplist; ++ip)
                pcoefficient[ip] = unscaledParameters[model->plist(ip)];
            if (sensi == AMICI_SENSI_ORDER_SECOND)
                if (o2mode == AMICI_O2MODE_FULL)
                    for (int ip = 0; ip < np; ++ip)
                        augcoefficient[ip] = unscaledParameters[ip];
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
                if (x)
                    if (sx)
                        for (int ip = 0; ip < nplist; ++ip)
                            for (int ix = 0; ix < nxtrue; ++ix)
                                for (int it = 0; it < nt; ++it)
                                    sx[(ip * nxtrue + ix) * nt + it] =
                                        x[(nxtrue + ip * nxtrue + ix) * nt +
                                          it];

                if (y)
                    if (sy)
                        for (int ip = 0; ip < nplist; ++ip)
                            for (int iy = 0; iy < nytrue; ++iy)
                                for (int it = 0; it < nt; ++it)
                                    sy[(ip * nytrue + iy) * nt + it] =
                                        y[(nytrue + ip * nytrue + iy) * nt +
                                          it];

                if (z)
                    if (sz)
                        for (int ip = 0; ip < nplist; ++ip)
                            for (int iz = 0; iz < nztrue; ++iz)
                                for (int it = 0; it < nt; ++it)
                                    sz[(ip * nztrue + iz) * nt + it] =
                                        z[(nztrue + ip * nztrue + iz) * nt +
                                          it];
            }
        }

        if (sllh)
            for (int ip = 0; ip < nplist; ++ip)
                sllh[ip] *= pcoefficient[ip];

#define chainRule(QUANT, IND1, N1T, N1, IND2, N2)                              \
    if (s##QUANT)                                                              \
        for (int ip = 0; ip < nplist; ++ip)                                    \
            for (int IND1 = 0; IND1 < N1T; ++IND1)                             \
                for (int IND2 = 0; IND2 < N2; ++IND2) {                        \
                    s##QUANT[(ip * N1 + IND1) * N2 + IND2] *=                  \
                        pcoefficient[ip];                                      \
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
        if (s2llh) {
            if (sllh) {
                for (int ip = 0; ip < nplist; ++ip) {
                    for (int iJ = 1; iJ < nJ; ++iJ) {
                        s2llh[ip * nplist + (iJ - 1)] *=
                            pcoefficient[ip] * augcoefficient[iJ - 1];
                        if (model->plist(ip) == iJ - 1)
                            s2llh[ip * nplist + (iJ - 1)] +=
                                sllh[ip] * coefficient;
                    }
                }
            }
        }

#define s2ChainRule(QUANT, IND1, N1T, N1, IND2, N2)                            \
    if (s##QUANT)                                                              \
        for (int ip = 0; ip < nplist; ++ip)                                    \
            for (int iJ = 1; iJ < nJ; ++iJ)                                    \
                for (int IND1 = 0; IND1 < N1T; ++IND1)                         \
                    for (int IND2 = 0; IND2 < N2; ++IND2) {                    \
                        s##QUANT[(ip * N1 + iJ * N1T + IND1) * N2 + IND2] *=   \
                            pcoefficient[ip] * augcoefficient[iJ - 1];         \
                        if (model->plist(ip) == iJ - 1)                        \
                            s##QUANT[(ip * N1 + iJ * N1T + IND1) * N2 +        \
                                     IND2] +=                                  \
                                s##QUANT[(ip * N1 + IND1) * N2 + IND2] *       \
                                coefficient;                                   \
                    }

        s2ChainRule(x, ix, nxtrue, nx, it, nt);
        s2ChainRule(y, iy, nytrue, ny, it, nt);
        s2ChainRule(sigmay, iy, nytrue, ny, it, nt);
        s2ChainRule(z, iz, nztrue, nz, ie, nmaxevent);
        s2ChainRule(sigmaz, iz, nztrue, nz, ie, nmaxevent);
        s2ChainRule(rz, iz, nztrue, nz, ie, nmaxevent);
    }

    if (o2mode == AMICI_O2MODE_DIR) { // directional
        if (s2llh) {
            if (sllh) {
                for (int ip = 0; ip < nplist; ++ip) {
                    s2llh[ip] *= pcoefficient[ip];
                    s2llh[ip] += model->k()[nk - nplist + ip] * sllh[ip] /
                                 unscaledParameters[model->plist(ip)];
                }
            }
        }

#define s2vecChainRule(QUANT, IND1, N1T, N1, IND2, N2)                         \
    if (s##QUANT)                                                              \
        for (int ip = 0; ip < nplist; ++ip)                                    \
            for (int IND1 = 0; IND1 < N1T; ++IND1)                             \
                for (int IND2 = 0; IND2 < N2; ++IND2) {                        \
                    s##QUANT[(ip * N1 + N1T + IND1) * N2 + IND2] *=            \
                        pcoefficient[ip];                                      \
                    s##QUANT[(ip * N1 + N1T + IND1) * N2 + IND2] +=            \
                        model->k()[nk - nplist + ip] *                           \
                        s##QUANT[(ip * N1 + IND1) * N2 + IND2] /               \
                        unscaledParameters[model->plist(ip)];                  \
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

ReturnData::~ReturnData() {
    /**
     * @brief default destructor
     */
    if (!freeFieldsOnDestruction)
        return;

    if (ts)
        delete[] ts;
    if (xdot)
        delete[] xdot;
    if (J)
        delete[] J;
    if (z)
        delete[] z;
    if (sigmaz)
        delete[] sigmaz;
    if (sz)
        delete[] sz;
    if (ssigmaz)
        delete[] ssigmaz;
    if (rz)
        delete[] rz;
    if (srz)
        delete[] srz;
    if (s2rz)
        delete[] s2rz;
    if (x)
        delete[] x;
    if (sx)
        delete[] sx;
    if (y)
        delete[] y;
    if (sigmay)
        delete[] sigmay;
    if (sy)
        delete[] sy;
    if (ssigmay)
        delete[] ssigmay;
    if (numsteps)
        delete[] numsteps;
    if (numrhsevals)
        delete[] numrhsevals;
    if (numerrtestfails)
        delete[] numerrtestfails;
    if (numnonlinsolvconvfails)
        delete[] numnonlinsolvconvfails;
    if (order)
        delete[] order;
    if (numstepsB)
        delete[] numstepsB;
    if (numrhsevalsB)
        delete[] numrhsevalsB;
    if (numerrtestfailsB)
        delete[] numerrtestfailsB;
    if (numnonlinsolvconvfailsB)
        delete[] numnonlinsolvconvfailsB;
    if (sx0)
        delete[] sx0;
    if (x0)
        delete[] x0;
    if (newton_status)
        delete[] newton_status;
    if (newton_numsteps)
        delete[] newton_numsteps;
    if (newton_numlinsteps)
        delete[] newton_numlinsteps;
    if (newton_time)
        delete[] newton_time;
    if (llh)
        delete[] llh;
    if (sllh)
        delete[] sllh;
    if (s2llh)
        delete[] s2llh;
    if (chi2)
        delete[] chi2;
    if (status)
        delete[] status;
}

void ReturnData::initFields() {
    /**
     * @brief initialises sol object with the corresponding fields
     */
    initField1(&status, "status", 1);

    initField1(&ts, "t", nt);
    initField1(&llh, "llh", 1);
    initField1(&chi2, "chi2", 1);
    initField2(&numsteps, "numsteps", nt, 1);
    initField2(&numrhsevals, "numrhsevals", nt, 1);
    initField2(&numerrtestfails, "numerrtestfails", nt, 1);
    initField2(&numnonlinsolvconvfails, "numnonlinsolvconvfails", nt, 1);
    initField2(&order, "order", nt, 1);

    if ((nz > 0) & (ne > 0)) {
        initField2(&z, "z", nmaxevent, nz);
        initField2(&rz, "rz", nmaxevent, nz);
        initField2(&sigmaz, "sigmaz", nmaxevent, nz);
    }
    if (nx > 0) {
        initField2(&x, "x", nt, nx);
        initField2(&xdot, "xdot", 1, nx);
        initField2(&J, "J", nx, nx);
        initField2(&x0, "x0", 1, nx);
        initField2(&sx0, "sx0", nx, nplist);
        initField2(&newton_status, "newton_status", 1, 1);
        initField2(&newton_numsteps, "newton_numsteps", 1, 2);
        initField2(&newton_numlinsteps, "newton_numlinsteps", newton_maxsteps,
                   2);
        initField2(&newton_time, "newton_time", 1, 1);
    }
    if (ny > 0) {
        initField2(&y, "y", nt, ny);
        initField2(&sigmay, "sigmay", nt, ny);
    }
    if (sensi >= AMICI_SENSI_ORDER_FIRST) {
        initField2(&sllh, "sllh", nplist, 1);

        if (sensi_meth == AMICI_SENSI_FSA) {
            initField3(&sx, "sx", nt, nx, nplist);
            if (ny > 0) {
                initField3(&sy, "sy", nt, ny, nplist);
                initField3(&ssigmay, "ssigmay", nt, ny, nplist);
            }
            if ((nz > 0) & (ne > 0)) {
                initField3(&srz, "srz", nmaxevent, nz, nplist);
                if (sensi >= AMICI_SENSI_ORDER_SECOND) {
                    initField4(&s2rz, "s2rz", nmaxevent, nztrue, nplist,
                               nplist);
                }
                initField3(&sz, "sz", nmaxevent, nz, nplist);
                initField3(&ssigmaz, "ssigmaz", nmaxevent, nz, nplist);
            }
        }

        if (sensi_meth == AMICI_SENSI_ASA) {
            if (ny > 0) {
                initField3(&ssigmay, "ssigmay", nt, ny, nplist);
            }
            if ((nz > 0) & (ne > 0)) {
                initField3(&ssigmaz, "ssigmaz", nmaxevent, nz, nplist);
            }
            initField2(&numstepsB, "numstepsB", nt, 1);
            initField2(&numrhsevalsB, "numrhsevalsB", nt, 1);
            initField2(&numerrtestfailsB, "numerrtestfailsB", nt, 1);
            initField2(&numnonlinsolvconvfailsB, "numnonlinsolvconvfailsB", nt,
                       1);
        }

        if (sensi >= AMICI_SENSI_ORDER_SECOND) {
            initField2(&s2llh, "s2llh", nJ - 1, nplist);
        }
    }
}

void ReturnData::initField1(double **fieldPointer, const char *fieldName,
                            int dim) {
    /**
     * @brief initialise vector and attach to the field
     * @param fieldPointer pointer of the field to which the vector will be
     * attached
     * @param fieldName Name of the field to which the vector will be attached
     * @param dim number of elements in the vector
     */
    *fieldPointer = new double[dim]();
}

void ReturnData::initField2(double **fieldPointer, const char *fieldName,
                            int dim1, int dim2) {
    /**
     * @brief initialise matrix and attach to the field
     * @param fieldPointer pointer of the field to which the matrix will be
     * attached
     * @param fieldName Name of the field to which the matrix will be attached
     * @param dim1 number of rows in the matrix
     * @param dim2 number of columns in the matrix
     */
    *fieldPointer = new double[(dim1) * (dim2)]();
}

void ReturnData::initField3(double **fieldPointer, const char *fieldName,
                            int dim1, int dim2, int dim3) {
    /**
     * @brief initialise 3D tensor and attach to the field
     * @param fieldPointer pointer of the field to which the tensor will be
     * attached
     * @param fieldName Name of the field to which the tensor will be attached
     * @param dim1 number of rows in the tensor
     * @param dim2 number of columns in the tensor
     * @param dim3 number of elements in the third dimension of the tensor
     */

    *fieldPointer = new double[(dim1) * (dim2) * (dim3)]();
}

void ReturnData::initField4(double **fieldPointer, const char *fieldName,
                            int dim1, int dim2, int dim3, int dim4) {
    /**
     * @brief initialise 4D tensor and attach to the field
     * @param fieldPointer pointer of the field to which the tensor will be
     * attached
     * @param fieldName Name of the field to which the tensor will be attached
     * @param dim1 number of rows in the tensor
     * @param dim2 number of columns in the tensor
     * @param dim3 number of elements in the third dimension of the tensor
     * @param dim4 number of elements in the fourth dimension of the tensor
     */
    *fieldPointer = new double[(dim1) * (dim2) * (dim3) * (dim4)]();
}

} // namespace amici
