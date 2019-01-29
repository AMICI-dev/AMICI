#include "amici/model.h"
#include "amici/amici.h"
#include "amici/misc.h"
#include "amici/exception.h"
#include "amici/symbolic_functions.h"

#include <numeric>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <typeinfo>
#include <utility>
#include <regex>

namespace amici {

void Model::fsy(const int it, const AmiVectorArray *sx, ReturnData *rdata) {
    if (!ny)
        return;

    // copy dydp for current time to sy
    std::copy(dydp.begin(), dydp.end(), &rdata->sy[it * nplist() * ny]);

    sx->flatten_to_vector(this->sx);

    // compute sy = 1.0*dydx*sx + 1.0*sy
    // dydx A[ny,nx_solver] * sx B[nx_solver,nplist] = sy C[ny,nplist]
    //        M  K                 K  N                     M  N
    //        lda                  ldb                      ldc
    amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans, ny, nplist(), nx_solver,
                1.0, dydx.data(), ny, this->sx.data(), nx_solver, 1.0,
                &rdata->sy[it*nplist()*ny], ny);

    if(alwaysCheckFinite)
        checkFinite(nplist()*ny, &rdata->sy[it*nplist()*ny], "sy");
}

void Model::fsz_tf(const int *nroots, const int ie, ReturnData *rdata) {
    for (int iz = 0; iz < nz; ++iz)
        if (z2event[iz] - 1 == ie)
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->sz.at((nroots[ie]*nplist()+ip)*nz + iz) = 0.0;
}

void Model::fsJy(const int it, const std::vector<realtype>& dJydx, const AmiVectorArray *sx, ReturnData *rdata) {

    // Compute dJydx*sx for current 'it'
    // dJydx        rdata->nt x nJ        x nx_solver
    // sx           rdata->nt x nx_solver x nplist()
    std::vector<realtype> multResult(nJ * nplist(), 0);
    sx->flatten_to_vector(this->sx);

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans, nJ,
                nplist(), nx_solver, 1.0, &dJydx.at(it*nJ*nx_solver), nJ, this->sx.data(), nx_solver, 0.0,
                multResult.data(), nJ);

    // multResult    nJ        x nplist()
    // dJydp         nJ        x nplist()
    // dJydxTmp      nJ        x nx_solver
    // sx            nx_solver x nplist()

    // sJy += multResult + dJydp
    for (int iJ = 0; iJ < nJ; ++iJ) {
        if (iJ == 0)
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->sllh.at(ip) -= multResult.at(ip * nJ) + dJydp.at(ip * nJ);
        else
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->s2llh.at((iJ - 1) + ip * (nJ - 1)) -=
                        multResult.at(iJ + ip * nJ) + dJydp.at(iJ + ip * nJ);
    }
}

void Model::fdJydp(const int it, ReturnData *rdata, const ExpData *edata) {

    // dJydy         nJ x nytrue x ny
    // dydp          nplist * ny
    // dJydp         nplist x nJ
    // dJydsigma

    getmy(it,edata);
    std::fill(dJydp.begin(), dJydp.end(), 0.0);

    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (isNaN(my.at(iyt)))
            continue;

        // dJydp = 1.0 * dJydp +  1.0 * dJydy * dydp
        amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
                    nJ, nplist(), ny,
                    1.0, &dJydy.at(iyt*nJ*ny), nJ,
                    dydp.data(), ny,
                    1.0, dJydp.data(), nJ);

        // dJydp = 1.0 * dJydp +  1.0 * dJydsigma * dsigmaydp
        amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
                    nJ, nplist(), ny,
                    1.0, &dJydsigma.at(iyt*nJ*ny), nJ,
                    dsigmaydp.data(), ny,
                    1.0, dJydp.data(), nJ);
    }

    if (rdata->sensi_meth != SensitivityMethod::adjoint)
        return;

    if(!ny)
        return;

    for (int iJ = 0; iJ < nJ; iJ++) {
        for (int ip = 0; ip < nplist(); ip++) {
            if (iJ == 0) {
                rdata->sllh.at(ip) -= dJydp[ip * nJ];
            } else {
                rdata->s2llh.at((iJ - 1) + ip * (nJ - 1)) -=
                        dJydp[iJ + ip * nJ];
            }
        }
    }
}

void Model::fdJydx(std::vector<realtype> *dJydx, const int it, const ExpData *edata) {

    // dJydy         nJ x ny        x nytrue
    // dydx          ny x nx_solver
    // dJydx         nJ x nx_solver x nt
    getmy(it,edata);
    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (isNaN(my.at(iyt)))
            continue;
    // dJydy A[nyt,nJ,ny] * dydx B[ny,nx_solver] = dJydx C[it,nJ,nx_solver]
    //         slice                                       slice
    //             M  K            K  N                       M  N
    //             lda             ldb                        ldc
        amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
                    nJ, nx_solver, ny, 1.0, &dJydy.at(iyt*ny*nJ), nJ, dydx.data(), ny, 1.0,
                    &dJydx->at(it*nx_solver*nJ), nJ);
    }
}

void Model::fsJz(const int nroots, const std::vector<realtype>& dJzdx, const AmiVectorArray *sx, ReturnData *rdata) {
    // sJz           nJ x nplist()
    // dJzdp         nJ x nplist()
    // dJzdx         nmaxevent x nJ        x nx_solver
    // sx            rdata->nt x nx_solver x nplist()

    // Compute dJzdx*sx for current 'ie'
    // dJzdx        rdata->nt x nJ        x nx_solver
    // sx           rdata->nt x nx_solver x nplist()

    std::vector<realtype> multResult(nJ * nplist(), 0);
    sx->flatten_to_vector(this->sx);

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans, nJ,
                nplist(), nx_solver, 1.0, &dJzdx.at(nroots*nx_solver*nJ), nJ, this->sx.data(), nx_solver, 1.0,
                multResult.data(), nJ);

    // sJy += multResult + dJydp
    for (int iJ = 0; iJ < nJ; ++iJ) {
        if (iJ == 0)
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->sllh.at(ip) -= multResult.at(ip * nJ) + dJzdp.at(ip * nJ);
        else
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->s2llh.at((iJ - 1) + ip * (nJ - 1)) -=
                        multResult.at(iJ + ip * nJ) + dJzdp.at(iJ + ip * nJ);
    }
}

void Model::fdJzdp(const int nroots, realtype t, const ExpData *edata,
                   const ReturnData *rdata) {
    // dJzdz         nJ x nz x nztrue
    // dJzdsigma     nJ x nz x nztrue
    // dzdp          nz x nplist()
    // dJzdp         nJ x nplist()

    getmz(nroots,edata);
    std::fill(dJzdp.begin(),dJzdp.end(),0.0);
    for (int izt = 0; izt < nztrue; ++izt) {
        if (isNaN(mz.at(izt)))
            continue;

        if (t < rdata->ts.at(rdata->ts.size() - 1)) {
            // with z
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
                        nJ, nplist(), nz, 1.0, &dJzdz.at(izt*nz*nJ), nJ, dzdp.data(), nz,
                        1.0, dJzdp.data(), nJ);
        } else {
            // with rz
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nplist(), nz, 1.0,
                        &dJrzdsigma.at(izt*nz*nJ), nJ, dsigmazdp.data(), nz, 1.0,
                        dJzdp.data(), nJ);

            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
                        nJ, nplist(), nz, 1.0, &dJrzdz.at(izt*nz*nJ), nJ, dzdp.data(), nz,
                        1.0, dJzdp.data(), nJ);
        }

        amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
                    nJ, nplist(), nz, 1.0, &dJzdsigma.at(izt*nz*nJ), nJ,
                    dsigmazdp.data(), nz, 1.0, dJzdp.data(), nJ);
    }
}

void Model::fdJzdx(std::vector<realtype> *dJzdx, const int nroots, realtype t, const ExpData *edata, const ReturnData *rdata) {
    // dJzdz         nJ x nz        x nztrue
    // dzdx          nz x nx_solver
    // dJzdx         nJ x nx_solver x nmaxevent
    getmz(nroots,edata);
    for (int izt = 0; izt < nztrue; ++izt) {
        if (isNaN(mz.at(izt)))
            continue;

        if (t < rdata->ts.at(rdata->ts.size() - 1)) {
            // z
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nx_solver, nz, 1.0, &dJzdz.at(izt*nz*nJ), nJ,
                        dzdx.data(), nz, 1.0, &dJzdx->at(nroots*nx_solver*nJ), nJ);
        } else {
            // rz
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nx_solver, nz, 1.0, &dJrzdz.at(izt*nz*nJ), nJ,
                        drzdx.data(), nz, 1.0, &dJzdx->at(nroots*nx_solver*nJ), nJ);
        }
    }
}

void Model::initialize(AmiVector *x, AmiVector *dx,
                       AmiVectorArray *sx, AmiVectorArray *sdx,
                       bool computeSensitivities) {
    initializeStates(x);
    if(computeSensitivities)
        initializeStateSensitivities(sx, x);

    fdx0(x, dx);
    if(computeSensitivities)
        fsdx0();

    if(ne)
        initHeaviside(x,dx);

}


void Model::initializeStates(AmiVector *x) {
    if (x0data.empty()) {
        fx0(x);
    } else {
        std::vector<realtype> x0_solver(nx_solver, 0.0);
        ftotal_cl(total_cl.data(), x0data.data());
        fx_solver(x0_solver.data(), x0data.data());
        for (int ix = 0; ix < nx_solver; ix++) {
            (*x)[ix] = (realtype) x0_solver.at(ix);
        }
    }
}

void Model::initializeStateSensitivities(AmiVectorArray *sx, AmiVector *x) {
    if (sx0data.empty()) {
        fsx0(sx, x);
    } else {
        realtype *stcl = nullptr;
        std::vector<realtype> sx0_solver_slice(nx_solver, 0.0);
        for (int ip = 0; ip < nplist(); ip++) {
            if (ncl() > 0)
                stcl = &stotal_cl.at(plist(ip) * ncl());
            fstotal_cl(stcl, &sx0data.at(ip*nx_rdata), plist(ip));
            fsx_solver(sx0_solver_slice.data(), &sx0data.at(ip*nx_rdata));
            for (int ix = 0; ix < nx_solver; ix++) {
                sx->at(ix,ip) = (realtype) sx0_solver_slice.at(ix);
            }
        }
    }
}

void Model::initHeaviside(AmiVector *x, AmiVector *dx) {
    std::vector<realtype> rootvals(ne,0.0);
    froot(tstart, x, dx, rootvals.data());
    for (int ie = 0; ie < ne; ie++) {
        if (rootvals.at(ie) < 0) {
            h.at(ie) = 0.0;
        } else if (rootvals.at(ie) == 0) {
            throw AmiException("Simulation started in an event. This could lead to "
                               "unexpected results, aborting simulation! Please "
                               "specify an earlier simulation start via "
                               "options.t0");
        } else {
            h.at(ie) = 1.0;
        }
    }
}

int Model::nplist() const {
    return plist_.size();
}

int Model::np() const {
    return originalParameters.size();
}

int Model::nk() const {
    return fixedParameters.size();
}

int Model::ncl() const {return nx_rdata-nx_solver;}

const double *Model::k() const{
    return fixedParameters.data();
}

int Model::nMaxEvent() const {
    return nmaxevent;
}

void Model::setNMaxEvent(int nmaxevent) {
    this->nmaxevent = nmaxevent;
}

int Model::nt() const {
    return ts.size();
}

const std::vector<ParameterScaling> &Model::getParameterScale() const {
    return pscale;
}

void Model::setParameterScale(ParameterScaling pscale) {
    this->pscale.assign(this->pscale.size(), pscale);
    scaleParameters(unscaledParameters, this->pscale, originalParameters);
    sx0data.clear();
}

void Model::setParameterScale(std::vector<ParameterScaling> const& pscaleVec) {
    if(pscale.size() != this->originalParameters.size())
        throw AmiException("Dimension mismatch. Size of parameter scaling does not match number of model parameters.");
    this->pscale = pscaleVec;
    scaleParameters(unscaledParameters, this->pscale, originalParameters);
    sx0data.clear();
}

std::vector<realtype> const& Model::getParameters() const {
    return originalParameters;
}

/**
 * @brief local helper function to get parameters
 * @param ids vector of name/ids of (fixed)Parameters
 * @param values values of the (fixed)Parameters
 * @param id name/id to look for in the vector
 * @param variable_name string indicating what variable we are lookin at
 * @param id_name string indicating whether name or id was specified
 * @return value of the selected parameter
 */
realtype getValueById(std::vector<std::string> const& ids, std::vector<realtype> const& values,
                      std::string const& id, const char* variable_name, const char* id_name) {
    auto it = std::find(ids.begin(), ids.end(), id);
    if(it != ids.end())
        return values.at(it - ids.begin());

    throw AmiException("Could not find %s with specified %s", variable_name, id_name);
}

/**
 * @brief local helper function to set parameters
 * @param ids vector of names/ids of (fixed)Parameters
 * @param values values of the (fixed)Parameters
 * @param value for the selected parameter
 * @param id name/id to look for in the vector
 * @param variable_name string indicating what variable we are lookin at
 * @param id_name string indicating whether name or id was specified
 */
void setValueById(std::vector<std::string> const& ids, std::vector<realtype> &values, realtype value,
                      std::string const& id, const char* variable_name, const char* id_name) {
    auto it = std::find(ids.begin(), ids.end(), id);
    if(it != ids.end())
        values.at(it - ids.begin()) = value;
    else
        throw AmiException("Could not find %s with specified %s", variable_name, id_name);
}

/**
 * @brief local helper function to set parameters via regex
 * @param ids vector of names/ids of (fixed)Parameters
 * @param values values of the (fixed)Parameters
 * @param value for the selected parameter
 * @param regex string according to which names/ids are to be matched
 * @param variable_name string indicating what variable we are lookin at
 * @param id_name string indicating whether name or id was specified
 * @return number of matched names/ids
 */
int setValueByIdRegex(std::vector<std::string> const& ids, std::vector<realtype> &values, realtype value,
                  std::string const& regex, const char* variable_name, const char* id_name) {
    try {
        std::regex pattern (regex);
        int n_found = 0;
        for(const auto &id: ids) {
            if(std::regex_match(id, pattern)) {
                values.at(&id - &ids[0]) = value;
                ++n_found;
            }
        }

        if(n_found == 0)
            throw AmiException("Could not find %s with specified %s", variable_name, id_name);

        return n_found;
    } catch (std::regex_error const& e) {
        throw AmiException("Specified regex pattern could not be compiled: %s", e.what());
    }
}

realtype Model::getParameterById(std::string const& par_id) const {
    if(!hasParameterIds())
        throw AmiException("Could not access parameters by id as they are not set");
    return getValueById(getParameterIds(), originalParameters, par_id, "parameters", "id");
}

realtype Model::getParameterByName(std::string const& par_name) const {
    if(!hasParameterNames())
        throw AmiException("Could not access parameters by name as they are not set");
    return getValueById(getParameterNames(),
                        originalParameters,
                        par_name,
                        "parameters",
                        "name");
}

void Model::setParameters(const std::vector<realtype> &p) {
    if(p.size() != (unsigned) np())
        throw AmiException("Dimension mismatch. Size of parameters does not match number of model parameters.");
    this->originalParameters = p;
    this->unscaledParameters.resize(originalParameters.size());
    unscaleParameters(originalParameters, pscale, unscaledParameters);
}

void Model::setParameterById(std::string const& par_id, realtype value) {
    if(!hasParameterIds())
        throw AmiException("Could not access parameters by id as they are not set");

    setValueById(getParameterIds(),
                 originalParameters,
                 value,
                 par_id,
                 "parameter",
                 "id");
    unscaleParameters(originalParameters, pscale, unscaledParameters);
}

int Model::setParametersByIdRegex(std::string const& par_id_regex, realtype value) {
    if(!hasParameterIds())
        throw AmiException("Could not access parameters by id as they are not set");
    int n_found = setValueByIdRegex(getParameterIds(),
                                    originalParameters,
                                    value,
                                    par_id_regex,
                                    "parameter",
                                    "id");
    unscaleParameters(originalParameters, pscale, unscaledParameters);
    return n_found;
}

void Model::setParameterByName(std::string const& par_name, realtype value) {
    if(!hasParameterNames())
        throw AmiException("Could not access parameters by name as they are not set");

    setValueById(getParameterNames(),
                 originalParameters,
                 value,
                 par_name,
                 "parameter",
                 "name");
    unscaleParameters(originalParameters, pscale, unscaledParameters);
}

int Model::setParametersByNameRegex(std::string const& par_name_regex, realtype value) {
    if(!hasParameterNames())
        throw AmiException("Could not access parameters by name as they are not set");

    int n_found = setValueByIdRegex(getParameterNames(),
                                    originalParameters,
                                    value,
                                    par_name_regex,
                                    "parameter",
                                    "name");

    unscaleParameters(originalParameters, pscale, unscaledParameters);
    return n_found;
}

bool Model::hasStateIds() const { return nx_rdata && !getStateIds().empty(); }

std::vector<std::string> Model::getStateIds() const {
    return std::vector<std::string>();
}

bool Model::hasFixedParameterIds() const { return nk() && !getFixedParameterIds().empty(); }

std::vector<std::string> Model::getFixedParameterIds() const {
    return std::vector<std::string>();
}

bool Model::hasObservableIds() const { return ny && !getObservableIds().empty(); }

std::vector<std::string> Model::getObservableIds() const {
    return std::vector<std::string>();
}

void Model::setSteadyStateSensitivityMode(const SteadyStateSensitivityMode mode) {
    steadyStateSensitivityMode = mode;
}

SteadyStateSensitivityMode Model::getSteadyStateSensitivityMode() const {
    return steadyStateSensitivityMode;
}

void Model::setReinitializeFixedParameterInitialStates(bool flag) {
    if (flag && !isFixedParameterStateReinitializationAllowed())
        throw AmiException("State reinitialization cannot be enabled for this model"
                           "as this feature was disabled at compile time. Most likely,"
                           " this was because some initial states depending on "
                           "fixedParameters also depended on parameters");
    reinitializeFixedParameterInitialStates = flag;
}

bool Model::getReinitializeFixedParameterInitialStates() const {
    return reinitializeFixedParameterInitialStates;
}

void Model::fx_rdata(realtype *x_rdata, const realtype *x_solver, const realtype *tcl) {
    if (nx_solver != nx_rdata)
        throw AmiException(
                "A model that has differing nx_solver and nx_rdata needs "
                "to implement its own fx_rdata");
    std::copy_n(x_solver, nx_solver, x_rdata);
}

void Model::fsx_rdata(realtype *sx_rdata, const realtype *sx_solver, const realtype *stcl, const int ip) {
    fx_rdata(sx_rdata, sx_solver, stcl);
}

void Model::fx_solver(realtype *x_solver, const realtype *x_rdata) {
    if (nx_solver != nx_rdata)
        throw AmiException(
                "A model that has differing nx_solver and nx_rdata needs "
                "to implement its own fx_solver");
    std::copy_n(x_rdata, nx_rdata, x_solver);
}

void Model::fsx_solver(realtype *sx_solver, const realtype *sx_rdata) {
    /* for the moment we do not need an implementation of fsx_solver as
     * we can simply reuse fx_solver and replace states by their
     * sensitivities */
    fx_solver(sx_solver, sx_rdata);
}

void Model::ftotal_cl(realtype *total_cl, const realtype *x_rdata) {
    if (nx_solver != nx_rdata)
        throw AmiException(
                "A model that has differing nx_solver and nx_rdata needs "
                "to implement its own ftotal_cl");
}

void Model::fstotal_cl(realtype *stotal_cl, const realtype *sx_rdata, const int ip) {
    /* for the moment we do not need an implementation of fstotal_cl as
     * we can simply reuse ftotal_cl and replace states by their
     * sensitivities */
    ftotal_cl(stotal_cl, sx_rdata);
}

const std::vector<realtype> &Model::getUnscaledParameters() const {
    return unscaledParameters;
}

const std::vector<realtype> &Model::getFixedParameters() const {
    return fixedParameters;
}

realtype Model::getFixedParameterById(std::string const& par_id) const {
    if(!hasFixedParameterIds())
        throw AmiException("Could not access fixed parameters by id as they are not set");

    return getValueById(getFixedParameterIds(),
                        fixedParameters,
                        par_id,
                        "fixedParameters",
                        "id");
}

realtype Model::getFixedParameterByName(std::string const& par_name) const {
    if(!hasFixedParameterNames())
        throw AmiException("Could not access fixed parameters by name as they are not set");

    return getValueById(getFixedParameterNames(),
                        fixedParameters,
                        par_name,
                        "fixedParameters",
                        "name");
}

void Model::setFixedParameters(const std::vector<realtype> &k) {
    if(k.size() != (unsigned) nk())
        throw AmiException("Dimension mismatch. Size of fixedParameters does not match number of fixed model parameters.");
    this->fixedParameters = k;
}

void Model::setFixedParameterById(std::string const& par_id, realtype value) {
    if(!hasFixedParameterIds())
        throw AmiException("Could not access fixed parameters by id as they are not set");

    setValueById(getFixedParameterIds(),
                 fixedParameters,
                 value,
                 par_id,
                 "fixedParameters",
                 "id");
}

int Model::setFixedParametersByIdRegex(std::string const& par_id_regex, realtype value) {
    if(!hasFixedParameterIds())
        throw AmiException("Could not access fixed parameters by id as they are not set");

    return  setValueByIdRegex(getFixedParameterIds(),
                              fixedParameters,
                              value,
                              par_id_regex,
                              "fixedParameters",
                              "id");
}

void Model::setFixedParameterByName(std::string const& par_name, realtype value) {
    if(!hasFixedParameterNames())
        throw AmiException("Could not access fixed parameters by name as they are not set");

    setValueById(getFixedParameterNames(),
                 fixedParameters,
                 value,
                 par_name,
                 "fixedParameters",
                 "name");
}

int Model::setFixedParametersByNameRegex(std::string const& par_name_regex, realtype value) {
    if(!hasFixedParameterNames())
        throw AmiException("Could not access fixed parameters by name as they are not set");

    return  setValueByIdRegex(getFixedParameterIds(),
                              fixedParameters,
                              value,
                              par_name_regex,
                              "fixedParameters",
                              "name");
}

std::vector<realtype> const& Model::getTimepoints() const {
    return ts;
}

void Model::setTimepoints(const std::vector<realtype> &ts) {
    if (!std::is_sorted(ts.begin(), ts.end()))
        throw AmiException("Encountered non-monotonic timepoints, please order"
                           " timepoints such that they are monotonically"
                           " increasing!");
    this->ts = std::move(ts);
}

std::vector<bool> const& Model::getStateIsNonNegative() const {
    return stateIsNonNegative;
}

void Model::setStateIsNonNegative(std::vector<bool> const& nonNegative) {
    if(nx_solver != nx_rdata) {
        throw AmiException("Nonnegative states are not supported whith"
                           " conservation laws enabled");
    }
    if (stateIsNonNegative.size() != static_cast<unsigned long>(nx_rdata)) {
        throw AmiException("Dimension of input stateIsNonNegative (%u) does "
                           "not agree with number of state variables (%d)",
                           stateIsNonNegative.size(), nx_rdata);
    }
    stateIsNonNegative=nonNegative;
    anyStateNonNegative=std::any_of(stateIsNonNegative.begin(),
                                    stateIsNonNegative.end(),
                                    [](bool x) { return x; });
}

void Model::setAllStatesNonNegative()
{
    setStateIsNonNegative(std::vector<bool>(nx_solver, true));
}

double Model::t(int idx) const {
    return ts.at(idx);
}

const std::vector<int> &Model::getParameterList() const {
    return plist_;
}

void Model::setParameterList(const std::vector<int> &plist) {
    int np = this->np(); // cannot capture 'this' in lambda expression
    if(std::any_of(plist.begin(), plist.end(),
                   [&np](int idx){return idx < 0 || idx >= np;})) {
        throw AmiException("Indices in plist must be in [0..np]");
    }
    this->plist_ = plist;

    initializeVectors();
}

std::vector<realtype> const& Model::getInitialStates() const {
    return x0data;
}

void Model::setInitialStates(const std::vector<realtype> &x0) {
    if (x0.size() != (unsigned)nx_rdata && x0.size() != 0)
        throw AmiException("Dimension mismatch. Size of x0 does not match "
                           "number of model states.");

    if (x0.size() == 0) {
        x0data.clear();
        return;
    }

    x0data = x0;
}

const std::vector<realtype> &Model::getInitialStateSensitivities() const {
    return sx0data;
}

void Model::setInitialStateSensitivities(const std::vector<realtype> &sx0) {
    if (sx0.size() != (unsigned)nx_rdata * nplist() && sx0.size() != 0)
        throw AmiException("Dimension mismatch. Size of sx0 does not match "
                           "number of model states * number of parameter "
                           "selected for sensitivities.");

    if (sx0.size() == 0) {
        sx0data.clear();
        return;
    }

    realtype chainrulefactor = 1.0;
    std::vector<realtype> sx0_rdata(nx_rdata * nplist(), 0.0);
    for (int ip = 0; ip < nplist(); ip++) {

        // revert chainrule
        switch (pscale.at(plist(ip))) {
            case ParameterScaling::log10:
                chainrulefactor = unscaledParameters.at(plist(ip))
                * log(10);
                break;
            case ParameterScaling::ln:
                chainrulefactor = unscaledParameters.at(plist(ip));
                break;
            case ParameterScaling::none:
                chainrulefactor = 1.0;
                break;
        }

        for(int ix=0; ix < nx_rdata; ++ix) {
            sx0_rdata.at(ip*nx_rdata + ix) =
            sx0.at(ip*nx_rdata + ix) / chainrulefactor;
        }
    }
    setUnscaledInitialStateSensitivities(sx0_rdata);
}

void Model::setUnscaledInitialStateSensitivities(
    const std::vector<realtype> &sx0) {
    if (sx0.size() != (unsigned)nx_rdata * nplist() && sx0.size() != 0)
        throw AmiException("Dimension mismatch. Size of sx0 does not match "
                           "number of model states * number of parameter "
                           "selected for sensitivities.");

    if (sx0.size() == 0) {
        sx0data.clear();
        return;
    }

    sx0data = sx0;
}

double Model::t0() const{
    return tstart;
}

void Model::setT0(double t0) {
    tstart = t0;
}

int Model::plist(int pos) const{
    return plist_.at(pos);
}

bool Model::hasParameterNames() const { return np() && !getParameterNames().empty(); }

std::vector<std::string> Model::getParameterNames() const { return std::vector<std::string>(); }

bool Model::hasStateNames() const { return nx_rdata && !getStateNames().empty(); }

std::vector<std::string> Model::getStateNames() const { return std::vector<std::string>(); }

bool Model::hasFixedParameterNames() const { return nk() && !getFixedParameterNames().empty(); }

std::vector<std::string> Model::getFixedParameterNames() const { return std::vector<std::string>(); }

bool Model::hasObservableNames() const { return ny && !getObservableNames().empty(); }

std::vector<std::string> Model::getObservableNames() const { return std::vector<std::string>(); }

bool Model::hasParameterIds() const { return np() && !getParameterIds().empty(); }

std::vector<std::string> Model::getParameterIds() const {
    return std::vector<std::string>();
}


Model::Model()
    : nx_rdata(0), nxtrue_rdata(0), nx_solver(0), nxtrue_solver(0), ny(0), nytrue(0), nz(0), nztrue(0),
      ne(0), nw(0), ndwdx(0), ndwdp(0), nnz(0), nJ(0), ubw(0), lbw(0),
      o2mode(SecondOrderMode::none), x_pos_tmp(0) {}

Model::Model(const int nx_rdata,
             const int nxtrue_rdata,
             const int nx_solver,
             const int nxtrue_solver,
             const int ny,
             const int nytrue,
             const int nz,
             const int nztrue,
             const int ne,
             const int nJ,
             const int nw,
             const int ndwdx,
             const int ndwdp,
             const int nnz,
             const int ubw,
             const int lbw,
             SecondOrderMode o2mode,
             const std::vector<realtype>& p,
             std::vector<realtype> k,
             const std::vector<int>& plist,
             std::vector<realtype> idlist,
             std::vector<int> z2event)
    : nx_rdata(nx_rdata), nxtrue_rdata(nxtrue_rdata),
      nx_solver(nx_solver), nxtrue_solver(nxtrue_solver),
      ny(ny), nytrue(nytrue),
      nz(nz), nztrue(nztrue),
      ne(ne),
      nw(nw),
      ndwdx(ndwdx),
      ndwdp(ndwdp),
      nnz(nnz),
      nJ(nJ),
      ubw(ubw),
      lbw(lbw),
      o2mode(o2mode),
      z2event(std::move(z2event)),
      idlist(std::move(idlist)),
      sigmay(ny, 0.0),
      dsigmaydp(ny*plist.size(), 0.0),
      sigmaz(nz, 0.0),
      dsigmazdp(nz*plist.size(), 0.0),
      dJydp(nJ*plist.size(), 0.0),
      dJzdp(nJ*plist.size(), 0.0),
      deltax(nx_solver, 0.0),
      deltasx(nx_solver*plist.size(), 0.0),
      deltaxB(nx_solver, 0.0),
      deltaqB(nJ*plist.size(), 0.0),
      dxdotdp(nx_solver*plist.size(), 0.0),
      J(nx_solver, nx_solver, nnz, CSC_MAT),
      my(nytrue, 0.0),
      mz(nztrue, 0.0),
      dJydy(nJ*nytrue*ny, 0.0),
      dJydsigma(nJ*nytrue*ny, 0.0),
      dJzdz(nJ*nztrue*nz, 0.0),
      dJzdsigma(nJ*nztrue*nz, 0.0),
      dJrzdz(nJ*nztrue*nz, 0.0),
      dJrzdsigma(nJ*nztrue*nz, 0.0),
      dzdx(nz*nx_solver, 0.0),
      dzdp(nz*plist.size(), 0.0),
      drzdx(nz*nx_solver, 0.0),
      drzdp(nz*plist.size(), 0.0),
      dydp(ny*plist.size(), 0.0),
      dydx(ny*nx_solver,0.0),
      w(nw, 0.0),
      dwdx(ndwdx, 0.0),
      dwdp(ndwdp, 0.0),
      M(nx_solver*nx_solver, 0.0),
      stau(plist.size(), 0.0),
      sx(nx_solver*plist.size(), 0.0),
      x_rdata(nx_rdata, 0.0),
      sx_rdata(nx_rdata, 0.0),
      h(ne,0.0),
      unscaledParameters(p),
      originalParameters(p),
      fixedParameters(std::move(k)),
      total_cl(nx_rdata-nx_solver),
      stotal_cl((nx_rdata-nx_solver) * np()),
      plist_(plist),
      stateIsNonNegative(nx_solver, false),
      x_pos_tmp(nx_solver),
      pscale(std::vector<ParameterScaling>(p.size(), ParameterScaling::none))
{
    requireSensitivitiesForAllParameters();
}

void Model::initializeVectors()
{
    dsigmaydp.resize(ny * nplist(), 0.0);
    dsigmazdp.resize(nz * nplist(), 0.0);
    dJydp.resize(nJ * nplist(), 0.0);
    dJzdp.resize(nJ * nplist(), 0.0);
    deltasx.resize(nx_solver * nplist(), 0.0);
    deltaqB.resize(nJ * nplist(), 0.0);
    dxdotdp.resize(nx_solver * nplist(), 0.0);
    dzdp.resize(nz * nplist(), 0.0);
    drzdp.resize(nz * nplist(), 0.0);
    dydp.resize(ny * nplist(), 0.0);
    stau.resize(nplist(), 0.0);
    sx.resize(nx_solver * nplist(), 0.0);
    sx0data.clear();
}

void Model::fx_rdata(AmiVector *x_rdata, const AmiVector *x) {
    fx_rdata(x_rdata->data(), x->data(), total_cl.data());
    if(alwaysCheckFinite)
        checkFinite(x_rdata->getLength(), x_rdata->data(), "x_rdata");
}

void Model::fx0(AmiVector *x) {
    std::fill(x_rdata.begin(), x_rdata.end(), 0.0);
    /* this function  also computes initial total abundances */
    fx0(x_rdata.data(), tstart, unscaledParameters.data(),
        fixedParameters.data());
    fx_solver(x->data(), x_rdata.data());
    ftotal_cl(total_cl.data(), x_rdata.data());

    if(alwaysCheckFinite) {
        checkFinite(x_rdata.size(), x_rdata.data(), "x0 x_rdata");
        checkFinite(x->getLength(), x->data(), "x0 x");
    }
}

void Model::fx0_fixedParameters(AmiVector *x) {
    if (!getReinitializeFixedParameterInitialStates())
        return;
    /* we transform to the unreduced states x_rdata and then apply
     x0_fixedparameters to (i) enable updates to states that were removed from
     conservation laws and (ii) be able to correctly compute total abundances
     after updating the state variables */
    fx_rdata(x_rdata.data(), x->data(), total_cl.data());
    fx0_fixedParameters(x_rdata.data(), tstart, unscaledParameters.data(),
                        fixedParameters.data());
    fx_solver(x->data(), x_rdata.data());
    /* update total abundances */
    ftotal_cl(total_cl.data(), x_rdata.data());
}



void Model::fsx_rdata(AmiVectorArray *sx_full, const AmiVectorArray *sx) {
    realtype *stcl = nullptr;
    for (int ip = 0; ip < nplist(); ip++) {
        if (ncl() > 0)
            stcl = &stotal_cl.at(plist(ip) * ncl());
        fsx_rdata(sx_full->data(ip), sx->data(ip), stcl, ip);
    }
}

void Model::fsx0(AmiVectorArray *sx, const AmiVector *x) {
    /* this function  also computes initial total abundance sensitivities */
    realtype *stcl = nullptr;
    for (int ip = 0; ip < nplist(); ip++) {
        if (ncl() > 0)
            stcl = &stotal_cl.at(plist(ip) * ncl());
        std::fill(sx_rdata.begin(), sx_rdata.end(), 0.0);
        fsx0(sx_rdata.data(), tstart, x->data(), unscaledParameters.data(),
             fixedParameters.data(), plist(ip));
        fsx_solver(sx->data(ip), sx_rdata.data());
        fstotal_cl(stcl, sx_rdata.data(), plist(ip));
    }
}

void Model::fsx0_fixedParameters(AmiVectorArray *sx, const AmiVector *x) {
    if (!getReinitializeFixedParameterInitialStates())
        return;
    realtype *stcl = nullptr;
    for (int ip = 0; ip < nplist(); ip++) {
        if (ncl() > 0)
            stcl = &stotal_cl.at(plist(ip) * ncl());
        fsx_rdata(sx_rdata.data(), sx->data(ip), stcl, plist(ip));
        fsx0_fixedParameters(sx_rdata.data(), tstart, x->data(),
                             unscaledParameters.data(), fixedParameters.data(),
                             plist(ip));
        fsx_solver(sx->data(ip), sx_rdata.data());
        fstotal_cl(stcl, sx_rdata.data(), plist(ip));
    }
}

void Model::fsdx0() {}

void Model::fstau(const realtype t, const int ie, const AmiVector *x, const AmiVectorArray *sx) {
    std::fill(stau.begin(),stau.end(),0.0);
    for(int ip = 0; ip < nplist(); ip++){
        fstau(&stau.at(ip),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),sx->data(ip),plist(ip),ie);
    }
}

void Model::fy(const realtype t, const int it, const AmiVector *x, ReturnData *rdata) {
    if (!ny)
        return;
    fw(t,x->data());
    fy(&rdata->y.at(it*ny),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),w.data());

    if(alwaysCheckFinite) {
        amici::checkFinite(ny, &rdata->y.at(it*ny), "y");
    }
}

void Model::fdydp(const realtype t, const AmiVector *x) {
    if (!ny)
        return;

    std::fill(dydp.begin(),dydp.end(),0.0);
    fw(t,x->data());
    fdwdp(t,x->data());
    for(int ip = 0; ip < nplist(); ip++){
        // get dydp slice (ny) for current time and parameter
        fdydp(&dydp.at(ip*ny),
              t,
              x->data(),
              unscaledParameters.data(),
              fixedParameters.data(),
              h.data(),
              plist(ip),
              w.data(),
              dwdp.data());
    }

    if(alwaysCheckFinite) {
        amici::checkFinite(dydp, "dydp");
    }
}

void Model::fdydx(const realtype t, const AmiVector *x) {
    if (!ny)
        return;

    std::fill(dydx.begin(),dydx.end(),0.0);
    fw(t,x->data());
    fdwdx(t,x->data());
    fdydx(dydx.data(),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),w.data(),dwdx.data());

    if(alwaysCheckFinite) {
        amici::checkFinite(dydx, "dydx");
    }
}

void Model::fz(const int nroots, const int ie, const realtype t, const AmiVector *x, ReturnData *rdata) {
    fz(&rdata->z.at(nroots*nz),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data());
}

void Model::fsz(const int nroots, const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata) {
    for(int ip = 0; ip < nplist();  ip++ ){
        fsz(&rdata->sz.at((nroots*nplist()+ip)*nz),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),sx->data(ip),plist(ip));
    }
}

void Model::frz(const int nroots, const int ie, const realtype t, const AmiVector *x, ReturnData *rdata) {
    frz(&rdata->rz.at(nroots*nz),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data());
}

void Model::fsrz(const int nroots, const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata) {
    for(int ip = 0; ip < nplist();  ip++ ){
        fsrz(&rdata->srz.at((nroots*nplist()+ip)*nz),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),sx->data(ip),plist(ip));
    }
}

void Model::fdzdp(const realtype t, const int ie, const AmiVector *x) {
    std::fill(dzdp.begin(),dzdp.end(),0.0);
    for(int ip = 0; ip < nplist(); ip++){
        fdzdp(dzdp.data(),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),plist(ip));
    }

    if(alwaysCheckFinite) {
        amici::checkFinite(dzdp, "dzdp");
    }
}

void Model::fdzdx(const realtype t, const int ie, const AmiVector *x) {
    std::fill(dzdx.begin(),dzdx.end(),0.0);
    fdzdx(dzdx.data(),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data());

    if(alwaysCheckFinite) {
        amici::checkFinite(dzdx, "dzdx");
    }
}

void Model::fdrzdp(const realtype t, const int ie, const AmiVector *x) {
    std::fill(drzdp.begin(),drzdp.end(),0.0);
    for(int ip = 0; ip < nplist(); ip++){
        fdrzdp(drzdp.data(),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),plist(ip));
    }

    if(alwaysCheckFinite) {
        amici::checkFinite(drzdp, "drzdp");
    }
}


void Model::fdrzdx(const realtype t, const int ie, const AmiVector *x) {
    std::fill(drzdx.begin(),drzdx.end(),0.0);
    fdrzdx(drzdx.data(),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data());

    if(alwaysCheckFinite) {
        amici::checkFinite(drzdx, "drzdx");
    }
}

void Model::fdeltax(const int ie, const realtype t, const AmiVector *x,
                    const AmiVector *xdot, const AmiVector *xdot_old) {
    std::fill(deltax.begin(),deltax.end(),0.0);
    fdeltax(deltax.data(),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),ie,xdot->data(),xdot_old->data());

    if(alwaysCheckFinite) {
        amici::checkFinite(deltax, "deltax");
    }
}

void Model::fdeltasx(const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx,
                     const AmiVector *xdot, const AmiVector *xdot_old) {
    fw(t,x->data());
    std::fill(deltasx.begin(),deltasx.end(),0.0);
    for(int ip = 0; ip < nplist(); ip++)
        fdeltasx(&deltasx.at(nx_solver*ip),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),w.data(),
                 plist(ip),ie,xdot->data(),xdot_old->data(),sx->data(ip),&stau.at(ip));

    if(alwaysCheckFinite) {
        amici::checkFinite(deltasx, "deltasx");
    }
}

void Model::fdeltaxB(const int ie, const realtype t, const AmiVector *x, const AmiVector *xB,
                     const AmiVector *xdot, const AmiVector *xdot_old) {
    std::fill(deltaxB.begin(),deltaxB.end(),0.0);
    fdeltaxB(deltaxB.data(),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),ie,xdot->data(),xdot_old->data(),xB->data());

    if(alwaysCheckFinite) {
        amici::checkFinite(deltaxB, "deltaxB");
    }
}

void Model::fdeltaqB(const int ie, const realtype t, const AmiVector *x, const AmiVector *xB,
                     const AmiVector *xdot, const AmiVector *xdot_old) {
    std::fill(deltaqB.begin(),deltaqB.end(),0.0);
    for(int ip = 0; ip < nplist(); ip++)
        fdeltaqB(deltaqB.data(),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                 plist(ip),ie,xdot->data(),xdot_old->data(),xB->data());

    if(alwaysCheckFinite) {
        amici::checkFinite(deltaqB, "deltaqB");
    }
}

void Model::fsigmay(const int it, ReturnData *rdata, const ExpData *edata) {
    if (!ny)
        return;

    std::fill(sigmay.begin(),sigmay.end(),0.0);

    fsigmay(sigmay.data(),rdata->ts.at(it), unscaledParameters.data(),fixedParameters.data());

    if(edata){
        auto sigmay_edata = edata->getObservedDataStdDevPtr(it);
        /* extract the value for the standard deviation from ExpData,
         * if the data value is NaN, use the parameter value */
        for (int iytrue = 0; iytrue < nytrue; iytrue++) {
            if (edata->isSetObservedDataStdDev(it, iytrue)) {
                sigmay.at(iytrue) = sigmay_edata[iytrue];
            }
        }
    }

    for(int i = 0; i < nytrue; ++i)
        checkSigmaPositivity(sigmay[i], "sigmay");

    std::copy_n(sigmay.data(), nytrue, &rdata->sigmay[it * rdata->ny]);
}

void Model::fdsigmaydp(const int it, ReturnData *rdata, const ExpData *edata) {
    if (!ny)
        return;

    std::fill(dsigmaydp.begin(), dsigmaydp.end(), 0.0);

    for(int ip = 0; ip < nplist(); ip++)
        // get dsigmaydp slice (ny) for current timepoint and parameter
        fdsigmaydp(&dsigmaydp.at(ip*ny),
                    rdata->ts.at(it),
                    unscaledParameters.data(),
                    fixedParameters.data(),
                    plist(ip));

    // sigmas in edata override model-sigma -> for those sigmas, set dsigmaydp to zero
    if(edata){
        for (int iy = 0; iy < nytrue; iy++) {
            if (!edata->isSetObservedDataStdDev(it, iy))
                continue;
            for (int ip = 0; ip < nplist(); ip++) {
                dsigmaydp[ip * ny + iy] = 0.0;
            }
        }
    }

    // copy dsigmaydp slice for current timepoint
    std::copy(dsigmaydp.begin(), dsigmaydp.end(), &rdata->ssigmay[it * nplist() * ny]);

    if(alwaysCheckFinite) {
        amici::checkFinite(dsigmaydp, "dsigmaydp");
    }
}

void Model::fsigmaz(const realtype t, const int ie, const int *nroots,
                    ReturnData *rdata, const ExpData *edata) {
    std::fill(sigmaz.begin(), sigmaz.end(), 0.0);
    fsigmaz(sigmaz.data(), t, unscaledParameters.data(),
            fixedParameters.data());
    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (z2event.at(iztrue) - 1 == ie) {
            if (edata) {
                if (edata->isSetObservedEventsStdDev(nroots[ie], iztrue)) {
                    auto sigmaz_edata =
                        edata->getObservedEventsStdDevPtr(nroots[ie]);
                    sigmaz.at(iztrue) = sigmaz_edata[iztrue];
                }
            }
            rdata->sigmaz[nroots[ie] * rdata->nz + iztrue] = sigmaz.at(iztrue);
        }
    }

    if(alwaysCheckFinite) {
        amici::checkFinite(sigmaz, "sigmaz");
    }
}

void Model::fdsigmazdp(const realtype t, const int ie, const int *nroots, ReturnData *rdata, const ExpData *edata) {
    std::fill(dsigmazdp.begin(),dsigmazdp.end(),0.0);
    for(int ip = 0; ip < nplist(); ip++) {
        // get dsigmazdp slice (nz) for current event and parameter
        fdsigmazdp(&dsigmazdp.at(ip*nz),
                   t,
                   unscaledParameters.data(),
                   fixedParameters.data(),
                   plist(ip));
    }

    // sigmas in edata override model-sigma -> for those sigmas, set dsigmazdp to zero
    if(edata) {
        for (int iz = 0; iz < nztrue; iz++) {
            if (z2event.at(iz) - 1 == ie && !edata->isSetObservedEventsStdDev(nroots[ie],iz)) {
                for(int ip = 0; ip < nplist(); ip++)
                    dsigmazdp.at(iz + nz * ip) = 0;
            }
        }
    }

    // copy dsigmazdp slice for current event
    std::copy(dsigmazdp.begin(), dsigmazdp.end(), &rdata->ssigmaz[nroots[ie] * nplist() * nz]);

    if(alwaysCheckFinite) {
        amici::checkFinite(dsigmazdp, "dsigmazdp");
    }
}

void Model::fJy(const int it, ReturnData *rdata, const ExpData *edata) {
    std::vector<realtype> nllh(nJ,0.0);
    getmy(it,edata);
    for(int iytrue = 0; iytrue < nytrue; iytrue++){
        if(!isNaN(my.at(iytrue))){
            std::fill(nllh.begin(),nllh.end(),0.0);
            fJy(nllh.data(),iytrue, unscaledParameters.data(),fixedParameters.data(),gety(it,rdata),sigmay.data(),my.data());
            rdata->llh -= nllh.at(0);
        }
    }
}

void Model::fJz(const int nroots, ReturnData *rdata, const ExpData *edata) {
    std::vector<realtype> nllh(nJ,0.0);
    getmz(nroots,edata);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            std::fill(nllh.begin(),nllh.end(),0.0);
            fJz(nllh.data(),iztrue, unscaledParameters.data(),fixedParameters.data(),getz(nroots,rdata),sigmaz.data(),mz.data());
            rdata->llh -= nllh.at(0);
        }
    }
}

void Model::fJrz(const int nroots, ReturnData *rdata, const ExpData *edata) {
    std::vector<realtype> nllh(nJ,0.0);
    getrz(nroots,rdata);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            std::fill(nllh.begin(),nllh.end(),0.0);
            fJrz(nllh.data(),iztrue, unscaledParameters.data(),fixedParameters.data(),getrz(nroots,rdata),sigmaz.data());
            rdata->llh -= nllh.at(0);
        }
    }
}

void Model::fdJydy(const int it, const ReturnData *rdata,
                   const ExpData *edata) {
    // load measurements to my
    getmy(it,edata);
    std::fill(dJydy.begin(),dJydy.end(),0.0);

    for(int iytrue = 0; iytrue < nytrue; iytrue++){
        if(!isNaN(my.at(iytrue))){
            // get dJydy slice (ny) for current timepoint and observable
            fdJydy(&dJydy.at(iytrue*ny*nJ),
                   iytrue,
                   unscaledParameters.data(),
                   fixedParameters.data(),
                   gety(it,rdata),
                   sigmay.data(),
                   my.data());
        }
    }

    if(alwaysCheckFinite) {
        amici::checkFinite(dJydy, "dJydy");
    }
}

void Model::fdJydsigma(const int it, const ReturnData *rdata,
                       const ExpData *edata) {
    // load measurements to my
    getmy(it, edata);
    std::fill(dJydsigma.begin(),dJydsigma.end(),0.0);

    for(int iytrue = 0; iytrue < nytrue; iytrue++){
        if(!isNaN(my.at(iytrue))){
            // get dJydsigma slice (ny) for current timepoint and observable
            fdJydsigma(&dJydsigma.at(iytrue*ny*nJ),
                       iytrue,
                       unscaledParameters.data(),
                       fixedParameters.data(),
                       gety(it,rdata),
                       sigmay.data(),
                       my.data());
        }
    }

    if(alwaysCheckFinite) {
        amici::checkFinite(dJydsigma, "dJydsigma");
    }
}

void Model::fdJzdz(const int nroots, const ReturnData *rdata,
                   const ExpData *edata) {
    getmz(nroots,edata);
    std::fill(dJzdz.begin(),dJzdz.end(),0.0);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            fdJzdz(&dJzdz.at(iztrue*nz*nJ),iztrue, unscaledParameters.data(),fixedParameters.data(),getz(nroots,rdata),sigmaz.data(),mz.data());
        }
    }

    if(alwaysCheckFinite) {
        amici::checkFinite(dJzdz, "dJzdz");
    }
}

void Model::fdJzdsigma(const int nroots, const ReturnData *rdata,
                       const ExpData *edata) {
    getmz(nroots,edata);
    std::fill(dJzdsigma.begin(),dJzdsigma.end(),0.0);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            fdJzdsigma(&dJzdsigma.at(iztrue*nz*nJ),iztrue, unscaledParameters.data(),fixedParameters.data(),getz(nroots,rdata),sigmaz.data(),mz.data());
        }
    }

    if(alwaysCheckFinite) {
        amici::checkFinite(dJzdsigma, "dJzdsigma");
    }
}

void Model::fdJrzdz(const int nroots, const ReturnData *rdata,
                    const ExpData *edata) {
    getmz(nroots,edata);
    std::fill(dJrzdz.begin(),dJrzdz.end(),0.0);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            fdJrzdz(&dJrzdz.at(iztrue*nz*nJ),iztrue, unscaledParameters.data(),fixedParameters.data(),getrz(nroots,rdata),sigmaz.data());
        }
    }

    if(alwaysCheckFinite) {
        amici::checkFinite(dJrzdz, "dJrzdz");
    }
}

void Model::fdJrzdsigma(const int nroots,const ReturnData *rdata,
                        const ExpData *edata) {
    std::fill(dJrzdsigma.begin(),dJrzdsigma.end(),0.0);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            fdJrzdsigma(&dJrzdsigma.at(iztrue*nz*nJ),iztrue, unscaledParameters.data(),fixedParameters.data(),getrz(nroots,rdata),sigmaz.data());
        }
    }

    if(alwaysCheckFinite) {
        amici::checkFinite(dJrzdsigma, "dJrzdsigma");
    }
}

void Model::fw(const realtype t, const realtype *x) {
    std::fill(w.begin(),w.end(),0.0);
    fw(w.data(), t, x, unscaledParameters.data(), fixedParameters.data(), h.data(), total_cl.data());

    if(alwaysCheckFinite) {
        amici::checkFinite(w, "w");
    }
}

void Model::fdwdp(const realtype t, const realtype *x) {
    fw(t,x);
    std::fill(dwdp.begin(),dwdp.end(),0.0);
    fdwdp(dwdp.data(), t, x, unscaledParameters.data(), fixedParameters.data(), h.data(), w.data(), total_cl.data(), stotal_cl.data());

    if(alwaysCheckFinite) {
        amici::checkFinite(dwdp, "dwdp");
    }
}

void Model::fdwdx(const realtype t, const realtype *x) {
    fw(t,x);
    std::fill(dwdx.begin(),dwdx.end(),0.0);
    fdwdx(dwdx.data(), t, x, unscaledParameters.data(), fixedParameters.data(), h.data(), w.data(), total_cl.data());

    if(alwaysCheckFinite) {
        amici::checkFinite(dwdx, "dwdx");
    }
}

void Model::fres(const int it, ReturnData *rdata, const ExpData *edata) {
    if (!edata || rdata->res.empty())
        return;

    auto observedData = edata->getObservedDataPtr(it);
    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * edata->nytrue();
        int iyt = iy + it * rdata->ny;
        if (!edata->isSetObservedData(it, iy))
            continue;
        rdata->res.at(iyt_true) = (rdata->y.at(iyt) - observedData[iy])/rdata->sigmay.at(iyt);
    }
}

void Model::fchi2(const int it, ReturnData *rdata) {
    if (rdata->res.empty())
        return;

    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * rdata->nytrue;
        rdata->chi2 += pow(rdata->res.at(iyt_true), 2);
    }
}

void Model::fsres(const int it, ReturnData *rdata, const ExpData *edata) {
    if (!edata || rdata->sres.empty())
        return;

    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * edata->nytrue();
        int iyt = iy + it * rdata->ny;
        if (!edata->isSetObservedData(it,iy))
            continue;
        for (int ip = 0; ip < nplist(); ++ip) {
            rdata->sres.at(iyt_true * nplist() + ip) =
            rdata->sy.at(iy + rdata->ny*(ip + it*nplist()))/rdata->sigmay.at(iyt);
        }
    }
}

void Model::fFIM(const int it, ReturnData *rdata) {
    if (rdata->sres.empty())
        return;

    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * rdata->nytrue;
        for (int ip = 0; ip < nplist(); ++ip) {
            for (int jp = 0; jp < nplist(); ++jp) {
                rdata->FIM.at(ip + nplist() * jp) +=
                rdata->sres.at(iyt_true * nplist() + ip)
                * rdata->sres.at(iyt_true * nplist() + jp);
            }
        }
    }
}

void Model::updateHeaviside(const std::vector<int>& rootsfound) {
    for (int ie = 0; ie < ne; ie++) {
        h.at(ie) += rootsfound.at(ie);
    }
}

void Model::updateHeavisideB(const int *rootsfound) {
    for (int ie = 0; ie < ne; ie++) {
        h.at(ie) -= rootsfound[ie];
    }
}


void Model::getmy(const int it, const ExpData *edata){
    if(edata) {
        std::copy_n(edata->getObservedDataPtr(it), nytrue, my.begin());
    } else {
        std::fill(my.begin(), my.end(), getNaN());
    }
}

const realtype *Model::gety(const int it, const ReturnData *rdata) const {
    return &rdata->y.at(it*ny);
}

realtype Model::gett(const int it, const ReturnData *rdata) const {
    return rdata->ts.at(it);
}

void Model::getmz(const int nroots, const ExpData *edata) {
    if(edata){
        std::copy_n(edata->getObservedEventsPtr(nroots), nztrue, mz.begin());
    } else {
        std::fill(mz.begin(), mz.end(), getNaN());
    }
}

const realtype *Model::getz(const int nroots, const ReturnData *rdata) const {
    return(&rdata->z.at(nroots*nz));
}

const realtype *Model::getrz(const int nroots, const ReturnData *rdata) const {
    return(&rdata->rz.at(nroots*nz));
}

const realtype *Model::getsz(const int nroots, const int ip, const ReturnData *rdata) const {
    return(&rdata->sz.at((nroots*nplist()+ip)*nz));
}

const realtype *Model::getsrz(const int nroots, const int ip, const ReturnData *rdata) const {
    return(&rdata->srz.at((nroots*nplist()+ip)*nz));
}

void Model::setAlwaysCheckFinite(bool alwaysCheck)
{
    alwaysCheckFinite = alwaysCheck;
}

bool Model::getAlwaysCheckFinite() const
{
    return alwaysCheckFinite;
}

int Model::checkFinite(const int N, const realtype *array,
                       const char *fun) const {
    auto result = amici::checkFinite(N, array, fun);

    if (result != AMICI_SUCCESS) {
        amici::checkFinite(ts.size(), ts.data(), "ts");
        amici::checkFinite(fixedParameters.size(), fixedParameters.data(), "k");
        amici::checkFinite(unscaledParameters.size(), unscaledParameters.data(),
                           "p");
        amici::checkFinite(w.size(), w.data(), "w");
    }

    return result;
}

void Model::requireSensitivitiesForAllParameters() {
    plist_.resize(np());
    std::iota(plist_.begin(), plist_.end(), 0);
    initializeVectors();
}

bool operator ==(const Model &a, const Model &b)
{
    if (typeid(a) != typeid(b))
        return false;

    return (a.nx_rdata == b.nx_rdata)
            && (a.nxtrue_rdata == b.nxtrue_rdata)
            && (a.nx_solver == b.nx_solver)
            && (a.nxtrue_solver == b.nxtrue_solver)
            && (a.ny == b.ny)
            && (a.nytrue == b.nytrue)
            && (a.nz == b.nz)
            && (a.nztrue == b.nztrue)
            && (a.ne == b.ne)
            && (a.nw == b.nw)
            && (a.ndwdx == b.ndwdx)
            && (a.ndwdp == b.ndwdp)
            && (a.nnz == b.nnz)
            && (a.nJ == b.nJ)
            && (a.ubw == b.ubw)
            && (a.lbw == b.lbw)
            && (a.o2mode == b.o2mode)
            && (a.z2event == b.z2event)
            && (a.idlist == b.idlist)
            && (a.h == b.h)
            && (a.unscaledParameters == b.unscaledParameters)
            && (a.originalParameters == b.originalParameters)
            && (a.fixedParameters == b.fixedParameters)
            && (a.plist_ == b.plist_)
            && (a.x0data == b.x0data)
            && (a.sx0data == b.sx0data)
            && (a.ts == b.ts)
            && (a.nmaxevent == b.nmaxevent)
            && (a.pscale == b.pscale)
            && (a.stateIsNonNegative == b.stateIsNonNegative)
            && (a.tstart == b.tstart);
}

N_Vector Model::computeX_pos(N_Vector x) {
    if (anyStateNonNegative){
        for (int ix = 0; ix < x_pos_tmp.getLength(); ++ix) {
            x_pos_tmp.at(ix) = (stateIsNonNegative.at(ix) && NV_Ith_S(x, ix) < 0) ? 0 : NV_Ith_S(x, ix);
        }
        return x_pos_tmp.getNVector();
    } else {
        return x;
    }
}

} // namespace amici
