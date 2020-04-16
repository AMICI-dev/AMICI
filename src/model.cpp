#include "amici/model.h"
#include "amici/amici.h"
#include "amici/exception.h"
#include "amici/misc.h"
#include "amici/symbolic_functions.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <regex>
#include <typeinfo>
#include <utility>

namespace amici {

/**
 * @brief local helper function to get parameters
 * @param ids vector of name/ids of (fixed)Parameters
 * @param values values of the (fixed)Parameters
 * @param id name/id to look for in the vector
 * @param variable_name string indicating what variable we are lookin at
 * @param id_name string indicating whether name or id was specified
 * @return value of the selected parameter
 */
static realtype getValueById(std::vector<std::string> const &ids,
                             std::vector<realtype> const &values,
                             std::string const &id, const char *variable_name,
                             const char *id_name) {
    auto it = std::find(ids.begin(), ids.end(), id);
    if (it != ids.end())
        return values.at(it - ids.begin());

    throw AmiException("Could not find %s with specified %s", variable_name,
                       id_name);
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
static void setValueById(std::vector<std::string> const &ids,
                         std::vector<realtype> &values, realtype value,
                         std::string const &id, const char *variable_name,
                         const char *id_name) {
    auto it = std::find(ids.begin(), ids.end(), id);
    if (it != ids.end())
        values.at(it - ids.begin()) = value;
    else
        throw AmiException("Could not find %s with specified %s", variable_name,
                           id_name);
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

static int setValueByIdRegex(std::vector<std::string> const &ids,
                             std::vector<realtype> &values, realtype value,
                             std::string const &regex,
                             const char *variable_name, const char *id_name) {
    try {
        /* For unknown reasons, the Intel compiler fails to compile patterns
         * such as p[\d]+, which work with g++ and clang.
         * Using std::regex_constants::extended fixes Intel issues, but g++
         * seems to not match the pattern correctly.
         * This is the best solution I was able to come up with...
         */
#ifdef __INTEL_COMPILER
        std::regex pattern(regex, std::regex_constants::extended);
#else
        std::regex pattern(regex);
#endif
        int n_found = 0;
        for (const auto &id : ids) {
            if (std::regex_match(id, pattern)) {
                values.at(&id - &ids[0]) = value;
                ++n_found;
            }
        }

        if (n_found == 0)
            throw AmiException("Could not find %s with specified %s (%s)",
                               variable_name, id_name, regex.c_str());

        return n_found;
    } catch (std::regex_error const &e) {
        auto err_string = regexErrorToString(e.code());
        throw AmiException("Specified regex pattern %s could not be compiled:"
                           " %s (%s)", regex.c_str(), e.what(),
                           err_string.c_str());
    }
}

Model::Model() : dxdotdp(0, 0), x_pos_tmp(0) {}

Model::Model(const int nx_rdata, const int nxtrue_rdata, const int nx_solver,
             const int nxtrue_solver, const int ny, const int nytrue,
             const int nz, const int nztrue, const int ne, const int nJ,
             const int nw, const int ndwdx, const int ndwdp, const int ndxdotdw,
             std::vector<int> ndJydy, const int nnz, const int ubw,
             const int lbw, SecondOrderMode o2mode,
             const std::vector<realtype> &p, std::vector<realtype> k,
             const std::vector<int> &plist, std::vector<realtype> idlist,
             std::vector<int> z2event, const bool pythonGenerated,
             const int ndxdotdp_explicit, const int ndxdotdp_implicit)
    : nx_rdata(nx_rdata), nxtrue_rdata(nxtrue_rdata), nx_solver(nx_solver),
      nxtrue_solver(nxtrue_solver), ny(ny), nytrue(nytrue), nz(nz),
      nztrue(nztrue), ne(ne), nw(nw), ndwdx(ndwdx), ndwdp(ndwdp),
      ndxdotdw(ndxdotdw), ndJydy(std::move(ndJydy)), nnz(nnz), nJ(nJ), ubw(ubw),
      lbw(lbw), pythonGenerated(pythonGenerated),
      ndxdotdp_explicit(ndxdotdp_explicit),
      ndxdotdp_implicit(ndxdotdp_implicit), o2mode(o2mode),
      idlist(std::move(idlist)), J(nx_solver, nx_solver, nnz, CSC_MAT),
      dxdotdw(nx_solver, nw, ndxdotdw, CSC_MAT),
      dwdp(nw, p.size(), ndwdp, CSC_MAT), dwdx(nw, nx_solver, ndwdx, CSC_MAT),
      M(nx_solver, nx_solver), w(nw), x_rdata(nx_rdata, 0.0),
      sx_rdata(nx_rdata, 0.0), x_pos_tmp(nx_solver), originalParameters(p),
      z2event(std::move(z2event)), stateIsNonNegative(nx_solver, false),
      pscale(std::vector<amici::ParameterScaling>(p.size(),
                                                  ParameterScaling::none)) {

    state.h.resize(ne, 0.0);
    state.total_cl.resize(nx_rdata - nx_solver, 0.0);
    state.stotal_cl.resize((nx_rdata - nx_solver) * p.size(), 0.0);
    state.unscaledParameters = p;
    state.fixedParameters = k;
    state.plist = plist;

    /* If Matlab wrapped: dxdotdp is a full AmiVector,
       if Python wrapped: dxdotdp_explicit and dxdotdp_implicit are CSC matrices
     */
    if (pythonGenerated) {
        dxdotdp_explicit =
            SUNMatrixWrapper(nx_solver, p.size(), ndxdotdp_explicit, CSC_MAT);
        dxdotdp_implicit =
            SUNMatrixWrapper(nx_solver, p.size(), ndxdotdp_implicit, CSC_MAT);

        // also dJydy depends on the way of wrapping
        if (static_cast<unsigned>(nytrue) != this->ndJydy.size())
            throw std::runtime_error(
                "Number of elements in ndJydy is not equal "
                " nytrue.");

        for (int iytrue = 0; iytrue < nytrue; ++iytrue)
            dJydy.emplace_back(
                SUNMatrixWrapper(nJ, ny, this->ndJydy[iytrue], CSC_MAT));
    } else {
        dJydy_matlab = std::vector<realtype>(nJ * nytrue * ny, 0.0);
    }
    requireSensitivitiesForAllParameters();
}

bool operator==(const Model &a, const Model &b) {
    if (typeid(a) != typeid(b))
        return false;

    bool bool_dxdotdp = true;
    if (a.pythonGenerated && b.pythonGenerated)
        bool_dxdotdp = (a.ndxdotdp_explicit == b.ndxdotdp_explicit) &&
            (a.ndxdotdp_implicit == b.ndxdotdp_implicit);
    if (a.pythonGenerated != b.pythonGenerated)
        bool_dxdotdp = false;

    return (a.nx_rdata == b.nx_rdata) && (a.nxtrue_rdata == b.nxtrue_rdata) &&
           (a.nx_solver == b.nx_solver) &&
           (a.nxtrue_solver == b.nxtrue_solver) && (a.ny == b.ny) &&
           (a.nytrue == b.nytrue) && (a.nz == b.nz) && (a.nztrue == b.nztrue) &&
           (a.ne == b.ne) && (a.nw == b.nw) && (a.ndwdx == b.ndwdx) &&
           (a.ndwdp == b.ndwdp) && (a.ndxdotdw == b.ndxdotdw) &&
           (a.nnz == b.nnz) && (a.nJ == b.nJ) && (a.ubw == b.ubw) &&
           (a.lbw == b.lbw) && (a.o2mode == b.o2mode) &&
           (a.z2event == b.z2event) && (a.idlist == b.idlist) &&
           (a.state.h == b.state.h) &&
           (a.state.unscaledParameters == b.state.unscaledParameters) &&
           (a.originalParameters == b.originalParameters) &&
           (a.state.fixedParameters == b.state.fixedParameters) &&
           (a.state.plist == b.state.plist) && (a.x0data == b.x0data) &&
           (a.sx0data == b.sx0data) && (a.ts == b.ts) &&
           (a.nmaxevent == b.nmaxevent) && (a.pscale == b.pscale) &&
           (a.stateIsNonNegative == b.stateIsNonNegative) &&
           (a.reinitializeFixedParameterInitialStates ==
            b.reinitializeFixedParameterInitialStates) &&
           (a.tstart == b.tstart) && bool_dxdotdp;
}

void Model::initialize(AmiVector &x, AmiVector &dx, AmiVectorArray &sx,
                       AmiVectorArray & /*sdx*/, bool computeSensitivities) {
    initializeStates(x);
    if (computeSensitivities)
        initializeStateSensitivities(sx, x);

    fdx0(x, dx);
    if (computeSensitivities)
        fsdx0();

    if (ne)
        initHeaviside(x, dx);
}

void Model::initializeB(AmiVector &xB, AmiVector &dxB, AmiVector &xQB) {
    xB.reset();
    dxB.reset();
    xQB.reset();
}

void Model::initializeStates(AmiVector &x) {
    if (x0data.empty()) {
        fx0(x);
    } else {
        std::vector<realtype> x0_solver(nx_solver, 0.0);
        ftotal_cl(state.total_cl.data(), x0data.data());
        fx_solver(x0_solver.data(), x0data.data());
        for (int ix = 0; ix < nx_solver; ix++) {
            x[ix] = (realtype)x0_solver.at(ix);
        }
    }
}

void Model::initializeStateSensitivities(AmiVectorArray &sx,
                                         AmiVector const &x) {
    if (sx0data.empty()) {
        fsx0(sx, x);
    } else {
        realtype *stcl = nullptr;
        std::vector<realtype> sx0_solver_slice(nx_solver, 0.0);
        for (int ip = 0; ip < nplist(); ip++) {
            if (ncl() > 0)
                stcl = &state.stotal_cl.at(plist(ip) * ncl());
            fstotal_cl(stcl, &sx0data.at(ip * nx_rdata), plist(ip));
            fsx_solver(sx0_solver_slice.data(), &sx0data.at(ip * nx_rdata));
            for (int ix = 0; ix < nx_solver; ix++) {
                sx.at(ix, ip) = (realtype)sx0_solver_slice.at(ix);
            }
        }
    }
}

void Model::initHeaviside(AmiVector const &x, AmiVector const &dx) {
    std::vector<realtype> rootvals(ne, 0.0);
    froot(tstart, x, dx, rootvals);
    for (int ie = 0; ie < ne; ie++) {
        if (rootvals.at(ie) < 0) {
            state.h.at(ie) = 0.0;
        } else if (rootvals.at(ie) == 0) {
            throw AmiException(
                "Simulation started in an event. This could lead to "
                "unexpected results, aborting simulation! Please "
                "specify an earlier simulation start via "
                "options.t0");
        } else {
            state.h.at(ie) = 1.0;
        }
    }
}

int Model::nplist() const { return state.plist.size(); }

int Model::np() const { return originalParameters.size(); }

int Model::nk() const { return state.fixedParameters.size(); }

int Model::ncl() const { return nx_rdata - nx_solver; }

const double *Model::k() const { return state.fixedParameters.data(); }

int Model::nMaxEvent() const { return nmaxevent; }

void Model::setNMaxEvent(int nmaxevent) { this->nmaxevent = nmaxevent; }

int Model::nt() const { return ts.size(); }

const std::vector<ParameterScaling> &Model::getParameterScale() const {
    return pscale;
}

void Model::setParameterScale(ParameterScaling pscale) {
    this->pscale.assign(this->pscale.size(), pscale);
    scaleParameters(state.unscaledParameters, this->pscale, originalParameters);
    sx0data.clear();
}

void Model::setParameterScale(std::vector<ParameterScaling> const &pscaleVec) {
    if (pscaleVec.size() != this->originalParameters.size())
        throw AmiException("Dimension mismatch. Size of parameter scaling does "
                           "not match number of model parameters.");
    this->pscale = pscaleVec;
    scaleParameters(state.unscaledParameters, this->pscale, originalParameters);
    sx0data.clear();
}

const std::vector<realtype> &Model::getUnscaledParameters() const {
    return state.unscaledParameters;
}

std::vector<realtype> const &Model::getParameters() const {
    return originalParameters;
}

realtype Model::getParameterById(std::string const &par_id) const {
    if (!hasParameterIds())
        throw AmiException(
            "Could not access parameters by id as they are not set");
    return getValueById(getParameterIds(), originalParameters, par_id,
                        "parameters", "id");
}

realtype Model::getParameterByName(std::string const &par_name) const {
    if (!hasParameterNames())
        throw AmiException(
            "Could not access parameters by name as they are not set");
    return getValueById(getParameterNames(), originalParameters, par_name,
                        "parameters", "name");
}

void Model::setParameters(const std::vector<realtype> &p) {
    if (p.size() != (unsigned)np())
        throw AmiException("Dimension mismatch. Size of parameters does not "
                           "match number of model parameters.");
    this->originalParameters = p;
    this->state.unscaledParameters.resize(originalParameters.size());
    unscaleParameters(originalParameters, pscale, state.unscaledParameters);
}

void Model::setParameterById(const std::map<std::string, realtype> &p,
                             bool ignoreErrors)
{
    for (auto& kv : p) {
        try {
            setParameterById(kv.first, kv.second);
        } catch (AmiException&) {
            if(!ignoreErrors)
                throw;
        }
    }
}

void Model::setParameterById(std::string const &par_id, realtype value) {
    if (!hasParameterIds())
        throw AmiException(
            "Could not access parameters by id as they are not set");

    setValueById(getParameterIds(), originalParameters, value, par_id,
                 "parameter", "id");
    unscaleParameters(originalParameters, pscale, state.unscaledParameters);
}

int Model::setParametersByIdRegex(std::string const &par_id_regex,
                                  realtype value) {
    if (!hasParameterIds())
        throw AmiException(
            "Could not access parameters by id as they are not set");
    int n_found = setValueByIdRegex(getParameterIds(), originalParameters,
                                    value, par_id_regex, "parameter", "id");
    unscaleParameters(originalParameters, pscale, state.unscaledParameters);
    return n_found;
}

void Model::setParameterByName(std::string const &par_name, realtype value) {
    if (!hasParameterNames())
        throw AmiException(
            "Could not access parameters by name as they are not set");

    setValueById(getParameterNames(), originalParameters, value, par_name,
                 "parameter", "name");
    unscaleParameters(originalParameters, pscale, state.unscaledParameters);
}

void Model::setParameterByName(const std::map<std::string, realtype> &p,
                               bool ignoreErrors)
{
    for (auto& kv : p) {
        try {
            setParameterByName(kv.first, kv.second);
        } catch (AmiException&) {
            if(!ignoreErrors)
                throw;
        }
    }
}

int Model::setParametersByNameRegex(std::string const &par_name_regex,
                                    realtype value) {
    if (!hasParameterNames())
        throw AmiException(
            "Could not access parameters by name as they are not set");

    int n_found = setValueByIdRegex(getParameterNames(),
                                    originalParameters,
                                    value, par_name_regex, "parameter", "name");

    unscaleParameters(originalParameters, pscale, state.unscaledParameters);
    return n_found;
}

const std::vector<realtype> &Model::getFixedParameters() const {
    return state.fixedParameters;
}

realtype Model::getFixedParameterById(std::string const &par_id) const {
    if (!hasFixedParameterIds())
        throw AmiException(
            "Could not access fixed parameters by id as they are not set");

    return getValueById(getFixedParameterIds(), state.fixedParameters,
                        par_id, "fixedParameters", "id");
}

realtype Model::getFixedParameterByName(std::string const &par_name) const {
    if (!hasFixedParameterNames())
        throw AmiException(
            "Could not access fixed parameters by name as they are not set");

    return getValueById(getFixedParameterNames(), state.fixedParameters,
                        par_name, "fixedParameters", "name");
}

void Model::setFixedParameters(const std::vector<realtype> &k) {
    if (k.size() != (unsigned)nk())
        throw AmiException("Dimension mismatch. Size of fixedParameters does "
                           "not match number of fixed model parameters.");
    this->state.fixedParameters = k;
}

void Model::setFixedParameterById(std::string const &par_id, realtype value) {
    if (!hasFixedParameterIds())
        throw AmiException(
            "Could not access fixed parameters by id as they are not set");

    setValueById(getFixedParameterIds(), state.fixedParameters, value, par_id,
                 "fixedParameters", "id");
}

int Model::setFixedParametersByIdRegex(std::string const &par_id_regex,
                                       realtype value) {
    if (!hasFixedParameterIds())
        throw AmiException(
            "Could not access fixed parameters by id as they are not set");

    return setValueByIdRegex(getFixedParameterIds(), state.fixedParameters,
                             value, par_id_regex, "fixedParameters", "id");
}

void Model::setFixedParameterByName(std::string const &par_name,
                                    realtype value) {
    if (!hasFixedParameterNames())
        throw AmiException(
            "Could not access fixed parameters by name as they are not set");

    setValueById(getFixedParameterNames(), state.fixedParameters, value,
                 par_name, "fixedParameters", "name");
}

int Model::setFixedParametersByNameRegex(std::string const &par_name_regex,
                                         realtype value) {
    if (!hasFixedParameterNames())
        throw AmiException(
            "Could not access fixed parameters by name as they are not set");

    return setValueByIdRegex(getFixedParameterIds(), state.fixedParameters,
                             value, par_name_regex, "fixedParameters", "name");
}

std::string Model::getName() const {
    return "";
}

bool Model::hasParameterNames() const {
    return np() == 0 || !getParameterNames().empty();
}

std::vector<std::string> Model::getParameterNames() const {
    return std::vector<std::string>();
}

bool Model::hasStateNames() const {
    return nx_rdata == 0 || !getStateNames().empty();
}

std::vector<std::string> Model::getStateNames() const {
    return std::vector<std::string>();
}

bool Model::hasFixedParameterNames() const {
    return nk() == 0 || !getFixedParameterNames().empty();
}

std::vector<std::string> Model::getFixedParameterNames() const {
    return std::vector<std::string>();
}

bool Model::hasObservableNames() const {
    return ny == 0 || !getObservableNames().empty();
}

std::vector<std::string> Model::getObservableNames() const {
    return std::vector<std::string>();
}

bool Model::hasParameterIds() const {
    return np() == 0 || !getParameterIds().empty();
}

std::vector<std::string> Model::getParameterIds() const {
    return std::vector<std::string>();
}

bool Model::hasStateIds() const {
    return nx_rdata == 0 || !getStateIds().empty();
}

std::vector<std::string> Model::getStateIds() const {
    return std::vector<std::string>();
}

bool Model::hasFixedParameterIds() const {
    return nk() == 0 || !getFixedParameterIds().empty();
}

std::vector<std::string> Model::getFixedParameterIds() const {
    return std::vector<std::string>();
}

bool Model::hasObservableIds() const {
    return ny == 0 || !getObservableIds().empty();
}

std::vector<std::string> Model::getObservableIds() const {
    return std::vector<std::string>();
}

std::vector<realtype> const &Model::getTimepoints() const { return ts; }

double Model::getTimepoint(const int it) const { return ts.at(it); }

void Model::setTimepoints(const std::vector<realtype> &ts) {
    if (!std::is_sorted(ts.begin(), ts.end()))
        throw AmiException("Encountered non-monotonic timepoints, please order"
                           " timepoints such that they are monotonically"
                           " increasing!");
    this->ts = ts;
}

double Model::t0() const { return tstart; }

void Model::setT0(double t0) { tstart = t0; }

std::vector<bool> const &Model::getStateIsNonNegative() const {
    return stateIsNonNegative;
}

void Model::setStateIsNonNegative(std::vector<bool> const &nonNegative) {
    if (nx_solver != nx_rdata) {
        throw AmiException("Nonnegative states are not supported whith"
                           " conservation laws enabled");
    }
    if (stateIsNonNegative.size() != static_cast<unsigned long>(nx_rdata)) {
        throw AmiException("Dimension of input stateIsNonNegative (%u) does "
                           "not agree with number of state variables (%d)",
                           stateIsNonNegative.size(), nx_rdata);
    }
    stateIsNonNegative = nonNegative;
    anyStateNonNegative =
        std::any_of(stateIsNonNegative.begin(), stateIsNonNegative.end(),
                    [](bool x) { return x; });
}

void Model::setAllStatesNonNegative() {
    setStateIsNonNegative(std::vector<bool>(nx_solver, true));
}

const std::vector<int> &Model::getParameterList() const { return state.plist; }

int Model::plist(int pos) const { return state.plist.at(pos); }

void Model::setParameterList(const std::vector<int> &plist) {
    int np = this->np(); // cannot capture 'this' in lambda expression
    if (std::any_of(plist.begin(), plist.end(),
                    [&np](int idx) { return idx < 0 || idx >= np; })) {
        throw AmiException("Indices in plist must be in [0..np]");
    }
    this->state.plist = plist;

    initializeVectors();
}

std::vector<realtype> Model::getInitialStates() {
    if(!x0data.empty()) {
        return x0data;
    }

    /* Initial states have not been set explicitly on this instance, so we
     * compute it, but don't save it, as this would have to be invalidated upon
     * changing parameters etc.
     */
    std::vector<realtype> x0(nx_rdata, 0.0);
    fx0(x0.data(), tstart, state.unscaledParameters.data(),
        state.fixedParameters.data());
    return x0;
}

void Model::setInitialStates(const std::vector<realtype> &x0) {
    if (x0.size() != (unsigned)nx_rdata && !x0.empty())
        throw AmiException("Dimension mismatch. Size of x0 does not match "
                           "number of model states.");

    if (x0.empty()) {
        x0data.clear();
        return;
    }

    x0data = x0;
}

bool Model::hasCustomInitialStates() const
{
    return !x0data.empty();
}

std::vector<realtype> Model::getInitialStateSensitivities() {
    if(!sx0data.empty()) {
        return sx0data;
    }

    /* Initial state sensitivities have not been set explicitly on this
     * instance, so we compute it, but don't save it, as this would have to be
     * invalidated upon changing parameters etc.
     */
    std::vector<realtype> sx0(nx_rdata * nplist(), 0.0);
    auto x0 = getInitialStates();
    for (int ip = 0; ip < nplist(); ip++) {
        fsx0(sx0.data(), tstart, x0.data(), state.unscaledParameters.data(),
             state.fixedParameters.data(), plist(ip));
    }
    return sx0;

}

void Model::setInitialStateSensitivities(const std::vector<realtype> &sx0) {
    if (sx0.size() != (unsigned)nx_rdata * nplist() && !sx0.empty())
        throw AmiException("Dimension mismatch. Size of sx0 does not match "
                           "number of model states * number of parameter "
                           "selected for sensitivities.");

    if (sx0.empty()) {
        sx0data.clear();
        return;
    }

    realtype chainrulefactor = 1.0;
    std::vector<realtype> sx0_rdata(nx_rdata * nplist(), 0.0);
    for (int ip = 0; ip < nplist(); ip++) {

        // revert chainrule
        switch (pscale.at(plist(ip))) {
        case ParameterScaling::log10:
            chainrulefactor = state.unscaledParameters.at(plist(ip)) * log(10);
            break;
        case ParameterScaling::ln:
            chainrulefactor = state.unscaledParameters.at(plist(ip));
            break;
        case ParameterScaling::none:
            chainrulefactor = 1.0;
            break;
        }

        for (int ix = 0; ix < nx_rdata; ++ix) {
            sx0_rdata.at(ip * nx_rdata + ix) =
                sx0.at(ip * nx_rdata + ix) / chainrulefactor;
        }
    }
    setUnscaledInitialStateSensitivities(sx0_rdata);
}

bool Model::hasCustomInitialStateSensitivities() const
{
    return !sx0data.empty();
}

void Model::setUnscaledInitialStateSensitivities(
    const std::vector<realtype> &sx0) {
    if (sx0.size() != (unsigned)nx_rdata * nplist() && !sx0.empty())
        throw AmiException("Dimension mismatch. Size of sx0 does not match "
                           "number of model states * number of parameter "
                           "selected for sensitivities.");

    if (sx0.empty()) {
        sx0data.clear();
        return;
    }

    sx0data = sx0;
}

void Model::setSteadyStateSensitivityMode(
    const SteadyStateSensitivityMode mode) {
    steadyStateSensitivityMode = mode;
}

SteadyStateSensitivityMode Model::getSteadyStateSensitivityMode() const {
    return steadyStateSensitivityMode;
}

void Model::setReinitializeFixedParameterInitialStates(bool flag) {
    if (flag && !isFixedParameterStateReinitializationAllowed())
        throw AmiException(
            "State reinitialization cannot be enabled for this model"
            "as this feature was disabled at compile time. Most likely,"
            " this was because some initial states depending on "
            "fixedParameters also depended on parameters");
    reinitializeFixedParameterInitialStates = flag;
}

bool Model::getReinitializeFixedParameterInitialStates() const {
    return reinitializeFixedParameterInitialStates;
}

void Model::requireSensitivitiesForAllParameters() {
    state.plist.resize(np());
    std::iota(state.plist.begin(), state.plist.end(), 0);
    initializeVectors();
}

void Model::getExpression(gsl::span<realtype> w, const realtype t, const AmiVector &x)
{
    fw(t, x.data());
    writeSlice(this->w, w);
}

void Model::getObservable(gsl::span<realtype> y, const realtype t,
                          const AmiVector &x) {
    fy(t, x);
    writeSlice(this->y, y);
}

void Model::getObservableSensitivity(gsl::span<realtype> sy, const realtype t,
                                     const AmiVector &x,
                                     const AmiVectorArray &sx) {
    if (!ny)
        return;

    fdydx(t, x);
    fdydp(t, x);

    this->sx.assign(nx_solver * nplist(), 0.0);
    sx.flatten_to_vector(this->sx);

    // compute sy = 1.0*dydx*sx + 1.0*sy
    // dydx A[ny,nx_solver] * sx B[nx_solver,nplist] = sy C[ny,nplist]
    //        M  K                 K  N                     M  N
    //        lda                  ldb                      ldc
    amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                BLASTranspose::noTrans, ny, nplist(), nx_solver, 1.0,
                dydx.data(), ny, this->sx.data(), nx_solver, 1.0, dydp.data(),
                ny);

    writeSlice(dydp, sy);

    if (alwaysCheckFinite)
        checkFinite(sy, "sy");
}

void Model::getObservableSigma(gsl::span<realtype> sigmay, const int it,
                               const ExpData *edata) {
    fsigmay(it, edata);
    writeSlice(this->sigmay, sigmay);
}

void Model::getObservableSigmaSensitivity(gsl::span<realtype> ssigmay,
                                          const int it, const ExpData *edata) {
    fdsigmaydp(it, edata);
    writeSlice(dsigmaydp, ssigmay);
}

void Model::addObservableObjective(realtype &Jy, const int it,
                                   const AmiVector &x, const ExpData &edata) {
    fy(edata.getTimepoint(it), x);
    fsigmay(it, &edata);

    std::vector<realtype> nllh(nJ, 0.0);
    for (int iyt = 0; iyt < nytrue; iyt++) {
        if (edata.isSetObservedData(it, iyt)) {
            std::fill(nllh.begin(), nllh.end(), 0.0);
            fJy(nllh.data(), iyt, state.unscaledParameters.data(),
                state.fixedParameters.data(), y.data(), sigmay.data(),
                edata.getObservedDataPtr(it));
            Jy -= nllh.at(0);
        }
    }
}

void Model::addObservableObjectiveSensitivity(std::vector<realtype> &sllh,
                                              std::vector<realtype> &s2llh,
                                              const int it, const AmiVector &x,
                                              const AmiVectorArray &sx,
                                              const ExpData &edata) {

    if (!ny)
        return;

    fdJydx(it, x, edata);
    fdJydp(it, x, edata);
    // Compute dJydx*sx for current 'it'
    // dJydx        rdata->nt x nJ        x nx_solver
    // sx           rdata->nt x nx_solver x nplist()
    sx.flatten_to_vector(this->sx);

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                BLASTranspose::noTrans, nJ, nplist(), nx_solver, 1.0,
                dJydx.data(), nJ, this->sx.data(), nx_solver, 1.0, dJydp.data(),
                nJ);

    writeLLHSensitivitySlice(dJydp, sllh, s2llh);
}

void Model::addPartialObservableObjectiveSensitivity(
    std::vector<realtype> &sllh, std::vector<realtype> &s2llh, const int it,
    const AmiVector &x, const ExpData &edata) {
    if (!ny)
        return;

    fdJydp(it, x, edata);

    writeLLHSensitivitySlice(dJydp, sllh, s2llh);
}

void Model::getAdjointStateObservableUpdate(gsl::span<realtype> dJydx,
                                            const int it, const AmiVector &x,
                                            const ExpData &edata) {
    fdJydx(it, x, edata);
    writeSlice(this->dJydx, dJydx);
}

void Model::getEvent(gsl::span<realtype> z, const int ie, const realtype t,
                     const AmiVector &x) {
    fz(ie, t, x);
    writeSliceEvent(this->z, z, ie);
}

void Model::getEventSensitivity(gsl::span<realtype> sz, const int ie,
                                const realtype t, const AmiVector &x,
                                const AmiVectorArray &sx) {
    for (int ip = 0; ip < nplist(); ip++) {
        fsz(&sz.at(ip * nz), ie, t, x.data(), state.unscaledParameters.data(),
            state.fixedParameters.data(), state.h.data(), sx.data(ip),
            plist(ip));
    }
}

void Model::getUnobservedEventSensitivity(gsl::span<realtype> sz,
                                          const int ie) {
    checkBufferSize(sz, nz * nplist());

    for (int iz = 0; iz < nz; ++iz)
        if (z2event[iz] - 1 == ie)
            for (int ip = 0; ip < nplist(); ++ip)
                sz.at(ip * nz + iz) = 0.0;
}

void Model::getEventRegularization(gsl::span<realtype> rz, const int ie,
                                   const realtype t, const AmiVector &x) {
    frz(ie, t, x);
    writeSliceEvent(this->rz, rz, ie);
}

void Model::getEventRegularizationSensitivity(gsl::span<realtype> srz,
                                              const int ie, const realtype t,
                                              const AmiVector &x,
                                              const AmiVectorArray &sx) {
    for (int ip = 0; ip < nplist(); ip++) {
        fsrz(&srz.at(ip * nz), ie, t, x.data(), state.unscaledParameters.data(),
             state.fixedParameters.data(), state.h.data(), sx.data(ip),
             plist(ip));
    }
}

void Model::getEventSigma(gsl::span<realtype> sigmaz, const int ie,
                          const int nroots, const realtype t,
                          const ExpData *edata) {
    fsigmaz(ie, nroots, t, edata);
    writeSliceEvent(this->sigmaz, sigmaz, ie);
}

void Model::getEventSigmaSensitivity(gsl::span<realtype> ssigmaz, const int ie,
                                     const int nroots, const realtype t,
                                     const ExpData *edata) {
    fdsigmazdp(ie, nroots, t, edata);
    writeSensitivitySliceEvent(dsigmazdp, ssigmaz, ie);
}

void Model::addEventObjective(realtype &Jz, const int ie, const int nroots,
                              const realtype t, const AmiVector &x,
                              const ExpData &edata) {
    fz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    std::vector<realtype> nllh(nJ, 0.0);
    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.isSetObservedEvents(nroots, iztrue)) {
            std::fill(nllh.begin(), nllh.end(), 0.0);
            fJz(nllh.data(), iztrue, state.unscaledParameters.data(),
                state.fixedParameters.data(), z.data(), sigmaz.data(),
                edata.getObservedEventsPtr(nroots));
            Jz -= nllh.at(0);
        }
    }
}

void Model::addEventObjectiveRegularization(realtype &Jrz, const int ie,
                                            const int nroots, const realtype t,
                                            const AmiVector &x,
                                            const ExpData &edata) {
    frz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    std::vector<realtype> nllh(nJ, 0.0);
    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.isSetObservedEvents(nroots, iztrue)) {
            std::fill(nllh.begin(), nllh.end(), 0.0);
            fJrz(nllh.data(), iztrue, state.unscaledParameters.data(),
                 state.fixedParameters.data(), rz.data(), sigmaz.data());
            Jrz -= nllh.at(0);
        }
    }
}

void Model::addEventObjectiveSensitivity(std::vector<realtype> &sllh,
                                         std::vector<realtype> &s2llh,
                                         const int ie, const int nroots,
                                         const realtype t, const AmiVector &x,
                                         const AmiVectorArray &sx,
                                         const ExpData &edata) {

    if (!nz)
        return;

    fdJzdx(ie, nroots, t, x, edata);
    fdJzdp(ie, nroots, t, x, edata);

    // sJz           nJ x nplist()
    // dJzdp         nJ x nplist()
    // dJzdx         nmaxevent x nJ        x nx_solver
    // sx            rdata->nt x nx_solver x nplist()

    // Compute dJzdx*sx for current 'ie'
    // dJzdx        rdata->nt x nJ        x nx_solver
    // sx           rdata->nt x nx_solver x nplist()
    sx.flatten_to_vector(this->sx);

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                BLASTranspose::noTrans, nJ, nplist(), nx_solver, 1.0,
                dJzdx.data(), nJ, this->sx.data(), nx_solver, 1.0, dJzdp.data(),
                nJ);

    // sJy += multResult + dJydp
    writeLLHSensitivitySlice(dJzdp, sllh, s2llh);
}

void Model::getAdjointStateEventUpdate(gsl::span<realtype> dJzdx, const int ie,
                                       const int nroots, const realtype t,
                                       const AmiVector &x,
                                       const ExpData &edata) {
    fdJzdx(ie, nroots, t, x, edata);
    writeSlice(this->dJzdx, dJzdx);
}

void Model::addPartialEventObjectiveSensitivity(std::vector<realtype> &sllh,
                                                std::vector<realtype> &s2llh,
                                                const int ie, const int nroots,
                                                const realtype t,
                                                const AmiVector &x,
                                                const ExpData &edata) {
    if (!nz)
        return;

    fdJzdp(ie, nroots, t, x, edata);

    writeLLHSensitivitySlice(dJzdp, sllh, s2llh);
}

void Model::getEventTimeSensitivity(std::vector<realtype> &stau,
                                    const realtype t, const int ie,
                                    const AmiVector &x,
                                    const AmiVectorArray &sx) {

    std::fill(stau.begin(), stau.end(), 0.0);

    for (int ip = 0; ip < nplist(); ip++) {
        fstau(&stau.at(ip), t, x.data(), state.unscaledParameters.data(),
              state.fixedParameters.data(), state.h.data(), sx.data(ip),
              plist(ip), ie);
    }
}

void Model::addStateEventUpdate(AmiVector &x, const int ie, const realtype t,
                                const AmiVector &xdot,
                                const AmiVector &xdot_old) {

    deltax.assign(nx_solver, 0.0);

    // compute update
    fdeltax(deltax.data(), t, x.data(), state.unscaledParameters.data(),
            state.fixedParameters.data(), state.h.data(), ie, xdot.data(),
            xdot_old.data());

    if (alwaysCheckFinite) {
        app->checkFinite(deltax, "deltax");
    }

    // update
    amici_daxpy(nx_solver, 1.0, deltax.data(), 1, x.data(), 1);
}

void Model::addStateSensitivityEventUpdate(AmiVectorArray &sx, const int ie,
                                           const realtype t,
                                           const AmiVector &x_old,
                                           const AmiVector &xdot,
                                           const AmiVector &xdot_old,
                                           const std::vector<realtype> &stau) {
    fw(t, x_old.data());

    for (int ip = 0; ip < nplist(); ip++) {

        deltasx.assign(nx_solver, 0.0);

        // compute update
        fdeltasx(deltasx.data(), t, x_old.data(),
                 state.unscaledParameters.data(), state.fixedParameters.data(),
                 state.h.data(), w.data(), plist(ip), ie, xdot.data(),
                 xdot_old.data(), sx.data(ip), &stau.at(ip));

        if (alwaysCheckFinite) {
            app->checkFinite(deltasx, "deltasx");
        }

        amici_daxpy(nx_solver, 1.0, deltasx.data(), 1, sx.data(ip), 1);
    }
}

void Model::addAdjointStateEventUpdate(AmiVector &xB, const int ie,
                                       const realtype t, const AmiVector &x,
                                       const AmiVector &xdot,
                                       const AmiVector &xdot_old) {

    deltaxB.assign(nx_solver, 0.0);

    // compute update
    fdeltaxB(deltaxB.data(), t, x.data(), state.unscaledParameters.data(),
             state.fixedParameters.data(), state.h.data(), ie, xdot.data(),
             xdot_old.data(), xB.data());

    if (alwaysCheckFinite) {
        app->checkFinite(deltaxB, "deltaxB");
    }

    // apply update
    for (int ix = 0; ix < nxtrue_solver; ++ix)
        for (int iJ = 0; iJ < nJ; ++iJ)
            xB.at(ix + iJ * nxtrue_solver) +=
                deltaxB.at(ix + iJ * nxtrue_solver);
}

void Model::addAdjointQuadratureEventUpdate(
    AmiVector xQB, const int ie, const realtype t, const AmiVector &x,
    const AmiVector &xB, const AmiVector &xdot, const AmiVector &xdot_old) {
    for (int ip = 0; ip < nplist(); ip++) {
        deltaqB.assign(nJ, 0.0);

        fdeltaqB(deltaqB.data(), t, x.data(), state.unscaledParameters.data(),
                 state.fixedParameters.data(), state.h.data(), plist(ip), ie,
                 xdot.data(), xdot_old.data(), xB.data());

        for (int iJ = 0; iJ < nJ; ++iJ)
            xQB.at(iJ) += deltaqB.at(iJ);
    }

    if (alwaysCheckFinite) {
        app->checkFinite(deltaqB, "deltaqB");
    }
}

void Model::updateHeaviside(const std::vector<int> &rootsfound) {
    for (int ie = 0; ie < ne; ie++) {
        state.h.at(ie) += rootsfound.at(ie);
    }
}

void Model::updateHeavisideB(const int *rootsfound) {
    for (int ie = 0; ie < ne; ie++) {
        state.h.at(ie) -= rootsfound[ie];
    }
}


int Model::checkFinite(gsl::span<const realtype> array, const char *fun) const {
    auto result = app->checkFinite(array, fun);

    if (result != AMICI_SUCCESS) {
        app->checkFinite(state.fixedParameters, "k");
        app->checkFinite(state.unscaledParameters, "p");
        app->checkFinite(w, "w");
    }

    return result;
}

void Model::setAlwaysCheckFinite(bool alwaysCheck) {
    alwaysCheckFinite = alwaysCheck;
}

bool Model::getAlwaysCheckFinite() const { return alwaysCheckFinite; }

void Model::fx0(AmiVector &x) {
    std::fill(x_rdata.begin(), x_rdata.end(), 0.0);
    /* this function  also computes initial total abundances */
    fx0(x_rdata.data(), tstart, state.unscaledParameters.data(),
        state.fixedParameters.data());
    fx_solver(x.data(), x_rdata.data());
    ftotal_cl(state.total_cl.data(), x_rdata.data());

    if (alwaysCheckFinite) {
        checkFinite(x_rdata, "x0 x_rdata");
        checkFinite(x.getVector(), "x0 x");
    }
}

void Model::fx0_fixedParameters(AmiVector &x) {
    if (!getReinitializeFixedParameterInitialStates())
        return;
    /* we transform to the unreduced states x_rdata and then apply
     x0_fixedparameters to (i) enable updates to states that were removed from
     conservation laws and (ii) be able to correctly compute total abundances
     after updating the state variables */
    fx_rdata(x_rdata.data(), x.data(), state.total_cl.data());
    fx0_fixedParameters(x_rdata.data(), tstart, state.unscaledParameters.data(),
                        state.fixedParameters.data());
    fx_solver(x.data(), x_rdata.data());
    /* update total abundances */
    ftotal_cl(state.total_cl.data(), x_rdata.data());
}

void Model::fsx0(AmiVectorArray &sx, const AmiVector &x) {
    /* this function  also computes initial total abundance sensitivities */
    realtype *stcl = nullptr;
    for (int ip = 0; ip < nplist(); ip++) {
        if (ncl() > 0)
            stcl = &state.stotal_cl.at(plist(ip) * ncl());
        std::fill(sx_rdata.begin(), sx_rdata.end(), 0.0);
        fsx0(sx_rdata.data(), tstart, x.data(), state.unscaledParameters.data(),
             state.fixedParameters.data(), plist(ip));
        fsx_solver(sx.data(ip), sx_rdata.data());
        fstotal_cl(stcl, sx_rdata.data(), plist(ip));
    }
}

void Model::fsx0_fixedParameters(AmiVectorArray &sx, const AmiVector &x) {
    if (!getReinitializeFixedParameterInitialStates())
        return;
    realtype *stcl = nullptr;
    for (int ip = 0; ip < nplist(); ip++) {
        if (ncl() > 0)
            stcl = &state.stotal_cl.at(plist(ip) * ncl());
        fsx_rdata(sx_rdata.data(), sx.data(ip), stcl, plist(ip));
        fsx0_fixedParameters(sx_rdata.data(), tstart, x.data(),
                             state.unscaledParameters.data(),
                             state.fixedParameters.data(),
                             plist(ip));
        fsx_solver(sx.data(ip), sx_rdata.data());
        fstotal_cl(stcl, sx_rdata.data(), plist(ip));
    }
}

void Model::fsdx0() {}

void Model::fx_rdata(AmiVector &x_rdata, const AmiVector &x) {
    fx_rdata(x_rdata.data(), x.data(), state.total_cl.data());
    if (alwaysCheckFinite)
        checkFinite(x_rdata.getVector(), "x_rdata");
}

void Model::fsx_rdata(AmiVectorArray &sx_rdata, const AmiVectorArray &sx) {
    realtype *stcl = nullptr;
    for (int ip = 0; ip < nplist(); ip++) {
        if (ncl() > 0)
            stcl = &state.stotal_cl.at(plist(ip) * ncl());
        fsx_rdata(sx_rdata.data(ip), sx.data(ip), stcl, ip);
    }
}

void Model::writeSliceEvent(gsl::span<const realtype> slice,
                            gsl::span<realtype> buffer, const int ie) {
    checkBufferSize(buffer, slice.size());
    checkBufferSize(buffer, z2event.size());
    for (unsigned izt = 0; izt < z2event.size(); ++izt)
        if (z2event.at(izt) - 1 == ie)
            buffer.at(izt) = slice.at(izt);
}

void Model::writeSensitivitySliceEvent(gsl::span<const realtype> slice,
                                       gsl::span<realtype> buffer,
                                       const int ie) {
    checkBufferSize(buffer, slice.size());
    checkBufferSize(buffer, z2event.size() * nplist());
    for (int ip = 0; ip < nplist(); ++ip)
        for (unsigned izt = 0; izt < z2event.size(); ++izt)
            if (z2event.at(izt) - 1 == ie)
                buffer.at(ip * nztrue + izt) = slice.at(ip * nztrue + izt);
}

void Model::writeLLHSensitivitySlice(const std::vector<realtype> &dLLhdp,
                                     std::vector<realtype> &sllh,
                                     std::vector<realtype> &s2llh) {
    checkLLHBufferSize(sllh, s2llh);

    amici_daxpy(nplist(), -1.0, dLLhdp.data(), nJ, sllh.data(), 1);
    for (int iJ = 1; iJ < nJ; ++iJ)
        amici_daxpy(nplist(), -1.0, &dLLhdp.at(iJ), nJ, &s2llh.at(iJ - 1),
                    nJ - 1);
}

void Model::checkLLHBufferSize(std::vector<realtype> &sllh,
                               std::vector<realtype> &s2llh) {
    if (sllh.size() != static_cast<unsigned>(nplist()))
        throw AmiException("Incorrect sllh buffer size! Was %u, expected %i.",
                           sllh.size(), nplist());

    if (s2llh.size() != static_cast<unsigned>((nJ - 1) * nplist()))
        throw AmiException("Incorrect s2llh buffer size! Was %u, expected %i.",
                           s2llh.size(), (nJ - 1) * nplist());
}

void Model::initializeVectors() {
    sx0data.clear();
    if (!pythonGenerated)
        dxdotdp = AmiVectorArray(nx_solver, nplist());
}

void Model::fy(const realtype t, const AmiVector &x) {
    if (!ny)
        return;

    y.assign(ny, 0.0);

    fw(t, x.data());
    fy(y.data(), t, x.data(), state.unscaledParameters.data(),
       state.fixedParameters.data(),
       state.h.data(), w.data());

    if (alwaysCheckFinite) {
        app->checkFinite(gsl::make_span(y.data(), ny), "y");
    }
}

void Model::fdydp(const realtype t, const AmiVector &x) {
    if (!ny)
        return;

    dydp.assign(ny * nplist(), 0.0);
    fw(t, x.data());
    fdwdp(t, x.data());

    /* get dydp slice (ny) for current time and parameter */
    for (int ip = 0; ip < nplist(); ip++)
        fdydp(&dydp.at(ip * ny), t, x.data(), state.unscaledParameters.data(),
              state.fixedParameters.data(), state.h.data(), plist(ip), w.data(),
              dwdp.data());

    if (alwaysCheckFinite) {
        app->checkFinite(dydp, "dydp");
    }
}

void Model::fdydx(const realtype t, const AmiVector &x) {
    if (!ny)
        return;

    dydx.assign(ny * nx_solver, 0.0);

    fw(t, x.data());
    fdwdx(t, x.data());
    fdydx(dydx.data(), t, x.data(), state.unscaledParameters.data(),
          state.fixedParameters.data(), state.h.data(), w.data(), dwdx.data());

    if (alwaysCheckFinite) {
        app->checkFinite(dydx, "dydx");
    }
}

void Model::fsigmay(const int it, const ExpData *edata) {
    if (!ny)
        return;

    sigmay.assign(ny, 0.0);

    fsigmay(sigmay.data(), getTimepoint(it), state.unscaledParameters.data(),
            state.fixedParameters.data());

    if (edata) {
        auto sigmay_edata = edata->getObservedDataStdDevPtr(it);
        /* extract the value for the standard deviation from ExpData,
         * if the data value is NaN, use the parameter value */
        for (int iytrue = 0; iytrue < nytrue; iytrue++) {
            if (edata->isSetObservedDataStdDev(it, iytrue))
                sigmay.at(iytrue) = sigmay_edata[iytrue];

            /* TODO: when moving second order code to cpp, verify
             * that this is actually what we want
             */
            for (int iJ = 1; iJ < nJ; iJ++)
                sigmay.at(iytrue + iJ*nytrue) = 0;

            if (edata->isSetObservedData(it, iytrue))
                checkSigmaPositivity(sigmay[iytrue], "sigmay");
        }
    }
}

void Model::fdsigmaydp(const int it, const ExpData *edata) {
    if (!ny)
        return;

    dsigmaydp.assign(ny * nplist(), 0.0);

    for (int ip = 0; ip < nplist(); ip++)
        // get dsigmaydp slice (ny) for current timepoint and parameter
        fdsigmaydp(&dsigmaydp.at(ip * ny), getTimepoint(it),
                   state.unscaledParameters.data(),
                   state.fixedParameters.data(),
                   plist(ip));

    // sigmas in edata override model-sigma -> for those sigmas, set dsigmaydp
    // to zero
    if (edata) {
        for (int iy = 0; iy < nytrue; iy++) {
            if (!edata->isSetObservedDataStdDev(it, iy))
                continue;
            for (int ip = 0; ip < nplist(); ip++) {
                dsigmaydp[ip * ny + iy] = 0.0;
            }
        }
    }

    if (alwaysCheckFinite) {
        app->checkFinite(dsigmaydp, "dsigmaydp");
    }
}

void Model::fdJydy_colptrs(sunindextype * /*indexptrs*/, int /*index*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__);
}

void Model::fdJydy_rowvals(sunindextype * /*indexptrs*/, int /*index*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__);
}

void Model::fdJydy(const int it, const AmiVector &x, const ExpData &edata) {

    fy(edata.getTimepoint(it), x);
    fsigmay(it, &edata);

    if (pythonGenerated) {
        for (int iyt = 0; iyt < nytrue; iyt++) {
            dJydy[iyt].zero();
            fdJydy_colptrs(dJydy[iyt].indexptrs(), iyt);
            fdJydy_rowvals(dJydy[iyt].indexvals(), iyt);

            if (!edata.isSetObservedData(it, iyt))
                continue;

            // get dJydy slice (ny) for current timepoint and observable
            fdJydy(dJydy[iyt].data(), iyt, state.unscaledParameters.data(),
                   state.fixedParameters.data(), y.data(), sigmay.data(),
                   edata.getObservedDataPtr(it));

            if (alwaysCheckFinite) {
                app->checkFinite(gsl::make_span(dJydy[iyt].get()), "dJydy");
            }
        }
    } else {
        std::fill(dJydy_matlab.begin(), dJydy_matlab.end(), 0.0);
        for (int iyt = 0; iyt < nytrue; iyt++) {
            if (!edata.isSetObservedData(it, iyt))
                continue;
            fdJydy(&dJydy_matlab.at(iyt * ny * nJ), iyt,
                   state.unscaledParameters.data(),
                   state.fixedParameters.data(), y.data(), sigmay.data(),
                   edata.getObservedDataPtr(it));
        }
        if (alwaysCheckFinite) {
            // get dJydy slice (ny) for current timepoint and observable
            app->checkFinite(dJydy_matlab, "dJydy");
        }
    }
}

void Model::fdJydsigma(const int it, const AmiVector &x, const ExpData &edata) {

    dJydsigma.assign(nytrue * ny * nJ, 0.0);

    fy(edata.getTimepoint(it), x);
    fsigmay(it, &edata);

    for (int iyt = 0; iyt < nytrue; iyt++) {
        if (edata.isSetObservedData(it, iyt))
            // get dJydsigma slice (ny) for current timepoint and observable
            fdJydsigma(&dJydsigma.at(iyt * ny * nJ), iyt,
                       state.unscaledParameters.data(),
                       state.fixedParameters.data(), y.data(), sigmay.data(),
                       edata.getObservedDataPtr(it));
    }

    if (alwaysCheckFinite) {
        app->checkFinite(dJydsigma, "dJydsigma");
    }
}

void Model::fdJydp(const int it, const AmiVector &x, const ExpData &edata) {
    // dJydy         nJ, nytrue x ny
    // dydp          nplist * ny
    // dJydp         nplist x nJ
    // dJydsigma

    dJydp.assign(nJ * nplist(), 0.0);

    fdJydy(it, x, edata);
    fdydp(edata.getTimepoint(it), x);

    fdJydsigma(it, x, edata);
    fdsigmaydp(it, &edata);

    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (!edata.isSetObservedData(it, iyt))
            continue;

        if (pythonGenerated) {
            // dJydp = 1.0 * dJydp +  1.0 * dJydy * dydp
            for (int iplist = 0; iplist < nplist(); ++iplist) {
                dJydy[iyt].multiply(
                    gsl::span<realtype>(&dJydp.at(iplist * nJ), nJ),
                    gsl::span<const realtype>(&dydp.at(iplist * ny), ny));
            }
        } else {
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nplist(), ny, 1.0,
                        &dJydy_matlab.at(iyt * nJ * ny), nJ, dydp.data(), ny,
                        1.0, dJydp.data(), nJ);
        }
        // dJydp = 1.0 * dJydp +  1.0 * dJydsigma * dsigmaydp
        amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                    BLASTranspose::noTrans, nJ, nplist(), ny, 1.0,
                    &dJydsigma.at(iyt * nJ * ny), nJ, dsigmaydp.data(), ny, 1.0,
                    dJydp.data(), nJ);
    }
}

void Model::fdJydx(const int it, const AmiVector &x, const ExpData &edata) {

    dJydx.assign(nJ * nx_solver, 0.0);

    fdydx(edata.getTimepoint(it), x);
    fdJydy(it, x, edata);

    // dJydy: nJ, ny x nytrue
    // dydx :     ny x nx_solver
    // dJydx:     nJ x nx_solver x nt
    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (!edata.isSetObservedData(it, iyt))
            continue;
        // dJydy A[nyt,nJ,ny] * dydx B[ny,nx_solver] = dJydx C[it,nJ,nx_solver]
        //         slice                                       slice
        //          M  K            K  N                       M  N
        //           lda             ldb                        ldc

        if (pythonGenerated) {
            for (int ix = 0; ix < nx_solver; ++ix) {
                dJydy[iyt].multiply(
                    gsl::span<realtype>(&dJydx.at(ix * nJ), nJ),
                    gsl::span<const realtype>(&dydx.at(ix * ny), ny));
            }
        } else {
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nx_solver, ny, 1.0,
                        &dJydy_matlab.at(iyt * ny * nJ), nJ, dydx.data(), ny,
                        1.0, dJydx.data(), nJ);
        }
    }

    if (alwaysCheckFinite) {
        app->checkFinite(dJydx, "dJydx");
    }
}

void Model::fz(const int ie, const realtype t, const AmiVector &x) {

    z.assign(nz, 0.0);

    fz(z.data(), ie, t, x.data(), state.unscaledParameters.data(),
       state.fixedParameters.data(), state.h.data());
}

void Model::fdzdp(const int ie, const realtype t, const AmiVector &x) {

    dzdp.assign(nz * nplist(), 0.0);

    for (int ip = 0; ip < nplist(); ip++) {
        fdzdp(dzdp.data(), ie, t, x.data(), state.unscaledParameters.data(),
              state.fixedParameters.data(), state.h.data(), plist(ip));
    }

    if (alwaysCheckFinite) {
        app->checkFinite(dzdp, "dzdp");
    }
}

void Model::fdzdx(const int ie, const realtype t, const AmiVector &x) {

    dzdx.assign(nz * nx_solver, 0.0);

    fdzdx(dzdx.data(), ie, t, x.data(), state.unscaledParameters.data(),
          state.fixedParameters.data(), state.h.data());

    if (alwaysCheckFinite) {
        app->checkFinite(dzdx, "dzdx");
    }
}

void Model::frz(const int ie, const realtype t, const AmiVector &x) {

    rz.assign(nz, 0.0);

    frz(rz.data(), ie, t, x.data(), state.unscaledParameters.data(),
        state.fixedParameters.data(), state.h.data());
}

void Model::fdrzdp(const int ie, const realtype t, const AmiVector &x) {

    drzdp.assign(nz * nplist(), 0.0);

    for (int ip = 0; ip < nplist(); ip++) {
        fdrzdp(drzdp.data(), ie, t, x.data(), state.unscaledParameters.data(),
               state.fixedParameters.data(), state.h.data(), plist(ip));
    }

    if (alwaysCheckFinite) {
        app->checkFinite(drzdp, "drzdp");
    }
}

void Model::fdrzdx(const int ie, const realtype t, const AmiVector &x) {

    drzdx.assign(nz * nx_solver, 0.0);

    fdrzdx(drzdx.data(), ie, t, x.data(), state.unscaledParameters.data(),
           state.fixedParameters.data(), state.h.data());

    if (alwaysCheckFinite) {
        app->checkFinite(drzdx, "drzdx");
    }
}

void Model::fsigmaz(const int ie, const int nroots, const realtype t,
                    const ExpData *edata) {
    if (!nz)
        return;

    sigmaz.assign(nz, 0.0);
    fsigmaz(sigmaz.data(), t, state.unscaledParameters.data(),
            state.fixedParameters.data());

    if (edata) {
        for (int iztrue = 0; iztrue < nztrue; iztrue++) {
            if (z2event.at(iztrue) - 1 == ie) {
                if (edata->isSetObservedEventsStdDev(nroots, iztrue)) {
                    auto sigmaz_edata =
                        edata->getObservedEventsStdDevPtr(nroots);
                    sigmaz.at(iztrue) = sigmaz_edata[iztrue];
                }

                /* TODO: when moving second order code to cpp, verify
                 * that this is actually what we want
                 */
                for (int iJ = 1; iJ < nJ; iJ++)
                    sigmaz.at(iztrue + iJ*nztrue) = 0;

                if (edata->isSetObservedEvents(nroots, iztrue))
                    checkSigmaPositivity(sigmaz[iztrue], "sigmaz");
            }
        }
    }
}

void Model::fdsigmazdp(const int ie, const int nroots, const realtype t,
                       const ExpData *edata) {

    dsigmazdp.assign(nz * nplist(), 0.0);

    for (int ip = 0; ip < nplist(); ip++) {
        // get dsigmazdp slice (nz) for current event and parameter
        fdsigmazdp(&dsigmazdp.at(ip * nz), t, state.unscaledParameters.data(),
                   state.fixedParameters.data(), plist(ip));
    }

    // sigmas in edata override model-sigma -> for those sigmas, set dsigmazdp
    // to zero
    if (edata) {
        for (int iz = 0; iz < nztrue; iz++) {
            if (z2event.at(iz) - 1 == ie &&
                !edata->isSetObservedEventsStdDev(nroots, iz)) {
                for (int ip = 0; ip < nplist(); ip++)
                    dsigmazdp.at(iz + nz * ip) = 0;
            }
        }
    }

    if (alwaysCheckFinite) {
        app->checkFinite(dsigmazdp, "dsigmazdp");
    }
}

void Model::fdJzdz(const int ie, const int nroots, const realtype t,
                   const AmiVector &x, const ExpData &edata) {

    dJzdz.assign(nztrue * nz * nJ, 0.0);

    fz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.isSetObservedEvents(nroots, iztrue)) {
            fdJzdz(&dJzdz.at(iztrue * nz * nJ), iztrue,
                   state.unscaledParameters.data(),
                   state.fixedParameters.data(), z.data(), sigmaz.data(),
                   edata.getObservedEventsPtr(nroots));
        }
    }

    if (alwaysCheckFinite) {
        app->checkFinite(dJzdz, "dJzdz");
    }
}

void Model::fdJzdsigma(const int ie, const int nroots, const realtype t,
                       const AmiVector &x, const ExpData &edata) {

    dJzdsigma.assign(nztrue * nz * nJ, 0.0);

    fz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.isSetObservedEvents(nroots, iztrue)) {
            fdJzdsigma(&dJzdsigma.at(iztrue * nz * nJ), iztrue,
                       state.unscaledParameters.data(),
                       state.fixedParameters.data(), z.data(), sigmaz.data(),
                       edata.getObservedEventsPtr(nroots));
        }
    }

    if (alwaysCheckFinite) {
        app->checkFinite(dJzdsigma, "dJzdsigma");
    }
}

void Model::fdJzdp(const int ie, const int nroots, realtype t,
                   const AmiVector &x, const ExpData &edata) {
    // dJzdz         nJ x nz x nztrue
    // dJzdsigma     nJ x nz x nztrue
    // dzdp          nz x nplist()
    // dJzdp         nJ x nplist()

    dJzdp.assign(nJ * nplist(), 0.0);

    fdJzdz(ie, nroots, t, x, edata);
    fdzdp(ie, t, x);

    fdJzdsigma(ie, nroots, t, x, edata);
    fdsigmazdp(ie, nroots, t, &edata);

    for (int izt = 0; izt < nztrue; ++izt) {
        if (!edata.isSetObservedEvents(nroots, izt))
            continue;

        if (t < edata.getTimepoint(edata.nt() - 1)) {
            // with z
            fdJzdz(ie, nroots, t, x, edata);
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nplist(), nz, 1.0,
                        &dJzdz.at(izt * nz * nJ), nJ, dzdp.data(), nz, 1.0,
                        dJzdp.data(), nJ);
        } else {
            // with rz
            fdJrzdz(ie, nroots, t, x, edata);
            fdJrzdsigma(ie, nroots, t, x, edata);

            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nplist(), nz, 1.0,
                        &dJrzdsigma.at(izt * nz * nJ), nJ, dsigmazdp.data(), nz,
                        1.0, dJzdp.data(), nJ);

            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nplist(), nz, 1.0,
                        &dJrzdz.at(izt * nz * nJ), nJ, dzdp.data(), nz, 1.0,
                        dJzdp.data(), nJ);
        }

        amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                    BLASTranspose::noTrans, nJ, nplist(), nz, 1.0,
                    &dJzdsigma.at(izt * nz * nJ), nJ, dsigmazdp.data(), nz, 1.0,
                    dJzdp.data(), nJ);
    }
}

void Model::fdJzdx(const int ie, const int nroots, const realtype t,
                   const AmiVector &x, const ExpData &edata) {
    // dJzdz         nJ x nz        x nztrue
    // dzdx          nz x nx_solver
    // dJzdx         nJ x nx_solver x nmaxevent

    dJzdx.assign(nJ * nx_solver, 0.0);

    fdJzdz(ie, nroots, t, x, edata);

    for (int izt = 0; izt < nztrue; ++izt) {
        if (!edata.isSetObservedEvents(nroots, izt))
            continue;

        if (t < edata.getTimepoint(edata.nt() - 1)) {
            // z
            fdzdx(ie, t, x);
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nx_solver, nz, 1.0,
                        &dJzdz.at(izt * nz * nJ), nJ, dzdx.data(), nz, 1.0,
                        dJzdx.data(), nJ);
        } else {
            // rz
            fdJrzdz(ie, nroots, t, x, edata);
            fdrzdx(ie, t, x);
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nx_solver, nz, 1.0,
                        &dJrzdz.at(izt * nz * nJ), nJ, drzdx.data(), nz, 1.0,
                        dJzdx.data(), nJ);
        }
    }
}

void Model::fdJrzdz(const int ie, const int nroots, const realtype t,
                    const AmiVector &x, const ExpData &edata) {

    dJrzdz.assign(nztrue * nz * nJ, 0.0);

    frz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.isSetObservedEvents(nroots, iztrue)) {
            fdJrzdz(&dJrzdz.at(iztrue * nz * nJ), iztrue,
                    state.unscaledParameters.data(),
                    state.fixedParameters.data(), rz.data(), sigmaz.data());
        }
    }

    if (alwaysCheckFinite) {
        app->checkFinite(dJrzdz, "dJrzdz");
    }
}

void Model::fdJrzdsigma(const int ie, const int nroots, const realtype t,
                        const AmiVector &x, const ExpData &edata) {

    dJrzdsigma.assign(nztrue * nz * nJ, 0.0);

    frz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.isSetObservedEvents(nroots, iztrue)) {
            fdJrzdsigma(&dJrzdsigma.at(iztrue * nz * nJ), iztrue,
                        state.unscaledParameters.data(),
                        state.fixedParameters.data(), rz.data(), sigmaz.data());
        }
    }

    if (alwaysCheckFinite) {
        app->checkFinite(dJrzdsigma, "dJrzdsigma");
    }
}

void Model::fw(const realtype t, const realtype *x) {
    std::fill(w.begin(), w.end(), 0.0);
    fw(w.data(), t, x, state.unscaledParameters.data(),
       state.fixedParameters.data(), state.h.data(), state.total_cl.data());

    if (alwaysCheckFinite) {
        app->checkFinite(w, "w");
    }
}

void Model::fdwdp(const realtype t, const realtype *x) {
    fw(t, x);
    if (pythonGenerated) {
        dwdp.reset();

        // avoid bad memory access when slicing
        if (!nw)
            return;

        fdwdp_colptrs(dwdp.indexptrs());
        fdwdp_rowvals(dwdp.indexvals());
        fdwdp(dwdp.data(), t, x, state.unscaledParameters.data(),
              state.fixedParameters.data(), state.h.data(), w.data(),
              state.total_cl.data(), state.stotal_cl.data());

    } else {
        // matlab generated
        fdwdp(dwdp.data(), t, x, state.unscaledParameters.data(),
              state.fixedParameters.data(), state.h.data(), w.data(),
              state.total_cl.data(), state.stotal_cl.data());
    }

    if (alwaysCheckFinite) {
        app->checkFinite(gsl::make_span(dwdp.get()), "dwdp");
    }
}

void Model::fdwdx(const realtype t, const realtype *x) {
    fw(t, x);
    dwdx.reset();

    fdwdx_colptrs(dwdx.indexptrs());
    fdwdx_rowvals(dwdx.indexvals());
    fdwdx(dwdx.data(), t, x, state.unscaledParameters.data(),
          state.fixedParameters.data(), state.h.data(), w.data(),
          state.total_cl.data());

    if (alwaysCheckFinite) {
        app->checkFinite(gsl::make_span(dwdx.get()), "dwdx");
    }
}

void Model::fx_rdata(realtype *x_rdata, const realtype *x_solver,
                     const realtype * /*tcl*/) {
    if (nx_solver != nx_rdata)
        throw AmiException(
            "A model that has differing nx_solver and nx_rdata needs "
            "to implement its own fx_rdata");
    std::copy_n(x_solver, nx_solver, x_rdata);
}

void Model::fsx_rdata(realtype *sx_rdata, const realtype *sx_solver,
                      const realtype *stcl, const int /*ip*/) {
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

void Model::ftotal_cl(realtype * /*total_cl*/, const realtype * /*x_rdata*/) {
    if (nx_solver != nx_rdata)
        throw AmiException(
            "A model that has differing nx_solver and nx_rdata needs "
            "to implement its own ftotal_cl");
}

void Model::fstotal_cl(realtype *stotal_cl, const realtype *sx_rdata,
                       const int /*ip*/) {
    /* for the moment we do not need an implementation of fstotal_cl as
     * we can simply reuse ftotal_cl and replace states by their
     * sensitivities */
    ftotal_cl(stotal_cl, sx_rdata);
}

N_Vector Model::computeX_pos(const_N_Vector x) {
    if (anyStateNonNegative) {
        for (int ix = 0; ix < x_pos_tmp.getLength(); ++ix) {
            x_pos_tmp.at(ix) =
                (stateIsNonNegative.at(ix) && NV_Ith_S(x, ix) < 0)
                    ? 0
                    : NV_Ith_S(x, ix);
        }
        return x_pos_tmp.getNVector();
    }

    return x;
}

} // namespace amici
