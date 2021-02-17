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
#include <assert.h>

namespace amici {

/**
 * @brief local helper function to get parameters
 * @param ids vector of name/ids of (fixed)Parameters
 * @param values values of the (fixed)Parameters
 * @param id name/id to look for in the vector
 * @param variable_name string indicating what variable we are looking at
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
 * @param variable_name string indicating what variable we are looking at
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
 * @param variable_name string indicating what variable we are looking at
 * @param id_name string indicating whether name or id was specified
 * @return number of matched names/ids
 */

static int setValueByIdRegex(std::vector<std::string> const &ids,
                             std::vector<realtype> &values, realtype value,
                             std::string const &regex,
                             const char *variable_name, const char *id_name) {
    try {
        std::regex pattern(regex);
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

Model::Model(ModelDimensions const & model_dimensions,
             SimulationParameters simulation_parameters,
             SecondOrderMode o2mode, std::vector<realtype> idlist, std::vector<int> z2event,
             const bool pythonGenerated, const int ndxdotdp_explicit,
             const int ndxdotdx_explicit, const int w_recursion_depth)
    : ModelDimensions(model_dimensions), pythonGenerated(pythonGenerated),
      o2mode(o2mode), idlist(std::move(idlist)),
      derived_state_(model_dimensions),
      z2event_(std::move(z2event)),
      state_is_non_negative_(nx_solver, false),
      w_recursion_depth_(w_recursion_depth),
      simulation_parameters_(std::move(simulation_parameters)) {
    Expects(model_dimensions.np == static_cast<int>(simulation_parameters_.parameters.size()));
    Expects(model_dimensions.nk == static_cast<int>(simulation_parameters_.fixedParameters.size()));

    simulation_parameters.pscale = std::vector<ParameterScaling>(model_dimensions.np, ParameterScaling::none);

    state_.h.resize(ne, 0.0);
    state_.total_cl.resize(nx_rdata - nx_solver, 0.0);
    state_.stotal_cl.resize((nx_rdata - nx_solver) * np(), 0.0);
    state_.unscaledParameters.resize(np());
    unscaleParameters(simulation_parameters_.parameters,
                      simulation_parameters_.pscale, state_.unscaledParameters);
    state_.fixedParameters = simulation_parameters_.fixedParameters;
    state_.plist = simulation_parameters_.plist;

    /* If Matlab wrapped: dxdotdp is a full AmiVector,
       if Python wrapped: dxdotdp_explicit and dxdotdp_implicit are CSC matrices
     */
    if (pythonGenerated) {

        dwdw_ = SUNMatrixWrapper(nw, nw, ndwdw, CSC_MAT);
        // size dynamically adapted for dwdx_ and dwdp_
        derived_state_.dwdx_ = SUNMatrixWrapper(nw, nx_solver, 0, CSC_MAT);
        derived_state_.dwdp_ = SUNMatrixWrapper(nw, np(), 0, CSC_MAT);

        for (int irec = 0; irec <= w_recursion_depth_; ++irec) {
            /* for the first element we know the exact size, while for all others we
               guess the size*/
            dwdp_hierarchical_.emplace_back(
                SUNMatrixWrapper(nw, np(), irec * ndwdw + ndwdp, CSC_MAT));
            dwdx_hierarchical_.emplace_back(
                SUNMatrixWrapper(nw, nx_solver, irec * ndwdw + ndwdx, CSC_MAT));
        }
        assert(static_cast<int>(dwdp_hierarchical_.size()) ==
               w_recursion_depth_ + 1);
        assert(static_cast<int>(dwdx_hierarchical_.size()) ==
               w_recursion_depth_ + 1);

        derived_state_.dxdotdp_explicit = SUNMatrixWrapper(
                    nx_solver, np(), ndxdotdp_explicit, CSC_MAT);
        // guess size, will be dynamically reallocated
        derived_state_.dxdotdp_implicit = SUNMatrixWrapper(
                    nx_solver, np(), ndwdp + ndxdotdw, CSC_MAT);
        derived_state_.dxdotdx_explicit = SUNMatrixWrapper(
                    nx_solver, nx_solver, ndxdotdx_explicit, CSC_MAT);
        // guess size, will be dynamically reallocated
        derived_state_.dxdotdx_implicit = SUNMatrixWrapper(
                    nx_solver, nx_solver, ndwdx + ndxdotdw, CSC_MAT);
        // dynamically allocate on first call
        derived_state_.dxdotdp_full = SUNMatrixWrapper(
                    nx_solver, np(), 0, CSC_MAT);

        for (int iytrue = 0; iytrue < nytrue; ++iytrue)
            derived_state_.dJydy_.emplace_back(
                SUNMatrixWrapper(nJ, ny, ndJydy.at(iytrue), CSC_MAT));
    } else {
        derived_state_.dwdx_ = SUNMatrixWrapper(nw, nx_solver, ndwdx, CSC_MAT);
        derived_state_.dwdp_ = SUNMatrixWrapper(nw, np(), ndwdp, CSC_MAT);
        derived_state_.dJydy_matlab_ = std::vector<realtype>(
                    nJ * nytrue * ny, 0.0);
    }
    requireSensitivitiesForAllParameters();
}

bool operator==(const Model &a, const Model &b) {
    if (typeid(a) != typeid(b))
        return false;

    return (static_cast<ModelDimensions const&>(a)
            == static_cast<ModelDimensions const&>(b))
            && (a.o2mode == b.o2mode) &&
           (a.z2event_ == b.z2event_) && (a.idlist == b.idlist) &&
           (a.state_.h == b.state_.h) &&
           (a.state_.unscaledParameters == b.state_.unscaledParameters) &&
           (a.simulation_parameters_ == b.simulation_parameters_) &&
           (a.state_.fixedParameters == b.state_.fixedParameters) &&
           (a.state_.plist == b.state_.plist) && (a.x0data_ == b.x0data_) &&
           (a.sx0data_ == b.sx0data_) &&
           (a.nmaxevent_ == b.nmaxevent_) &&
           (a.state_is_non_negative_ == b.state_is_non_negative_);
}

bool operator==(const ModelDimensions &a, const ModelDimensions &b) {
    if (typeid(a) != typeid(b))
        return false;
    return (a.nx_rdata == b.nx_rdata) && (a.nxtrue_rdata == b.nxtrue_rdata) &&
           (a.nx_solver == b.nx_solver) &&
           (a.nxtrue_solver == b.nxtrue_solver) &&
           (a.nx_solver_reinit == b.nx_solver_reinit) &&
           (a.np == b.np) && (a.nk == b.nk) && (a.ny == b.ny) &&
           (a.nytrue == b.nytrue) && (a.nz == b.nz) && (a.nztrue == b.nztrue) &&
           (a.ne == b.ne) && (a.nw == b.nw) && (a.ndwdx == b.ndwdx) &&
           (a.ndwdp == b.ndwdp) && (a.ndwdw == b.ndwdw) &&
           (a.ndxdotdw == b.ndxdotdw) && (a.ndJydy == b.ndJydy) &&
           (a.nnz == b.nnz) && (a.nJ == b.nJ) && (a.ubw == b.ubw) &&
           (a.lbw == b.lbw);
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

void Model::initializeB(AmiVector &xB, AmiVector &dxB, AmiVector &xQB,
                        bool posteq) const {
    xB.zero();
    dxB.zero();
    if (!posteq)
        xQB.zero();
}

void Model::initializeStates(AmiVector &x) {
    if (x0data_.empty()) {
        fx0(x);
    } else {
        std::vector<realtype> x0_solver(nx_solver, 0.0);
        ftotal_cl(state_.total_cl.data(), x0data_.data());
        fx_solver(x0_solver.data(), x0data_.data());
        std::copy(x0_solver.cbegin(), x0_solver.cend(), x.data());
    }
}

void Model::initializeStateSensitivities(AmiVectorArray &sx,
                                         AmiVector const &x) {
    if (sx0data_.empty()) {
        fsx0(sx, x);
    } else {
        realtype *stcl = nullptr;
        std::vector<realtype> sx0_solver_slice(nx_solver, 0.0);
        for (int ip = 0; ip < nplist(); ip++) {
            if (ncl() > 0)
                stcl = &state_.stotal_cl.at(plist(ip) * ncl());
            fstotal_cl(stcl, &sx0data_.at(ip * nx_rdata), plist(ip));
            fsx_solver(sx0_solver_slice.data(), &sx0data_.at(ip * nx_rdata));
            for (int ix = 0; ix < nx_solver; ix++) {
                sx.at(ix, ip) = sx0_solver_slice.at(ix);
            }
        }
    }
}

void Model::initHeaviside(AmiVector const &x, AmiVector const &dx) {
    std::vector<realtype> rootvals(ne, 0.0);
    froot(simulation_parameters_.tstart_, x, dx, rootvals);
    for (int ie = 0; ie < ne; ie++) {
        if (rootvals.at(ie) < 0) {
            state_.h.at(ie) = 0.0;
        } else {
            state_.h.at(ie) = 1.0;
        }
    }
}

int Model::nplist() const { return static_cast<int>(state_.plist.size()); }

int Model::np() const { return static_cast<int>(static_cast<ModelDimensions const&>(*this).np); }

int Model::nk() const { return static_cast<int>(state_.fixedParameters.size()); }

int Model::ncl() const { return nx_rdata - nx_solver; }

int Model::nx_reinit() const { return nx_solver_reinit; }

const double *Model::k() const { return state_.fixedParameters.data(); }

int Model::nMaxEvent() const { return nmaxevent_; }

void Model::setNMaxEvent(int nmaxevent) { nmaxevent_ = nmaxevent; }

int Model::nt() const { return static_cast<int>(simulation_parameters_.ts_.size()); }

const std::vector<ParameterScaling> &Model::getParameterScale() const {
    return simulation_parameters_.pscale;
}

void Model::setParameterScale(ParameterScaling pscale) {
    simulation_parameters_.pscale.assign(simulation_parameters_.pscale.size(), pscale);
    scaleParameters(state_.unscaledParameters, simulation_parameters_.pscale,
                    simulation_parameters_.parameters);
    sx0data_.clear();
}

void Model::setParameterScale(std::vector<ParameterScaling> const &pscaleVec) {
    if (pscaleVec.size() != simulation_parameters_.parameters.size())
        throw AmiException("Dimension mismatch. Size of parameter scaling does "
                           "not match number of model parameters.");
    simulation_parameters_.pscale = pscaleVec;
    scaleParameters(state_.unscaledParameters, simulation_parameters_.pscale,
                    simulation_parameters_.parameters);
    sx0data_.clear();
}

const std::vector<realtype> &Model::getUnscaledParameters() const {
    return state_.unscaledParameters;
}

std::vector<realtype> const &Model::getParameters() const {
    return simulation_parameters_.parameters;
}

realtype Model::getParameterById(std::string const &par_id) const {
    if (!hasParameterIds())
        throw AmiException(
            "Could not access parameters by id as they are not set");
    return getValueById(getParameterIds(), simulation_parameters_.parameters,
                        par_id, "parameters", "id");
}

realtype Model::getParameterByName(std::string const &par_name) const {
    if (!hasParameterNames())
        throw AmiException(
            "Could not access parameters by name as they are not set");
    return getValueById(getParameterNames(), simulation_parameters_.parameters,
                        par_name, "parameters", "name");
}

void Model::setParameters(const std::vector<realtype> &p) {
    if (p.size() != (unsigned)np())
        throw AmiException("Dimension mismatch. Size of parameters does not "
                           "match number of model parameters.");
    simulation_parameters_.parameters = p;
    state_.unscaledParameters.resize(simulation_parameters_.parameters.size());
    unscaleParameters(simulation_parameters_.parameters,
                      simulation_parameters_.pscale,
                      state_.unscaledParameters);
}

void Model::setParameterById(const std::map<std::string, realtype> &p,
                             bool ignoreErrors)
{
    for (auto& kv : p) {
        try {
            setParameterById(kv.first, kv.second);
        } catch (AmiException const&) {
            if(!ignoreErrors)
                throw;
        }
    }
}

void Model::setParameterById(std::string const &par_id, realtype value) {
    if (!hasParameterIds())
        throw AmiException(
            "Could not access parameters by id as they are not set");

    setValueById(getParameterIds(), simulation_parameters_.parameters,
                 value, par_id, "parameter", "id");
    unscaleParameters(simulation_parameters_.parameters,
                      simulation_parameters_.pscale,
                      state_.unscaledParameters);
}

int Model::setParametersByIdRegex(std::string const &par_id_regex,
                                  realtype value) {
    if (!hasParameterIds())
        throw AmiException(
            "Could not access parameters by id as they are not set");
    int n_found = setValueByIdRegex(getParameterIds(),
                                    simulation_parameters_.parameters,
                                    value, par_id_regex, "parameter", "id");
    unscaleParameters(simulation_parameters_.parameters,
                      simulation_parameters_.pscale, state_.unscaledParameters);
    return n_found;
}

void Model::setParameterByName(std::string const &par_name, realtype value) {
    if (!hasParameterNames())
        throw AmiException(
            "Could not access parameters by name as they are not set");

    setValueById(getParameterNames(), simulation_parameters_.parameters,
                 value, par_name, "parameter", "name");
    unscaleParameters(simulation_parameters_.parameters,
                      simulation_parameters_.pscale, state_.unscaledParameters);
}

void Model::setParameterByName(const std::map<std::string, realtype> &p,
                               bool ignoreErrors)
{
    for (auto& kv : p) {
        try {
            setParameterByName(kv.first, kv.second);
        } catch (AmiException const&) {
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
                                    simulation_parameters_.parameters,
                                    value, par_name_regex, "parameter", "name");

    unscaleParameters(simulation_parameters_.parameters,
                      simulation_parameters_.pscale, state_.unscaledParameters);
    return n_found;
}

const std::vector<realtype> &Model::getFixedParameters() const {
    return state_.fixedParameters;
}

realtype Model::getFixedParameterById(std::string const &par_id) const {
    if (!hasFixedParameterIds())
        throw AmiException(
            "Could not access fixed parameters by id as they are not set");

    return getValueById(getFixedParameterIds(), state_.fixedParameters,
                        par_id, "fixedParameters", "id");
}

realtype Model::getFixedParameterByName(std::string const &par_name) const {
    if (!hasFixedParameterNames())
        throw AmiException(
            "Could not access fixed parameters by name as they are not set");

    return getValueById(getFixedParameterNames(), state_.fixedParameters,
                        par_name, "fixedParameters", "name");
}

void Model::setFixedParameters(const std::vector<realtype> &k) {
    if (k.size() != (unsigned)nk())
        throw AmiException("Dimension mismatch. Size of fixedParameters does "
                           "not match number of fixed model parameters.");
    state_.fixedParameters = k;
}

void Model::setFixedParameterById(std::string const &par_id, realtype value) {
    if (!hasFixedParameterIds())
        throw AmiException(
            "Could not access fixed parameters by id as they are not set");

    setValueById(getFixedParameterIds(), state_.fixedParameters, value, par_id,
                 "fixedParameters", "id");
}

int Model::setFixedParametersByIdRegex(std::string const &par_id_regex,
                                       realtype value) {
    if (!hasFixedParameterIds())
        throw AmiException(
            "Could not access fixed parameters by id as they are not set");

    return setValueByIdRegex(getFixedParameterIds(), state_.fixedParameters,
                             value, par_id_regex, "fixedParameters", "id");
}

void Model::setFixedParameterByName(std::string const &par_name,
                                    realtype value) {
    if (!hasFixedParameterNames())
        throw AmiException(
            "Could not access fixed parameters by name as they are not set");

    setValueById(getFixedParameterNames(), state_.fixedParameters, value,
                 par_name, "fixedParameters", "name");
}

int Model::setFixedParametersByNameRegex(std::string const &par_name_regex,
                                         realtype value) {
    if (!hasFixedParameterNames())
        throw AmiException(
            "Could not access fixed parameters by name as they are not set");

    return setValueByIdRegex(getFixedParameterIds(), state_.fixedParameters,
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

bool Model::hasExpressionNames() const {
    return ny == 0 || !getExpressionNames().empty();
}

std::vector<std::string> Model::getExpressionNames() const {
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

bool Model::hasExpressionIds() const {
    return ny == 0 || !getExpressionIds().empty();
}

std::vector<std::string> Model::getExpressionIds() const {
    return std::vector<std::string>();
}


bool Model::hasQuadraticLLH() const {
    return true;
}

std::vector<realtype> const &Model::getTimepoints() const { return simulation_parameters_.ts_; }

double Model::getTimepoint(const int it) const { return simulation_parameters_.ts_.at(it); }

void Model::setTimepoints(const std::vector<realtype> &ts) {
    if (!std::is_sorted(ts.begin(), ts.end()))
        throw AmiException("Encountered non-monotonic timepoints, please order"
                           " timepoints such that they are monotonically"
                           " increasing!");
    simulation_parameters_.ts_ = ts;
}

double Model::t0() const { return simulation_parameters_.tstart_; }

void Model::setT0(double t0) { simulation_parameters_.tstart_ = t0; }

std::vector<bool> const &Model::getStateIsNonNegative() const {
    return state_is_non_negative_;
}

void Model::setStateIsNonNegative(std::vector<bool> const &nonNegative) {
    if (nx_solver != nx_rdata) {
        throw AmiException("Non-negative states are not supported with"
                           " conservation laws enabled");
    }
    if (state_is_non_negative_.size() != static_cast<unsigned long>(nx_rdata)) {
        throw AmiException("Dimension of input stateIsNonNegative (%u) does "
                           "not agree with number of state variables (%d)",
                           state_is_non_negative_.size(), nx_rdata);
    }
    state_is_non_negative_ = nonNegative;
    any_state_non_negative_ =
        std::any_of(state_is_non_negative_.begin(), state_is_non_negative_.end(),
                    [](bool x) { return x; });
}

void Model::setAllStatesNonNegative() {
    setStateIsNonNegative(std::vector<bool>(nx_solver, true));
}

const std::vector<int> &Model::getParameterList() const { return state_.plist; }

int Model::plist(int pos) const { return state_.plist.at(pos); }

void Model::setParameterList(const std::vector<int> &plist) {
    int np = this->np(); // cannot capture 'this' in lambda expression
    if (std::any_of(plist.begin(), plist.end(),
                    [&np](int idx) { return idx < 0 || idx >= np; })) {
        throw AmiException("Indices in plist must be in [0..np]");
    }
    state_.plist = plist;

    initializeVectors();
}

std::vector<realtype> Model::getInitialStates() {
    if(!x0data_.empty()) {
        return x0data_;
    }

    /* Initial states have not been set explicitly on this instance, so we
     * compute it, but don't save it, as this would have to be invalidated upon
     * changing parameters etc.
     */
    std::vector<realtype> x0(nx_rdata, 0.0);
    fx0(x0.data(), simulation_parameters_.tstart_, state_.unscaledParameters.data(),
        state_.fixedParameters.data());
    return x0;
}

void Model::setInitialStates(const std::vector<realtype> &x0) {
    if (x0.size() != (unsigned)nx_rdata && !x0.empty())
        throw AmiException("Dimension mismatch. Size of x0 does not match "
                           "number of model states.");

    if (x0.empty()) {
        x0data_.clear();
        return;
    }

    x0data_ = x0;
}

bool Model::hasCustomInitialStates() const
{
    return !x0data_.empty();
}

std::vector<realtype> Model::getInitialStateSensitivities() {
    if(!sx0data_.empty()) {
        return sx0data_;
    }

    /* Initial state sensitivities have not been set explicitly on this
     * instance, so we compute it, but don't save it, as this would have to be
     * invalidated upon changing parameters etc.
     */
    std::vector<realtype> sx0(nx_rdata * nplist(), 0.0);
    auto x0 = getInitialStates();
    for (int ip = 0; ip < nplist(); ip++) {
        fsx0(sx0.data(), simulation_parameters_.tstart_, x0.data(),
             state_.unscaledParameters.data(),
             state_.fixedParameters.data(), plist(ip));
    }
    return sx0;

}

void Model::setInitialStateSensitivities(const std::vector<realtype> &sx0) {
    if (sx0.size() != (unsigned)nx_rdata * nplist() && !sx0.empty())
        throw AmiException("Dimension mismatch. Size of sx0 does not match "
                           "number of model states * number of parameter "
                           "selected for sensitivities.");

    if (sx0.empty()) {
        sx0data_.clear();
        return;
    }

    realtype chainrulefactor = 1.0;
    std::vector<realtype> sx0_rdata(nx_rdata * nplist(), 0.0);
    for (int ip = 0; ip < nplist(); ip++) {

        // revert chainrule
        switch (simulation_parameters_.pscale.at(plist(ip))) {
        case ParameterScaling::log10:
            chainrulefactor = state_.unscaledParameters.at(plist(ip)) * log(10);
            break;
        case ParameterScaling::ln:
            chainrulefactor = state_.unscaledParameters.at(plist(ip));
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
    return !sx0data_.empty();
}

void Model::setUnscaledInitialStateSensitivities(
    const std::vector<realtype> &sx0) {
    if (sx0.size() != (unsigned)nx_rdata * nplist() && !sx0.empty())
        throw AmiException("Dimension mismatch. Size of sx0 does not match "
                           "number of model states * number of parameter "
                           "selected for sensitivities.");

    if (sx0.empty()) {
        sx0data_.clear();
        return;
    }

    sx0data_ = sx0;
}

void Model::setSteadyStateSensitivityMode(
    const SteadyStateSensitivityMode mode) {
    steadystate_sensitivity_mode_ = mode;
}

SteadyStateSensitivityMode Model::getSteadyStateSensitivityMode() const {
    return steadystate_sensitivity_mode_;
}

void Model::setReinitializeFixedParameterInitialStates(bool flag) {
    if (flag && !isFixedParameterStateReinitializationAllowed())
        throw AmiException(
            "State reinitialization cannot be enabled for this model "
            "as this feature was disabled at compile time. Most likely,"
            " this was because some initial states depending on "
            "fixedParameters also depended on parameters.");
    simulation_parameters_.reinitializeFixedParameterInitialStates = flag;

    if(flag) {
        simulation_parameters_.reinitializeAllFixedParameterDependentInitialStatesForSimulation(nx_rdata);
    } else {
        simulation_parameters_.reinitialization_state_idxs_sim.clear();
    }
}

bool Model::getReinitializeFixedParameterInitialStates() const {
    return simulation_parameters_.reinitializeFixedParameterInitialStates
            || !simulation_parameters_.reinitialization_state_idxs_sim.empty();
}

void Model::requireSensitivitiesForAllParameters() {
    state_.plist.resize(np());
    std::iota(state_.plist.begin(), state_.plist.end(), 0);
    initializeVectors();
}

void Model::getExpression(gsl::span<realtype> w, const realtype t, const AmiVector &x)
{
    fw(t, x.data());
    writeSlice(derived_state_.w_, w);
}

void Model::getObservable(gsl::span<realtype> y, const realtype t,
                          const AmiVector &x) {
    fy(t, x);
    writeSlice(derived_state_.y_, y);
}

void Model::getObservableSensitivity(gsl::span<realtype> sy, const realtype t,
                                     const AmiVector &x,
                                     const AmiVectorArray &sx) {
    if (!ny)
        return;

    fdydx(t, x);
    fdydp(t, x);

    derived_state_.sx_.resize(nx_solver * nplist());
    sx.flatten_to_vector(derived_state_.sx_);

    // compute sy = 1.0*dydx*sx + 1.0*sy
    // dydx A[ny,nx_solver] * sx B[nx_solver,nplist] = sy C[ny,nplist]
    //        M  K                 K  N                     M  N
    //        lda                  ldb                      ldc
    amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                BLASTranspose::noTrans, ny, nplist(), nx_solver, 1.0,
                derived_state_.dydx_.data(), ny,
                derived_state_.sx_.data(), nx_solver, 1.0,
                derived_state_.dydp_.data(),
                ny);

    writeSlice(derived_state_.dydp_, sy);

    if (always_check_finite_)
        checkFinite(sy, "sy");
}

void Model::getObservableSigma(gsl::span<realtype> sigmay, const int it,
                               const ExpData *edata) {
    fsigmay(it, edata);
    writeSlice(derived_state_.sigmay_, sigmay);
}

void Model::getObservableSigmaSensitivity(gsl::span<realtype> ssigmay,
                                          const int it, const ExpData *edata) {
    fdsigmaydp(it, edata);
    writeSlice(derived_state_.dsigmaydp_, ssigmay);
}

void Model::addObservableObjective(realtype &Jy, const int it,
                                   const AmiVector &x, const ExpData &edata) {
    fy(edata.getTimepoint(it), x);
    fsigmay(it, &edata);

    std::vector<realtype> nllh(nJ, 0.0);
    for (int iyt = 0; iyt < nytrue; iyt++) {
        if (edata.isSetObservedData(it, iyt)) {
            std::fill(nllh.begin(), nllh.end(), 0.0);
            fJy(nllh.data(), iyt, state_.unscaledParameters.data(),
                state_.fixedParameters.data(),
                derived_state_.y_.data(),
                derived_state_.sigmay_.data(),
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
    derived_state_.sx_.resize(nx_solver * nplist());
    sx.flatten_to_vector(derived_state_.sx_);

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                BLASTranspose::noTrans, nJ, nplist(), nx_solver, 1.0,
                derived_state_.dJydx_.data(), nJ,
                derived_state_.sx_.data(), nx_solver, 1.0,
                derived_state_.dJydp_.data(),
                nJ);

    writeLLHSensitivitySlice(derived_state_.dJydp_, sllh, s2llh);
}

void Model::addPartialObservableObjectiveSensitivity(
    std::vector<realtype> &sllh, std::vector<realtype> &s2llh, const int it,
    const AmiVector &x, const ExpData &edata) {
    if (!ny)
        return;

    fdJydp(it, x, edata);

    writeLLHSensitivitySlice(derived_state_.dJydp_, sllh, s2llh);
}

void Model::getAdjointStateObservableUpdate(gsl::span<realtype> dJydx,
                                            const int it, const AmiVector &x,
                                            const ExpData &edata) {
    fdJydx(it, x, edata);
    writeSlice(derived_state_.dJydx_, dJydx);
}

void Model::getEvent(gsl::span<realtype> z, const int ie, const realtype t,
                     const AmiVector &x) {
    fz(ie, t, x);
    writeSliceEvent(derived_state_.z_, z, ie);
}

void Model::getEventSensitivity(gsl::span<realtype> sz, const int ie,
                                const realtype t, const AmiVector &x,
                                const AmiVectorArray &sx) {
    for (int ip = 0; ip < nplist(); ip++) {
        fsz(&sz.at(ip * nz), ie, t, x.data(), state_.unscaledParameters.data(),
            state_.fixedParameters.data(), state_.h.data(), sx.data(ip),
            plist(ip));
    }
}

void Model::getUnobservedEventSensitivity(gsl::span<realtype> sz,
                                          const int ie) {
    checkBufferSize(sz, nz * nplist());

    for (int iz = 0; iz < nz; ++iz)
        if (z2event_.at(iz) - 1 == ie)
            for (int ip = 0; ip < nplist(); ++ip)
                sz.at(ip * nz + iz) = 0.0;
}

void Model::getEventRegularization(gsl::span<realtype> rz, const int ie,
                                   const realtype t, const AmiVector &x) {
    frz(ie, t, x);
    writeSliceEvent(derived_state_.rz_, rz, ie);
}

void Model::getEventRegularizationSensitivity(gsl::span<realtype> srz,
                                              const int ie, const realtype t,
                                              const AmiVector &x,
                                              const AmiVectorArray &sx) {
    for (int ip = 0; ip < nplist(); ip++) {
        fsrz(&srz.at(ip * nz), ie, t, x.data(), state_.unscaledParameters.data(),
             state_.fixedParameters.data(), state_.h.data(), sx.data(ip),
             plist(ip));
    }
}

void Model::getEventSigma(gsl::span<realtype> sigmaz, const int ie,
                          const int nroots, const realtype t,
                          const ExpData *edata) {
    fsigmaz(ie, nroots, t, edata);
    writeSliceEvent(derived_state_.sigmaz_, sigmaz, ie);
}

void Model::getEventSigmaSensitivity(gsl::span<realtype> ssigmaz, const int ie,
                                     const int nroots, const realtype t,
                                     const ExpData *edata) {
    fdsigmazdp(ie, nroots, t, edata);
    writeSensitivitySliceEvent(derived_state_.dsigmazdp_, ssigmaz, ie);
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
            fJz(nllh.data(), iztrue, state_.unscaledParameters.data(),
                state_.fixedParameters.data(),
                derived_state_.z_.data(), derived_state_.sigmaz_.data(),
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
            fJrz(nllh.data(), iztrue, state_.unscaledParameters.data(),
                 state_.fixedParameters.data(),
                 derived_state_.rz_.data(), derived_state_.sigmaz_.data());
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
    sx.flatten_to_vector(derived_state_.sx_);

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                BLASTranspose::noTrans, nJ, nplist(), nx_solver, 1.0,
                derived_state_.dJzdx_.data(), nJ,
                derived_state_.sx_.data(), nx_solver, 1.0,
                derived_state_.dJzdp_.data(),
                nJ);

    // sJy += multResult + dJydp
    writeLLHSensitivitySlice(derived_state_.dJzdp_, sllh, s2llh);
}

void Model::getAdjointStateEventUpdate(gsl::span<realtype> dJzdx, const int ie,
                                       const int nroots, const realtype t,
                                       const AmiVector &x,
                                       const ExpData &edata) {
    fdJzdx(ie, nroots, t, x, edata);
    writeSlice(derived_state_.dJzdx_, dJzdx);
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

    writeLLHSensitivitySlice(derived_state_.dJzdp_, sllh, s2llh);
}

void Model::getEventTimeSensitivity(std::vector<realtype> &stau,
                                    const realtype t, const int ie,
                                    const AmiVector &x,
                                    const AmiVectorArray &sx) {

    std::fill(stau.begin(), stau.end(), 0.0);

    for (int ip = 0; ip < nplist(); ip++) {
        fstau(&stau.at(ip), t, x.data(), state_.unscaledParameters.data(),
              state_.fixedParameters.data(), state_.h.data(), sx.data(ip),
              plist(ip), ie);
    }
}

void Model::addStateEventUpdate(AmiVector &x, const int ie, const realtype t,
                                const AmiVector &xdot,
                                const AmiVector &xdot_old) {

    derived_state_.deltax_.assign(nx_solver, 0.0);

    // compute update
    fdeltax(derived_state_.deltax_.data(), t, x.data(), state_.unscaledParameters.data(),
            state_.fixedParameters.data(), state_.h.data(), ie, xdot.data(),
            xdot_old.data());

    if (always_check_finite_) {
        app->checkFinite(derived_state_.deltax_, "deltax");
    }

    // update
    amici_daxpy(nx_solver, 1.0, derived_state_.deltax_.data(), 1, x.data(), 1);
}

void Model::addStateSensitivityEventUpdate(AmiVectorArray &sx, const int ie,
                                           const realtype t,
                                           const AmiVector &x_old,
                                           const AmiVector &xdot,
                                           const AmiVector &xdot_old,
                                           const std::vector<realtype> &stau) {
    fw(t, x_old.data());

    for (int ip = 0; ip < nplist(); ip++) {

        derived_state_.deltasx_.assign(nx_solver, 0.0);

        // compute update
        fdeltasx(derived_state_.deltasx_.data(), t, x_old.data(),
                 state_.unscaledParameters.data(), state_.fixedParameters.data(),
                 state_.h.data(), derived_state_.w_.data(), plist(ip), ie,
                 xdot.data(), xdot_old.data(), sx.data(ip), &stau.at(ip));

        if (always_check_finite_) {
            app->checkFinite(derived_state_.deltasx_, "deltasx");
        }

        amici_daxpy(nx_solver, 1.0, derived_state_.deltasx_.data(), 1,
                    sx.data(ip), 1);
    }
}

void Model::addAdjointStateEventUpdate(AmiVector &xB, const int ie,
                                       const realtype t, const AmiVector &x,
                                       const AmiVector &xdot,
                                       const AmiVector &xdot_old) {

    derived_state_.deltaxB_.assign(nx_solver, 0.0);

    // compute update
    fdeltaxB(derived_state_.deltaxB_.data(), t, x.data(),
             state_.unscaledParameters.data(),
             state_.fixedParameters.data(), state_.h.data(), ie, xdot.data(),
             xdot_old.data(), xB.data());

    if (always_check_finite_) {
        app->checkFinite(derived_state_.deltaxB_, "deltaxB");
    }

    // apply update
    for (int ix = 0; ix < nxtrue_solver; ++ix)
        for (int iJ = 0; iJ < nJ; ++iJ)
            xB.at(ix + iJ * nxtrue_solver) +=
                derived_state_.deltaxB_.at(ix + iJ * nxtrue_solver);
}

void Model::addAdjointQuadratureEventUpdate(
    AmiVector xQB, const int ie, const realtype t, const AmiVector &x,
    const AmiVector &xB, const AmiVector &xdot, const AmiVector &xdot_old) {
    for (int ip = 0; ip < nplist(); ip++) {
        derived_state_.deltaqB_.assign(nJ, 0.0);

        fdeltaqB(derived_state_.deltaqB_.data(), t, x.data(), state_.unscaledParameters.data(),
                 state_.fixedParameters.data(), state_.h.data(), plist(ip), ie,
                 xdot.data(), xdot_old.data(), xB.data());

        for (int iJ = 0; iJ < nJ; ++iJ)
            xQB.at(iJ) += derived_state_.deltaqB_.at(iJ);
    }

    if (always_check_finite_) {
        app->checkFinite(derived_state_.deltaqB_, "deltaqB");
    }
}

void Model::updateHeaviside(const std::vector<int> &rootsfound) {
    for (int ie = 0; ie < ne; ie++) {
        state_.h.at(ie) += rootsfound.at(ie);
    }
}

void Model::updateHeavisideB(const int *rootsfound) {
    for (int ie = 0; ie < ne; ie++) {
        state_.h.at(ie) -= rootsfound[ie];
    }
}


int Model::checkFinite(gsl::span<const realtype> array, const char *fun) const {
    auto result = app->checkFinite(array, fun);

    if (result != AMICI_SUCCESS) {
        app->checkFinite(state_.fixedParameters, "k");
        app->checkFinite(state_.unscaledParameters, "p");
        app->checkFinite(derived_state_.w_, "w");
        app->checkFinite(simulation_parameters_.ts_, "t");
    }

    return result;
}

void Model::setAlwaysCheckFinite(bool alwaysCheck) {
    always_check_finite_ = alwaysCheck;
}

bool Model::getAlwaysCheckFinite() const { return always_check_finite_; }

void Model::fx0(AmiVector &x) {
    std::fill(derived_state_.x_rdata_.begin(), derived_state_.x_rdata_.end(), 0.0);
    /* this function  also computes initial total abundances */
    fx0(derived_state_.x_rdata_.data(), simulation_parameters_.tstart_,
        state_.unscaledParameters.data(),
        state_.fixedParameters.data());
    fx_solver(x.data(), derived_state_.x_rdata_.data());
    ftotal_cl(state_.total_cl.data(), derived_state_.x_rdata_.data());

    if (always_check_finite_) {
        checkFinite(derived_state_.x_rdata_, "x0 x_rdata");
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
    fx_rdata(derived_state_.x_rdata_.data(), x.data(), state_.total_cl.data());
    fx0_fixedParameters(derived_state_.x_rdata_.data(),
                        simulation_parameters_.tstart_,
                        state_.unscaledParameters.data(),
                        state_.fixedParameters.data(),
                        simulation_parameters_.reinitialization_state_idxs_sim
                        );
    fx_solver(x.data(), derived_state_.x_rdata_.data());
    /* update total abundances */
    ftotal_cl(state_.total_cl.data(), derived_state_.x_rdata_.data());
}

void Model::fsx0(AmiVectorArray &sx, const AmiVector &x) {
    /* this function  also computes initial total abundance sensitivities */
    realtype *stcl = nullptr;
    for (int ip = 0; ip < nplist(); ip++) {
        if (ncl() > 0)
            stcl = &state_.stotal_cl.at(plist(ip) * ncl());
        std::fill(derived_state_.sx_rdata_.begin(),
                  derived_state_.sx_rdata_.end(), 0.0);
        fsx0(derived_state_.sx_rdata_.data(), simulation_parameters_.tstart_,
             x.data(), state_.unscaledParameters.data(),
             state_.fixedParameters.data(), plist(ip));
        fsx_solver(sx.data(ip), derived_state_.sx_rdata_.data());
        fstotal_cl(stcl, derived_state_.sx_rdata_.data(), plist(ip));
    }
}

void Model::fsx0_fixedParameters(AmiVectorArray &sx, const AmiVector &x) {
    if (!getReinitializeFixedParameterInitialStates())
        return;
    realtype *stcl = nullptr;
    for (int ip = 0; ip < nplist(); ip++) {
        if (ncl() > 0)
            stcl = &state_.stotal_cl.at(plist(ip) * ncl());
        fsx_rdata(derived_state_.sx_rdata_.data(), sx.data(ip), stcl, plist(ip));
        fsx0_fixedParameters(derived_state_.sx_rdata_.data(),
                             simulation_parameters_.tstart_, x.data(),
                             state_.unscaledParameters.data(),
                             state_.fixedParameters.data(),
                             plist(ip),
                             simulation_parameters_.reinitialization_state_idxs_sim);
        fsx_solver(sx.data(ip), derived_state_.sx_rdata_.data());
        fstotal_cl(stcl, derived_state_.sx_rdata_.data(), plist(ip));
    }
}

void Model::fsdx0() {}

void Model::fx_rdata(AmiVector &x_rdata, const AmiVector &x) {
    fx_rdata(x_rdata.data(), x.data(), state_.total_cl.data());
    if (always_check_finite_)
        checkFinite(x_rdata.getVector(), "x_rdata");
}

void Model::fsx_rdata(AmiVectorArray &sx_rdata, const AmiVectorArray &sx) {
    realtype *stcl = nullptr;
    for (int ip = 0; ip < nplist(); ip++) {
        if (ncl() > 0)
            stcl = &state_.stotal_cl.at(plist(ip) * ncl());
        fsx_rdata(sx_rdata.data(ip), sx.data(ip), stcl, ip);
    }
}

void Model::writeSliceEvent(gsl::span<const realtype> slice,
                            gsl::span<realtype> buffer, const int ie) {
    checkBufferSize(buffer, slice.size());
    checkBufferSize(buffer, z2event_.size());
    for (unsigned izt = 0; izt < z2event_.size(); ++izt)
        if (z2event_.at(izt) - 1 == ie)
            buffer.at(izt) = slice.at(izt);
}

void Model::writeSensitivitySliceEvent(gsl::span<const realtype> slice,
                                       gsl::span<realtype> buffer,
                                       const int ie) {
    checkBufferSize(buffer, slice.size());
    checkBufferSize(buffer, z2event_.size() * nplist());
    for (int ip = 0; ip < nplist(); ++ip)
        for (unsigned izt = 0; izt < z2event_.size(); ++izt)
            if (z2event_.at(izt) - 1 == ie)
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

void Model::checkLLHBufferSize(std::vector<realtype> const &sllh,
                               std::vector<realtype> const &s2llh) const {
    if (sllh.size() != static_cast<unsigned>(nplist()))
        throw AmiException("Incorrect sllh buffer size! Was %u, expected %i.",
                           sllh.size(), nplist());

    if (s2llh.size() != static_cast<unsigned>((nJ - 1) * nplist()))
        throw AmiException("Incorrect s2llh buffer size! Was %u, expected %i.",
                           s2llh.size(), (nJ - 1) * nplist());
}

void Model::initializeVectors() {
    sx0data_.clear();
    if (!pythonGenerated)
        derived_state_.dxdotdp = AmiVectorArray(nx_solver, nplist());
}

void Model::fy(const realtype t, const AmiVector &x) {
    if (!ny)
        return;

    derived_state_.y_.assign(ny, 0.0);

    fw(t, x.data());
    fy(derived_state_.y_.data(), t, x.data(), state_.unscaledParameters.data(),
       state_.fixedParameters.data(),
       state_.h.data(), derived_state_.w_.data());

    if (always_check_finite_) {
        app->checkFinite(gsl::make_span(derived_state_.y_.data(), ny), "y");
    }
}

void Model::fdydp(const realtype t, const AmiVector &x) {
    if (!ny)
        return;

    derived_state_.dydp_.assign(ny * nplist(), 0.0);
    fw(t, x.data());
    fdwdp(t, x.data());

    /* get dydp slice (ny) for current time and parameter */
    for (int ip = 0; ip < nplist(); ip++)
        if (pythonGenerated) {
            fdydp(&derived_state_.dydp_.at(ip * ny), t, x.data(),
                  state_.unscaledParameters.data(),
                  state_.fixedParameters.data(), state_.h.data(), plist(ip),
                  derived_state_.w_.data(), state_.stotal_cl.data());
        } else {
            fdydp(&derived_state_.dydp_.at(ip * ny), t, x.data(),
                  state_.unscaledParameters.data(),
                  state_.fixedParameters.data(), state_.h.data(), plist(ip),
                  derived_state_.w_.data(), derived_state_.dwdp_.data());
        }

    if (always_check_finite_) {
        app->checkFinite(derived_state_.dydp_, "dydp");
    }
}

void Model::fdydx(const realtype t, const AmiVector &x) {
    if (!ny)
        return;

    derived_state_.dydx_.assign(ny * nx_solver, 0.0);

    fw(t, x.data());
    fdwdx(t, x.data());
    fdydx(derived_state_.dydx_.data(), t, x.data(), state_.unscaledParameters.data(),
          state_.fixedParameters.data(), state_.h.data(),
          derived_state_.w_.data(), derived_state_.dwdx_.data());

    if (always_check_finite_) {
        app->checkFinite(derived_state_.dydx_, "dydx");
    }
}

void Model::fsigmay(const int it, const ExpData *edata) {
    if (!ny)
        return;

    derived_state_.sigmay_.assign(ny, 0.0);

    fsigmay(derived_state_.sigmay_.data(), getTimepoint(it), state_.unscaledParameters.data(),
            state_.fixedParameters.data());

    if (edata) {
        auto sigmay_edata = edata->getObservedDataStdDevPtr(it);
        /* extract the value for the standard deviation from ExpData,
         * if the data value is NaN, use the parameter value */
        for (int iytrue = 0; iytrue < nytrue; iytrue++) {
            if (edata->isSetObservedDataStdDev(it, iytrue))
                derived_state_.sigmay_.at(iytrue) = sigmay_edata[iytrue];

            /* TODO: when moving second order code to cpp, verify
             * that this is actually what we want
             */
            for (int iJ = 1; iJ < nJ; iJ++)
                derived_state_.sigmay_.at(iytrue + iJ*nytrue) = 0;

            if (edata->isSetObservedData(it, iytrue))
                checkSigmaPositivity(derived_state_.sigmay_.at(iytrue), "sigmay");
        }
    }
}

void Model::fdsigmaydp(const int it, const ExpData *edata) {
    if (!ny)
        return;

    derived_state_.dsigmaydp_.assign(ny * nplist(), 0.0);

    for (int ip = 0; ip < nplist(); ip++)
        // get dsigmaydp slice (ny) for current timepoint and parameter
        fdsigmaydp(&derived_state_.dsigmaydp_.at(ip * ny), getTimepoint(it),
                   state_.unscaledParameters.data(),
                   state_.fixedParameters.data(),
                   plist(ip));

    // sigmas in edata override model-sigma -> for those sigmas, set dsigmaydp
    // to zero
    if (edata) {
        for (int iy = 0; iy < nytrue; iy++) {
            if (!edata->isSetObservedDataStdDev(it, iy))
                continue;
            for (int ip = 0; ip < nplist(); ip++) {
                derived_state_.dsigmaydp_.at(ip * ny + iy) = 0.0;
            }
        }
    }

    if (always_check_finite_) {
        app->checkFinite(derived_state_.dsigmaydp_, "dsigmaydp");
    }
}

void Model::fdJydy(const int it, const AmiVector &x, const ExpData &edata) {
    if (!ny)
        return;

    fy(edata.getTimepoint(it), x);
    fsigmay(it, &edata);

    if (pythonGenerated) {
        for (int iyt = 0; iyt < nytrue; iyt++) {
            if (!derived_state_.dJydy_.at(iyt).capacity())
                continue;
            derived_state_.dJydy_.at(iyt).zero();
            fdJydy_colptrs(derived_state_.dJydy_.at(iyt), iyt);
            fdJydy_rowvals(derived_state_.dJydy_.at(iyt), iyt);

            if (!edata.isSetObservedData(it, iyt))
                continue;

            // get dJydy slice (ny) for current timepoint and observable
            fdJydy(derived_state_.dJydy_.at(iyt).data(), iyt, state_.unscaledParameters.data(),
                   state_.fixedParameters.data(), derived_state_.y_.data(),
                   derived_state_.sigmay_.data(),
                   edata.getObservedDataPtr(it));

            if (always_check_finite_) {
                app->checkFinite(
                            gsl::make_span(derived_state_.dJydy_.at(iyt).get()),
                            "dJydy");
            }
        }
    } else {
        std::fill(derived_state_.dJydy_matlab_.begin(),
                  derived_state_.dJydy_matlab_.end(), 0.0);
        for (int iyt = 0; iyt < nytrue; iyt++) {
            if (!edata.isSetObservedData(it, iyt))
                continue;
            fdJydy(&derived_state_.dJydy_matlab_.at(iyt * ny * nJ), iyt,
                   state_.unscaledParameters.data(),
                   state_.fixedParameters.data(), derived_state_.y_.data(),
                   derived_state_.sigmay_.data(),
                   edata.getObservedDataPtr(it));
        }
        if (always_check_finite_) {
            // get dJydy slice (ny) for current timepoint and observable
            app->checkFinite(derived_state_.dJydy_matlab_, "dJydy");
        }
    }
}

void Model::fdJydsigma(const int it, const AmiVector &x, const ExpData &edata) {
    if (!ny)
        return;

    derived_state_.dJydsigma_.assign(nytrue * ny * nJ, 0.0);

    fy(edata.getTimepoint(it), x);
    fsigmay(it, &edata);

    for (int iyt = 0; iyt < nytrue; iyt++) {
        if (edata.isSetObservedData(it, iyt))
            // get dJydsigma slice (ny) for current timepoint and observable
            fdJydsigma(&derived_state_.dJydsigma_.at(iyt * ny * nJ), iyt,
                       state_.unscaledParameters.data(),
                       state_.fixedParameters.data(), derived_state_.y_.data(),
                       derived_state_.sigmay_.data(),
                       edata.getObservedDataPtr(it));
    }

    if (always_check_finite_) {
        app->checkFinite(derived_state_.dJydsigma_, "dJydsigma");
    }
}

void Model::fdJydp(const int it, const AmiVector &x, const ExpData &edata) {
    // dJydy         nJ, nytrue x ny
    // dydp          nplist * ny
    // dJydp         nplist x nJ
    // dJydsigma
    if (!ny)
        return;

    derived_state_.dJydp_.assign(nJ * nplist(), 0.0);

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
                derived_state_.dJydy_.at(iyt).multiply(
                    gsl::span<realtype>(&derived_state_.dJydp_.at(iplist * nJ), nJ),
                    gsl::span<const realtype>(&derived_state_.dydp_.at(iplist * ny), ny));
            }
        } else {
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nplist(), ny, 1.0,
                        &derived_state_.dJydy_matlab_.at(iyt * nJ * ny), nJ,
                        derived_state_.dydp_.data(), ny,
                        1.0, derived_state_.dJydp_.data(), nJ);
        }
        // dJydp = 1.0 * dJydp +  1.0 * dJydsigma * dsigmaydp
        amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                    BLASTranspose::noTrans, nJ, nplist(), ny, 1.0,
                    &derived_state_.dJydsigma_.at(iyt * nJ * ny), nJ,
                    derived_state_.dsigmaydp_.data(), ny, 1.0,
                    derived_state_.dJydp_.data(), nJ);
    }
}

void Model::fdJydx(const int it, const AmiVector &x, const ExpData &edata) {
    if (!ny)
        return;

    derived_state_.dJydx_.assign(nJ * nx_solver, 0.0);

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
                derived_state_.dJydy_.at(iyt).multiply(
                    gsl::span<realtype>(&derived_state_.dJydx_.at(ix * nJ), nJ),
                    gsl::span<const realtype>(&derived_state_.dydx_.at(ix * ny), ny));
            }
        } else {
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nx_solver, ny, 1.0,
                        &derived_state_.dJydy_matlab_.at(iyt * ny * nJ), nJ,
                        derived_state_.dydx_.data(), ny,
                        1.0, derived_state_.dJydx_.data(), nJ);
        }
    }

    if (always_check_finite_) {
        app->checkFinite(derived_state_.dJydx_, "dJydx");
    }
}

void Model::fz(const int ie, const realtype t, const AmiVector &x) {

    derived_state_.z_.assign(nz, 0.0);

    fz(derived_state_.z_.data(), ie, t, x.data(), state_.unscaledParameters.data(),
       state_.fixedParameters.data(), state_.h.data());
}

void Model::fdzdp(const int ie, const realtype t, const AmiVector &x) {
    if (!nz)
        return;

    derived_state_.dzdp_.assign(nz * nplist(), 0.0);

    for (int ip = 0; ip < nplist(); ip++) {
        fdzdp(derived_state_.dzdp_.data(), ie, t, x.data(),
              state_.unscaledParameters.data(),
              state_.fixedParameters.data(), state_.h.data(), plist(ip));
    }

    if (always_check_finite_) {
        app->checkFinite(derived_state_.dzdp_, "dzdp");
    }
}

void Model::fdzdx(const int ie, const realtype t, const AmiVector &x) {
    if (!nz)
        return;

    derived_state_.dzdx_.assign(nz * nx_solver, 0.0);

    fdzdx(derived_state_.dzdx_.data(), ie, t, x.data(), state_.unscaledParameters.data(),
          state_.fixedParameters.data(), state_.h.data());

    if (always_check_finite_) {
        app->checkFinite(derived_state_.dzdx_, "dzdx");
    }
}

void Model::frz(const int ie, const realtype t, const AmiVector &x) {

    derived_state_.rz_.assign(nz, 0.0);

    frz(derived_state_.rz_.data(), ie, t, x.data(),
        state_.unscaledParameters.data(),
        state_.fixedParameters.data(), state_.h.data());
}

void Model::fdrzdp(const int ie, const realtype t, const AmiVector &x) {
    if (!nz)
        return;

    derived_state_.drzdp_.assign(nz * nplist(), 0.0);

    for (int ip = 0; ip < nplist(); ip++) {
        fdrzdp(derived_state_.drzdp_.data(), ie, t, x.data(),
               state_.unscaledParameters.data(),
               state_.fixedParameters.data(), state_.h.data(), plist(ip));
    }

    if (always_check_finite_) {
        app->checkFinite(derived_state_.drzdp_, "drzdp");
    }
}

void Model::fdrzdx(const int ie, const realtype t, const AmiVector &x) {
    if (!nz)
        return;

    derived_state_.drzdx_.assign(nz * nx_solver, 0.0);

    fdrzdx(derived_state_.drzdx_.data(), ie, t, x.data(), state_.unscaledParameters.data(),
           state_.fixedParameters.data(), state_.h.data());

    if (always_check_finite_) {
        app->checkFinite(derived_state_.drzdx_, "drzdx");
    }
}

void Model::fsigmaz(const int ie, const int nroots, const realtype t,
                    const ExpData *edata) {
    if (!nz)
        return;

    derived_state_.sigmaz_.assign(nz, 0.0);
    fsigmaz(derived_state_.sigmaz_.data(), t, state_.unscaledParameters.data(),
            state_.fixedParameters.data());

    if (edata) {
        for (int iztrue = 0; iztrue < nztrue; iztrue++) {
            if (z2event_.at(iztrue) - 1 == ie) {
                if (edata->isSetObservedEventsStdDev(nroots, iztrue)) {
                    auto sigmaz_edata =
                        edata->getObservedEventsStdDevPtr(nroots);
                    derived_state_.sigmaz_.at(iztrue) = sigmaz_edata[iztrue];
                }

                /* TODO: when moving second order code to cpp, verify
                 * that this is actually what we want
                 */
                for (int iJ = 1; iJ < nJ; iJ++)
                    derived_state_.sigmaz_.at(iztrue + iJ*nztrue) = 0;

                if (edata->isSetObservedEvents(nroots, iztrue))
                    checkSigmaPositivity(derived_state_.sigmaz_.at(iztrue),
                                         "sigmaz");
            }
        }
    }
}

void Model::fdsigmazdp(const int ie, const int nroots, const realtype t,
                       const ExpData *edata) {
    if (!nz)
        return;

    derived_state_.dsigmazdp_.assign(nz * nplist(), 0.0);

    for (int ip = 0; ip < nplist(); ip++) {
        // get dsigmazdp slice (nz) for current event and parameter
        fdsigmazdp(&derived_state_.dsigmazdp_.at(ip * nz), t, state_.unscaledParameters.data(),
                   state_.fixedParameters.data(), plist(ip));
    }

    // sigmas in edata override model-sigma -> for those sigmas, set dsigmazdp
    // to zero
    if (edata) {
        for (int iz = 0; iz < nztrue; iz++) {
            if (z2event_.at(iz) - 1 == ie &&
                !edata->isSetObservedEventsStdDev(nroots, iz)) {
                for (int ip = 0; ip < nplist(); ip++)
                    derived_state_.dsigmazdp_.at(iz + nz * ip) = 0;
            }
        }
    }

    if (always_check_finite_) {
        app->checkFinite(derived_state_.dsigmazdp_, "dsigmazdp");
    }
}

void Model::fdJzdz(const int ie, const int nroots, const realtype t,
                   const AmiVector &x, const ExpData &edata) {
    if (!nz)
        return;

    derived_state_.dJzdz_.assign(nztrue * nz * nJ, 0.0);

    fz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.isSetObservedEvents(nroots, iztrue)) {
            fdJzdz(&derived_state_.dJzdz_.at(iztrue * nz * nJ), iztrue,
                   state_.unscaledParameters.data(),
                   state_.fixedParameters.data(),
                   derived_state_.z_.data(), derived_state_.sigmaz_.data(),
                   edata.getObservedEventsPtr(nroots));
        }
    }

    if (always_check_finite_) {
        app->checkFinite(derived_state_.dJzdz_, "dJzdz");
    }
}

void Model::fdJzdsigma(const int ie, const int nroots, const realtype t,
                       const AmiVector &x, const ExpData &edata) {
    if (!nz)
        return;

    derived_state_.dJzdsigma_.assign(nztrue * nz * nJ, 0.0);

    fz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.isSetObservedEvents(nroots, iztrue)) {
            fdJzdsigma(&derived_state_.dJzdsigma_.at(iztrue * nz * nJ), iztrue,
                       state_.unscaledParameters.data(),
                       state_.fixedParameters.data(), derived_state_.z_.data(),
                       derived_state_.sigmaz_.data(),
                       edata.getObservedEventsPtr(nroots));
        }
    }

    if (always_check_finite_) {
        app->checkFinite(derived_state_.dJzdsigma_, "dJzdsigma");
    }
}

void Model::fdJzdp(const int ie, const int nroots, realtype t,
                   const AmiVector &x, const ExpData &edata) {
    if (!nz)
        return;
    // dJzdz         nJ x nz x nztrue
    // dJzdsigma     nJ x nz x nztrue
    // dzdp          nz x nplist()
    // dJzdp         nJ x nplist()

    derived_state_.dJzdp_.assign(nJ * nplist(), 0.0);

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
                        &derived_state_.dJzdz_.at(izt * nz * nJ), nJ,
                        derived_state_.dzdp_.data(), nz, 1.0,
                        derived_state_.dJzdp_.data(), nJ);
        } else {
            // with rz
            fdJrzdz(ie, nroots, t, x, edata);
            fdJrzdsigma(ie, nroots, t, x, edata);

            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nplist(), nz, 1.0,
                        &derived_state_.dJrzdsigma_.at(izt * nz * nJ), nJ,
                        derived_state_.dsigmazdp_.data(), nz,
                        1.0, derived_state_.dJzdp_.data(), nJ);

            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nplist(), nz, 1.0,
                        &derived_state_.dJrzdz_.at(izt * nz * nJ), nJ,
                        derived_state_.dzdp_.data(), nz, 1.0,
                        derived_state_.dJzdp_.data(), nJ);
        }

        amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                    BLASTranspose::noTrans, nJ, nplist(), nz, 1.0,
                    &derived_state_.dJzdsigma_.at(izt * nz * nJ), nJ,
                    derived_state_.dsigmazdp_.data(), nz, 1.0,
                    derived_state_.dJzdp_.data(), nJ);
    }
}

void Model::fdJzdx(const int ie, const int nroots, const realtype t,
                   const AmiVector &x, const ExpData &edata) {
    // dJzdz         nJ x nz        x nztrue
    // dzdx          nz x nx_solver
    // dJzdx         nJ x nx_solver x nmaxevent
    if(!nz)
        return;

    derived_state_.dJzdx_.assign(nJ * nx_solver, 0.0);

    fdJzdz(ie, nroots, t, x, edata);

    for (int izt = 0; izt < nztrue; ++izt) {
        if (!edata.isSetObservedEvents(nroots, izt))
            continue;

        if (t < edata.getTimepoint(edata.nt() - 1)) {
            // z
            fdzdx(ie, t, x);
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nx_solver, nz, 1.0,
                        &derived_state_.dJzdz_.at(izt * nz * nJ), nJ,
                        derived_state_.dzdx_.data(), nz, 1.0,
                        derived_state_.dJzdx_.data(), nJ);
        } else {
            // rz
            fdJrzdz(ie, nroots, t, x, edata);
            fdrzdx(ie, t, x);
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nx_solver, nz, 1.0,
                        &derived_state_.dJrzdz_.at(izt * nz * nJ), nJ,
                        derived_state_.drzdx_.data(), nz, 1.0,
                        derived_state_.dJzdx_.data(), nJ);
        }
    }
}

void Model::fdJrzdz(const int ie, const int nroots, const realtype t,
                    const AmiVector &x, const ExpData &edata) {
    if (!nz)
        return;

    derived_state_.dJrzdz_.assign(nztrue * nz * nJ, 0.0);

    frz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.isSetObservedEvents(nroots, iztrue)) {
            fdJrzdz(&derived_state_.dJrzdz_.at(iztrue * nz * nJ), iztrue,
                    state_.unscaledParameters.data(),
                    state_.fixedParameters.data(), derived_state_.rz_.data(),
                    derived_state_.sigmaz_.data());
        }
    }

    if (always_check_finite_) {
        app->checkFinite(derived_state_.dJrzdz_, "dJrzdz");
    }
}

void Model::fdJrzdsigma(const int ie, const int nroots, const realtype t,
                        const AmiVector &x, const ExpData &edata) {
    if (!nz)
        return;

    derived_state_.dJrzdsigma_.assign(nztrue * nz * nJ, 0.0);

    frz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.isSetObservedEvents(nroots, iztrue)) {
            fdJrzdsigma(&derived_state_.dJrzdsigma_.at(iztrue * nz * nJ), iztrue,
                        state_.unscaledParameters.data(),
                        state_.fixedParameters.data(), derived_state_.rz_.data(),
                        derived_state_.sigmaz_.data());
        }
    }

    if (always_check_finite_) {
        app->checkFinite(derived_state_.dJrzdsigma_, "dJrzdsigma");
    }
}

void Model::fw(const realtype t, const realtype *x) {
    std::fill(derived_state_.w_.begin(), derived_state_.w_.end(), 0.0);
    fw(derived_state_.w_.data(), t, x, state_.unscaledParameters.data(),
       state_.fixedParameters.data(), state_.h.data(), state_.total_cl.data());

    if (always_check_finite_) {
        app->checkFinite(derived_state_.w_, "w");
    }
}

void Model::fdwdp(const realtype t, const realtype *x) {
    if (!nw)
        return;

    fw(t, x);
    derived_state_.dwdp_.zero();
    if (pythonGenerated) {
        if (!dwdp_hierarchical_.at(0).capacity())
            return;
        fdwdw(t,x);
        dwdp_hierarchical_.at(0).zero();
        fdwdp_colptrs(dwdp_hierarchical_.at(0));
        fdwdp_rowvals(dwdp_hierarchical_.at(0));
        fdwdp(dwdp_hierarchical_.at(0).data(), t, x,
              state_.unscaledParameters.data(), state_.fixedParameters.data(),
              state_.h.data(), derived_state_.w_.data(), state_.total_cl.data(),
              state_.stotal_cl.data());

        for (int irecursion = 1; irecursion <= w_recursion_depth_;
             irecursion++) {
            dwdw_.sparse_multiply(dwdp_hierarchical_.at(irecursion),
                                  dwdp_hierarchical_.at(irecursion - 1));
        }
        derived_state_.dwdp_.sparse_sum(dwdp_hierarchical_);

    } else {
        if (!derived_state_.dwdp_.capacity())
            return;
        // matlab generated
        fdwdp(derived_state_.dwdp_.data(), t, x,
              state_.unscaledParameters.data(), state_.fixedParameters.data(),
              state_.h.data(), derived_state_.w_.data(),
              state_.total_cl.data(), state_.stotal_cl.data());
    }

    if (always_check_finite_) {
        app->checkFinite(gsl::make_span(derived_state_.dwdp_.get()), "dwdp");
    }
}

void Model::fdwdx(const realtype t, const realtype *x) {
    if (!nw)
        return;

    fw(t, x);

    derived_state_.dwdx_.zero();
    if (pythonGenerated) {
        if (!dwdx_hierarchical_.at(0).capacity())
                return;
        fdwdw(t,x);
        dwdx_hierarchical_.at(0).zero();
        fdwdx_colptrs(dwdx_hierarchical_.at(0));
        fdwdx_rowvals(dwdx_hierarchical_.at(0));
        fdwdx(dwdx_hierarchical_.at(0).data(), t, x,
              state_.unscaledParameters.data(), state_.fixedParameters.data(),
              state_.h.data(), derived_state_.w_.data(), state_.total_cl.data());

        for (int irecursion = 1; irecursion <= w_recursion_depth_;
             irecursion++) {
            dwdw_.sparse_multiply(dwdx_hierarchical_.at(irecursion),
                                  dwdx_hierarchical_.at(irecursion - 1));
        }
        derived_state_.dwdx_.sparse_sum(dwdx_hierarchical_);

    } else {
        if (!derived_state_.dwdx_.capacity())
            return;
        derived_state_.dwdx_.zero();
        fdwdx(derived_state_.dwdx_.data(), t, x,
              state_.unscaledParameters.data(),
              state_.fixedParameters.data(), state_.h.data(),
              derived_state_.w_.data(),
              state_.total_cl.data());
    }

    if (always_check_finite_) {
        app->checkFinite(gsl::make_span(derived_state_.dwdx_.get()), "dwdx");
    }
}

void Model::fdwdw(const realtype t, const realtype *x) {
    if (!nw || !dwdw_.capacity())
        return;
    dwdw_.zero();
    fdwdw_colptrs(dwdw_);
    fdwdw_rowvals(dwdw_);
    fdwdw(dwdw_.data(), t, x, state_.unscaledParameters.data(),
          state_.fixedParameters.data(), state_.h.data(),
          derived_state_.w_.data(), state_.total_cl.data());

    if (always_check_finite_) {
        app->checkFinite(gsl::make_span(dwdw_.get()), "dwdw");
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

const_N_Vector Model::computeX_pos(const_N_Vector x) {
    if (any_state_non_negative_) {
        for (int ix = 0; ix < derived_state_.x_pos_tmp_.getLength(); ++ix) {
            derived_state_.x_pos_tmp_.at(ix) =
                (state_is_non_negative_.at(ix) && NV_Ith_S(x, ix) < 0)
                    ? 0
                    : NV_Ith_S(x, ix);
        }
        return derived_state_.x_pos_tmp_.getNVector();
    }

    return x;
}

void Model::setReinitializationStateIdxs(std::vector<int> const& idxs)
{
    for(auto idx: idxs) {
        if (idx < 0 || idx >= nx_rdata)
            throw AmiException("Invalid state index given: %d", idx);
    }

    simulation_parameters_.reinitialization_state_idxs_sim = idxs;
}

const std::vector<int> &Model::getReinitializationStateIdxs() const
{
    return simulation_parameters_.reinitialization_state_idxs_sim;
}

const AmiVectorArray &Model::get_dxdotdp() const{
    assert(!pythonGenerated);
    return derived_state_.dxdotdp;
}

const SUNMatrixWrapper &Model::get_dxdotdp_full() const{
    assert(pythonGenerated);
    return derived_state_.dxdotdp_full;
}

} // namespace amici
