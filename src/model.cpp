#include <amici/amici.h>
#include <amici/cblas.h>
#include <amici/exception.h>
#include <amici/misc.h>
#include <amici/model.h>
#include <amici/symbolic_functions.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <numeric>
#include <regex>
#include <sstream>
#include <utility>

namespace amici {

/**
 * @brief Maps ModelQuantity items to their string value
 */
std::map<ModelQuantity, std::string> const model_quantity_to_str{
    {ModelQuantity::J, "J"},
    {ModelQuantity::JB, "JB"},
    {ModelQuantity::Jv, "Jv"},
    {ModelQuantity::JvB, "JvB"},
    {ModelQuantity::JDiag, "JDiag"},
    {ModelQuantity::sx, "sx"},
    {ModelQuantity::sy, "sy"},
    {ModelQuantity::sz, "sz"},
    {ModelQuantity::srz, "srz"},
    {ModelQuantity::ssigmay, "ssigmay"},
    {ModelQuantity::ssigmaz, "ssigmaz"},
    {ModelQuantity::xdot, "xdot"},
    {ModelQuantity::sxdot, "sxdot"},
    {ModelQuantity::xBdot, "xBdot"},
    {ModelQuantity::x0, "x0"},
    {ModelQuantity::x0_rdata, "x0_rdata"},
    {ModelQuantity::x, "x"},
    {ModelQuantity::x_rdata, "x_rdata"},
    {ModelQuantity::dwdw, "dwdw"},
    {ModelQuantity::dwdx, "dwdx"},
    {ModelQuantity::dwdp, "dwdp"},
    {ModelQuantity::y, "y"},
    {ModelQuantity::dydp, "dydp"},
    {ModelQuantity::dydx, "dydx"},
    {ModelQuantity::w, "w"},
    {ModelQuantity::root, "root"},
    {ModelQuantity::qBdot, "qBdot"},
    {ModelQuantity::qBdot_ss, "qBdot_ss"},
    {ModelQuantity::xBdot_ss, "xBdot_ss"},
    {ModelQuantity::JSparseB_ss, "JSparseB_ss"},
    {ModelQuantity::deltax, "deltax"},
    {ModelQuantity::deltasx, "deltasx"},
    {ModelQuantity::deltaxB, "deltaxB"},
    {ModelQuantity::k, "k"},
    {ModelQuantity::p, "p"},
    {ModelQuantity::ts, "ts"},
    {ModelQuantity::dJydy, "dJydy"},
    {ModelQuantity::deltaqB, "deltaqB"},
    {ModelQuantity::dsigmaydp, "dsigmaydp"},
    {ModelQuantity::dsigmaydy, "dsigmaydy"},
    {ModelQuantity::dsigmazdp, "dsigmazdp"},
    {ModelQuantity::dJydsigma, "dJydsigma"},
    {ModelQuantity::dJydx, "dJydx"},
    {ModelQuantity::dJrzdx, "dJrzdx"},
    {ModelQuantity::dJzdx, "dJzdx"},
    {ModelQuantity::dzdp, "dzdp"},
    {ModelQuantity::dzdx, "dzdx"},
    {ModelQuantity::dJrzdsigma, "dJrzdsigma"},
    {ModelQuantity::dJrzdz, "dJrzdz"},
    {ModelQuantity::dJzdsigma, "dJzdsigma"},
    {ModelQuantity::dJzdz, "dJzdz"},
    {ModelQuantity::drzdp, "drzdp"},
    {ModelQuantity::drzdx, "drzdx"},

};

static void set_nan_to_zero(std::vector<realtype>& vec) {
    std::ranges::for_each(vec, [](double& val) {
        if (std::isnan(val)) {
            val = 0.0;
        }
    });
}

/**
 * @brief local helper function to get parameters
 * @param ids vector of name/ids of (fixed)Parameters
 * @param values values of the (fixed)Parameters
 * @param id name/id to look for in the vector
 * @param variable_name string indicating what variable we are looking at
 * @param id_name string indicating whether name or id was specified
 * @return value of the selected parameter
 */
static realtype get_value_by_id(
    std::vector<std::string> const& ids, std::vector<realtype> const& values,
    std::string const& id, char const* variable_name, char const* id_name
) {
    auto it = std::ranges::find(ids, id);
    if (it != ids.end())
        return values.at(it - ids.begin());

    throw AmiException(
        "Could not find %s with specified %s", variable_name, id_name
    );
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
static void set_value_by_id(
    std::vector<std::string> const& ids, std::vector<realtype>& values,
    realtype const value, std::string const& id, char const* variable_name,
    char const* id_name
) {
    auto it = std::ranges::find(ids, id);
    if (it != ids.end())
        values.at(it - ids.begin()) = value;
    else
        throw AmiException(
            "Could not find %s with specified %s", variable_name, id_name
        );
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

static int set_value_by_id_regex(
    std::vector<std::string> const& ids, std::vector<realtype>& values,
    realtype value, std::string const& regex, char const* variable_name,
    char const* id_name
) {
    try {
        std::regex pattern(regex);
        int n_found = 0;
        for (auto const& id : ids) {
            if (std::regex_match(id, pattern)) {
                values.at(&id - &ids[0]) = value;
                ++n_found;
            }
        }

        if (n_found == 0)
            throw AmiException(
                "Could not find %s with specified %s (%s)", variable_name,
                id_name, regex.c_str()
            );

        return n_found;
    } catch (std::regex_error const& e) {
        auto err_string = regex_error_to_string(e.code());
        throw AmiException(
            "Specified regex pattern %s could not be compiled:"
            " %s (%s)",
            regex.c_str(), e.what(), err_string.c_str()
        );
    }
}

Model::Model(
    ModelDimensions const& model_dimensions,
    SimulationParameters simulation_parameters, SecondOrderMode const o2mode,
    std::vector<realtype> idlist, std::vector<int> z2event,
    std::vector<Event> events
)
    : ModelDimensions(model_dimensions)
    , state_(*this)
    , derived_state_(*this)
    , z2event_(std::move(z2event))
    , state_is_non_negative_(nx_solver, false)
    , o2_mode_(o2mode)
    , id_list_(std::move(idlist))
    , simulation_parameters_(std::move(simulation_parameters))
    , events_(std::move(events)) {
    model_dimensions.validate();
    Expects(
        model_dimensions.np
        == gsl::narrow<int>(simulation_parameters_.free_parameters.size())
    );
    Expects(
        model_dimensions.nk
        == gsl::narrow<int>(simulation_parameters_.fixed_parameters.size())
    );

    Expects((events_.size() == (unsigned long)ne));

    simulation_parameters_.pscale = std::vector<ParameterScaling>(
        model_dimensions.np, ParameterScaling::none
    );

    unscale_parameters(
        simulation_parameters_.free_parameters, simulation_parameters_.pscale,
        state_.unscaled_parameters
    );
    state_.fixed_parameters = simulation_parameters_.fixed_parameters;
    state_.plist = simulation_parameters_.plist;

    require_sensitivities_for_all_parameters();
}

bool operator==(Model const& a, Model const& b) {
    if (typeid(a) != typeid(b))
        return false;

    return (static_cast<ModelDimensions const&>(a)
            == static_cast<ModelDimensions const&>(b))
           && (a.o2_mode_ == b.o2_mode_) && (a.z2event_ == b.z2event_)
           && (a.id_list_ == b.id_list_)
           && (a.simulation_parameters_ == b.simulation_parameters_)
           && (a.x0data_ == b.x0data_) && (a.sx0data_ == b.sx0data_)
           && (a.nmaxevent_ == b.nmaxevent_)
           && (a.state_is_non_negative_ == b.state_is_non_negative_)
           && (a.sigma_res_ == b.sigma_res_) && (a.min_sigma_ == b.min_sigma_)
           && (a.state_ == b.state_)
           && (a.steadystate_mask_ == b.steadystate_mask_);
}

bool operator==(ModelDimensions const& a, ModelDimensions const& b) {
    if (typeid(a) != typeid(b))
        return false;
    return (a.nx_rdata == b.nx_rdata) && (a.nxtrue_rdata == b.nxtrue_rdata)
           && (a.nx_solver == b.nx_solver)
           && (a.nxtrue_solver == b.nxtrue_solver)
           && (a.nx_solver_reinit == b.nx_solver_reinit) && (a.np == b.np)
           && (a.nk == b.nk) && (a.ny == b.ny) && (a.nytrue == b.nytrue)
           && (a.nz == b.nz) && (a.nztrue == b.nztrue) && (a.ne == b.ne)
           && (a.ne_solver == b.ne_solver) && (a.nspl == b.nspl)
           && (a.nw == b.nw) && (a.ndwdx == b.ndwdx) && (a.ndwdp == b.ndwdp)
           && (a.ndwdw == b.ndwdw) && (a.ndxdotdw == b.ndxdotdw)
           && (a.ndJydy == b.ndJydy) && (a.nnz == b.nnz) && (a.nJ == b.nJ)
           && (a.ubw == b.ubw) && (a.lbw == b.lbw);
}

void Model::initialize(
    realtype t, AmiVector& x, AmiVector& dx, AmiVectorArray& sx,
    AmiVectorArray& /*sdx*/, bool const computeSensitivities,
    std::vector<int>& roots_found
) {
    initialize_state(t, x);
    initialize_splines();
    if (computeSensitivities) {
        initialize_state_sensitivities(t, sx, x);
        initialize_spline_sensitivities();
    }

    fdx0(x, dx);
    if (computeSensitivities)
        fsdx0();

    if (ne)
        initialize_events(t, x, dx, roots_found);

    // evaluate static expressions once
    auto x_pos = compute_x_pos(x);
    fw(t, x_pos, true);
    fdwdw(t, x_pos, true);
    fdwdx(t, x_pos, true);
    if (computeSensitivities) {
        fdwdp(t, x_pos, true);
    }
}

void Model::reinitialize(
    realtype const t, AmiVector& x, AmiVectorArray& sx,
    bool const computeSensitivities
) {
    fx0_fixedParameters(t, x);

    // re-evaluate static expressions once
    auto x_pos = compute_x_pos(x);
    fw(t, x_pos, true);
    fdwdw(t, x_pos, true);
    fdwdx(t, x_pos, true);
    if (computeSensitivities) {
        fsx0_fixedParameters(t, sx, x);
        fdwdp(t, x_pos, true);
    }
}

void Model::initialize_b(
    AmiVector& xB, AmiVector& dxB, AmiVector& xQB, bool posteq
) const {
    xB.zero();
    dxB.zero();
    if (!posteq)
        xQB.zero();
}

void Model::initialize_state(realtype t, AmiVector& x) {
    if (x0data_.empty()) {
        fx0(t, x);
    } else {
        std::vector<realtype> x0_solver(nx_solver, 0.0);
        ftotal_cl(
            state_.total_cl.data(), x0data_.data(),
            state_.unscaled_parameters.data(), state_.fixed_parameters.data()
        );
        fx_solver(x0_solver.data(), x0data_.data());
        std::copy(x0_solver.cbegin(), x0_solver.cend(), x.data());
    }

    check_finite(x.get_vector(), ModelQuantity::x0, t);
}

void Model::initialize_splines() {
    splines_ = fcreate_splines(
        state_.unscaled_parameters.data(), state_.fixed_parameters.data()
    );
    derived_state_.spl_.resize(splines_.size(), 0.0);
    for (auto& spline : splines_) {
        spline.compute_coefficients();
    }
}

void Model::initialize_spline_sensitivities() {
    derived_state_.sspl_
        = SUNMatrixWrapper(splines_.size(), np(), derived_state_.sunctx_);
    int allnodes = 0;
    for (auto const& spline : splines_) {
        allnodes += spline.n_nodes();
    }

    std::vector<realtype> dspline_valuesdp(allnodes * nplist(), 0.0);
    std::vector<realtype> dspline_slopesdp(allnodes * nplist(), 0.0);
    std::vector<realtype> tmp_dvalues(allnodes, 0.0);
    std::vector<realtype> tmp_dslopes(allnodes, 0.0);
    for (int ip = 0; ip < nplist(); ip++) {
        std::ranges::fill(tmp_dvalues, 0.0);
        std::ranges::fill(tmp_dslopes, 0.0);
        fdspline_valuesdp(
            tmp_dvalues.data(), state_.unscaled_parameters.data(),
            state_.fixed_parameters.data(), plist(ip)
        );
        fdspline_slopesdp(
            tmp_dslopes.data(), state_.unscaled_parameters.data(),
            state_.fixed_parameters.data(), plist(ip)
        );
        /* NB dspline_valuesdp/dspline_slopesdp must be filled
         * using the following order for the indices
         * (from slower to faster): spline, node, parameter.
         * That is what the current spline implementation expects.
         */
        int k = 0;
        int offset = ip;
        for (auto const& spline : splines_) {
            for (int n = 0; n < spline.n_nodes(); n++) {
                dspline_valuesdp[offset] = tmp_dvalues[k];
                dspline_slopesdp[offset] = tmp_dslopes[k];
                offset += nplist();
                k += 1;
            }
        }
        assert(k == allnodes);
    }

    int spline_offset = 0;
    for (auto& spline : splines_) {
        spline.compute_coefficients_sensi(
            nplist(), spline_offset, dspline_valuesdp, dspline_slopesdp
        );
        spline_offset += spline.n_nodes() * nplist();
    }
}

void Model::initialize_state_sensitivities(
    realtype t, AmiVectorArray& sx, AmiVector const& x
) {
    if (sx0data_.empty()) {
        fsx0(t, sx, x);
    } else {
        realtype* stcl = nullptr;
        std::vector<realtype> sx0_solver_slice(nx_solver, 0.0);
        for (int ip = 0; ip < nplist(); ip++) {
            if (ncl() > 0)
                stcl = &state_.stotal_cl.at(plist(ip) * ncl());
            fstotal_cl(
                stcl, &sx0data_.at(ip * nx_rdata), plist(ip),
                derived_state_.x_rdata_.data(),
                state_.unscaled_parameters.data(),
                state_.fixed_parameters.data(), state_.total_cl.data()
            );
            fsx_solver(sx0_solver_slice.data(), &sx0data_.at(ip * nx_rdata));
            for (int ix = 0; ix < nx_solver; ix++) {
                sx.at(ix, ip) = sx0_solver_slice.at(ix);
            }
        }
    }
}

void Model::initialize_events(
    realtype t, AmiVector const& x, AmiVector const& dx,
    std::vector<int>& roots_found
) {
    // construct "old" Heaviside vector from initial values
    std::vector<realtype> h_old(ne, 1.0);
    for (int ie = 0; ie < ne; ie++) {
        if (events_.at(ie).get_initial_value() == false)
            h_old.at(ie) = 0.0;
    }

    reinit_events(t, x, dx, h_old, roots_found);
}

void Model::reinit_events(
    realtype t, AmiVector const& x, AmiVector const& dx,
    std::vector<realtype> const& h_old, std::vector<int>& roots_found
) {
    std::vector<realtype> rootvals(ne, 0.0);
    froot(t, x, dx, rootvals);
    std::ranges::fill(roots_found, 0);
    for (int ie = 0; ie < ne; ie++) {
        if (rootvals.at(ie) < 0.0) {
            state_.h.at(ie) = 0.0;
        } else {
            state_.h.at(ie) = 1.0;
            if (h_old.at(ie) <= 0.0) {
                // only false->true triggers event
                roots_found.at(ie) = 1;
            }
        }
    }

    // re-compute parameter-dependent but state-independent roots
    reinit_explicit_roots();
}

void Model::reinit_explicit_roots() {
    explicit_roots_.clear();

    // Evaluate event trigger timepoints.
    // This assumes that `fw` was called before.
    // That happens during model initialization or pre-equilibration.
    auto const exp_roots = fexplicit_roots(
        state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
        derived_state_.w_.data()
    );
    Expects(exp_roots.size() == gsl::narrow<size_t>(ne - ne_solver));

    // group events by timepoints
    for (decltype(exp_roots)::size_type iee = 0; iee < exp_roots.size();
         ++iee) {
        // index within all events / root functions (not just explicit ones)
        int const ie = ne_solver + gsl::narrow<int>(iee);
        auto const& cur_roots = exp_roots[iee];
        Expects(!cur_roots.empty());
        for (auto const& root : cur_roots) {
            auto it = explicit_roots_.find(root);
            if (it != explicit_roots_.end()) {
                it->second.push_back(ie);
            } else {
                explicit_roots_[root] = {ie};
            }
        }
    }
}

int Model::nplist() const { return gsl::narrow<int>(state_.plist.size()); }

int Model::np() const {
    return gsl::narrow<int>(static_cast<ModelDimensions const&>(*this).np);
}

int Model::nk() const {
    return gsl::narrow<int>(state_.fixed_parameters.size());
}

int Model::ncl() const { return nx_rdata - nx_solver; }

int Model::nx_reinit() const { return nx_solver_reinit; }

double const* Model::k() const { return state_.fixed_parameters.data(); }

int Model::n_max_event() const { return nmaxevent_; }

void Model::set_n_max_event(int nmaxevent) { nmaxevent_ = nmaxevent; }

int Model::nt() const {
    return gsl::narrow<int>(simulation_parameters_.timepoints.size());
}

std::vector<ParameterScaling> const& Model::get_parameter_scale() const {
    return simulation_parameters_.pscale;
}

void Model::set_parameter_scale(ParameterScaling pscale) {
    simulation_parameters_.pscale.assign(
        simulation_parameters_.pscale.size(), pscale
    );
    scale_parameters(
        state_.unscaled_parameters, simulation_parameters_.pscale,
        simulation_parameters_.free_parameters
    );
    sx0data_.clear();
}

void Model::set_parameter_scale(
    std::vector<ParameterScaling> const& pscaleVec
) {
    if (pscaleVec.size() != simulation_parameters_.free_parameters.size())
        throw AmiException(
            "Dimension mismatch. Size of parameter scaling does "
            "not match number of model parameters."
        );
    simulation_parameters_.pscale = pscaleVec;
    scale_parameters(
        state_.unscaled_parameters, simulation_parameters_.pscale,
        simulation_parameters_.free_parameters
    );
    sx0data_.clear();
}

std::vector<realtype> const& Model::get_unscaled_parameters() const {
    return state_.unscaled_parameters;
}

std::vector<realtype> const& Model::get_free_parameters() const {
    return simulation_parameters_.free_parameters;
}

realtype Model::get_free_parameter_by_id(std::string const& par_id) const {
    if (!has_free_parameter_ids())
        throw AmiException(
            "Could not access parameters by id as they are not set"
        );
    return get_value_by_id(
        get_free_parameter_ids(), simulation_parameters_.free_parameters,
        par_id, "parameters", "id"
    );
}

realtype Model::get_free_parameter_by_name(std::string const& par_name) const {
    if (!has_free_parameter_names())
        throw AmiException(
            "Could not access parameters by name as they are not set"
        );
    return get_value_by_id(
        get_free_parameter_names(), simulation_parameters_.free_parameters,
        par_name, "parameters", "name"
    );
}

void Model::set_free_parameters(std::vector<realtype> const& p) {
    if (p.size() != (unsigned)np())
        throw AmiException(
            "Dimension mismatch. Size of parameters does not "
            "match number of model parameters."
        );
    simulation_parameters_.free_parameters = p;
    state_.unscaled_parameters.resize(
        simulation_parameters_.free_parameters.size()
    );
    unscale_parameters(
        simulation_parameters_.free_parameters, simulation_parameters_.pscale,
        state_.unscaled_parameters
    );
}

void Model::set_free_parameter_by_id(
    std::map<std::string, realtype> const& p, bool const ignoreErrors
) {
    for (auto const& [parameter_id, value] : p) {
        try {
            set_free_parameter_by_id(parameter_id, value);
        } catch (AmiException const&) {
            if (!ignoreErrors)
                throw;
        }
    }
}

void Model::set_free_parameter_by_id(
    std::string const& par_id, realtype const value
) {
    if (!has_free_parameter_ids())
        throw AmiException(
            "Could not access parameters by id as they are not set"
        );

    set_value_by_id(
        get_free_parameter_ids(), simulation_parameters_.free_parameters, value,
        par_id, "parameter", "id"
    );
    unscale_parameters(
        simulation_parameters_.free_parameters, simulation_parameters_.pscale,
        state_.unscaled_parameters
    );
}

int Model::set_free_parameters_by_id_regex(
    std::string const& par_id_regex, realtype const value
) {
    if (!has_free_parameter_ids())
        throw AmiException(
            "Could not access parameters by id as they are not set"
        );
    int n_found = set_value_by_id_regex(
        get_free_parameter_ids(), simulation_parameters_.free_parameters, value,
        par_id_regex, "parameter", "id"
    );
    unscale_parameters(
        simulation_parameters_.free_parameters, simulation_parameters_.pscale,
        state_.unscaled_parameters
    );
    return n_found;
}

void Model::set_free_parameter_by_name(
    std::string const& par_name, realtype const value
) {
    if (!has_free_parameter_names())
        throw AmiException(
            "Could not access parameters by name as they are not set"
        );

    set_value_by_id(
        get_free_parameter_names(), simulation_parameters_.free_parameters,
        value, par_name, "parameter", "name"
    );
    unscale_parameters(
        simulation_parameters_.free_parameters, simulation_parameters_.pscale,
        state_.unscaled_parameters
    );
}

void Model::set_free_parameter_by_name(
    std::map<std::string, realtype> const& p, bool ignoreErrors
) {
    for (auto const& [name, value] : p) {
        try {
            set_free_parameter_by_name(name, value);
        } catch (AmiException const&) {
            if (!ignoreErrors)
                throw;
        }
    }
}

int Model::set_free_parameters_by_name_regex(
    std::string const& par_name_regex, realtype value
) {
    if (!has_free_parameter_names())
        throw AmiException(
            "Could not access parameters by name as they are not set"
        );

    int n_found = set_value_by_id_regex(
        get_free_parameter_names(), simulation_parameters_.free_parameters,
        value, par_name_regex, "parameter", "name"
    );

    unscale_parameters(
        simulation_parameters_.free_parameters, simulation_parameters_.pscale,
        state_.unscaled_parameters
    );
    return n_found;
}

std::vector<realtype> const& Model::get_fixed_parameters() const {
    return state_.fixed_parameters;
}

realtype Model::get_fixed_parameter_by_id(std::string const& par_id) const {
    if (!has_fixed_parameter_ids())
        throw AmiException(
            "Could not access fixed parameters by id as they are not set"
        );

    return get_value_by_id(
        get_fixed_parameter_ids(), state_.fixed_parameters, par_id,
        "fixedParameters", "id"
    );
}

realtype Model::get_fixed_parameter_by_name(std::string const& par_name) const {
    if (!has_fixed_parameter_names())
        throw AmiException(
            "Could not access fixed parameters by name as they are not set"
        );

    return get_value_by_id(
        get_fixed_parameter_names(), state_.fixed_parameters, par_name,
        "fixedParameters", "name"
    );
}

void Model::set_fixed_parameters(std::vector<realtype> const& k) {
    if (k.size() != (unsigned)nk())
        throw AmiException(
            "Dimension mismatch. Size of fixedParameters does "
            "not match number of fixed model parameters."
        );
    state_.fixed_parameters = k;
}

void Model::set_fixed_parameter_by_id(
    std::string const& par_id, realtype value
) {
    if (!has_fixed_parameter_ids())
        throw AmiException(
            "Could not access fixed parameters by id as they are not set"
        );

    set_value_by_id(
        get_fixed_parameter_ids(), state_.fixed_parameters, value, par_id,
        "fixedParameters", "id"
    );
}

int Model::set_fixed_parameters_by_id_regex(
    std::string const& par_id_regex, realtype value
) {
    if (!has_fixed_parameter_ids())
        throw AmiException(
            "Could not access fixed parameters by id as they are not set"
        );

    return set_value_by_id_regex(
        get_fixed_parameter_ids(), state_.fixed_parameters, value, par_id_regex,
        "fixedParameters", "id"
    );
}

void Model::set_fixed_parameter_by_name(
    std::string const& par_name, realtype value
) {
    if (!has_fixed_parameter_names())
        throw AmiException(
            "Could not access fixed parameters by name as they are not set"
        );

    set_value_by_id(
        get_fixed_parameter_names(), state_.fixed_parameters, value, par_name,
        "fixedParameters", "name"
    );
}

int Model::set_fixed_parameters_by_name_regex(
    std::string const& par_name_regex, realtype value
) {
    if (!has_fixed_parameter_names())
        throw AmiException(
            "Could not access fixed parameters by name as they are not set"
        );

    return set_value_by_id_regex(
        get_fixed_parameter_ids(), state_.fixed_parameters, value,
        par_name_regex, "fixedParameters", "name"
    );
}

std::string Model::get_name() const { return ""; }

bool Model::has_free_parameter_names() const {
    return np() == 0 || !get_free_parameter_names().empty();
}

std::vector<std::string> Model::get_free_parameter_names() const { return {}; }

bool Model::has_state_names() const {
    return nx_rdata == 0 || !get_state_names().empty();
}

std::vector<std::string> Model::get_state_names() const { return {}; }

std::vector<std::string> Model::get_state_names_solver() const { return {}; }

bool Model::has_fixed_parameter_names() const {
    return nk() == 0 || !get_fixed_parameter_names().empty();
}

std::vector<std::string> Model::get_fixed_parameter_names() const { return {}; }

bool Model::has_observable_names() const {
    return ny == 0 || !get_observable_names().empty();
}

std::vector<std::string> Model::get_observable_names() const { return {}; }

bool Model::has_expression_names() const {
    return ny == 0 || !get_expression_names().empty();
}

std::vector<std::string> Model::get_expression_names() const { return {}; }

bool Model::has_free_parameter_ids() const {
    return np() == 0 || !get_free_parameter_ids().empty();
}

std::vector<std::string> Model::get_free_parameter_ids() const { return {}; }

bool Model::has_state_ids() const {
    return nx_rdata == 0 || !get_state_ids().empty();
}

std::vector<std::string> Model::get_state_ids() const { return {}; }

std::vector<std::string> Model::get_state_ids_solver() const { return {}; }

bool Model::has_fixed_parameter_ids() const {
    return nk() == 0 || !get_fixed_parameter_ids().empty();
}

std::vector<std::string> Model::get_fixed_parameter_ids() const { return {}; }

bool Model::has_observable_ids() const {
    return ny == 0 || !get_observable_ids().empty();
}

std::vector<std::string> Model::get_observable_ids() const { return {}; }

bool Model::has_expression_ids() const {
    return ny == 0 || !get_expression_ids().empty();
}

std::vector<std::string> Model::get_expression_ids() const { return {}; }

bool Model::has_quadratic_llh() const { return true; }

std::vector<realtype> const& Model::get_timepoints() const {
    return simulation_parameters_.timepoints;
}

double Model::get_timepoint(int const it) const {
    return simulation_parameters_.timepoints.at(it);
}

void Model::set_timepoints(std::vector<realtype> const& ts) {
    if (!std::ranges::is_sorted(ts))
        throw AmiException(
            "Encountered non-monotonic timepoints, please order"
            " timepoints such that they are monotonically"
            " increasing!"
        );
    simulation_parameters_.timepoints = ts;
}

double Model::t0() const { return simulation_parameters_.t_start; }

void Model::set_t0(double t0) { simulation_parameters_.t_start = t0; }

double Model::t0_preeq() const { return simulation_parameters_.t_start_preeq; }

void Model::set_t0_preeq(double t0_preeq) {
    simulation_parameters_.t_start_preeq = t0_preeq;
}

std::vector<bool> const& Model::get_state_is_non_negative() const {
    return state_is_non_negative_;
}

void Model::set_state_is_non_negative(std::vector<bool> const& nonNegative) {
    auto any_state_non_negative
        = std::ranges::any_of(nonNegative, [](bool x) { return x; });
    if (nx_solver != nx_rdata) {
        if (any_state_non_negative)
            throw AmiException(
                "Non-negative states are not supported with"
                " conservation laws enabled."
            );
        // nothing to do, as `state_is_non_negative_` will always be all-false
        // in case of conservation laws
        return;
    }
    if (nonNegative.size() != gsl::narrow<unsigned long>(nx_rdata)) {
        throw AmiException(
            "Dimension of input stateIsNonNegative (%u) does "
            "not agree with number of state variables (%d)",
            nonNegative.size(), nx_rdata
        );
    }
    state_is_non_negative_ = nonNegative;
    any_state_non_negative_ = any_state_non_negative;
}

void Model::set_all_states_non_negative() {
    set_state_is_non_negative(std::vector<bool>(nx_solver, true));
}

std::vector<int> const& Model::get_parameter_list() const {
    return state_.plist;
}

int Model::plist(int pos) const { return state_.plist.at(pos); }

void Model::set_parameter_list(std::vector<int> const& plist) {
    int np = this->np(); // cannot capture 'this' in lambda expression
    if (std::ranges::any_of(plist, [&np](int idx) {
            return idx < 0 || idx >= np;
        })) {
        throw AmiException("Indices in plist must be in [0..np]");
    }
    state_.plist = plist;

    initialize_vectors();
}

std::vector<realtype> Model::get_initial_state(realtype t0) {
    if (!x0data_.empty()) {
        return x0data_;
    }

    /* Initial states have not been set explicitly on this instance, so we
     * compute it, but don't save it, as this would have to be invalidated upon
     * changing parameters etc.
     */
    std::vector<realtype> x0(nx_rdata, 0.0);
    fx0(x0.data(), t0, state_.unscaled_parameters.data(),
        state_.fixed_parameters.data());
    return x0;
}

void Model::set_initial_state(std::vector<realtype> const& x0) {
    if (x0.size() != (unsigned)nx_rdata && !x0.empty())
        throw AmiException(
            "Dimension mismatch. Size of x0 does not match "
            "number of model states."
        );

    if (x0.empty()) {
        x0data_.clear();
        return;
    }

    x0data_ = x0;
}

bool Model::has_custom_initial_state() const { return !x0data_.empty(); }

std::vector<realtype> Model::get_initial_state_sensitivities(realtype t0) {
    if (!sx0data_.empty()) {
        return sx0data_;
    }

    /* Initial state sensitivities have not been set explicitly on this
     * instance, so we compute it, but don't save it, as this would have to be
     * invalidated upon changing parameters etc.
     */
    std::vector<realtype> sx0(nx_rdata * nplist(), 0.0);
    auto x0 = get_initial_state(t0);
    for (int ip = 0; ip < nplist(); ip++) {
        fsx0(
            sx0.data(), t0, x0.data(), state_.unscaled_parameters.data(),
            state_.fixed_parameters.data(), plist(ip)
        );
    }
    return sx0;
}

void Model::set_initial_state_sensitivities(std::vector<realtype> const& sx0) {
    if (sx0.size() != (unsigned)nx_rdata * nplist() && !sx0.empty())
        throw AmiException(
            "Dimension mismatch. Size of sx0 does not match "
            "number of model states * number of parameter "
            "selected for sensitivities."
        );

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
            chainrulefactor
                = state_.unscaled_parameters.at(plist(ip)) * log(10);
            break;
        case ParameterScaling::ln:
            chainrulefactor = state_.unscaled_parameters.at(plist(ip));
            break;
        case ParameterScaling::none:
            chainrulefactor = 1.0;
            break;
        }

        for (int ix = 0; ix < nx_rdata; ++ix) {
            sx0_rdata.at(ip * nx_rdata + ix)
                = sx0.at(ip * nx_rdata + ix) / chainrulefactor;
        }
    }
    set_unscaled_initial_state_sensitivities(sx0_rdata);
}

bool Model::has_custom_initial_state_sensitivities() const {
    return !sx0data_.empty();
}

void Model::set_unscaled_initial_state_sensitivities(
    std::vector<realtype> const& sx0
) {
    if (sx0.size() != (unsigned)nx_rdata * nplist() && !sx0.empty())
        throw AmiException(
            "Dimension mismatch. Size of sx0 does not match "
            "number of model states * number of parameter "
            "selected for sensitivities."
        );

    if (sx0.empty()) {
        sx0data_.clear();
        return;
    }

    sx0data_ = sx0;
}

void Model::set_steady_state_computation_mode(
    SteadyStateComputationMode const mode
) {
    if (mode != SteadyStateComputationMode::integrationOnly && ne && logger_) {
        logger_->log(
            LogSeverity::warning, "WARNING",
            "Non-initial events will not be handled if Newton's method is used "
            "for steady state computation."
        );
    }
    steadystate_computation_mode_ = mode;
}

SteadyStateComputationMode Model::get_steady_state_computation_mode() const {
    return steadystate_computation_mode_;
}

void Model::set_steady_state_sensitivity_mode(
    SteadyStateSensitivityMode const mode
) {
    steadystate_sensitivity_mode_ = mode;
}

SteadyStateSensitivityMode Model::get_steady_state_sensitivity_mode() const {
    return steadystate_sensitivity_mode_;
}

void Model::set_reinitialize_fixed_parameter_initial_states(bool flag) {
    if (flag && !is_fixed_parameter_state_reinitialization_allowed())
        throw AmiException(
            "State reinitialization cannot be enabled for this model "
            "as this feature was disabled at compile time. Most likely,"
            " this was because some initial states depending on "
            "fixedParameters also depended on parameters."
        );
    simulation_parameters_.reinitialize_fixed_parameter_initial_states = flag;

    if (flag) {
        simulation_parameters_
            .reinitialize_all_fixed_parameter_dependent_initial_states_for_simulation(
                nx_rdata
            );
    } else {
        simulation_parameters_.reinitialization_state_idxs_sim.clear();
    }
}

bool Model::get_reinitialize_fixed_parameter_initial_states() const {
    return simulation_parameters_.reinitialize_fixed_parameter_initial_states
           || !simulation_parameters_.reinitialization_state_idxs_sim.empty();
}

void Model::require_sensitivities_for_all_parameters() {
    state_.plist.resize(np());
    std::iota(state_.plist.begin(), state_.plist.end(), 0);
    initialize_vectors();
}

void Model::get_expression(
    gsl::span<realtype> w, realtype const t, AmiVector const& x
) {
    fw(t, compute_x_pos(x), false);
    write_slice(derived_state_.w_, w);
}

void Model::get_observable(
    gsl::span<realtype> y, realtype const t, AmiVector const& x
) {
    fy(t, x);
    write_slice(derived_state_.y_, y);
}

ObservableScaling Model::get_observable_scaling(int /*iy*/) const {
    return ObservableScaling::lin;
}

void Model::get_observable_sensitivity(
    gsl::span<realtype> sy, realtype const t, AmiVector const& x,
    AmiVectorArray const& sx
) {
    if (!ny)
        return;

    fdydx(t, x);
    fsspl(t);
    fdydp(t, x);

    derived_state_.sx_.resize(nx_solver * nplist());
    sx.flatten_to_vector(derived_state_.sx_);

    // compute sy = 1.0*dydx*sx + 1.0*dydp
    // dydx A[ny,nx_solver] * sx B[nx_solver,nplist] = sy C[ny,nplist]
    //        M  K                 K  N                     M  N
    //        lda                  ldb                      ldc

    if (nx_solver) {
        set_nan_to_zero(derived_state_.dydx_);
        set_nan_to_zero(derived_state_.sx_);
        amici_dgemm(
            BLASLayout::colMajor, BLASTranspose::noTrans,
            BLASTranspose::noTrans, ny, nplist(), nx_solver, 1.0,
            derived_state_.dydx_.data(), ny, derived_state_.sx_.data(),
            nx_solver, 1.0, derived_state_.dydp_.data(), ny
        );
    }
    write_slice(derived_state_.dydp_, sy);

    if (always_check_finite_)
        check_finite(sy, ModelQuantity::sy, nplist());
}

void Model::get_observable_sigma(
    gsl::span<realtype> sigmay, int const it, ExpData const* edata
) {
    fsigmay(it, edata);
    write_slice(derived_state_.sigmay_, sigmay);
}

void Model::get_observable_sigma_sensitivity(
    gsl::span<realtype> ssigmay, gsl::span<realtype const> sy, int const it,
    ExpData const* edata
) {
    fdsigmaydp(it, edata);
    write_slice(derived_state_.dsigmaydp_, ssigmay);

    // ssigmay = dsigmaydy*(dydx_solver*sx+dydp)+dsigmaydp
    //         = dsigmaydy*sy+dsigmaydp

    fdsigmaydy(it, edata);

    // compute ssigmay = 1.0 * dsigmaydp + 1.0 * dsigmaydy * sy
    // dsigmaydp C[ny,nplist] += dsigmaydy A[ny,ny] * sy B[ny,nplist]
    //             M  N                      M  K          K  N
    //             ldc                       lda           ldb
    set_nan_to_zero(derived_state_.dsigmaydy_);
    derived_state_.sy_.assign(sy.begin(), sy.end());
    set_nan_to_zero(derived_state_.sy_);
    amici_dgemm(
        BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
        ny, nplist(), ny, 1.0, derived_state_.dsigmaydy_.data(), ny,
        derived_state_.sy_.data(), ny, 1.0, ssigmay.data(), ny
    );

    if (always_check_finite_)
        check_finite(ssigmay, ModelQuantity::ssigmay, nplist());
}

void Model::add_observable_objective(
    realtype& Jy, int const it, AmiVector const& x, ExpData const& edata
) {
    fy(edata.get_timepoint(it), x);
    fsigmay(it, &edata);

    std::vector<realtype> nllh(nJ, 0.0);
    for (int iyt = 0; iyt < nytrue; iyt++) {
        if (edata.is_set_measurement(it, iyt)) {
            std::ranges::fill(nllh, 0.0);
            fJy(nllh.data(), iyt, state_.unscaled_parameters.data(),
                state_.fixed_parameters.data(), derived_state_.y_.data(),
                derived_state_.sigmay_.data(), edata.get_measurements_ptr(it));
            Jy -= nllh.at(0);
        }
    }
}

void Model::add_observable_objective_sensitivity(
    std::vector<realtype>& sllh, std::vector<realtype>& s2llh, int const it,
    AmiVector const& x, AmiVectorArray const& sx, ExpData const& edata
) {

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
    set_nan_to_zero(derived_state_.dJydx_);
    set_nan_to_zero(derived_state_.sx_);
    amici_dgemm(
        BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
        nJ, nplist(), nx_solver, 1.0, derived_state_.dJydx_.data(), nJ,
        derived_state_.sx_.data(), nx_solver, 1.0, derived_state_.dJydp_.data(),
        nJ
    );

    write_llh_sensitivity_slice(derived_state_.dJydp_, sllh, s2llh);
}

void Model::add_partial_observable_objective_sensitivity(
    std::vector<realtype>& sllh, std::vector<realtype>& s2llh, int const it,
    AmiVector const& x, ExpData const& edata
) {
    if (!ny)
        return;

    fdJydp(it, x, edata);

    write_llh_sensitivity_slice(derived_state_.dJydp_, sllh, s2llh);
}

void Model::get_adjoint_state_observable_update(
    gsl::span<realtype> dJydx, int const it, AmiVector const& x,
    ExpData const& edata
) {
    fdJydx(it, x, edata);
    write_slice(derived_state_.dJydx_, dJydx);
}

void Model::get_event(
    gsl::span<realtype> z, int const ie, realtype const t, AmiVector const& x
) {
    fz(ie, t, x);
    write_slice_event(derived_state_.z_, z, ie);
}

void Model::get_event_sensitivity(
    gsl::span<realtype> sz, int const ie, realtype const t, AmiVector const& x,
    AmiVectorArray const& sx
) {
    if (!nz)
        return;

    fdzdx(ie, t, x);
    fdzdp(ie, t, x);

    derived_state_.sx_.resize(nx_solver * nplist());
    sx.flatten_to_vector(derived_state_.sx_);

    // compute sz = 1.0*dzdx*sx + 1.0*dzdp
    // dzdx A[nz,nx_solver] * sx B[nx_solver,nplist] = sz C[nz,nplist]
    //        M  K                 K  N                     M  N
    //        lda                  ldb                      ldc
    set_nan_to_zero(derived_state_.dzdx_);
    set_nan_to_zero(derived_state_.sx_);
    amici_dgemm(
        BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
        nz, nplist(), nx_solver, 1.0, derived_state_.dzdx_.data(), nz,
        derived_state_.sx_.data(), nx_solver, 1.0, derived_state_.dzdp_.data(),
        nz
    );

    add_slice(derived_state_.dzdp_, sz);

    if (always_check_finite_)
        check_finite(sz, ModelQuantity::sz, nplist());
}

void Model::get_unobserved_event_sensitivity(
    gsl::span<realtype> sz, int const ie
) {
    check_buffer_size(sz, nz * nplist());

    for (int iz = 0; iz < nz; ++iz)
        if (z2event_.at(iz) == ie)
            for (int ip = 0; ip < nplist(); ++ip)
                sz[ip * nz + iz] = 0.0;
}

void Model::get_event_regularization(
    gsl::span<realtype> rz, int const ie, realtype const t, AmiVector const& x
) {
    frz(ie, t, x);
    write_slice_event(derived_state_.rz_, rz, ie);
}

void Model::get_event_regularization_sensitivity(
    gsl::span<realtype> srz, int const ie, realtype const t, AmiVector const& x,
    AmiVectorArray const& sx
) {

    if (!nz)
        return;

    fdrzdx(ie, t, x);
    fdrzdp(ie, t, x);

    derived_state_.sx_.resize(nx_solver * nplist());
    sx.flatten_to_vector(derived_state_.sx_);

    // compute srz = 1.0*drzdx*sx + 1.0*drzdp
    // drzdx A[nz,nx_solver] * sx B[nx_solver,nplist] = srz C[nz,nplist]
    //         M  K                 K  N                      M  N
    //         lda                  ldb                       ldc
    set_nan_to_zero(derived_state_.drzdx_);
    set_nan_to_zero(derived_state_.sx_);
    amici_dgemm(
        BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
        nz, nplist(), nx_solver, 1.0, derived_state_.drzdx_.data(), nz,
        derived_state_.sx_.data(), nx_solver, 1.0, derived_state_.drzdp_.data(),
        nz
    );

    add_slice(derived_state_.drzdp_, srz);

    if (always_check_finite_)
        check_finite(srz, ModelQuantity::srz, nplist());
}

void Model::get_event_sigma(
    gsl::span<realtype> sigmaz, int const ie, int const nroots,
    realtype const t, ExpData const* edata
) {
    fsigmaz(ie, nroots, t, edata);
    write_slice_event(derived_state_.sigmaz_, sigmaz, ie);
}

void Model::get_event_sigma_sensitivity(
    gsl::span<realtype> ssigmaz, int const ie, int const nroots,
    realtype const t, ExpData const* edata
) {
    fdsigmazdp(ie, nroots, t, edata);
    write_sensitivity_slice_event(derived_state_.dsigmazdp_, ssigmaz, ie);
}

void Model::add_event_objective(
    realtype& Jz, int const ie, int const nroots, realtype const t,
    AmiVector const& x, ExpData const& edata
) {
    fz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    std::vector<realtype> nllh(nJ, 0.0);
    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.is_set_event_measurement(nroots, iztrue)) {
            std::ranges::fill(nllh, 0.0);
            fJz(nllh.data(), iztrue, state_.unscaled_parameters.data(),
                state_.fixed_parameters.data(), derived_state_.z_.data(),
                derived_state_.sigmaz_.data(),
                edata.get_event_measurements_ptr(nroots));
            Jz -= nllh.at(0);
        }
    }
}

void Model::add_event_objective_regularization(
    realtype& Jrz, int const ie, int const nroots, realtype const t,
    AmiVector const& x, ExpData const& edata
) {
    frz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    std::vector<realtype> nllh(nJ, 0.0);
    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.is_set_event_measurement(nroots, iztrue)) {
            std::ranges::fill(nllh, 0.0);
            fJrz(
                nllh.data(), iztrue, state_.unscaled_parameters.data(),
                state_.fixed_parameters.data(), derived_state_.rz_.data(),
                derived_state_.sigmaz_.data()
            );
            Jrz -= nllh.at(0);
        }
    }
}

void Model::add_event_objective_sensitivity(
    std::vector<realtype>& sllh, std::vector<realtype>& s2llh, int const ie,
    int const nroots, realtype const t, AmiVector const& x,
    AmiVectorArray const& sx, ExpData const& edata
) {

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
    set_nan_to_zero(derived_state_.dJzdx_);
    set_nan_to_zero(derived_state_.sx_);
    amici_dgemm(
        BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
        nJ, nplist(), nx_solver, 1.0, derived_state_.dJzdx_.data(), nJ,
        derived_state_.sx_.data(), nx_solver, 1.0, derived_state_.dJzdp_.data(),
        nJ
    );

    // sJy += multResult + dJydp
    write_llh_sensitivity_slice(derived_state_.dJzdp_, sllh, s2llh);
}

void Model::get_adjoint_state_event_update(
    gsl::span<realtype> dJzdx, int const ie, int const nroots, realtype const t,
    AmiVector const& x, ExpData const& edata
) {
    fdJzdx(ie, nroots, t, x, edata);
    write_slice(derived_state_.dJzdx_, dJzdx);
}

void Model::add_partial_event_objective_sensitivity(
    std::vector<realtype>& sllh, std::vector<realtype>& s2llh, int const ie,
    int const nroots, realtype const t, AmiVector const& x, ExpData const& edata
) {
    if (!nz)
        return;

    fdJzdp(ie, nroots, t, x, edata);

    write_llh_sensitivity_slice(derived_state_.dJzdp_, sllh, s2llh);
}

void Model::get_event_time_sensitivity(
    std::vector<realtype>& stau, realtype const t, int const ie,
    AmiVector const& x, AmiVectorArray const& sx, AmiVector const& dx
) {

    std::ranges::fill(stau, 0.0);

    for (int ip = 0; ip < nplist(); ip++) {
        fstau(
            &stau.at(ip), t, compute_x_pos(x),
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            state_.h.data(), derived_state_.w_.data(), dx.data(),
            state_.total_cl.data(), sx.data(ip), plist(ip), ie
        );
    }
}

void Model::add_state_event_update(
    AmiVector& x, int const ie, realtype const t, AmiVector const& xdot,
    AmiVector const& xdot_old, AmiVector const& x_old, ModelState const& state
) {

    derived_state_.deltax_.assign(nx_solver, 0.0);

    std::copy_n(compute_x_pos(x), nx_solver, x.data());

    // compute update
    fdeltax(
        derived_state_.deltax_.data(), t, x.data(),
        state.unscaled_parameters.data(), state.fixed_parameters.data(),
        state.h.data(), ie, xdot.data(), xdot_old.data(), x_old.data()
    );

    if (always_check_finite_) {
        check_finite(derived_state_.deltax_, ModelQuantity::deltax, t);
    }

    // update
    amici_daxpy(nx_solver, 1.0, derived_state_.deltax_.data(), 1, x.data(), 1);
}

void Model::add_state_sensitivity_event_update(
    AmiVectorArray& sx, int const ie, realtype const t, AmiVector const& x,
    AmiVector const& x_old, AmiVector const& xdot, AmiVector const& xdot_old,
    AmiVectorArray const& sx_old, std::vector<realtype> const& stau
) {
    fw(t, x_old.data(), false);

    for (int ip = 0; ip < nplist(); ip++) {

        derived_state_.deltasx_.assign(nx_solver, 0.0);

        // compute update
        fdeltasx(
            derived_state_.deltasx_.data(), t, x.data(),
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            state_.h.data(), derived_state_.w_.data(), plist(ip), ie,
            xdot.data(), xdot_old.data(), sx_old.data(ip), &stau.at(ip),
            state_.total_cl.data(), x_old.data()
        );
        if (always_check_finite_) {
            check_finite(
                derived_state_.deltasx_, ModelQuantity::deltasx, nplist()
            );
        }

        amici_daxpy(
            nx_solver, 1.0, derived_state_.deltasx_.data(), 1, sx.data(ip), 1
        );
    }
}

void Model::add_adjoint_state_event_update(
    AmiVector& xB, int const ie, realtype const t, AmiVector const& x,
    AmiVector const& xdot, AmiVector const& xdot_old, AmiVector const& x_old,
    AmiVector const& dx
) {

    derived_state_.deltaxB_.assign(nxtrue_solver * nJ, 0.0);

    // compute update
    fdeltaxB(
        derived_state_.deltaxB_.data(), t, compute_x_pos(x),
        state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
        state_.h.data(), derived_state_.w_.data(), dx.data(), ie, xdot.data(),
        xdot_old.data(), x_old.data(), xB.data(), state_.total_cl.data()
    );

    if (always_check_finite_) {
        check_finite(derived_state_.deltaxB_, ModelQuantity::deltaxB, t);
    }

    // apply update
    for (int ix = 0; ix < nxtrue_solver; ++ix)
        for (int iJ = 0; iJ < nJ; ++iJ)
            xB.at(ix + iJ * nxtrue_solver)
                += derived_state_.deltaxB_.at(ix + iJ * nxtrue_solver);
}

void Model::add_adjoint_quadrature_event_update(
    AmiVector& xQB, int const ie, realtype const t, AmiVector const& x,
    AmiVector const& xB, AmiVector const& xdot, AmiVector const& xdot_old,
    AmiVector const& x_old, AmiVector const& dx
) {
    for (int ip = 0; ip < nplist(); ip++) {
        derived_state_.deltaqB_.assign(nJ, 0.0);

        fdeltaqB(
            derived_state_.deltaqB_.data(), t, compute_x_pos(x),
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            state_.h.data(), derived_state_.w_.data(), dx.data(), plist(ip), ie,
            xdot.data(),
            xdot_old.data(), x_old.data(), xB.data()
        );

        for (int iJ = 0; iJ < nJ; ++iJ)
            xQB.at(ip * nJ + iJ) += derived_state_.deltaqB_.at(iJ);
    }

    if (always_check_finite_) {
        check_finite(derived_state_.deltaqB_, ModelQuantity::deltaqB, t);
    }
}

void Model::update_heaviside(std::vector<int> const& rootsfound) {
    for (int ie = 0; ie < ne; ie++) {
        state_.h.at(ie) += rootsfound.at(ie);
    }
}

int Model::check_finite(
    gsl::span<realtype const> array, ModelQuantity model_quantity, realtype t
) const {
    auto it = std::ranges::find_if(array, [](realtype x) {
        return !std::isfinite(x);
    });
    if (it == array.end()) {
        return AMICI_SUCCESS;
    }

    // there is some issue - produce a meaningful message
    auto flat_index = it - array.begin();

    std::string msg_id;
    std::string non_finite_type;
    if (std::isnan(array[flat_index])) {
        msg_id = "AMICI:NaN";
        non_finite_type = "NaN";
    } else if (std::isinf(array[flat_index])) {
        msg_id = "AMICI:Inf";
        non_finite_type = "Inf";
    }
    std::string element_id = std::to_string(flat_index);

    switch (model_quantity) {
    case ModelQuantity::xdot:
    case ModelQuantity::xBdot:
    case ModelQuantity::x0:
    case ModelQuantity::x:
    case ModelQuantity::x_rdata:
    case ModelQuantity::x0_rdata:
    case ModelQuantity::Jv:
    case ModelQuantity::JvB:
    case ModelQuantity::JDiag:
    case ModelQuantity::deltax:
    case ModelQuantity::deltaxB:
        if (has_state_ids()) {
            element_id = get_state_ids_solver()[flat_index];
        }
        break;
    case ModelQuantity::y:
        if (has_observable_ids()) {
            element_id = get_observable_ids()[flat_index];
        }
        break;
    case ModelQuantity::w:
        if (has_expression_ids()) {
            element_id = get_expression_ids()[flat_index];
        }
        break;
    case ModelQuantity::k:
        if (has_fixed_parameter_ids()) {
            element_id = get_fixed_parameter_ids()[flat_index];
        }
        break;
    case ModelQuantity::p:
        if (has_free_parameter_ids()) {
            element_id = get_free_parameter_ids()[flat_index];
        }
        break;
    default:
        break;
    }

    std::string model_quantity_str;
    try {
        model_quantity_str = model_quantity_to_str.at(model_quantity);
    } catch (std::out_of_range const&) {
        // Missing model quantity string - terminate if this is a debug build,
        // but show the quantity number if non-debug.
        gsl_ExpectsDebug(false);
        model_quantity_str = std::to_string(static_cast<int>(model_quantity));
    }
    if (logger_) {
        auto t_msg = std::isfinite(t)
                         ? std::string(" at t=" + std::to_string(t) + " ")
                         : std::string();

        logger_->log(
            LogSeverity::warning, msg_id,
            "AMICI encountered a %s value for %s[%i] (%s)%s",
            non_finite_type.c_str(), model_quantity_str.c_str(),
            gsl::narrow<int>(flat_index), element_id.c_str(), t_msg.c_str()
        );
    }
    // check upstream, without infinite recursion
    if (model_quantity != ModelQuantity::k && model_quantity != ModelQuantity::p
        && model_quantity != ModelQuantity::ts) {
        check_finite(state_.fixed_parameters, ModelQuantity::k, t);
        check_finite(state_.unscaled_parameters, ModelQuantity::p, t);
        if (!always_check_finite_ && model_quantity != ModelQuantity::w) {
            // don't check twice if always_check_finite_ is true
            check_finite(derived_state_.w_, ModelQuantity::w, t);
        }
    }
    return AMICI_RECOVERABLE_ERROR;
}

int Model::check_finite(
    gsl::span<realtype const> array, ModelQuantity model_quantity,
    size_t num_cols, realtype t
) const {
    auto it = std::ranges::find_if(array, [](realtype x) {
        return !std::isfinite(x);
    });
    if (it == array.end()) {
        return AMICI_SUCCESS;
    }

    // there is some issue - produce a meaningful message
    auto flat_index = it - array.begin();
    sunindextype row, col;
    std::tie(row, col) = unravel_index(flat_index, num_cols);
    std::string msg_id;
    std::string non_finite_type;
    if (std::isnan(array[flat_index])) {
        msg_id = "AMICI:NaN";
        non_finite_type = "NaN";
    } else if (std::isinf(array[flat_index])) {
        msg_id = "AMICI:Inf";
        non_finite_type = "Inf";
    }
    std::string row_id = std::to_string(row);
    std::string col_id = std::to_string(col);

    switch (model_quantity) {
    case ModelQuantity::sy:
    case ModelQuantity::ssigmay:
    case ModelQuantity::dydp:
    case ModelQuantity::dsigmaydp:
        if (has_observable_ids())
            row_id += " " + get_observable_ids()[row];
        if (has_free_parameter_ids())
            col_id
                += " " + get_free_parameter_ids()[plist(gsl::narrow<int>(col))];
        break;
    case ModelQuantity::dydx:
        if (has_observable_ids())
            row_id += " " + get_observable_ids()[row];
        if (has_state_ids())
            col_id += " " + get_state_ids_solver()[col];
        break;
    case ModelQuantity::deltasx:
        if (has_state_ids())
            row_id += " " + get_state_ids_solver()[row];
        if (has_free_parameter_ids())
            col_id
                += " " + get_free_parameter_ids()[plist(gsl::narrow<int>(col))];
        break;
    case ModelQuantity::dJydy:
    case ModelQuantity::dJydsigma:
        if (has_observable_ids())
            col_id += " " + get_observable_ids()[col];
        break;
    case ModelQuantity::dJydx:
    case ModelQuantity::dJzdx:
    case ModelQuantity::dJrzdx:
    case ModelQuantity::dzdx:
    case ModelQuantity::drzdx:
        if (has_state_ids())
            col_id += " " + get_state_ids_solver()[col];
        break;
    case ModelQuantity::deltaqB:
    case ModelQuantity::sz:
    case ModelQuantity::dzdp:
    case ModelQuantity::drzdp:
    case ModelQuantity::dsigmazdp:
        if (has_free_parameter_ids())
            col_id
                += " " + get_free_parameter_ids()[plist(gsl::narrow<int>(col))];
        break;
    case ModelQuantity::dsigmaydy:
        if (has_observable_ids()) {
            auto obs_ids = get_observable_ids();
            row_id += " " + obs_ids[row];
            col_id += " " + obs_ids[col];
        }
        break;
    default:
        break;
    }

    std::string model_quantity_str;
    try {
        model_quantity_str = model_quantity_to_str.at(model_quantity);
    } catch (std::out_of_range const&) {
        // Missing model quantity string - terminate if this is a debug build,
        // but show the quantity number if non-debug.
        gsl_ExpectsDebug(false);
        model_quantity_str = std::to_string(static_cast<int>(model_quantity));
    }

    if (logger_) {
        auto t_msg = std::isfinite(t)
                         ? std::string(" at t=" + std::to_string(t) + " ")
                         : std::string();

        logger_->log(
            LogSeverity::warning, msg_id,
            "AMICI encountered a %s value for %s[%i] (%s, %s)%s",
            non_finite_type.c_str(), model_quantity_str.c_str(),
            gsl::narrow<int>(flat_index), row_id.c_str(), col_id.c_str(),
            t_msg.c_str()
        );
    }

    // check upstream
    check_finite(state_.fixed_parameters, ModelQuantity::k, t);
    check_finite(state_.unscaled_parameters, ModelQuantity::p, t);
    check_finite(derived_state_.w_, ModelQuantity::w, t);

    return AMICI_RECOVERABLE_ERROR;
}

int Model::check_finite(
    SUNMatrix m, ModelQuantity model_quantity, realtype t
) const {
    // check flat array, to see if there are any issues
    // (faster, in particular for sparse arrays)
    auto m_flat = gsl::make_span(m);
    auto it = std::ranges::find_if(m_flat, [](realtype x) {
        return !std::isfinite(x);
    });
    if (it == m_flat.end()) {
        return AMICI_SUCCESS;
    }

    // there is some issue - produce a meaningful message
    auto flat_index = it - m_flat.begin();
    auto [row, col] = unravel_index(flat_index, m);
    std::string msg_id;
    std::string non_finite_type;
    if (std::isnan(m_flat[flat_index])) {
        msg_id = "AMICI:NaN";
        non_finite_type = "NaN";
    } else if (std::isinf(m_flat[flat_index])) {
        msg_id = "AMICI:Inf";
        non_finite_type = "Inf";
    } else {
        throw std::runtime_error(
            "Value is not finite, but neither infinite nor NaN."
        );
    }
    std::string row_id = std::to_string(row);
    std::string col_id = std::to_string(col);

    switch (model_quantity) {
    case ModelQuantity::J:
    case ModelQuantity::JB:
        if (has_state_ids()) {
            auto state_ids = get_state_ids_solver();
            row_id += " " + state_ids[row];
            col_id += " " + state_ids[col];
        }
        break;
    case ModelQuantity::dwdx:
        if (has_expression_ids())
            row_id += " " + get_expression_ids()[row];
        if (has_state_ids())
            col_id += " " + get_state_ids_solver()[col];
        break;
    case ModelQuantity::dwdw:
        if (has_expression_ids()) {
            auto expr_ids = get_expression_ids();
            row_id += " " + expr_ids[row];
            col_id += " " + expr_ids[col];
        }
        break;
    case ModelQuantity::dwdp:
        if (has_expression_ids())
            row_id += " " + get_expression_ids()[row];
        if (has_free_parameter_ids())
            col_id += " " + get_free_parameter_ids()[col];
        break;
    default:
        break;
    }

    std::string model_quantity_str;
    try {
        model_quantity_str = model_quantity_to_str.at(model_quantity);
    } catch (std::out_of_range const&) {
        // Missing model quantity string - terminate if this is a debug build,
        // but show the quantity number if non-debug.
        gsl_ExpectsDebug(false);
        model_quantity_str = std::to_string(static_cast<int>(model_quantity));
    }

    if (logger_)
        logger_->log(
            LogSeverity::warning, msg_id,
            "AMICI encountered a %s value for %s[%i] (%s, %s) at t=%g",
            non_finite_type.c_str(), model_quantity_str.c_str(),
            gsl::narrow<int>(flat_index), row_id.c_str(), col_id.c_str(), t
        );

    // check upstream
    check_finite(state_.fixed_parameters, ModelQuantity::k, t);
    check_finite(state_.unscaled_parameters, ModelQuantity::p, t);
    check_finite(derived_state_.w_, ModelQuantity::w, t);

    return AMICI_RECOVERABLE_ERROR;
}

void Model::set_always_check_finite(bool alwaysCheck) {
    always_check_finite_ = alwaysCheck;
}

bool Model::get_always_check_finite() const { return always_check_finite_; }

void Model::fx0(realtype t, AmiVector& x) {
    std::ranges::fill(derived_state_.x_rdata_, 0.0);
    /* this function  also computes initial total abundances */
    fx0(derived_state_.x_rdata_.data(), t, state_.unscaled_parameters.data(),
        state_.fixed_parameters.data());
    fx_solver(x.data(), derived_state_.x_rdata_.data());
    ftotal_cl(
        state_.total_cl.data(), derived_state_.x_rdata_.data(),
        state_.unscaled_parameters.data(), state_.fixed_parameters.data()
    );

    check_finite(derived_state_.x_rdata_, ModelQuantity::x0_rdata, t);
}

void Model::fx0_fixedParameters(realtype t, AmiVector& x) {
    if (!get_reinitialize_fixed_parameter_initial_states())
        return;

    /* we transform to the unreduced states x_rdata and then apply
     x0_fixedparameters to (i) enable updates to states that were removed from
     conservation laws and (ii) be able to correctly compute total abundances
     after updating the state variables */
    fx_rdata(
        derived_state_.x_rdata_.data(), compute_x_pos(x),
        state_.total_cl.data(), state_.unscaled_parameters.data(),
        state_.fixed_parameters.data()
    );
    fx0_fixedParameters(
        derived_state_.x_rdata_.data(), t, state_.unscaled_parameters.data(),
        state_.fixed_parameters.data(),
        simulation_parameters_.reinitialization_state_idxs_sim
    );
    fx_solver(x.data(), derived_state_.x_rdata_.data());
    /* update total abundances */
    ftotal_cl(
        state_.total_cl.data(), derived_state_.x_rdata_.data(),
        state_.unscaled_parameters.data(), state_.fixed_parameters.data()
    );
}

void Model::fsx0(realtype t, AmiVectorArray& sx, AmiVector const& x) {
    /* this function  also computes initial total abundance sensitivities */
    realtype* stcl = nullptr;
    for (int ip = 0; ip < nplist(); ip++) {
        if (ncl() > 0)
            stcl = &state_.stotal_cl.at(plist(ip) * ncl());
        std::ranges::fill(derived_state_.sx_rdata_, 0.0);
        fsx0(
            derived_state_.sx_rdata_.data(), t, compute_x_pos(x),
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            plist(ip)
        );
        fsx_solver(sx.data(ip), derived_state_.sx_rdata_.data());
        fstotal_cl(
            stcl, derived_state_.sx_rdata_.data(), plist(ip),
            derived_state_.x_rdata_.data(), state_.unscaled_parameters.data(),
            state_.fixed_parameters.data(), state_.total_cl.data()
        );
    }
}

void Model::fsx0_fixedParameters(
    realtype t, AmiVectorArray& sx, AmiVector const& x
) {
    if (!get_reinitialize_fixed_parameter_initial_states())
        return;
    realtype* stcl = nullptr;
    for (int ip = 0; ip < nplist(); ip++) {
        if (ncl() > 0)
            stcl = &state_.stotal_cl.at(plist(ip) * ncl());
        fsx_rdata(
            derived_state_.sx_rdata_.data(), sx.data(ip), stcl,
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            x.data(), state_.total_cl.data(), plist(ip)
        );
        fsx0_fixedParameters(
            derived_state_.sx_rdata_.data(), t, compute_x_pos(x),
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            plist(ip), simulation_parameters_.reinitialization_state_idxs_sim
        );
        fsx_solver(sx.data(ip), derived_state_.sx_rdata_.data());
        fstotal_cl(
            stcl, derived_state_.sx_rdata_.data(), plist(ip),
            derived_state_.x_rdata_.data(), state_.unscaled_parameters.data(),
            state_.fixed_parameters.data(), state_.total_cl.data()
        );
    }
}

void Model::fsdx0() {}

void Model::fx_rdata(gsl::span<realtype> x_rdata, AmiVector const& x) {
    fx_rdata(
        x_rdata.data(), compute_x_pos(x), state_.total_cl.data(),
        state_.unscaled_parameters.data(), state_.fixed_parameters.data()
    );
    if (always_check_finite_)
        check_finite(
            x_rdata, ModelQuantity::x_rdata,
            std::numeric_limits<realtype>::quiet_NaN()
        );
}

void Model::fsx_rdata(
    gsl::span<realtype> sx_rdata, AmiVectorArray const& sx,
    AmiVector const& x_solver
) {
    realtype* stcl = nullptr;
    for (int ip = 0; ip < nplist(); ip++) {
        if (ncl() > 0)
            stcl = &state_.stotal_cl.at(plist(ip) * ncl());
        fsx_rdata(
            &sx_rdata[ip * nx_rdata], sx.data(ip), stcl,
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            x_solver.data(), state_.total_cl.data(), plist(ip)
        );
    }
}

void Model::write_slice_event(
    gsl::span<realtype const> slice, gsl::span<realtype> buffer, int const ie
) {
    check_buffer_size(buffer, slice.size());
    check_buffer_size(buffer, z2event_.size());
    for (unsigned izt = 0; izt < z2event_.size(); ++izt)
        if (z2event_.at(izt) == ie)
            buffer[izt] = slice[izt];
}

void Model::write_sensitivity_slice_event(
    gsl::span<realtype const> slice, gsl::span<realtype> buffer, int const ie
) {
    check_buffer_size(buffer, slice.size());
    check_buffer_size(buffer, z2event_.size() * nplist());
    for (int ip = 0; ip < nplist(); ++ip)
        for (unsigned izt = 0; izt < z2event_.size(); ++izt)
            if (z2event_.at(izt) == ie)
                buffer[ip * nztrue + izt] = slice[ip * nztrue + izt];
}

void Model::write_llh_sensitivity_slice(
    std::vector<realtype> const& dLLhdp, std::vector<realtype>& sllh,
    std::vector<realtype>& s2llh
) {
    check_llh_buffer_size(sllh, s2llh);

    amici_daxpy(nplist(), -1.0, dLLhdp.data(), nJ, sllh.data(), 1);
    for (int iJ = 1; iJ < nJ; ++iJ)
        amici_daxpy(
            nplist(), -1.0, &dLLhdp.at(iJ), nJ, &s2llh.at(iJ - 1), nJ - 1
        );
}

void Model::check_llh_buffer_size(
    std::vector<realtype> const& sllh, std::vector<realtype> const& s2llh
) const {
    if (sllh.size() != gsl::narrow<unsigned>(nplist()))
        throw AmiException(
            "Incorrect sllh buffer size! Was %u, expected %i.", sllh.size(),
            nplist()
        );

    if (s2llh.size() != gsl::narrow<unsigned>((nJ - 1) * nplist()))
        throw AmiException(
            "Incorrect s2llh buffer size! Was %u, expected %i.", s2llh.size(),
            (nJ - 1) * nplist()
        );
}

void Model::initialize_vectors() { sx0data_.clear(); }

void Model::fy(realtype const t, AmiVector const& x) {
    if (!ny)
        return;

    auto x_pos = compute_x_pos(x);

    derived_state_.y_.assign(ny, 0.0);

    fw(t, x_pos, false);
    fy(derived_state_.y_.data(), t, x_pos, state_.unscaled_parameters.data(),
       state_.fixed_parameters.data(), state_.h.data(),
       derived_state_.w_.data());

    if (always_check_finite_) {
        check_finite(
            gsl::make_span(derived_state_.y_.data(), ny), ModelQuantity::y, t
        );
    }
}

void Model::fdydp(realtype const t, AmiVector const& x) {
    if (!ny)
        return;

    auto x_pos = compute_x_pos(x);

    derived_state_.dydp_.assign(ny * nplist(), 0.0);
    fw(t, x_pos, false);
    fdwdp(t, x_pos, false);

    /* get dydp slice (ny) for current time and parameter */
    for (int ip = 0; ip < nplist(); ip++)
        fdydp(
            &derived_state_.dydp_.at(ip * ny), t, x_pos,
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            state_.h.data(), plist(ip), derived_state_.w_.data(),
            state_.total_cl.data(), state_.stotal_cl.data(),
            derived_state_.spl_.data(), derived_state_.sspl_.data()
        );

    if (always_check_finite_) {
        check_finite(derived_state_.dydp_, ModelQuantity::dydp, nplist());
    }
}

void Model::fdydx(realtype const t, AmiVector const& x) {
    if (!ny)
        return;

    auto x_pos = compute_x_pos(x);

    derived_state_.dydx_.assign(ny * nx_solver, 0.0);

    fw(t, x_pos, false);
    fdydx(
        derived_state_.dydx_.data(), t, x_pos,
        state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
        state_.h.data(), derived_state_.w_.data()
    );

    if (always_check_finite_) {
        check_finite(derived_state_.dydx_, ModelQuantity::dydx, ny);
    }
}

void Model::fsigmay(int const it, ExpData const* edata) {
    if (!ny)
        return;

    derived_state_.sigmay_.assign(ny, 0.0);

    fsigmay(
        derived_state_.sigmay_.data(), get_timepoint(it),
        state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
        derived_state_.y_.data()
    );

    if (edata) {
        auto sigmay_edata = edata->get_noise_scales_ptr(it);
        /* extract the value for the standard deviation from ExpData,
         * if the data value is NaN, use the parameter value */
        for (int iytrue = 0; iytrue < nytrue; iytrue++) {
            if (edata->is_set_noise_scale(it, iytrue))
                derived_state_.sigmay_.at(iytrue) = sigmay_edata[iytrue];

            /* TODO: when moving second order code to cpp, verify
             * that this is actually what we want
             */
            for (int iJ = 1; iJ < nJ; iJ++)
                derived_state_.sigmay_.at(iytrue + iJ * nytrue) = 0;

            if (edata->is_set_measurement(it, iytrue)) {
                std::string obs_id = has_observable_ids()
                                         ? get_observable_ids().at(iytrue)
                                         : std::to_string(iytrue);
                std::stringstream ss;
                ss << "sigmay (" << obs_id << ", ExpData::id=" << edata->id
                   << ", t=" << get_timepoint(it) << ")";
                check_sigma_positivity(
                    derived_state_.sigmay_.at(iytrue), ss.str().c_str()
                );
            }
        }
    }
}

void Model::fdsigmaydp(int const it, ExpData const* edata) {
    if (!ny)
        return;

    derived_state_.dsigmaydp_.assign(ny * nplist(), 0.0);

    for (int ip = 0; ip < nplist(); ip++)
        // get dsigmaydp slice (ny) for current timepoint and parameter
        fdsigmaydp(
            &derived_state_.dsigmaydp_.at(ip * ny), get_timepoint(it),
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            derived_state_.y_.data(), plist(ip)
        );

    // sigmas in edata override model-sigma -> for those sigmas, set dsigmaydp
    // to zero
    if (edata) {
        for (int iy = 0; iy < nytrue; iy++) {
            if (!edata->is_set_noise_scale(it, iy))
                continue;
            for (int ip = 0; ip < nplist(); ip++) {
                derived_state_.dsigmaydp_.at(ip * ny + iy) = 0.0;
            }
        }
    }

    if (always_check_finite_) {
        check_finite(
            derived_state_.dsigmaydp_, ModelQuantity::dsigmaydp, nplist()
        );
    }
}

void Model::fdsigmaydy(int const it, ExpData const* edata) {
    if (!ny)
        return;

    derived_state_.dsigmaydy_.assign(ny * ny, 0.0);

    // get dsigmaydy slice (ny) for current timepoint
    fdsigmaydy(
        derived_state_.dsigmaydy_.data(), get_timepoint(it),
        state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
        derived_state_.y_.data()
    );

    // sigmas in edata override model-sigma -> for those sigmas, set dsigmaydy
    // to zero
    if (edata) {
        for (int isigmay = 0; isigmay < nytrue; ++isigmay) {
            if (!edata->is_set_noise_scale(it, isigmay))
                continue;
            for (int iy = 0; iy < nytrue; ++iy) {
                derived_state_.dsigmaydy_.at(isigmay * ny + iy) = 0.0;
            }
        }
    }

    if (always_check_finite_) {
        check_finite(derived_state_.dsigmaydy_, ModelQuantity::dsigmaydy, ny);
    }
}

void Model::fdJydy(int const it, AmiVector const& x, ExpData const& edata) {
    if (!ny)
        return;

    fy(edata.get_timepoint(it), x);
    fsigmay(it, &edata);

    fdJydsigma(it, x, edata);
    fdsigmaydy(it, &edata);

    set_nan_to_zero(derived_state_.dJydsigma_);
    set_nan_to_zero(derived_state_.dsigmaydy_);
    for (int iyt = 0; iyt < nytrue; iyt++) {
        if (!derived_state_.dJydy_.at(iyt).capacity())
            continue;
        derived_state_.dJydy_.at(iyt).zero();
        fdJydy_colptrs(derived_state_.dJydy_.at(iyt), iyt);
        fdJydy_rowvals(derived_state_.dJydy_.at(iyt), iyt);

        if (!edata.is_set_measurement(it, iyt))
            continue;

        // get dJydy slice (ny) for current timepoint and observable
        fdJydy(
            derived_state_.dJydy_.at(iyt).data(), iyt,
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            derived_state_.y_.data(), derived_state_.sigmay_.data(),
            edata.get_measurements_ptr(it)
        );

        // dJydy += dJydsigma * dsigmaydy
        // C(nJ,ny)  A(nJ,ny)  * B(ny,ny)
        // sparse    dense       dense
        derived_state_.dJydy_dense_.zero();
        amici_dgemm(
            BLASLayout::colMajor, BLASTranspose::noTrans,
            BLASTranspose::noTrans, nJ, ny, ny, 1.0,
            &derived_state_.dJydsigma_.at(iyt * nJ * ny), nJ,
            derived_state_.dsigmaydy_.data(), ny, 1.0,
            derived_state_.dJydy_dense_.data(), nJ
        );

        auto tmp_sparse
            = SUNMatrixWrapper(derived_state_.dJydy_dense_, 0.0, CSC_MAT);
        auto ret
            = SUNMatScaleAdd(1.0, derived_state_.dJydy_.at(iyt), tmp_sparse);
        if (ret != SUN_SUCCESS) {
            throw AmiException(
                "SUNMatScaleAdd failed with status %d in %s", ret, __func__
            );
        }
        derived_state_.dJydy_.at(iyt).refresh();

        if (always_check_finite_) {
            check_finite(
                gsl::make_span(derived_state_.dJydy_.at(iyt).get()),
                ModelQuantity::dJydy, ny
            );
        }
    }
}

void Model::fdJydsigma(int const it, AmiVector const& x, ExpData const& edata) {
    if (!ny)
        return;

    derived_state_.dJydsigma_.assign(nytrue * ny * nJ, 0.0);

    fy(edata.get_timepoint(it), x);
    fsigmay(it, &edata);

    for (int iyt = 0; iyt < nytrue; iyt++) {
        if (edata.is_set_measurement(it, iyt)) {
            // get dJydsigma slice (ny) for current timepoint and observable
            fdJydsigma(
                &derived_state_.dJydsigma_.at(iyt * ny * nJ), iyt,
                state_.unscaled_parameters.data(),
                state_.fixed_parameters.data(), derived_state_.y_.data(),
                derived_state_.sigmay_.data(), edata.get_measurements_ptr(it)
            );
            if (always_check_finite_) {
                check_finite(
                    gsl::span<realtype>(
                        &derived_state_.dJydsigma_.at(iyt * ny * nJ), ny * nJ
                    ),
                    ModelQuantity::dJydsigma, ny
                );
            }
        }
    }
}

void Model::fdJydp(int const it, AmiVector const& x, ExpData const& edata) {
    // dJydy         nJ, nytrue x ny
    // dydp          nplist * ny
    // dJydp         nplist x nJ
    if (!ny)
        return;

    derived_state_.dJydp_.assign(nJ * nplist(), 0.0);

    fdJydy(it, x, edata);
    fdydp(edata.get_timepoint(it), x);

    fdJydsigma(it, x, edata);
    fdsigmaydp(it, &edata);

    set_nan_to_zero(derived_state_.dJydsigma_);
    set_nan_to_zero(derived_state_.dsigmaydp_);
    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (!edata.is_set_measurement(it, iyt))
            continue;

        // dJydp = 1.0 * dJydp +  1.0 * dJydy * dydp
        for (int iplist = 0; iplist < nplist(); ++iplist) {
            derived_state_.dJydy_.at(iyt).multiply(
                gsl::span<realtype>(&derived_state_.dJydp_.at(iplist * nJ), nJ),
                gsl::span<realtype const>(
                    &derived_state_.dydp_.at(iplist * ny), ny
                )
            );
        }

        // dJydp = 1.0 * dJydp +  1.0 * dJydsigma * dsigmaydp
        amici_dgemm(
            BLASLayout::colMajor, BLASTranspose::noTrans,
            BLASTranspose::noTrans, nJ, nplist(), ny, 1.0,
            &derived_state_.dJydsigma_.at(iyt * nJ * ny), nJ,
            derived_state_.dsigmaydp_.data(), ny, 1.0,
            derived_state_.dJydp_.data(), nJ
        );
    }
}

void Model::fdJydx(int const it, AmiVector const& x, ExpData const& edata) {
    if (!ny)
        return;

    derived_state_.dJydx_.assign(nJ * nx_solver, 0.0);

    fdydx(edata.get_timepoint(it), x);
    fdJydy(it, x, edata);

    // dJydy: nJ, ny x nytrue
    // dydx :     ny x nx_solver
    // dJydx:     nJ x nx_solver x nt
    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (!edata.is_set_measurement(it, iyt))
            continue;
        // dJydy A[nyt,nJ,ny] * dydx B[ny,nx_solver] = dJydx C[it,nJ,nx_solver]
        //         slice                                       slice
        //          M  K            K  N                       M  N
        //           lda             ldb                        ldc

        for (int ix = 0; ix < nx_solver; ++ix) {
            derived_state_.dJydy_.at(iyt).multiply(
                gsl::span<realtype>(&derived_state_.dJydx_.at(ix * nJ), nJ),
                gsl::span<realtype const>(&derived_state_.dydx_.at(ix * ny), ny)
            );
        }
    }

    if (always_check_finite_) {
        check_finite(derived_state_.dJydx_, ModelQuantity::dJydx, nx_solver);
    }
}

void Model::fz(int const ie, realtype const t, AmiVector const& x) {

    derived_state_.z_.assign(nz, 0.0);

    fz(derived_state_.z_.data(), ie, t, compute_x_pos(x),
       state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
       state_.h.data());
}

void Model::fdzdp(int const ie, realtype const t, AmiVector const& x) {
    if (!nz)
        return;

    derived_state_.dzdp_.assign(nz * nplist(), 0.0);

    for (int ip = 0; ip < nplist(); ip++) {
        fdzdp(
            derived_state_.dzdp_.data(), ie, t, compute_x_pos(x),
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            state_.h.data(), plist(ip)
        );
    }

    if (always_check_finite_) {
        check_finite(derived_state_.dzdp_, ModelQuantity::dzdp, nplist());
    }
}

void Model::fdzdx(int const ie, realtype const t, AmiVector const& x) {
    if (!nz)
        return;

    derived_state_.dzdx_.assign(nz * nx_solver, 0.0);

    fdzdx(
        derived_state_.dzdx_.data(), ie, t, compute_x_pos(x),
        state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
        state_.h.data()
    );

    if (always_check_finite_) {
        check_finite(derived_state_.dzdx_, ModelQuantity::dzdx, nx_solver);
    }
}

void Model::frz(int const ie, realtype const t, AmiVector const& x) {

    derived_state_.rz_.assign(nz, 0.0);

    frz(derived_state_.rz_.data(), ie, t, compute_x_pos(x),
        state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
        state_.h.data());
}

void Model::fdrzdp(int const ie, realtype const t, AmiVector const& x) {
    if (!nz)
        return;

    derived_state_.drzdp_.assign(nz * nplist(), 0.0);

    for (int ip = 0; ip < nplist(); ip++) {
        fdrzdp(
            derived_state_.drzdp_.data(), ie, t, compute_x_pos(x),
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            state_.h.data(), plist(ip)
        );
    }

    if (always_check_finite_) {
        check_finite(derived_state_.drzdp_, ModelQuantity::drzdp, nplist());
    }
}

void Model::fdrzdx(int const ie, realtype const t, AmiVector const& x) {
    if (!nz)
        return;

    derived_state_.drzdx_.assign(nz * nx_solver, 0.0);

    fdrzdx(
        derived_state_.drzdx_.data(), ie, t, compute_x_pos(x),
        state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
        state_.h.data()
    );

    if (always_check_finite_) {
        check_finite(derived_state_.drzdx_, ModelQuantity::drzdx, nx_solver);
    }
}

void Model::fsigmaz(
    int const ie, int const nroots, realtype const t, ExpData const* edata
) {
    if (!nz)
        return;

    derived_state_.sigmaz_.assign(nz, 0.0);
    fsigmaz(
        derived_state_.sigmaz_.data(), t, state_.unscaled_parameters.data(),
        state_.fixed_parameters.data()
    );

    if (edata) {
        for (int iztrue = 0; iztrue < nztrue; iztrue++) {
            if (z2event_.at(iztrue) == ie) {
                if (edata->is_set_event_noise_scale(nroots, iztrue)) {
                    auto sigmaz_edata
                        = edata->get_event_noise_scales_ptr(nroots);
                    derived_state_.sigmaz_.at(iztrue) = sigmaz_edata[iztrue];
                }

                /* TODO: when moving second order code to cpp, verify
                 * that this is actually what we want
                 */
                for (int iJ = 1; iJ < nJ; iJ++)
                    derived_state_.sigmaz_.at(iztrue + iJ * nztrue) = 0;

                if (edata->is_set_event_measurement(nroots, iztrue))
                    check_sigma_positivity(
                        derived_state_.sigmaz_.at(iztrue), "sigmaz"
                    );
            }
        }
    }
}

void Model::fdsigmazdp(
    int const ie, int const nroots, realtype const t, ExpData const* edata
) {
    if (!nz)
        return;

    derived_state_.dsigmazdp_.assign(nz * nplist(), 0.0);

    for (int ip = 0; ip < nplist(); ip++) {
        // get dsigmazdp slice (nz) for current event and parameter
        fdsigmazdp(
            &derived_state_.dsigmazdp_.at(ip * nz), t,
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            plist(ip)
        );
    }

    // sigmas in edata override model-sigma -> for those sigmas, set dsigmazdp
    // to zero
    if (edata) {
        for (int iz = 0; iz < nztrue; iz++) {
            if (z2event_.at(iz) == ie
                && !edata->is_set_event_noise_scale(nroots, iz)) {
                for (int ip = 0; ip < nplist(); ip++)
                    derived_state_.dsigmazdp_.at(iz + nz * ip) = 0;
            }
        }
    }

    if (always_check_finite_) {
        check_finite(
            derived_state_.dsigmazdp_, ModelQuantity::dsigmazdp, nplist()
        );
    }
}

void Model::fdJzdz(
    int const ie, int const nroots, realtype const t, AmiVector const& x,
    ExpData const& edata
) {
    if (!nz)
        return;

    derived_state_.dJzdz_.assign(nztrue * nz * nJ, 0.0);

    fz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.is_set_event_measurement(nroots, iztrue)) {
            fdJzdz(
                &derived_state_.dJzdz_.at(iztrue * nz * nJ), iztrue,
                state_.unscaled_parameters.data(),
                state_.fixed_parameters.data(), derived_state_.z_.data(),
                derived_state_.sigmaz_.data(),
                edata.get_event_measurements_ptr(nroots)
            );
            if (always_check_finite_) {
                check_finite(
                    gsl::span<realtype>(
                        &derived_state_.dJzdz_.at(iztrue * nz * nJ), nz * nJ
                    ),
                    ModelQuantity::dJzdz, nz
                );
            }
        }
    }
}

void Model::fdJzdsigma(
    int const ie, int const nroots, realtype const t, AmiVector const& x,
    ExpData const& edata
) {
    if (!nz)
        return;

    derived_state_.dJzdsigma_.assign(nztrue * nz * nJ, 0.0);

    fz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.is_set_event_measurement(nroots, iztrue)) {
            fdJzdsigma(
                &derived_state_.dJzdsigma_.at(iztrue * nz * nJ), iztrue,
                state_.unscaled_parameters.data(),
                state_.fixed_parameters.data(), derived_state_.z_.data(),
                derived_state_.sigmaz_.data(),
                edata.get_event_measurements_ptr(nroots)
            );
            if (always_check_finite_) {
                check_finite(
                    gsl::span<realtype>(
                        &derived_state_.dJzdsigma_.at(iztrue * nz * nJ), nz * nJ
                    ),
                    ModelQuantity::dJzdsigma, nz
                );
            }
        }
    }
}

void Model::fdJzdp(
    int const ie, int const nroots, realtype t, AmiVector const& x,
    ExpData const& edata
) {
    if (!nz)
        return;
    // dJzdz         nJ x nz x nztrue
    // dJzdsigma     nJ x nz x nztrue
    // dzdp          nz x nplist()
    // dJzdp         nJ x nplist()

    derived_state_.dJzdp_.assign(nJ * nplist(), 0.0);

    fdzdp(ie, t, x);
    fdsigmazdp(ie, nroots, t, &edata);
    fdJzdz(ie, nroots, t, x, edata);
    fdJrzdz(ie, nroots, t, x, edata);
    fdJzdsigma(ie, nroots, t, x, edata);
    fdJrzdsigma(ie, nroots, t, x, edata);

    set_nan_to_zero(derived_state_.dzdp_);
    set_nan_to_zero(derived_state_.dsigmazdp_);
    set_nan_to_zero(derived_state_.dJzdz_);
    set_nan_to_zero(derived_state_.dJrzdz_);
    set_nan_to_zero(derived_state_.dJzdsigma_);
    set_nan_to_zero(derived_state_.dJrzdsigma_);
    for (int izt = 0; izt < nztrue; ++izt) {
        if (!edata.is_set_event_measurement(nroots, izt))
            continue;

        if (t < edata.get_timepoint(edata.nt() - 1)) {
            // with z
            amici_dgemm(
                BLASLayout::colMajor, BLASTranspose::noTrans,
                BLASTranspose::noTrans, nJ, nplist(), nz, 1.0,
                &derived_state_.dJzdz_.at(izt * nz * nJ), nJ,
                derived_state_.dzdp_.data(), nz, 1.0,
                derived_state_.dJzdp_.data(), nJ
            );
        } else {
            // with rz
            amici_dgemm(
                BLASLayout::colMajor, BLASTranspose::noTrans,
                BLASTranspose::noTrans, nJ, nplist(), nz, 1.0,
                &derived_state_.dJrzdsigma_.at(izt * nz * nJ), nJ,
                derived_state_.dsigmazdp_.data(), nz, 1.0,
                derived_state_.dJzdp_.data(), nJ
            );

            amici_dgemm(
                BLASLayout::colMajor, BLASTranspose::noTrans,
                BLASTranspose::noTrans, nJ, nplist(), nz, 1.0,
                &derived_state_.dJrzdz_.at(izt * nz * nJ), nJ,
                derived_state_.dzdp_.data(), nz, 1.0,
                derived_state_.dJzdp_.data(), nJ
            );
        }

        amici_dgemm(
            BLASLayout::colMajor, BLASTranspose::noTrans,
            BLASTranspose::noTrans, nJ, nplist(), nz, 1.0,
            &derived_state_.dJzdsigma_.at(izt * nz * nJ), nJ,
            derived_state_.dsigmazdp_.data(), nz, 1.0,
            derived_state_.dJzdp_.data(), nJ
        );
    }
}

void Model::fdJzdx(
    int const ie, int const nroots, realtype const t, AmiVector const& x,
    ExpData const& edata
) {
    // dJzdz         nJ x nz        x nztrue
    // dzdx          nz x nx_solver
    // dJzdx         nJ x nx_solver x nmaxevent
    if (!nz)
        return;

    derived_state_.dJzdx_.assign(nJ * nx_solver, 0.0);

    fdJzdz(ie, nroots, t, x, edata);
    fdJrzdz(ie, nroots, t, x, edata);
    fdzdx(ie, t, x);
    fdrzdx(ie, t, x);

    set_nan_to_zero(derived_state_.dJzdz_);
    set_nan_to_zero(derived_state_.dJrzdz_);
    set_nan_to_zero(derived_state_.dzdx_);
    set_nan_to_zero(derived_state_.drzdx_);

    for (int izt = 0; izt < nztrue; ++izt) {
        if (!edata.is_set_event_measurement(nroots, izt))
            continue;

        if (t < edata.get_timepoint(edata.nt() - 1)) {
            // z
            amici_dgemm(
                BLASLayout::colMajor, BLASTranspose::noTrans,
                BLASTranspose::noTrans, nJ, nx_solver, nz, 1.0,
                &derived_state_.dJzdz_.at(izt * nz * nJ), nJ,
                derived_state_.dzdx_.data(), nz, 1.0,
                derived_state_.dJzdx_.data(), nJ
            );
        } else {
            // rz
            amici_dgemm(
                BLASLayout::colMajor, BLASTranspose::noTrans,
                BLASTranspose::noTrans, nJ, nx_solver, nz, 1.0,
                &derived_state_.dJrzdz_.at(izt * nz * nJ), nJ,
                derived_state_.drzdx_.data(), nz, 1.0,
                derived_state_.dJzdx_.data(), nJ
            );
        }
    }
}

void Model::fdJrzdz(
    int const ie, int const nroots, realtype const t, AmiVector const& x,
    ExpData const& edata
) {
    if (!nz)
        return;

    derived_state_.dJrzdz_.assign(nztrue * nz * nJ, 0.0);

    frz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.is_set_event_measurement(nroots, iztrue)) {
            fdJrzdz(
                &derived_state_.dJrzdz_.at(iztrue * nz * nJ), iztrue,
                state_.unscaled_parameters.data(),
                state_.fixed_parameters.data(), derived_state_.rz_.data(),
                derived_state_.sigmaz_.data()
            );
            if (always_check_finite_) {
                check_finite(
                    gsl::span<realtype>(
                        &derived_state_.dJrzdz_.at(iztrue * nz * nJ), nz * nJ
                    ),
                    ModelQuantity::dJrzdz, nz
                );
            }
        }
    }
}

void Model::fdJrzdsigma(
    int const ie, int const nroots, realtype const t, AmiVector const& x,
    ExpData const& edata
) {
    if (!nz)
        return;

    derived_state_.dJrzdsigma_.assign(nztrue * nz * nJ, 0.0);

    frz(ie, t, x);
    fsigmaz(ie, nroots, t, &edata);

    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (edata.is_set_event_measurement(nroots, iztrue)) {
            fdJrzdsigma(
                &derived_state_.dJrzdsigma_.at(iztrue * nz * nJ), iztrue,
                state_.unscaled_parameters.data(),
                state_.fixed_parameters.data(), derived_state_.rz_.data(),
                derived_state_.sigmaz_.data()
            );
            if (always_check_finite_) {
                check_finite(
                    gsl::span<realtype>(
                        &derived_state_.dJrzdsigma_.at(iztrue * nz * nJ),
                        nz * nJ
                    ),
                    ModelQuantity::dJrzdsigma, nz
                );
            }
        }
    }
}

void Model::fspl(realtype const t) {
    for (int ispl = 0; ispl < nspl; ispl++)
        derived_state_.spl_[ispl] = splines_[ispl].get_value(t);
}

void Model::fsspl(realtype const t) {
    derived_state_.sspl_.zero();
    realtype* sspl_data = derived_state_.sspl_.data();
    for (int ip = 0; ip < nplist(); ip++) {
        for (int ispl = 0; ispl < nspl; ispl++)
            sspl_data[ispl + nspl * plist(ip)] = splines_[ispl].get_sensitivity(
                t, ip, derived_state_.spl_[ispl]
            );
    }
}

void Model::fw(realtype const t, realtype const* x, bool include_static) {
    if (include_static) {
        std::ranges::fill(derived_state_.w_, 0.0);
    }
    fspl(t);
    fw(derived_state_.w_.data(), t, x, state_.unscaled_parameters.data(),
       state_.fixed_parameters.data(), state_.h.data(), state_.total_cl.data(),
       derived_state_.spl_.data(), include_static);

    if (always_check_finite_) {
        check_finite(derived_state_.w_, ModelQuantity::w, t);
    }
}

void Model::fdwdp(realtype const t, realtype const* x, bool include_static) {
    if (!nw)
        return;

    fw(t, x, include_static);
    derived_state_.dwdp_.zero();
    if (!ndwdp)
        return;
    fsspl(t);
    fdwdw(t, x, include_static);
    if (include_static) {
        derived_state_.dwdp_hierarchical_.at(0).zero();
        fdwdp_colptrs(derived_state_.dwdp_hierarchical_.at(0));
        fdwdp_rowvals(derived_state_.dwdp_hierarchical_.at(0));
    }
    fdwdp(
        derived_state_.dwdp_hierarchical_.at(0).data(), t, x,
        state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
        state_.h.data(), derived_state_.w_.data(), state_.total_cl.data(),
        state_.stotal_cl.data(), derived_state_.spl_.data(),
        derived_state_.sspl_.data(), include_static
    );

    for (int irecursion = 1; irecursion <= w_recursion_depth; irecursion++) {
        derived_state_.dwdw_.sparse_multiply(
            derived_state_.dwdp_hierarchical_.at(irecursion),
            derived_state_.dwdp_hierarchical_.at(irecursion - 1)
        );
    }
    derived_state_.dwdp_.sparse_sum(derived_state_.dwdp_hierarchical_);

    if (always_check_finite_) {
        check_finite(derived_state_.dwdp_, ModelQuantity::dwdp, t);
    }
}

void Model::fdwdx(realtype const t, realtype const* x, bool include_static) {
    // NOTE: (at least) Model_{ODE,DAE}::fJSparse rely on `fw` and `fdwdw`
    //  being called from here. They need to be executed even if nx_solver==0.
    if (!nw)
        return;

    fw(t, x, include_static);

    derived_state_.dwdx_.zero();
    fdwdw(t, x, include_static);

    auto&& dwdx_hierarchical_0 = derived_state_.dwdx_hierarchical_.at(0);
    if (!dwdx_hierarchical_0.data() || !dwdx_hierarchical_0.capacity())
        return;

    if (include_static) {
        dwdx_hierarchical_0.zero();
        fdwdx_colptrs(dwdx_hierarchical_0);
        fdwdx_rowvals(dwdx_hierarchical_0);
    }
    fdwdx(
        dwdx_hierarchical_0.data(), t, x, state_.unscaled_parameters.data(),
        state_.fixed_parameters.data(), state_.h.data(),
        derived_state_.w_.data(), state_.total_cl.data(),
        derived_state_.spl_.data(), include_static
    );

    for (int irecursion = 1; irecursion <= w_recursion_depth; irecursion++) {
        derived_state_.dwdw_.sparse_multiply(
            derived_state_.dwdx_hierarchical_.at(irecursion),
            derived_state_.dwdx_hierarchical_.at(irecursion - 1)
        );
    }
    derived_state_.dwdx_.sparse_sum(derived_state_.dwdx_hierarchical_);

    if (always_check_finite_) {
        check_finite(derived_state_.dwdx_, ModelQuantity::dwdx, t);
    }
}

void Model::fdwdw(realtype const t, realtype const* x, bool include_static) {
    if (!ndwdw)
        return;

    if (include_static) {
        derived_state_.dwdw_.zero();
        fdwdw_colptrs(derived_state_.dwdw_);
        fdwdw_rowvals(derived_state_.dwdw_);
    }

    fdwdw(
        derived_state_.dwdw_.data(), t, x, state_.unscaled_parameters.data(),
        state_.fixed_parameters.data(), state_.h.data(),
        derived_state_.w_.data(), state_.total_cl.data(), include_static
    );

    if (always_check_finite_) {
        check_finite(derived_state_.dwdw_, ModelQuantity::dwdw, t);
    }
}

void Model::fx_rdata(
    realtype* x_rdata, realtype const* x_solver, realtype const* /*tcl*/,
    realtype const* /*p*/, realtype const* /*k*/
) {
    if (nx_solver != nx_rdata)
        throw AmiException(
            "A model that has differing nx_solver and nx_rdata needs "
            "to implement its own fx_rdata"
        );
    std::copy_n(x_solver, nx_solver, x_rdata);
}

void Model::fsx_rdata(
    realtype* sx_rdata, realtype const* sx_solver, realtype const* stcl,
    realtype const* p, realtype const* k, realtype const* x_solver,
    realtype const* tcl, int const ip
) {
    if (nx_solver == nx_rdata) {
        std::copy_n(sx_solver, nx_solver, sx_rdata);
        return;
    }

    // sx_rdata = dx_rdata/dx_solver * sx_solver
    //             + dx_rdata/d_tcl * stcl + dxrdata/dp

    // 1) sx_rdata(nx_rdata, 1) = dx_rdatadp
    std::fill_n(sx_rdata, nx_rdata, 0.0);
    fdx_rdatadp(sx_rdata, x_solver, tcl, p, k, ip);

    // the following could be moved to the calling function, as it's independent
    //  of `ip`

    // 2) sx_rdata(nx_rdata, 1) +=
    //          dx_rdata/dx_solver(nx_rdata,nx_solver) * sx_solver(nx_solver, 1)
    derived_state_.dx_rdatadx_solver.zero();
    fdx_rdatadx_solver(
        derived_state_.dx_rdatadx_solver.data(), x_solver, tcl, p, k
    );
    fdx_rdatadx_solver_colptrs(derived_state_.dx_rdatadx_solver);
    fdx_rdatadx_solver_rowvals(derived_state_.dx_rdatadx_solver);
    derived_state_.dx_rdatadx_solver.multiply(
        gsl::make_span(sx_rdata, nx_rdata), gsl::make_span(sx_solver, nx_solver)
    );

    // 3) sx_rdata(nx_rdata, 1) += dx_rdata/d_tcl(nx_rdata,ntcl) * stcl
    derived_state_.dx_rdatadtcl.zero();
    fdx_rdatadtcl(derived_state_.dx_rdatadtcl.data(), x_solver, tcl, p, k);
    fdx_rdatadtcl_colptrs(derived_state_.dx_rdatadtcl);
    fdx_rdatadtcl_rowvals(derived_state_.dx_rdatadtcl);
    derived_state_.dx_rdatadtcl.multiply(
        gsl::make_span(sx_rdata, nx_rdata),
        gsl::make_span(stcl, (nx_rdata - nx_solver))
    );
}

void Model::fx_solver(realtype* x_solver, realtype const* x_rdata) {
    if (nx_solver != nx_rdata)
        throw AmiException(
            "A model that has differing nx_solver and nx_rdata needs "
            "to implement its own fx_solver"
        );
    std::copy_n(x_rdata, nx_rdata, x_solver);
}

void Model::fsx_solver(realtype* sx_solver, realtype const* sx_rdata) {
    /* for the moment we do not need an implementation of fsx_solver as
     * we can simply reuse fx_solver and replace states by their
     * sensitivities */
    fx_solver(sx_solver, sx_rdata);
}

void Model::ftotal_cl(
    realtype* /*total_cl*/, realtype const* /*x_rdata*/, realtype const* /*p*/,
    realtype const* /*k*/
) {
    if (nx_solver != nx_rdata)
        throw AmiException(
            "A model that has differing nx_solver and nx_rdata needs "
            "to implement its own ftotal_cl"
        );
}

void Model::fstotal_cl(
    realtype* stotal_cl, realtype const* sx_rdata, int const ip,
    realtype const* x_rdata, realtype const* p, realtype const* k,
    realtype const* tcl
) {
    if (nx_solver == nx_rdata)
        return;

    // stotal_cl(ncl,1) =
    //            dtotal_cl/dp(ncl,1)
    //          + dtotal_cl/dx_rdata(ncl,nx_rdata) * sx_rdata(nx_rdata,1)

    // 1) stotal_cl = dtotal_cl/dp
    std::fill_n(stotal_cl, ncl(), 0.0);
    fdtotal_cldp(stotal_cl, x_rdata, p, k, ip);

    // 2) stotal_cl += dtotal_cl/dx_rdata(ncl,nx_rdata) * sx_rdata(nx_rdata,1)
    derived_state_.dtotal_cldx_rdata.zero();
    fdtotal_cldx_rdata(
        derived_state_.dtotal_cldx_rdata.data(), x_rdata, p, k, tcl
    );
    fdtotal_cldx_rdata_colptrs(derived_state_.dtotal_cldx_rdata);
    fdtotal_cldx_rdata_rowvals(derived_state_.dtotal_cldx_rdata);
    derived_state_.dtotal_cldx_rdata.multiply(
        gsl::make_span(stotal_cl, ncl()), gsl::make_span(sx_rdata, nx_rdata)
    );
}

std::vector<double> Model::get_trigger_timepoints() const {
    std::vector<double> trigger_timepoints(explicit_roots_.size(), 0.0);
    // collect keys from state_independent_events_ which are the trigger
    // timepoints
    auto it = trigger_timepoints.begin();
    for (auto const& kv : explicit_roots_) {
        *(it++) = kv.first;
    }
    std::ranges::sort(trigger_timepoints);
    return trigger_timepoints;
}

void Model::set_steadystate_mask(std::vector<realtype> const& mask) {
    if (mask.empty()) {
        steadystate_mask_.clear();
        return;
    }

    if (gsl::narrow<int>(mask.size()) != nx_solver)
        throw AmiException(
            "Steadystate mask has wrong size: %d, expected %d",
            gsl::narrow<int>(mask.size()), nx_solver
        );

    steadystate_mask_ = mask;
}

const_N_Vector Model::compute_x_pos(const_N_Vector x) {
    if (any_state_non_negative_) {
        for (int ix = 0; ix < derived_state_.x_pos_tmp_.size(); ++ix) {
            derived_state_.x_pos_tmp_.at(ix)
                = (state_is_non_negative_.at(ix) && NV_Ith_S(x, ix) < 0)
                      ? 0
                      : NV_Ith_S(x, ix);
        }
        return derived_state_.x_pos_tmp_.get_nvector();
    }

    return x;
}

realtype const* Model::compute_x_pos(AmiVector const& x) {
    if (any_state_non_negative_) {
        compute_x_pos(x.get_nvector());
        return derived_state_.x_pos_tmp_.data();
    }
    return x.data();
}

void Model::set_reinitialization_state_idxs(std::vector<int> const& idxs) {
    for (auto idx : idxs) {
        if (idx < 0 || idx >= nx_rdata)
            throw AmiException("Invalid state index given: %d", idx);
    }

    simulation_parameters_.reinitialization_state_idxs_sim = idxs;
}

std::vector<int> const& Model::get_reinitialization_state_idxs() const {
    return simulation_parameters_.reinitialization_state_idxs_sim;
}

SUNMatrixWrapper const& Model::get_dxdotdp_full() const {
    return derived_state_.dxdotdp_full;
}

} // namespace amici
