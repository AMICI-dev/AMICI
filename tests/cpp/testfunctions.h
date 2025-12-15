#ifndef TESTFUNCTIONS_H
#define TESTFUNCTIONS_H

#include <amici/amici.h>
#include <amici/hdf5.h>

#include <H5Cpp.h>

#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace amici {

class ReturnData;
class ExpData;

#if !defined(NEW_OPTION_FILE) || !defined(HDFFILE) || !defined(HDFFILEWRITE)
#error "Must define NEW_OPTION_FILE HDFFILE HDFFILEWRITE"
#endif

#define TEST_ATOL 1e-10
#define TEST_RTOL 1e-05

/**
 * @brief helper function to initialize default names/ids
 * @param name name of variables
 * @param length number of variables
 * @return default names/ids
 */
std::span<std::string_view const>
getVariableNames(std::string const& name, int length);

/**
 * @brief The Model_Test class is a model-unspecific implementation
 of model designed for unit testing.
 */
class Model_Test : public Model {
  public:
    /** constructor with model dimensions
     * @param model_dimensions ModelDimensions model_dimensions,
     * @param o2mode second order sensitivity mode
     * @param simulation_parameters Simulation parameters
     * @param idlist indexes indicating algebraic components (DAE only)
     * @param z2event mapping of event outputs to events
     */
    Model_Test(
        ModelDimensions const& model_dimensions,
        SimulationParameters simulation_parameters,
        SecondOrderMode const o2mode, std::vector<realtype> const idlist,
        std::vector<int> const z2event, std::vector<Event> const events
    )
        : Model(
              model_dimensions, simulation_parameters, o2mode, idlist, z2event,
              events
          ) {}

    /** default constructor */
    Model_Test()
        : Model(
              ModelDimensions(), SimulationParameters(), SecondOrderMode::none,
              std::vector<realtype>(), std::vector<int>()
          ) {}

    Model* clone() const override { return new Model_Test(*this); }

    std::unique_ptr<Solver> create_solver() override {
        throw AmiException("not implemented");
    }
    void froot(
        realtype const /*t*/, AmiVector const& /*x*/, AmiVector const& /*dx*/,
        gsl::span<realtype> /*root*/
    ) override {
        throw AmiException("not implemented");
    }
    void fxdot(
        realtype const /*t*/, AmiVector const& /*x*/, AmiVector const& /*dx*/,
        AmiVector& /*xdot*/
    ) override {
        throw AmiException("not implemented");
    }
    void fsxdot(
        realtype const /*t*/, AmiVector const& /*x*/, AmiVector const& /*dx*/,
        int const /*ip*/, AmiVector const& /*sx*/, AmiVector const& /*sdx*/,
        AmiVector& /*sxdot*/
    ) override {
        throw AmiException("not implemented");
    }
    void
    fJ(realtype const /*t*/, realtype const /*cj*/, AmiVector const& /*x*/,
       AmiVector const& /*dx*/, AmiVector const& /*xdot*/,
       SUNMatrix /*J*/) override {
        throw AmiException("not implemented");
    }
    void fJSparse(
        realtype const /*t*/, realtype const /*cj*/, AmiVector const& /*x*/,
        AmiVector const& /*dx*/, AmiVector const& /*xdot*/, SUNMatrix /*J*/
    ) override {
        throw AmiException("not implemented");
    }
    void
    fJB(realtype const /*t*/, realtype /*cj*/, AmiVector const& /*x*/,
        AmiVector const& /*dx*/, AmiVector const& /*xB*/,
        AmiVector const& /*dxB*/, AmiVector const& /*xBdot*/,
        SUNMatrix /*JB*/) override {
        throw AmiException("not implemented");
    }
    void fJSparseB(
        realtype const /*t*/, realtype /*cj*/, AmiVector const& /*x*/,
        AmiVector const& /*dx*/, AmiVector const& /*xB*/,
        AmiVector const& /*dxB*/, AmiVector const& /*xBdot*/, SUNMatrix /*JB*/
    ) override {
        throw AmiException("not implemented");
    }
    void fJDiag(
        realtype const /*t*/, AmiVector& /*Jdiag*/, realtype const /*cj*/,
        AmiVector const& /*x*/, AmiVector const& /*dx*/
    ) override {
        throw AmiException("not implemented");
    }
    void fdxdotdp(
        realtype const /*t*/, AmiVector const& /*x*/, AmiVector const& /*dx*/
    ) override {
        throw AmiException("not implemented");
    }
    void
    fJv(realtype const /*t*/, AmiVector const& /*x*/, AmiVector const& /*dx*/,
        AmiVector const& /*xdot*/, AmiVector const& /*v*/, AmiVector& /*nJv*/,
        realtype const /*cj*/) override {
        throw AmiException("not implemented");
    }

    void fxBdot_ss(
        realtype const /*t*/, AmiVector const& /*xB*/, AmiVector const& /*dxB*/,
        AmiVector& /*xBdot*/
    ) override {
        throw AmiException("not implemented");
    }

    void fJSparseB_ss(SUNMatrix /*JB*/) override {
        throw AmiException("not implemented");
    }

    void write_steady_state_JB(
        realtype const /*t*/, realtype /*cj*/, AmiVector const& /*x*/,
        AmiVector const& /*dx*/, AmiVector const& /*xB*/,
        AmiVector const& /*dxB*/, AmiVector const& /*xBdot*/
    ) override {
        throw AmiException("not implemented");
    }

    std::span<std::string_view const>
    get_free_parameter_names() const override {
        return getVariableNames("p", np());
    }

    std::span<std::string_view const> get_state_names() const override {
        return getVariableNames("x", nx_rdata);
    }

    std::span<std::string_view const>
    get_fixed_parameter_names() const override {
        return getVariableNames("k", nk());
    }

    std::span<std::string_view const> get_observable_names() const override {
        return getVariableNames("y", ny);
    }

    std::span<std::string_view const> get_free_parameter_ids() const override {
        return getVariableNames("p", np());
    }

    std::span<std::string_view const> get_state_ids() const override {
        return getVariableNames("x", nx_rdata);
    }

    std::span<std::string_view const> get_fixed_parameter_ids() const override {
        return getVariableNames("k", nk());
    }

    std::span<std::string_view const> get_observable_ids() const override {
        return getVariableNames("y", ny);
    }
};

void simulateWithDefaultOptions();

void simulateVerifyWrite(std::string const& path);

void simulateVerifyWrite(std::string const& path, double atol, double rtol);

void simulateVerifyWrite(
    std::string const& hdffileOptions, std::string const& hdffileResults,
    std::string const& hdffilewrite, std::string const& path, double atol,
    double rtol
);

std::unique_ptr<ExpData> getTestExpData(Model const& model);

bool withinTolerance(
    double expected, double actual, double atol, double rtol, int index,
    std::string_view name
);

void checkEqualArray(
    double const* expected, double const* actual, int length, double atol,
    double rtol, std::string_view name
);

void checkEqualArray(
    std::vector<double> const& expected, std::vector<double> const& actual,
    double atol, double rtol, std::string const& name
);

void verifyReturnData(
    std::string const& hdffile, std::string const& resultPath,
    ReturnData const* rdata, Model const* model, double atol, double rtol
);

void verifyReturnDataSensitivities(
    const H5::H5File& file, std::string const& resultPath,
    ReturnData const* rdata, Model const* model, double atol, double rtol
);

void printBacktrace(int depth);

} // namespace amici

#endif
