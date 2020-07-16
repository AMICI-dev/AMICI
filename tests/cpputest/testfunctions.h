#ifndef TESTFUNCTIONS_H
#define TESTFUNCTIONS_H

#include <amici/hdf5.h>
#include <amici/amici.h>

#include <H5Cpp.h>

#ifndef __APPLE__
#include <iostream>
#endif
#include <string>
#include <sstream>    // make std::ostringstream available (needs to come before TestHarness.h)
#include <CppUTest/TestHarness.h>
#include <CppUTestExt/MockSupport.h>

namespace amici {

class ReturnData;
class ExpData;

#if !defined(NEW_OPTION_FILE) || !defined(HDFFILE) || !defined(HDFFILEWRITE)
# error "Must define NEW_OPTION_FILE HDFFILE HDFFILEWRITE"
#endif

#define TEST_ATOL 1e-10
#define TEST_RTOL 1e-05

/**
 * @brief helper function to initialize default names/ids
 * @param name name of variables
 * @param length number of variables
 * @return default names/ids
 */
std::vector<std::string> getVariableNames(const char* name, int length);

/**
 * @brief The Model_Test class is a model-unspecific implementation
 of model designed for unit testing.
 */
class Model_Test : public Model {
public:
  /** constructor with model dimensions
   * @param nx number of state variables
   * @param nxtrue number of state variables of the non-augmented model
   * @param ny number of observables
   * @param nytrue number of observables of the non-augmented model
   * @param nz number of event observables
   * @param nztrue number of event observables of the non-augmented model
   * @param ne number of events
   * @param nJ number of objective functions
   * @param nw number of repeating elements
   * @param ndwdx number of nonzero elements in the x derivative of the
   * repeating elements
   * @param ndwdp number of nonzero elements in the p derivative of the
   * repeating elements
   * @param nnz number of nonzero elements in Jacobian
   * @param ubw upper matrix bandwidth in the Jacobian
   * @param lbw lower matrix bandwidth in the Jacobian
   * @param o2mode second order sensitivity mode
   * @param p parameters
   * @param k constants
   * @param plist indexes wrt to which sensitivities are to be computed
   * @param idlist indexes indicating algebraic components (DAE only)
   * @param z2event mapping of event outputs to events
   */
  Model_Test(const int nx_rdata, const int nxtrue_rdata, const int nx_solver,
             const int nxtrue_solver, const int nx_solver_reinit, const int ny, 
             const int nytrue, const int nz, const int nztrue, const int ne, 
             const int nJ, const int nw, const int ndwdx, const int ndwdp, 
             const int ndxdotdw, const int nnz, const int ubw, const int lbw,
             const SecondOrderMode o2mode, const std::vector<realtype> p,
             const std::vector<realtype> k, const std::vector<int> plist,
             const std::vector<realtype> idlist, const std::vector<int> z2event)
      : Model(nx_rdata, nxtrue_rdata, nx_solver, nxtrue_solver, nx_solver_reinit, ny, nytrue, nz,
              nztrue, ne, nJ, nw, ndwdx, ndwdp, ndxdotdw, {}, nnz, ubw, lbw, o2mode,
              p, k, plist, idlist, z2event) {}

  /** default constructor */
  Model_Test()
      : Model(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {}, 0, 0, 0,
              SecondOrderMode::none, std::vector<realtype>(),
              std::vector<realtype>(), std::vector<int>(),
              std::vector<realtype>(), std::vector<int>()) {}

  virtual Model *clone() const override { return new Model_Test(*this); }

  virtual std::unique_ptr<Solver> getSolver() override {
      throw AmiException("not implemented");
    }
    virtual void froot(const realtype t, const AmiVector &x,
                       const AmiVector &dx, gsl::span<realtype> root) override {
        throw AmiException("not implemented");
    }
    virtual void fxdot(const realtype t, const AmiVector &x,
                       const AmiVector &dx, AmiVector &xdot) override {
        throw AmiException("not implemented");
    }
    virtual void fsxdot(const realtype t, const AmiVector &x,
                        const AmiVector &dx, const int ip, const AmiVector &sx,
                        const AmiVector &sdx, AmiVector &sxdot) override {
        throw AmiException("not implemented");
    }
    virtual void fJ(const realtype t, const realtype cj, const AmiVector &x,
                    const AmiVector &dx, const AmiVector &xdot, SUNMatrix J)
    override {
        throw AmiException("not implemented");
    }
    virtual void fJSparse(const realtype t, const realtype cj,
                          const AmiVector &x, const AmiVector &dx,
                          const AmiVector &xdot, SUNMatrix J) override {
        throw AmiException("not implemented");
    }
    virtual void fJB(const realtype t, realtype cj, const AmiVector &x,
                     const AmiVector &dx, const AmiVector &xB, 
                     const AmiVector &dxB, const AmiVector &xBdot, 
                     SUNMatrix JB)
    override {
        throw AmiException("not implemented");
    }
    virtual void fJSparseB(const realtype t, realtype cj, const AmiVector &x,
                           const AmiVector &dx, const AmiVector &xB,
                           const AmiVector &dxB, const AmiVector &xBdot,
                           SUNMatrix JB) 
    override {
        throw AmiException("not implemented");
    }
    virtual void fJDiag(const realtype t, AmiVector &Jdiag,
                        const realtype cj, const AmiVector &x,
                        const AmiVector &dx) override {
        throw AmiException("not implemented");
    }
    virtual void fdxdotdp(const realtype t, const AmiVector &x,
                          const AmiVector &dx) override {
        throw AmiException("not implemented");
    }
    virtual void fJv(const realtype t, const AmiVector &x, const AmiVector &dx,
                     const AmiVector &xdot,const AmiVector &v, AmiVector &nJv,
                     const realtype cj) override {
        throw AmiException("not implemented");
    }

    virtual void fxBdot_ss(const realtype t, const AmiVector &xB,
                           const AmiVector &dxB, AmiVector &xBdot) override {
        throw AmiException("not implemented");
    }

    virtual void fJSparseB_ss(SUNMatrix JB) override {
        throw AmiException("not implemented");
    }

    virtual void writeSteadystateJB(const realtype t, realtype cj, 
                                    const AmiVector &x, const AmiVector &dx, 
                                    const AmiVector &xB, const AmiVector &dxB, 
                                    const AmiVector &xBdot) override {
        throw AmiException("not implemented");
    }

    virtual std::vector<std::string> getParameterNames() const override
    {
        return getVariableNames("p", np());
    }

    virtual std::vector<std::string> getStateNames() const override
    {
        return getVariableNames("x", nx_rdata);
    }

    virtual std::vector<std::string> getFixedParameterNames() const override
    {
        return getVariableNames("k", nk());
    }

    virtual std::vector<std::string> getObservableNames() const override
    {
        return getVariableNames("y", ny);
    }

    virtual std::vector<std::string> getParameterIds() const override
    {
        return getVariableNames("p", np());
    }

    virtual std::vector<std::string> getStateIds() const override
    {
        return getVariableNames("x", nx_rdata);
    }

    virtual std::vector<std::string> getFixedParameterIds() const override
    {
        return getVariableNames("k", nk());
    }

    virtual std::vector<std::string> getObservableIds() const override
    {
        return getVariableNames("y", ny);
    }


};

void simulateWithDefaultOptions();

void simulateVerifyWrite(const std::string& path);

void simulateVerifyWrite(std::string const& path, double atol, double rtol);

void simulateVerifyWrite(const std::string& hdffileOptions, const std::string& hdffileResults,
                         const std::string& hdffilewrite, const std::string& path,
                         double atol, double rtol);

std::unique_ptr<ExpData> getTestExpData(const Model &model);

bool withinTolerance(double expected, double actual, double atol, double rtol, int index, const char *name);

void checkEqualArray(const double *expected, const double *actual, int length, double atol, double rtol, const char *name);

void checkEqualArray(std::vector<double> const& expected, std::vector<double> const& actual,
                     double atol, double rtol, std::string const& name);

// TODO: delete after transitioning to C++-written test results
void verifyReturnDataMatlab(const std::string &hdffile, const std::string &resultPath, const ReturnData *rdata, const Model *model, double atol, double rtol);

// TODO: delete after transitioning to C++-written test results
void verifyReturnDataSensitivitiesMatlab(const H5::H5File &file_id, const std::string &resultPath, const ReturnData *rdata, const Model *model, double atol, double rtol);

void verifyReturnData(const std::string &hdffile, const std::string &resultPath, const ReturnData *rdata, const Model *model, double atol, double rtol);

void verifyReturnDataSensitivities(const H5::H5File &file_id, const std::string &resultPath, const ReturnData *rdata, const Model *model, double atol, double rtol);

void printBacktrace(int depth);

} // namespace amici

#endif
