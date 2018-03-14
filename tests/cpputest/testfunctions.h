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
#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

namespace amici {

class ReturnData;
class ExpData;

#define NEW_OPTION_FILE "../../testOptions.h5"
#define HDFFILE "../../expectedResults.h5"
#define HDFFILEWRITE "../../writeResults.h5"
#define TEST_ATOL 1e-10
#define TEST_RTOL 1e-05


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
    Model_Test(const int nx, const int nxtrue,
               const int ny, const int nytrue, const int nz, const int nztrue,
               const int ne, const int nJ, const int nw, const int ndwdx,
               const int ndwdp, const int nnz, const int ubw, const int lbw,
               const AMICI_o2mode o2mode, const std::vector<realtype> p,
               const std::vector<realtype> k, const std::vector<int> plist,
               const std::vector<realtype> idlist, const std::vector<int> z2event)
        : Model(nx,nxtrue,ny,nytrue,nz,nztrue,ne,nJ,nw,ndwdx,ndwdp,nnz,ubw,lbw,o2mode,p,k,plist,idlist,z2event) {}

    /** default constructor */
    Model_Test()
        : Model(0,0,0,0,0,0,0,0,0,0,0,0,0,0,AMICI_O2MODE_NONE,std::vector<realtype>(),std::vector<realtype>(),std::vector<int>(),std::vector<realtype>(),std::vector<int>()) {}

    virtual Model* clone() const override { return new Model_Test(*this); }

    virtual std::unique_ptr<Solver> getSolver() override {
        throw AmiException("not implemented");
    }
    virtual void froot(realtype t, AmiVector *x, AmiVector *dx, realtype *root) override {
        throw AmiException("not implemented");
    }
    virtual void fxdot(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot) override {
        throw AmiException("not implemented");
    }
    virtual void fJ(realtype t, realtype cj, AmiVector *x, AmiVector *dx, AmiVector *xdot, DlsMat J) override {
        throw AmiException("not implemented");
    }
    virtual void fJSparse(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                          AmiVector *xdot, SlsMat J) override {
        throw AmiException("not implemented");
    }
    virtual void fJDiag(realtype t, AmiVector *Jdiag, realtype cj, AmiVector *x,
                        AmiVector *dx) override {
        throw AmiException("not implemented");
    }
    virtual void fdxdotdp(realtype t, AmiVector *x, AmiVector *dx) override {
        throw AmiException("not implemented");
    }
    virtual void fJv(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot,
                     AmiVector *v, AmiVector *nJv, realtype cj) override {
        throw AmiException("not implemented");
    }
};


void simulateVerifyWrite(const std::string path);

void simulateVerifyWrite(std::string path, double atol, double rtol);

void simulateVerifyWrite(const std::string hdffileOptions, const std::string hdffileResults,
                         const std::string hdffilewrite, std::string path,
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
