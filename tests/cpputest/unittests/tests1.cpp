#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

#include <include/amici_solver_idas.h>
#include <include/symbolic_functions.h>
#include <include/amici_model.h>
#include <cstring>
#include <cmath>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

using namespace amici;

TEST_GROUP(amici)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(amici, testRunAmiciSimulationStateMismatch) {
    UserData udata(1, 2, 3);
    Model model(4, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, AMICI_O2MODE_NONE);
    bool exceptionThrown = false;
    try {
        runAmiciSimulation(&udata, nullptr, nullptr, &model);
    } catch (amici::AmiException& ex) {
        exceptionThrown = true;
    }
    CHECK_TRUE(exceptionThrown)
}

TEST(amici, testRunAmiciSimulationRdataMissing) {
    UserData udata(1, 2, 3);
    Model model(1, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, AMICI_O2MODE_NONE);
    bool exceptionThrown = false;
    try {
        runAmiciSimulation(&udata, nullptr, nullptr, &model);
    } catch (amici::AmiException& ex) {
        exceptionThrown = true;
    }
    CHECK_TRUE(exceptionThrown)
}

TEST_GROUP(userData)
{
    void setup() {

    }

    void teardown() {

    }
};


TEST(userData, testConstructionDestruction) {
    UserData udata;
}

TEST(userData, testPrint) {
    UserData udata(1, 2, 3);
    udata.print();
}

TEST(userData, setPbar) {
    UserData udata(1, 2, 3);
    std::vector<double> pbar = {1.0};
    udata.setPbar(pbar.data());
    udata.setPbar(nullptr);
    udata.setPbar(pbar.data());
}

TEST(userData, setXbar) {
    UserData udata(1, 2, 3);
    std::vector<double> xbar = {1.0, 2.0, 3.0};
    udata.setXbar(xbar.data());
    udata.setXbar(nullptr);
    udata.setXbar(xbar.data());
}

TEST(userData, setPlist) {
    UserData udata(1, 2, 3);
    std::vector<double> plistdouble = {1.0};
    std::vector<double> plistint = {1};
    udata.setPlist(plistdouble.data(),plistdouble.size());
    double* doubleptr = nullptr;
    udata.setPlist(doubleptr,plistdouble.size());
    int* intptr = nullptr;
    udata.setPlist(intptr,plistdouble.size());
    udata.setPlist(plistint.data(),plistint.size());
    udata.setPlist(plistdouble.data(),plistdouble.size());
}

TEST(userData, testCopy) {
    int np = 1;
    int nk = 2;
    int nx = 3;
    int nt = 4;
    UserData udata1(np, nk, nx);
    std::vector<double> ts = {1.0, 2.0, 3.0, 4.0};
    udata1.setTimepoints(ts.data(),ts.size());
    
    std::vector<double> plist = {1.0};
    udata1.setPlist(plist.data(),plist.size());
    
    std::vector<double> pbar = {1.0};
    udata1.setPbar(pbar.data());
    
    std::vector<double> xbar = {1.0, 2.0, 3.0};
    udata1.setXbar(xbar.data());
    
    std::vector<double> x0 = {1.0, 2.0, 3.0};
    udata1.setStateInitialization(x0.data());
    
    std::vector<double> sx0 = {1.0, 2.0, 3.0};
    udata1.setSensitivityInitialization(sx0.data());
    
    UserData udata2(udata1);
}


TEST(userData, testScalingLin) {
    UserData udata(1, 0, 0);
    const double p[1] = {1};
    udata.setParameters(p);

    udata.pscale = AMICI_SCALING_NONE;
    double unscaled[1];
    udata.unscaleParameters(unscaled);

    CHECK_EQUAL(p[0], unscaled[0]);
}

TEST(userData, testScalingLog) {
    UserData udata(1, 0, 0);
    const double p[1] = {1};
    udata.setParameters(p);

    udata.pscale = AMICI_SCALING_LN;
    double unscaled[1];
    udata.unscaleParameters(unscaled);

    DOUBLES_EQUAL(exp(p[0]), unscaled[0], 1e-16);
}

TEST(userData, testScalingLog10) {
    UserData udata(1, 0, 0);
    const double p[1] = {1};
    udata.setParameters(p);

    udata.pscale = AMICI_SCALING_LOG10;
    double unscaled[1];
    udata.unscaleParameters(unscaled);

    DOUBLES_EQUAL(pow(10, p[0]), unscaled[0], 1e-16);
}


TEST_GROUP(symbolicFunctions)
{
    void setup() {

    }

    void teardown() {

    }
};


TEST(symbolicFunctions, testSign) {
    CHECK_EQUAL(-1, sign(-2));
    CHECK_EQUAL(0, sign(0));
    CHECK_EQUAL(1, sign(2));
}

TEST(symbolicFunctions, testHeaviside) {
    CHECK_EQUAL(0, heaviside(-1));
    CHECK_EQUAL(0, heaviside(0));
    CHECK_EQUAL(1, heaviside(1));
}


TEST(symbolicFunctions, testMin) {
    CHECK_EQUAL(-1, am_min(-1, 2, 0));
    CHECK_EQUAL(-2, am_min(1, -2, 0));
    CHECK_TRUE(amiIsNaN(am_min(amiGetNaN(), amiGetNaN(), 0)));
    CHECK_EQUAL(-1, am_min(-1, amiGetNaN(), 0));
    CHECK_EQUAL(-1, am_min(amiGetNaN(), -1, 0));
}

TEST(symbolicFunctions, testMax) {
    CHECK_EQUAL(2, am_max(-1, 2, 0));
    CHECK_EQUAL(1, am_max(1, -2, 0));
    CHECK_TRUE(amiIsNaN(am_max(amiGetNaN(), amiGetNaN(), 0)));
    CHECK_EQUAL(-1, am_max(-1, amiGetNaN(), 0));
    CHECK_EQUAL(-1, am_max(amiGetNaN(), -1, 0));
}


TEST(symbolicFunctions, testDMin) {
    CHECK_EQUAL(0, Dam_min(1, -1, -2, 0));
    CHECK_EQUAL(1, Dam_min(1, -1, 2, 0));
    CHECK_EQUAL(1, Dam_min(2, -1, -2, 0));
    CHECK_EQUAL(0, Dam_min(2, -1, 2, 0));
}

TEST(symbolicFunctions, testDMax) {
    CHECK_EQUAL(1, Dam_max(1, -1, -2, 0));
    CHECK_EQUAL(0, Dam_max(1, -1, 2, 0));
    CHECK_EQUAL(0, Dam_max(2, -1, -2, 0));
    CHECK_EQUAL(1, Dam_max(2, -1, 2, 0));
}

TEST_GROUP(amiciSolverIdas)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(amiciSolverIdas, testConstructionDestruction) {
    IDASolver solver;
}
