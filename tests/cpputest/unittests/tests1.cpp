#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

#include <include/amici_solver_idas.h>
#include <include/symbolic_functions.h>
#include <include/amici_model.h>
#include <cstring>
#include <cmath>

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

    CHECK_EQUAL(AMICI_ERROR_UDATA, runAmiciSimulation(&udata, nullptr, nullptr, &model));
}

TEST(amici, testRunAmiciSimulationRdataMissing) {
    UserData udata(1, 2, 3);
    Model model(1, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, AMICI_O2MODE_NONE);

    CHECK_EQUAL(AMICI_ERROR_RDATA, runAmiciSimulation(&udata, nullptr, nullptr, &model));
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

TEST(userData, testCopy) {
    UserData udata1(1, 2, 3);
    udata1.k = new double[udata1.nk];
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
