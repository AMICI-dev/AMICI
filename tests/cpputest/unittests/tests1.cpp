#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

#include <include/amici_solver_idas.h>
#include <include/symbolic_functions.h>
#include <cstring>
#include <cmath>

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

    DOUBLES_EQUAL(pow10(p[0]), unscaled[0], 1e-16);
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
