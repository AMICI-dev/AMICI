#include "testfunctions.h"

#include <amici/amici.h>
#include <amici/amici_solver_idas.h>
#include <amici/amici_solver_cvodes.h>
#include <amici/symbolic_functions.h>
#include <amici/amici_model_ode.h>

#include <cstring>
#include <cmath>
#include <vector>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"


std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(new amici::Model_Test());
}

using namespace amici;

TEST_GROUP(amici)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST_GROUP(model)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(model, testScalingLin) {
    Model_Test model(3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, AMICI_O2MODE_NONE,
                     std::vector<realtype>(1,0.0),std::vector<realtype>(3,0),std::vector<int>(2,1),
                     std::vector<realtype>(0,0.0),std::vector<int>(0,1));

    std::vector<double> p {1};
    model.setParameters(p);

    model.setParameterScale(AMICI_SCALING_NONE);
    double unscaled[1];
    model.unscaleParameters(unscaled);

    CHECK_EQUAL(p[0], unscaled[0]);
}

TEST(model, testScalingLog) {
    Model_Test model(3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, AMICI_O2MODE_NONE,
                     std::vector<realtype>(1,0.0),std::vector<realtype>(3,0),std::vector<int>(2,1),
                     std::vector<realtype>(0,0.0),std::vector<int>(0,1));

    std::vector<double> p {1};
    model.setParameters(p);

    model.setParameterScale(AMICI_SCALING_LN);
    double unscaled[1];
    model.unscaleParameters(unscaled);

    DOUBLES_EQUAL(exp(p[0]), unscaled[0], 1e-16);
}

TEST(model, testScalingLog10) {
    Model_Test model(3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, AMICI_O2MODE_NONE,
                     std::vector<realtype>(1,0.0),std::vector<realtype>(3,0),std::vector<int>(2,1),
                     std::vector<realtype>(0,0.0),std::vector<int>(0,1));

    std::vector<double> p {1};
    model.setParameters(p);

    model.setParameterScale(AMICI_SCALING_LOG10);
    double unscaled[1];
    model.unscaleParameters(unscaled);

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
    CHECK_EQUAL(-1, amici::min(-1, 2, 0));
    CHECK_EQUAL(-2, amici::min(1, -2, 0));
    CHECK_TRUE(amici::isNaN(amici::min(amici::getNaN(), amici::getNaN(), 0)));
    CHECK_EQUAL(-1, amici::min(-1, amici::getNaN(), 0));
    CHECK_EQUAL(-1, amici::min(amici::getNaN(), -1, 0));
}

TEST(symbolicFunctions, testMax) {
    CHECK_EQUAL(2, amici::max(-1, 2, 0));
    CHECK_EQUAL(1, amici::max(1, -2, 0));
    CHECK_TRUE(amici::isNaN(amici::max(amici::getNaN(), amici::getNaN(), 0)));
    CHECK_EQUAL(-1, amici::max(-1, amici::getNaN(), 0));
    CHECK_EQUAL(-1, amici::max(amici::getNaN(), -1, 0));
}


TEST(symbolicFunctions, testDMin) {
    CHECK_EQUAL(0, amici::Dmin(1, -1, -2, 0));
    CHECK_EQUAL(1, amici::Dmin(1, -1, 2, 0));
    CHECK_EQUAL(1, amici::Dmin(2, -1, -2, 0));
    CHECK_EQUAL(0, amici::Dmin(2, -1, 2, 0));
}

TEST(symbolicFunctions, testDMax) {
    CHECK_EQUAL(1, amici::Dmax(1, -1, -2, 0));
    CHECK_EQUAL(0, amici::Dmax(1, -1, 2, 0));
    CHECK_EQUAL(0, amici::Dmax(2, -1, -2, 0));
    CHECK_EQUAL(1, amici::Dmax(2, -1, 2, 0));
}

TEST_GROUP(amiciSolver)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(amiciSolver, testEquality) {
    IDASolver i1, i2;
    CVodeSolver c1, c2;

    CHECK_TRUE(i1 == i2);
    CHECK_TRUE(c1 == c2);
    CHECK_FALSE(i1 == c1);
}

TEST(amiciSolver, testClone) {
    IDASolver i1;
    auto i2 = std::unique_ptr<Solver>(i1.clone());
    CHECK_TRUE(i1 == *i2);

    CVodeSolver c1;
    auto c2 = std::unique_ptr<Solver>(c1.clone());
    CHECK_TRUE(c1 == *c2);
    CHECK_FALSE(*i2 == *c2);
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

