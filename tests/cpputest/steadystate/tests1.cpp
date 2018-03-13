#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici.h>

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupSteadystate)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(groupSteadystate, testModelFromHDF5) {
    // Test reading some python-written options
    std::vector<double> pExp {1, 0.5, 0.4, 2, 0.1};
    std::vector<double> kExp {0.1, 0.4, 0.7, 1.0};

    auto model = getModel();
    amici::hdf5::readModelDataFromHDF5(NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");

    amici::checkEqualArray(kExp, model->getFixedParameters(), TEST_ATOL, TEST_RTOL, "k");
    CHECK_EQUAL(51, model->nt());
    CHECK_EQUAL(0.0, model->t(0));
    CHECK_EQUAL(100.0, model->t(model->nt() - 2));
    CHECK_EQUAL(INFINITY, model->t(model->nt() - 1));

    for(int i = 0; i < model->np(); ++i) {
        CHECK_EQUAL(pExp[i], model->getUnscaledParameters()[i]);
        CHECK_EQUAL(log10(pExp[i]), model->getParameters()[i]);
    }
}

TEST(groupSteadystate, testInequality) {
    auto modelA = getModel();
    auto modelB = std::unique_ptr<amici::Model>(new amici::Model_Test());

    CHECK_FALSE(*modelA == *modelB);
}


TEST(groupSteadystate, testCopyModel) {
    auto modelA = getModel();
    auto modelB = std::unique_ptr<amici::Model>(modelA->clone());

    CHECK_TRUE(*modelA == *modelB);
}


TEST(groupSteadystate, testCloneModel) {
    auto modelA = getModel();
    auto modelB = std::unique_ptr<amici::Model>(new Model_model_steadystate());

    CHECK_TRUE(*modelA == *modelB);
}


TEST(groupSteadystate, testReuseSolver) {
    auto model = getModel();
    auto solver = model->getSolver();

    amici::hdf5::readModelDataFromHDF5(NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");
    amici::hdf5::readSolverSettingsFromHDF5(NEW_OPTION_FILE, *solver, "/model_steadystate/nosensi/options");

    runAmiciSimulation(*solver, nullptr, *model);
    runAmiciSimulation(*solver, nullptr, *model);
}



TEST(groupSteadystate, testSimulation) {
    amici::simulateVerifyWrite("/model_steadystate/nosensi/");
}

TEST(groupSteadystate, testSensitivityForward) {
    amici::simulateVerifyWrite("/model_steadystate/sensiforward/");
}

TEST(groupSteadystate, testSensitivityForwardPlist) {
    amici::simulateVerifyWrite("/model_steadystate/sensiforwardplist/");
}


TEST(groupSteadystate, testSensitivityForwardErrorInt) {
    amici::simulateVerifyWrite("/model_steadystate/sensiforwarderrorint/");
}

TEST(groupSteadystate, testSensitivityForwardErrorNewt) {
    amici::simulateVerifyWrite("/model_steadystate/sensiforwarderrornewt/");
}


TEST(groupSteadystate, testSensitivityForwardDense) {
    amici::simulateVerifyWrite("/model_steadystate/sensiforwarddense/");
}

TEST(groupSteadystate, testSensitivityForwardSPBCG) {
    amici::simulateVerifyWrite("/model_steadystate/nosensiSPBCG/", 10*TEST_ATOL, 10*TEST_RTOL);
}


