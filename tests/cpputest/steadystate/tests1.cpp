#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupSteadystate)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(groupSteadystate, testDefault) {
    amici::simulateWithDefaultOptions();
}


TEST(groupSteadystate, testModelFromHDF5) {
    // Test reading some python-written options
    std::vector<double> pExp {1, 0.5, 0.4, 2, 0.1};
    std::vector<double> kExp {0.1, 0.4, 0.7, 1.0};

    auto model = amici::generic_model::getModel();
    amici::hdf5::readModelDataFromHDF5(NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");

    amici::checkEqualArray(kExp, model->getFixedParameters(), TEST_ATOL, TEST_RTOL, "k");
    CHECK_EQUAL(51, model->nt());
    CHECK_EQUAL(0.0, model->getTimepoint(0));
    CHECK_EQUAL(100.0, model->getTimepoint(model->nt() - 2));
    CHECK_EQUAL(INFINITY, model->getTimepoint(model->nt() - 1));

    for(int i = 0; i < model->np(); ++i) {
        CHECK_EQUAL(pExp[i], model->getUnscaledParameters()[i]);
        CHECK_EQUAL(log10(pExp[i]), model->getParameters()[i]);
    }
}

TEST(groupSteadystate, testInequality) {
    auto modelA = amici::generic_model::getModel();
    auto modelB = std::unique_ptr<amici::Model>(new amici::Model_Test());

    CHECK_FALSE(*modelA == *modelB);
}


TEST(groupSteadystate, testCopyModel) {
    auto modelA = amici::generic_model::getModel();
    auto modelB = std::unique_ptr<amici::Model>(modelA->clone());

    CHECK_TRUE(*modelA == *modelB);
}


TEST(groupSteadystate, testCloneModel) {
    auto modelA = amici::generic_model::getModel();
    auto modelB = std::unique_ptr<amici::Model>(
        new amici::model_model_steadystate::Model_model_steadystate());

    CHECK_TRUE(*modelA == *modelB);
}

TEST(groupSteadystate, testExpDataFromReturnData) {
    auto model = amici::generic_model::getModel();
    auto solver = model->getSolver();

    amici::hdf5::readModelDataFromHDF5(NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");
    amici::hdf5::readSolverSettingsFromHDF5(NEW_OPTION_FILE, *solver, "/model_steadystate/nosensi/options");

    auto rdata = runAmiciSimulation(*solver, nullptr, *model);
    auto edata = amici::ExpData(*rdata, 0.1, 0.1);
    runAmiciSimulation(*solver, &edata, *model);
}

TEST(groupSteadystate, testReuseSolver) {
    auto model = amici::generic_model::getModel();
    auto solver = model->getSolver();

    amici::hdf5::readModelDataFromHDF5(NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");
    amici::hdf5::readSolverSettingsFromHDF5(NEW_OPTION_FILE, *solver, "/model_steadystate/nosensi/options");

    runAmiciSimulation(*solver, nullptr, *model);
    runAmiciSimulation(*solver, nullptr, *model);
}

TEST(groupSteadystate, testRethrow) {
    auto model = amici::generic_model::getModel();
    auto solver = model->getSolver();

    amici::hdf5::readModelDataFromHDF5(NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");
    amici::hdf5::readSolverSettingsFromHDF5(NEW_OPTION_FILE, *solver, "/model_steadystate/nosensi/options");

    // p = NaN will raise amici::IntegrationFailure
    auto p = model->getParameters();
    std::fill(p.begin(), p.end(), std::nan(""));
    model->setParameters(p);

    // must not throw
    runAmiciSimulation(*solver, nullptr, *model);

    // must throw
    CHECK_THROWS(amici::IntegrationFailure, runAmiciSimulation(*solver, nullptr, *model, true));
}

TEST(groupSteadystate, testInitialStatesNonEmpty) {
    auto model = amici::generic_model::getModel();
    CHECK_FALSE(model->getInitialStates().empty());
}

TEST(groupSteadystate, testInitialStateSensitivitiesNonEmpty) {
    auto model = amici::generic_model::getModel();
    CHECK_FALSE(model->getInitialStateSensitivities().empty());
}

TEST(groupSteadystate, testSimulation) {
    amici::simulateVerifyWrite("/model_steadystate/nosensi/",
                               100*TEST_ATOL, 100*TEST_RTOL);
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

TEST(groupSteadystate, testSensiFwdNewtonPreeq) {
    amici::simulateVerifyWrite("/model_steadystate/sensifwdnewtonpreeq/");
}

TEST(groupSteadystate, testSensiAdjNewtonPreeq) {
    amici::simulateVerifyWrite("/model_steadystate/sensiadjnewtonpreeq/");
}

TEST(groupSteadystate, testSensiFwdSimPreeq) {
    amici::simulateVerifyWrite("/model_steadystate/sensifwdsimpreeq/");
}

TEST(groupSteadystate, testSensiFwdSimPreeqFSA) {
    amici::simulateVerifyWrite("/model_steadystate/sensifwdsimpreeqFSA/");
}

TEST(groupSteadystate, testSensiAdjSimPreeq) {
    amici::simulateVerifyWrite("/model_steadystate/sensiadjsimpreeq/");
}

TEST(groupSteadystate, testSensiAdjSimPreeqFSA) {
    amici::simulateVerifyWrite("/model_steadystate/sensiadjsimpreeqFSA/");
}

TEST(groupSteadystate, testSensiFwdByhandPreeq) {
    amici::simulateVerifyWrite("/model_steadystate/sensifwdbyhandpreeq/");
}

TEST(groupSteadystate, testSensiAdjByhandPreeq) {
    amici::simulateVerifyWrite("/model_steadystate/sensiadjbyhandpreeq/");
}
