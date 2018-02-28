#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupSteadystate)
{
    void setup() {

    }

    void teardown() {

    }
};

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

    amici::readModelDataFromHDF5(HDFFILE, *model, "/model_steadystate/nosensi/options");
    amici::readSolverSettingsFromHDF5(HDFFILE, *solver, "/model_steadystate/nosensi/options");

    std::unique_ptr<amici::ReturnData>(amici::getSimulationResults(*model, nullptr, *solver));
    std::unique_ptr<amici::ReturnData>(amici::getSimulationResults(*model, nullptr, *solver));
}



TEST(groupSteadystate, testSimulation) {
    amici::simulateAndVerifyFromFile("/model_steadystate/nosensi/");
    amici::simulateAndWriteToFile("/model_steadystate/nosensi/");
}

TEST(groupSteadystate, testSensitivityForward) {
    amici::simulateAndVerifyFromFile("/model_steadystate/sensiforward/");
    amici::simulateAndWriteToFile("/model_steadystate/sensiforward/");
}

TEST(groupSteadystate, testSensitivityForwardPlist) {
    amici::simulateAndVerifyFromFile("/model_steadystate/sensiforwardplist/");
    amici::simulateAndWriteToFile("/model_steadystate/sensiforwardplist/");
}


TEST(groupSteadystate, testSensitivityForwardErrorInt) {
    amici::simulateAndVerifyFromFile("/model_steadystate/sensiforwarderrorint/");
    amici::simulateAndWriteToFile("/model_steadystate/sensiforwarderrorint/");
}

TEST(groupSteadystate, testSensitivityForwardErrorNewt) {
    amici::simulateAndVerifyFromFile("/model_steadystate/sensiforwarderrornewt/");
    amici::simulateAndWriteToFile("/model_steadystate/sensiforwarderrornewt/");
}


TEST(groupSteadystate, testSensitivityForwardDense) {
    amici::simulateAndVerifyFromFile("/model_steadystate/sensiforwarddense/");
}

TEST(groupSteadystate, testSensitivityForwardSPBCG) {
    amici::simulateAndVerifyFromFile("/model_steadystate/nosensiSPBCG/",10*TEST_ATOL, 10*TEST_RTOL);
}


