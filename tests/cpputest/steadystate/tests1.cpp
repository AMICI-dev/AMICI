#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

#include <cstring>
#include "wrapfunctions.h"

#define NEW_OPTION_FILE "../../testOptions.h5"
TEST_GROUP(groupSteadystate)
{
    void setup() {

    }

    void teardown() {

    }
};

void simulateAndVerifyFromFile(const std::string hdffileOptions, const std::string hdffileResults, std::string path, double atol, double rtol)
{
    using namespace amici;
    // read options from file
    std::string optionsPath = path + "/options";
    auto model = getModel();
    auto solver = model->getSolver();
    hdf5::readModelDataFromHDF5(hdffileOptions, *model, optionsPath);
    hdf5::readSolverSettingsFromHDF5(hdffileOptions, *solver, optionsPath);

    // read measurements from file
    std::string measurementPath = path + "/data";

    std::unique_ptr<const ExpData> edata;
    if(hdf5::locationExists(hdffileOptions, measurementPath))
        edata = hdf5::readSimulationExpData(hdffileResults, measurementPath, *model);

    // simulate & verify
    auto rdata = std::unique_ptr<ReturnData>(getSimulationResults(*model, edata.get(), *solver));
    std::string resultPath = path + "/results";
    verifyReturnData(hdffileResults.c_str(), resultPath.c_str(), rdata.get(), model.get(), atol, rtol);
}
void simulateAndVerifyFromFile(std::string path) {
    simulateAndVerifyFromFile(NEW_OPTION_FILE, HDFFILE, path, TEST_ATOL, TEST_RTOL);

}


TEST(groupSteadystate, testModelFromHDF5) {
    // Test reading some python-written options
    std::vector<double> pExp {1, 0.5, 0.4, 2, 0.1};
    std::vector<double> kExp {0.1, 0.4, 0.7, 1.0};

    auto model = getModel();
    amici::hdf5::readModelDataFromHDF5(NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");

    amici::checkEqualArray(kExp.data(), model->k(), kExp.size(), TEST_ATOL, TEST_RTOL, "k");
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

    amici::hdf5::readModelDataFromHDF5(HDFFILE, *model, "/model_steadystate/nosensi/options");
    amici::hdf5::readSolverSettingsFromHDF5(HDFFILE, *solver, "/model_steadystate/nosensi/options");

    std::unique_ptr<amici::ReturnData>(amici::getSimulationResults(*model, nullptr, *solver));
    std::unique_ptr<amici::ReturnData>(amici::getSimulationResults(*model, nullptr, *solver));
}



TEST(groupSteadystate, testSimulation) {
    simulateAndVerifyFromFile("/model_steadystate/nosensi/");
    amici::simulateAndWriteToFile("/model_steadystate/nosensi/");
}

TEST(groupSteadystate, testSensitivityForward) {
    simulateAndVerifyFromFile("/model_steadystate/sensiforward/");
    amici::simulateAndWriteToFile("/model_steadystate/sensiforward/");
}

TEST(groupSteadystate, testSensitivityForwardPlist) {
    simulateAndVerifyFromFile("/model_steadystate/sensiforwardplist/");
    amici::simulateAndWriteToFile("/model_steadystate/sensiforwardplist/");
}


TEST(groupSteadystate, testSensitivityForwardErrorInt) {
    simulateAndVerifyFromFile("/model_steadystate/sensiforwarderrorint/");
    amici::simulateAndWriteToFile("/model_steadystate/sensiforwarderrorint/");
}

TEST(groupSteadystate, testSensitivityForwardErrorNewt) {
    simulateAndVerifyFromFile("/model_steadystate/sensiforwarderrornewt/");
    amici::simulateAndWriteToFile("/model_steadystate/sensiforwarderrornewt/");
}


TEST(groupSteadystate, testSensitivityForwardDense) {
    simulateAndVerifyFromFile("/model_steadystate/sensiforwarddense/");
}

TEST(groupSteadystate, testSensitivityForwardSPBCG) {
    simulateAndVerifyFromFile(NEW_OPTION_FILE, HDFFILE, "/model_steadystate/nosensiSPBCG/",10*TEST_ATOL, 10*TEST_RTOL);
}


