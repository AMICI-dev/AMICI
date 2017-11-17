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
    amici::simulateAndVerifyFromFile("/model_steadystate/nosensiSPBCG/");
}


