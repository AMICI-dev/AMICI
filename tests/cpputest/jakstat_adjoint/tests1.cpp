#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupJakstatAdjoint)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(groupJakstatAdjoint, testSimulation) {
    amici::simulateAndVerifyFromFile("/model_jakstat_adjoint/nosensi/");
}

TEST(groupJakstatAdjoint, testSensitivityForward) {
    amici::simulateAndVerifyFromFile("/model_jakstat_adjoint/sensiforward/");
}

TEST(groupJakstatAdjoint, testSensitivityForwardLogParam) {
    amici::simulateAndVerifyFromFile("/model_jakstat_adjoint/sensiforwardlogparam/");
}

TEST(groupJakstatAdjoint, testSensitivityAdjoint) {
    amici::simulateAndVerifyFromFile("/model_jakstat_adjoint/sensiadjoint/");
    amici::simulateAndWriteToFile("/model_jakstat_adjoint/sensiadjoint/");
}


