#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

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
    amici::simulateVerifyWrite("/model_jakstat_adjoint/nosensi/");
}

TEST(groupJakstatAdjoint, testSensitivityForward) {
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensiforward/");
}

TEST(groupJakstatAdjoint, testSensitivityForwardLogParam) {
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensiforwardlogparam/");
}

TEST(groupJakstatAdjoint, testSensitivityAdjoint) {
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensiadjoint/");
}

TEST(groupJakstatAdjoint, testSensitivityForwardEmptySensInd) {
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensiforwardemptysensind/");
}

TEST(groupJakstatAdjoint, testSensitivityAdjointEmptySensInd) {
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensiadjointemptysensind/");
}

