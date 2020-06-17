#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

TEST_GROUP(groupJakstatAdjointO2){};

TEST(groupJakstatAdjointO2, testSensitivityForward2)
{
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensi2forward/");
}

TEST(groupJakstatAdjointO2, testSensitivityForward2LogParam)
{
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensi2forwardlogparam/");
}

TEST(groupJakstatAdjointO2, testSensitivityAdjoint2)
{
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensi2adjoint/");
}
