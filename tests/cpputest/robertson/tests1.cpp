#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "wrapfunctions.h"
#include <cstring>
#include <testfunctions.h>

TEST_GROUP(groupRobertson){};

TEST(groupRobertson, testSimulation)
{
    amici::simulateVerifyWrite("/model_robertson/nosensi/");
}

TEST(groupRobertson, testSensitivityForward)
{
    amici::simulateVerifyWrite(
      "/model_robertson/sensiforward/", 1e6 * TEST_ATOL, 1e2 * TEST_RTOL);
}

TEST(groupRobertson, testSensitivityForwardDense)
{
    amici::simulateVerifyWrite(
      "/model_robertson/sensiforwarddense/", 1e6 * TEST_ATOL, 1e2 * TEST_RTOL);
}

TEST(groupRobertson, testSensitivityForwardSPBCG)
{
    amici::simulateVerifyWrite(
      "/model_robertson/sensiforwardSPBCG/", 1e7 * TEST_ATOL, 1e2 * TEST_RTOL);
}
