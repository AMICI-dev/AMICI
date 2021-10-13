#include "wrapfunctions.h"
#include <cstring>
#include <testfunctions.h>

#include <gtest/gtest.h>

TEST(ExampleRobertson, Simulation)
{
    amici::simulateVerifyWrite("/model_robertson/nosensi/");
}

TEST(ExampleRobertson, SensitivityForward)
{
    amici::simulateVerifyWrite(
      "/model_robertson/sensiforward/", 1e6 * TEST_ATOL, 1e2 * TEST_RTOL);
}

TEST(ExampleRobertson, SensitivityForwardDense)
{
    amici::simulateVerifyWrite(
      "/model_robertson/sensiforwarddense/", 1e6 * TEST_ATOL, 1e2 * TEST_RTOL);
}

TEST(ExampleRobertson, SensitivityForwardSPBCG)
{
    amici::simulateVerifyWrite(
      "/model_robertson/sensiforwardSPBCG/", 1e7 * TEST_ATOL, 1e2 * TEST_RTOL);
}
