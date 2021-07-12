#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

#include <gtest/gtest.h>

TEST(groupNeuron, testSimulation)
{
    amici::simulateVerifyWrite(
      "/model_neuron/nosensi/", 100 * TEST_ATOL, 100 * TEST_RTOL);
}

TEST(groupNeuron, testSensitivityForward)
{
    amici::simulateVerifyWrite(
      "/model_neuron/sensiforward/", 10 * TEST_ATOL, 10 * TEST_RTOL);
}
