#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

#include <gtest/gtest.h>

TEST(ExampleNeuron, Simulation)
{
    amici::simulateVerifyWrite(
      "/model_neuron/nosensi/", 100 * TEST_ATOL, 100 * TEST_RTOL);
}

TEST(ExampleNeuron, SensitivityForward)
{
    amici::simulateVerifyWrite(
      "/model_neuron/sensiforward/", 10 * TEST_ATOL, 10 * TEST_RTOL);
}
