#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

#include <gtest/gtest.h>

TEST(ExampleNeuronO2, Sensitivity2)
{
    amici::simulateVerifyWrite(
      "/model_neuron/sensi2forward/", 10 * TEST_ATOL, 10 * TEST_RTOL);
}
