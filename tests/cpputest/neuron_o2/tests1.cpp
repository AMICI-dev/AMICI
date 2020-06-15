#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

TEST_GROUP(groupNeuronO2){};

TEST(groupNeuronO2, testSensitivity2)
{
    amici::simulateVerifyWrite(
      "/model_neuron/sensi2forward/", 10 * TEST_ATOL, 10 * TEST_RTOL);
}
