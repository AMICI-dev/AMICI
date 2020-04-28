#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupNeuron)
{
    void setup() {

    }

    void teardown() {

    }
};


TEST(groupNeuron, testSimulation) {
    amici::simulateVerifyWrite("/model_neuron/nosensi/",
                               100*TEST_ATOL, 100*TEST_RTOL);
}

TEST(groupNeuron, testSensitivityForward) {
    amici::simulateVerifyWrite("/model_neuron/sensiforward/",
                               10*TEST_ATOL, 10*TEST_RTOL);
}
