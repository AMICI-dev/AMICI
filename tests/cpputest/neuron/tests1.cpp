#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

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
    amici::simulateVerifyWrite("/model_neuron/nosensi/", 10*TEST_ATOL, 10*TEST_RTOL);
}

TEST(groupNeuron, testSensitivityForward) {
    amici::simulateVerifyWrite("/model_neuron/sensiforward/", 10*TEST_ATOL, 10*TEST_RTOL);
}
