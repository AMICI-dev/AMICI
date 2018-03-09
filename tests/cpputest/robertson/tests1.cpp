#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <testfunctions.h>
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>
#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupRobertson)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(groupRobertson, testSimulation) {
    amici::simulateVerifyWrite("/model_robertson/nosensi/");
}

TEST(groupRobertson, testSensitivityForward) {
    amici::simulateVerifyWrite("/model_robertson/sensiforward/");
}




