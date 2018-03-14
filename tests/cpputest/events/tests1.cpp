#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupEvents)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(groupEvents, testSimulation) {
    amici::simulateVerifyWrite("/model_events/nosensi/");
}

TEST(groupEvents, testSensitivityForward) {
    amici::simulateVerifyWrite("/model_events/sensiforward/");
}


