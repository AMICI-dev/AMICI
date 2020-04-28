#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.hpp"

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
    amici::simulateVerifyWrite("/model_nested_events/nosensi/");
}

TEST(groupEvents, testSensitivityForward) {
    amici::simulateVerifyWrite("/model_nested_events/sensiforward/");
}


