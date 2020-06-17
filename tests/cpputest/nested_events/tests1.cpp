#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

TEST_GROUP(groupEvents){};

TEST(groupEvents, testSimulation)
{
    amici::simulateVerifyWrite("/model_nested_events/nosensi/");
}

TEST(groupEvents, testSensitivityForward)
{
    amici::simulateVerifyWrite("/model_nested_events/sensiforward/");
}
