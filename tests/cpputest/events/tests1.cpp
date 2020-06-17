#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

TEST_GROUP(groupEvents){};

TEST(groupEvents, testDefault)
{
    amici::simulateWithDefaultOptions();
}

TEST(groupEvents, testSimulation)
{
    amici::simulateVerifyWrite("/model_events/nosensi/");
}

TEST(groupEvents, testSensitivityForward)
{
    amici::simulateVerifyWrite("/model_events/sensiforward/");
}
