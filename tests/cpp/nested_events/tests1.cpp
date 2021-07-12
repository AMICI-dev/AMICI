#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

#include <gtest/gtest.h>

TEST(groupEvents, testSimulation)
{
    amici::simulateVerifyWrite("/model_nested_events/nosensi/");
}

TEST(groupEvents, testSensitivityForward)
{
    amici::simulateVerifyWrite("/model_nested_events/sensiforward/");
}
