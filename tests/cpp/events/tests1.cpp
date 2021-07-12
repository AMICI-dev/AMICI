#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

#include <gtest/gtest.h>

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
