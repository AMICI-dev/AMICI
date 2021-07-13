#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

#include <gtest/gtest.h>

TEST(ExampleEvents, Default)
{
    amici::simulateWithDefaultOptions();
}

TEST(ExampleEvents, Simulation)
{
    amici::simulateVerifyWrite("/model_events/nosensi/");
}

TEST(ExampleEvents, SensitivityForward)
{
    amici::simulateVerifyWrite("/model_events/sensiforward/");
}
