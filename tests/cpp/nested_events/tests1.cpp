#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

#include <gtest/gtest.h>

TEST(ExampleNestedEvents, Simulation)
{
    amici::simulateVerifyWrite("/model_nested_events/nosensi/");
}

TEST(ExampleNestedEvents, SensitivityForward)
{
    amici::simulateVerifyWrite("/model_nested_events/sensiforward/");
}
