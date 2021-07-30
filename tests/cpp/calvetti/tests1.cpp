#include "wrapfunctions.h"
#include <cstring>
#include <testfunctions.h>

#include <gtest/gtest.h>

TEST(ExampleCalvetti, Simulation)
{
    amici::simulateVerifyWrite("/model_calvetti/nosensi/");
}
