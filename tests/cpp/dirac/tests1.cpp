#include <testfunctions.h>

#include "wrapfunctions.h"
#include <cstring>

#include <gtest/gtest.h>


TEST(ExampleDirac, Simulation)
{
    amici::simulateVerifyWrite("/model_dirac/nosensi/");
}

TEST(ExampleDirac, SensitivityForward)
{
    amici::simulateVerifyWrite("/model_dirac/sensiforward/");
}
