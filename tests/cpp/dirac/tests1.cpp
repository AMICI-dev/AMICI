#include <testfunctions.h>

#include "wrapfunctions.h"
#include <cstring>

#include <gtest/gtest.h>


TEST(groupDirac, testSimulation)
{
    amici::simulateVerifyWrite("/model_dirac/nosensi/");
}

TEST(groupDirac, testSensitivityForward)
{
    amici::simulateVerifyWrite("/model_dirac/sensiforward/");
}
