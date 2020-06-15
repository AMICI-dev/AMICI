#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <testfunctions.h>

#include "wrapfunctions.h"
#include <cstring>

TEST_GROUP(groupDirac){};

TEST(groupDirac, testSimulation)
{
    amici::simulateVerifyWrite("/model_dirac/nosensi/");
}

TEST(groupDirac, testSensitivityForward)
{
    amici::simulateVerifyWrite("/model_dirac/sensiforward/");
}
