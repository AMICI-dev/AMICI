#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "wrapfunctions.h"
#include <cstring>
#include <testfunctions.h>

TEST_GROUP(groupCalvetti){};

TEST(groupCalvetti, testSimulation)
{
    amici::simulateVerifyWrite("/model_calvetti/nosensi/");
}
