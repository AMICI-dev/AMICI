#include "wrapfunctions.h"
#include <cstring>
#include <testfunctions.h>

#include <gtest/gtest.h>

TEST(groupCalvetti, testSimulation)
{
    amici::simulateVerifyWrite("/model_calvetti/nosensi/");
}
