#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <testfunctions.hpp>
#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupCalvetti)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(groupCalvetti, testSimulation) {
    amici::simulateVerifyWrite("/model_calvetti/nosensi/");
}





