#include <time.h>
#include <stdlib.h>

#include "CppUTest/CommandLineTestRunner.h"
#include "CppUTest/TestHarness.h"

int main(int argc, char** argv)
{
    srand(time(NULL));
    return CommandLineTestRunner::RunAllTests(argc, argv);
}
