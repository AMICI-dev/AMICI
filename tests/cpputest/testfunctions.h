#ifndef TESTFUNCTIONS_H
#define TESTFUNCTIONS_H

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "include/ami_hdf5.h"

#define HDFFILE "../expectedResults.h5"
#define TEST_EPSILON 1e-12

#ifndef __APPLE__
#include <iostream>
#endif

ExpData *getTestExpData();

void checkEqualArray(const double *expected, const double *actual, int length, double epsilon);

void verifyReturnData(const char* resultPath, const ReturnData *rdata, const UserData*udata, double epsilon);

#endif
