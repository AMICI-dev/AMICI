#ifndef TESTFUNCTIONS_H
#define TESTFUNCTIONS_H

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include <hdf5.h>

class UserData;
class ReturnData;
class ExpData;

#define HDFFILE "../../expectedResults.h5"
#define TEST_ATOL 1e-10
#define TEST_RTOL 1e-05

#ifndef __APPLE__
#include <iostream>
#endif

ExpData *getTestExpData(const UserData *udata);

bool withinTolerance(double expected, double actual, double atol, double rtol);

void checkEqualArray(const double *expected, const double *actual, int length, double atol, double rtol);

void verifyReturnData(const char* resultPath, const ReturnData *rdata, const UserData*udata, double atol, double rtol);

void verifyReturnDataSensitivities(hid_t file_id, const char* resultPath, const ReturnData *rdata, const UserData*udata, double atol, double rtol);

void printBacktrace(int depth);
#endif
