#ifndef TESTFUNCTIONS_H
#define TESTFUNCTIONS_H

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include <hdf5.h>
#include <string>

#ifndef __APPLE__
#include <iostream>
#endif

namespace amici {

class UserData;
class ReturnData;
class ExpData;
class Model;

#define HDFFILE "../../expectedResults.h5"
#define TEST_ATOL 1e-10
#define TEST_RTOL 1e-05


void simulateAndVerifyFromFile(Model *model, const std::string path);

void simulateAndVerifyFromFile(Model *model, std::string path, double atol, double rtol);

void simulateAndVerifyFromFile(Model *model, const std::string hdffile, std::string path, double atol, double rtol);

ExpData *getTestExpData(const UserData *udata, Model *model);

bool withinTolerance(double expected, double actual, double atol, double rtol, int index);

void checkEqualArray(const double *expected, const double *actual, int length, double atol, double rtol);

void verifyReturnData(const char *hdffile, const char* resultPath, const ReturnData *rdata, const UserData*udata, const Model *model, double atol, double rtol);

void verifyReturnDataSensitivities(hid_t file_id, const char* resultPath, const ReturnData *rdata, const UserData*udata, const Model *model, double atol, double rtol);

void printBacktrace(int depth);

} // namespace amici

#endif
