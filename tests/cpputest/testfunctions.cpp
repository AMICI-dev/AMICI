#include "testfunctions.h"

#include <include/amici_hdf5.h>
#include <amici.h>

#include <cstring>
#include <execinfo.h>
#include <dlfcn.h>    // for dladdr
#include <cxxabi.h>   // for __cxa_demangle
#include <cstdio>
#include <cmath>
#include <unistd.h>

extern std::unique_ptr<amici::Model> getModel();

namespace amici {

void simulateAndVerifyFromFile(const std::string path)
{
    simulateAndVerifyFromFile(NEW_OPTION_FILE, HDFFILE, path, TEST_ATOL, TEST_RTOL);
}

void simulateAndVerifyFromFile(std::string path, double atol, double rtol)
{
    simulateAndVerifyFromFile(NEW_OPTION_FILE, HDFFILE, path, atol, rtol);
}

void simulateAndWriteToFile(const std::string path)
{
    simulateAndWriteToFile(HDFFILE, HDFFILEWRITE, path, TEST_ATOL, TEST_RTOL);
}

void simulateAndWriteToFile(std::string path, double atol, double rtol)
{
    simulateAndWriteToFile(HDFFILE, HDFFILEWRITE, path, atol, rtol);
}


void simulateAndVerifyFromFile(const std::string hdffileOptions, const std::string hdffileResults, std::string path, double atol, double rtol)
{
    using namespace amici;
    // read options from file
    std::string optionsPath = path + "/options";
    auto model = getModel();
    auto solver = model->getSolver();
    hdf5::readModelDataFromHDF5(hdffileOptions, *model, optionsPath);
    hdf5::readSolverSettingsFromHDF5(hdffileOptions, *solver, optionsPath);

    // read measurements from file
    std::string measurementPath = path + "/data";

    std::unique_ptr<const ExpData> edata;
    if(hdf5::locationExists(hdffileOptions, measurementPath))
        edata = hdf5::readSimulationExpData(hdffileResults, measurementPath, *model);

    // simulate & verify
    auto rdata = std::unique_ptr<ReturnData>(getSimulationResults(*model, edata.get(), *solver));
    std::string resultPath = path + "/results";
    verifyReturnData(hdffileResults.c_str(), resultPath.c_str(), rdata.get(), model.get(), atol, rtol);
}

void simulateAndWriteToFile(const std::string hdffile, const std::string hdffilewrite, std::string path, double atol, double rtol)
{
    // read simulation options
    std::string optionsPath = path + "/options";
    auto model = getModel();
    auto solver = model->getSolver();

    hdf5::readModelDataFromHDF5(hdffile, *model, optionsPath);
    hdf5::readSolverSettingsFromHDF5(hdffile, *solver, optionsPath);

    std::string measurementPath = path + "/data";
    std::unique_ptr<const ExpData> edata;
    if(hdf5::locationExists(hdffile, measurementPath))
        edata = hdf5::readSimulationExpData(hdffile, measurementPath, *model);

    auto rdata = std::unique_ptr<ReturnData>(getSimulationResults(*model, edata.get(), *solver));

    std::string writePath = path + "/write";
    hdf5::writeReturnData(*rdata, hdffilewrite, writePath);
    verifyReturnData(hdffilewrite, writePath, rdata.get(), model.get(), atol, rtol);
    remove(hdffilewrite.c_str());
}


std::unique_ptr<ExpData> getTestExpData(Model const& model) {
    return std::unique_ptr<ExpData>(new ExpData(model));
}

bool withinTolerance(double expected, double actual, double atol, double rtol, int index, const char *name) {
    bool withinTol =  fabs(expected - actual) <= atol || fabs((expected - actual) / (rtol + expected)) <= rtol;

    if(!withinTol && std::isnan(expected) && std::isnan(actual))
        withinTol = true;

    if(!withinTol && std::isinf(expected) && std::isinf(actual))
        withinTol = true;

    if(!withinTol) {
        fprintf(stderr, "ERROR: Expected value %e, but was %e in %s at index %d.\n",expected, actual, name, index);
        fprintf(stderr, "       Relative error: %e (tolerance was %e)\n", fabs((expected - actual) / (rtol + expected)), rtol);
        fprintf(stderr, "       Absolute error: %e (tolerance was %e)\n", fabs(expected - actual), atol);
        printBacktrace(12);
    }

    return withinTol;
}

void checkEqualArray(std::vector<double> const& expected, std::vector<double> const& actual,
                     double atol, double rtol, std::string const& name) {
    CHECK_EQUAL(expected.size(), actual.size());

    for(int i = 0; i < expected.size(); ++i)
    {
        bool withinTol = withinTolerance(expected[i], actual[i], atol, rtol, i, name.c_str());
        CHECK_TRUE(withinTol);
    }
}

void checkEqualArray(const double *expected, const double *actual, const int length, double atol, double rtol, const char *name) {
    if(!expected && !actual)
        return;
    if(!length)
        return;

    CHECK_TRUE(expected && actual);

    for(int i = 0; i < length; ++i)
    {
        bool withinTol = withinTolerance(expected[i], actual[i], atol, rtol, i, name);
        CHECK_TRUE(withinTol);
    }
}

void checkEqualArrayStrided(const double *expected, const double *actual, int length, int strideExpected, int strideActual, double atol, double rtol, const char *name) {
    if(!expected && !actual)
        return;

    CHECK_TRUE(expected && actual);

    for(int i = 0; i < length; ++i)
    {
        bool withinTol = withinTolerance(expected[i * strideExpected], actual[i * strideActual], atol, rtol, i, name);
        CHECK_TRUE(withinTol);
    }
}

void verifyReturnData(std::string const& hdffile, std::string const& resultPath,
                      const ReturnData *rdata, const Model *model, double atol, double rtol) {
    CHECK_FALSE(rdata == nullptr);

    // compare to saved data in hdf file
    H5::H5File file(hdffile, H5F_ACC_RDONLY);

    hsize_t m, n;


    double statusExp = hdf5::getDoubleScalarAttribute(file, resultPath, "status");
    CHECK_EQUAL((int) statusExp, *rdata->status);

    double llhExp = hdf5::getDoubleScalarAttribute(file, resultPath, "llh");
    CHECK_TRUE(withinTolerance(llhExp, *rdata->llh, atol, rtol, 1, "llh"));

    auto expected = hdf5::getDoubleArrayAttribute2D(file, resultPath, "x", m, n);
    checkEqualArray(expected.data(), rdata->x, model->nt() * model->nxtrue, atol, rtol, "x");

    //    CHECK_EQUAL(AMICI_O2MODE_FULL, udata->o2mode);

    if(hdf5::attributeExists(file, resultPath, "J")) {
        expected = hdf5::getDoubleArrayAttribute2D(file, resultPath, "J", m, n);
        checkEqualArray(expected.data(), rdata->J, model->nx * model->nx, atol, rtol, "J");
    }

    expected = hdf5::getDoubleArrayAttribute2D(file, resultPath, "y", m, n);
    checkEqualArray(expected.data(), rdata->y, model->nt() * model->nytrue, atol, rtol, "y");

    if(model->nz>0) {
        expected = hdf5::getDoubleArrayAttribute2D(file, resultPath, "z", m, n);
        checkEqualArray(expected.data(), rdata->z, model->nMaxEvent() * model->nztrue, atol, rtol, "z");

        expected = hdf5::getDoubleArrayAttribute2D(file, resultPath, "rz", m, n);
        checkEqualArray(expected.data(), rdata->rz, model->nMaxEvent() * model->nztrue, atol, rtol, "rz");

        expected = hdf5::getDoubleArrayAttribute2D(file, resultPath, "sigmaz", m, n);
        checkEqualArray(expected.data(), rdata->sigmaz, model->nMaxEvent() * model->nztrue, atol, rtol, "sigmaz");
    }

    expected = hdf5::getDoubleArrayAttribute2D(file, resultPath, "xdot", m, n);
    checkEqualArray(expected.data(), rdata->xdot, model->nxtrue, atol, rtol, "xdot");

    if(rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        verifyReturnDataSensitivities(file, resultPath, rdata, model, atol, rtol);
    } else {
        POINTERS_EQUAL(NULL, rdata->sllh);
        POINTERS_EQUAL(NULL, rdata->s2llh);
    }
}

void verifyReturnDataSensitivities(H5::H5File const& file, std::string const& resultPath,
                                   const ReturnData *rdata, const Model *model, double atol, double rtol) {
    hsize_t m, n, o;
    int status;

    auto expected = hdf5::getDoubleArrayAttribute2D(file, resultPath, "sllh", m, n);
    checkEqualArray(expected.data(), rdata->sllh, rdata->nplist, atol, rtol, "sllh");

    if(rdata->sensi_meth == AMICI_SENSI_FSA) {
        expected = hdf5::getDoubleArrayAttribute3D(file, resultPath, "sx", m, n, o);
        for(int ip = 0; ip < model->nplist(); ++ip)
            checkEqualArray(&expected[ip * model->nt() * model->nxtrue],
                    &rdata->sx[ip * model->nt() * model->nx],
                    model->nt() * model->nxtrue, atol, rtol, "sx");

        expected = hdf5::getDoubleArrayAttribute3D(file, resultPath, "sy", m, n, o);
        for(int ip = 0; ip < model->nplist(); ++ip)
            checkEqualArray(&expected[ip * model->nt() * model->nytrue],
                    &rdata->sy[ip * model->nt() * model->ny],
                    model->nt() * model->nytrue, atol, rtol, "sy");


        if(model->nz>0) {
            expected = hdf5::getDoubleArrayAttribute3D(file, resultPath, "sz", m, n, o);
            for(int ip = 0; ip < model->nplist(); ++ip)
                checkEqualArray(&expected[ip * model->nMaxEvent() * model->nztrue],
                        &rdata->sz[ip * model->nMaxEvent() * model->nz],
                        model->nMaxEvent() * model->nztrue, atol, rtol, "sz");

            expected = hdf5::getDoubleArrayAttribute3D(file, resultPath, "srz", m, n, o);
            for(int ip = 0; ip < model->nplist(); ++ip)
                checkEqualArray(&expected[ip * model->nMaxEvent() * model->nztrue],
                        &rdata->srz[ip * model->nMaxEvent() * model->nz],
                        model->nMaxEvent() * model->nztrue, atol, rtol, "srz");
        }

        expected = hdf5::getDoubleArrayAttribute3D(file, resultPath, "ssigmay", m, n, o);
        for(int ip = 0; ip < model->nplist(); ++ip)
            checkEqualArray(&expected[ip * model->nt() * model->nytrue],
                    &rdata->ssigmay[ip * model->nt() * model->ny],
                    model->nt() * model->nytrue, atol, rtol, "ssigmay");

        if(model->nz>0) {
            expected = hdf5::getDoubleArrayAttribute3D(file, resultPath, "ssigmaz", m, n, o);
            for(int ip = 0; ip < model->nplist(); ++ip)
                checkEqualArray(&expected[ip * model->nMaxEvent() * model->nztrue],
                        &rdata->ssigmaz[ip * model->nMaxEvent() * model->nz],
                        model->nMaxEvent() * model->nztrue, atol, rtol, "ssigmaz");
        }
    }

    if(rdata->sensi >= AMICI_SENSI_ORDER_SECOND) {
        expected = hdf5::getDoubleArrayAttribute2D(file, resultPath, "s2llh", m, n);
        checkEqualArray(expected.data(), rdata->s2llh, (model->nJ-1) * model->nplist(), atol, rtol, "s2llh");
    } else {
        POINTERS_EQUAL(nullptr, rdata->s2llh);
        POINTERS_EQUAL(nullptr, rdata->s2rz);
    }

}


void printBacktrace(const int nMaxFrames) {
    void *callstack[nMaxFrames];
    char buf[1024];
    int nFrames = backtrace(callstack, nMaxFrames);
    char **symbols = backtrace_symbols(callstack, nFrames);

    std::ostringstream trace_buf;
    for (int i = 0; i < nFrames; i++) {
        Dl_info info;
        if (dladdr(callstack[i], &info) && info.dli_sname) {
            char *demangled = NULL;
            int status = -1;
            if (info.dli_sname[0] == '_')
                demangled = abi::__cxa_demangle(info.dli_sname, NULL, 0, &status);
            snprintf(buf, sizeof(buf), "%-3d %*p %s + %zd\n",
                     i, int(2 + sizeof(void*) * 2), callstack[i],
                     status == 0 ? demangled :
                                   info.dli_sname == 0 ? symbols[i] : info.dli_sname,
                     (char *)callstack[i] - (char *)info.dli_saddr);
            free(demangled);
        } else {
            snprintf(buf, sizeof(buf), "%-3d %*p %s\n",
                     i, int(2 + sizeof(void*) * 2), callstack[i], symbols[i]);
        }
        trace_buf << buf;
    }
    free(symbols);
    if (nFrames == nMaxFrames)
        trace_buf << "[truncated]\n";
    printf("%s\n",trace_buf.str().c_str());
}

} // namespace amici
