#include "testfunctions.h"

#include <cstring>
#include <execinfo.h>
#include <dlfcn.h>    // for dladdr
#include <cxxabi.h>   // for __cxa_demangle
#include <cstdio>
#include <cmath>
#include <utility>
#include <unistd.h>

#include <CppUTest/TestHarness.h>
#include <CppUTestExt/MockSupport.h>

namespace amici {

namespace generic_model {

extern std::unique_ptr<amici::Model> getModel();

} // namespace generic_model

std::vector<std::string> getVariableNames(const char* name, int length)
{
    std::vector<std::string> names;
    names.resize(length);
    for (auto& it: names) {
        auto index = &it - &names[0];
        it += name + std::to_string(index);
    }
    return names;
}

void simulateVerifyWrite(const std::string& path)
{
    simulateVerifyWrite(NEW_OPTION_FILE, HDFFILE, HDFFILEWRITE, path, TEST_ATOL, TEST_RTOL);
}

void simulateVerifyWrite(std::string path, double atol, double rtol)
{
    simulateVerifyWrite(NEW_OPTION_FILE, HDFFILE, HDFFILEWRITE, std::move(path), atol, rtol);
}

void simulateWithDefaultOptions() {
    using namespace amici;
    auto model = amici::generic_model::getModel();
    auto solver = model->getSolver();
    std::unique_ptr<const ExpData> edata;
    auto rdata = runAmiciSimulation(*solver, edata.get(), *model);
}

void simulateVerifyWrite(const std::string& hdffileOptions, const std::string& hdffileResults, const std::string& hdffilewrite, const std::string& path, double atol, double rtol)
{
    using namespace amici;
    // read options from file
    std::string optionsPath = path + "/options";
    auto model = amici::generic_model::getModel();
    auto solver = model->getSolver();
    hdf5::readModelDataFromHDF5(hdffileOptions, *model, optionsPath);
    hdf5::readSolverSettingsFromHDF5(hdffileOptions, *solver, optionsPath);

    // read measurements from file
    std::string measurementPath = path + "/data";
    std::unique_ptr<const ExpData> edata;
    if(hdf5::locationExists(hdffileOptions, measurementPath))
        edata = hdf5::readSimulationExpData(hdffileOptions, measurementPath, *model);

    // simulate & verify
    auto rdata = runAmiciSimulation(*solver, edata.get(), *model);
    // perform second simulation to check reuse of solver and model object
    auto rdata_resimulation = runAmiciSimulation(*solver, edata.get(), *model);
    std::string resultPath = path + "/results";

    // write
    // delete destination group
    H5::H5File in(hdffileOptions, H5F_ACC_RDONLY);
    auto out = amici::hdf5::createOrOpenForWriting(hdffilewrite);
    if(hdf5::locationExists(out, path))
        H5Ldelete(out.getId(), path.c_str(), H5P_DEFAULT);
    hdf5::createGroup(out, path);
    std::string writePath = path + "/options";
    H5Ocopy(in.getId(), writePath.c_str(), out.getId(), writePath.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    writePath = path + "/data";
    if(edata)
        hdf5::writeSimulationExpData(*edata, out, writePath);

    writePath = path + "/results";
    hdf5::writeReturnData(*rdata, out, writePath);

    // verify simulated results
    verifyReturnData(hdffileResults, resultPath, rdata.get(), model.get(), atol, rtol);
    // verify resimulated results
    verifyReturnData(hdffileResults, resultPath, rdata_resimulation.get(), model.get(), atol, rtol);
    // verify written results
    verifyReturnData(hdffilewrite, writePath, rdata.get(), model.get(), atol, rtol);
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

    for(int i = 0; (unsigned) i < expected.size(); ++i)
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

    if(!hdf5::locationExists(hdffile, resultPath)) {
        fprintf(stderr, "ERROR: No results available for %s!\n",
                resultPath.c_str());
        return;
    }

    // compare to saved data in hdf file
    H5::H5File file(hdffile, H5F_ACC_RDONLY);

    hsize_t m, n;

    std::vector<realtype> expected;

    auto statusExp = hdf5::getIntScalarAttribute(file, resultPath, "status");
    CHECK_EQUAL(statusExp, rdata->status);

    double llhExp = hdf5::getDoubleScalarAttribute(file, resultPath, "llh");
    CHECK_TRUE(withinTolerance(llhExp, rdata->llh, atol, rtol, 1, "llh"));

    double chi2Exp = hdf5::getDoubleScalarAttribute(file, resultPath, "chi2");
    CHECK_TRUE(withinTolerance(chi2Exp, rdata->chi2, atol, rtol, 1, "chi2"));

    if(hdf5::locationExists(file, resultPath + "/x")) {
        expected = hdf5::getDoubleDataset2D(file, resultPath + "/x", m, n);
        checkEqualArray(expected, rdata->x, atol, rtol, "x");
    } else {
        CHECK_TRUE(rdata->x.empty());
    }

    //    CHECK_EQUAL(AMICI_O2MODE_FULL, udata->o2mode);

    if(hdf5::locationExists(file, resultPath + "/diagnosis/J")) {
        expected = hdf5::getDoubleDataset2D(file, resultPath + "/diagnosis/J", m, n);
        checkEqualArray(expected, rdata->J, atol, rtol, "J");
    } else {
        CHECK_TRUE(rdata->J.empty());
    }

    if(hdf5::locationExists(file, resultPath + "/y")) {
        expected = hdf5::getDoubleDataset2D(file, resultPath + "/y", m, n);
        checkEqualArray(expected, rdata->y, atol, rtol, "y");
    } else {
        CHECK_TRUE(rdata->y.empty());
    }

    if(hdf5::locationExists(file, resultPath + "/res")) {
        expected = hdf5::getDoubleDataset1D(file, resultPath + "/res");
        checkEqualArray(expected, rdata->res, atol, rtol, "res");
    } else {
        CHECK_TRUE(rdata->res.empty());
    }

    if(hdf5::locationExists(file, resultPath + "/z")) {
        expected = hdf5::getDoubleDataset2D(file, resultPath + "/z", m, n);
        checkEqualArray(expected, rdata->z, atol, rtol, "z");
    } else {
        CHECK_TRUE(rdata->z.empty());
    }

    if(hdf5::locationExists(file, resultPath + "/rz")) {
        expected = hdf5::getDoubleDataset2D(file, resultPath + "/rz", m, n);
        checkEqualArray(expected, rdata->rz, atol, rtol, "rz");
    } else {
        CHECK_TRUE(rdata->rz.empty());
    }

    if(hdf5::locationExists(file, resultPath + "/sigmaz")) {
        expected = hdf5::getDoubleDataset2D(file, resultPath + "/sigmaz", m, n);
        checkEqualArray(expected, rdata->sigmaz, atol, rtol, "sigmaz");
    } else {
        CHECK_TRUE(rdata->sigmaz.empty());
    }

    expected = hdf5::getDoubleDataset1D(file, resultPath + "/diagnosis/xdot");
    checkEqualArray(expected, rdata->xdot, atol, rtol, "xdot");

    expected = hdf5::getDoubleDataset1D(file, resultPath + "/x0");
    checkEqualArray(expected, rdata->x0, atol, rtol, "x0");

    if(rdata->sensi >= SensitivityOrder::first) {
        verifyReturnDataSensitivities(file, resultPath, rdata, model, atol, rtol);
    } else {
        CHECK_EQUAL(0, rdata->sllh.size());
        CHECK_EQUAL(0, rdata->s2llh.size());
    }
}

void verifyReturnDataSensitivities(H5::H5File const& file, std::string const& resultPath,
                                   const ReturnData *rdata, const Model *model, double atol, double rtol) {
    hsize_t m, n, o;
    std::vector<double> expected;
    if(hdf5::locationExists(file, resultPath + "/sllh")) {
        expected = hdf5::getDoubleDataset1D(file, resultPath + "/sllh");
        checkEqualArray(expected, rdata->sllh, atol, rtol, "sllh");
    } else {
        CHECK_TRUE(rdata->sllh.empty());
    }

    if(rdata->sensi_meth == SensitivityMethod::forward) {

        if(hdf5::locationExists(file, resultPath + "/sx0")) {
            expected = hdf5::getDoubleDataset2D(file, resultPath + "/sx0", m, n);
            checkEqualArray(expected, rdata->sx0, atol, rtol, "sx0");
        } else {
            CHECK_TRUE(rdata->sx0.empty());
        }

        if(hdf5::locationExists(file, resultPath + "/sres")) {
            expected = hdf5::getDoubleDataset2D(file, resultPath + "/sres", m, n);
            checkEqualArray(expected, rdata->sres, atol, rtol, "sres");
        } else {
            CHECK_TRUE(rdata->sres.empty());
        }

        if(hdf5::locationExists(file, resultPath + "/FIM")) {
            expected = hdf5::getDoubleDataset2D(file, resultPath + "/FIM", m, n);
            checkEqualArray(expected, rdata->FIM, atol, rtol, "FIM");
        } else {
            CHECK_TRUE(rdata->FIM.empty());
        }


        /* TODO REMOVE ASAP */
        if(rdata->sensi < SensitivityOrder::second) {
        /* /TODO REMOVE ASAP */

            if(hdf5::locationExists(file, resultPath + "/sx")) {
                expected = hdf5::getDoubleDataset3D(file, resultPath + "/sx", m, n, o);
                for(int ip = 0; ip < model->nplist(); ++ip)
                    checkEqualArray(&expected[ip * model->nt() * model->nxtrue_rdata],
                            &rdata->sx[ip * model->nt() * rdata->nx],
                            model->nt() * model->nxtrue_rdata, atol, rtol, "sx");
            } else {
                CHECK_TRUE(rdata->sx.empty());
            }

            if(hdf5::locationExists(file, resultPath + "/sy")) {
                expected = hdf5::getDoubleDataset3D(file, resultPath + "/sy", m, n, o);
                for(int ip = 0; ip < model->nplist(); ++ip)
                    checkEqualArray(&expected[ip * model->nt() * model->nytrue],
                            &rdata->sy[ip * model->nt() * model->ny],
                            model->nt() * model->nytrue, atol, rtol, "sy");
            } else {
                CHECK_TRUE(rdata->sy.empty());
            }


            if(model->nz>0) {
                expected = hdf5::getDoubleDataset3D(file, resultPath + "/sz", m, n, o);
                for(int ip = 0; ip < model->nplist(); ++ip)
                    checkEqualArray(&expected[ip * model->nMaxEvent() * model->nztrue],
                            &rdata->sz[ip * model->nMaxEvent() * model->nz],
                            model->nMaxEvent() * model->nztrue, atol, rtol, "sz");

                expected = hdf5::getDoubleDataset3D(file, resultPath + "/srz", m, n, o);
                for(int ip = 0; ip < model->nplist(); ++ip)
                    checkEqualArray(&expected[ip * model->nMaxEvent() * model->nztrue],
                            &rdata->srz[ip * model->nMaxEvent() * model->nz],
                            model->nMaxEvent() * model->nztrue, atol, rtol, "srz");
            }

            if(hdf5::locationExists(file, resultPath + "/ssigmay") || !rdata->ssigmay.empty()) {
                expected = hdf5::getDoubleDataset3D(file, resultPath + "/ssigmay", m, n, o);
                for(int ip = 0; ip < model->nplist(); ++ip)
                    checkEqualArray(&expected[ip * model->nt() * model->nytrue],
                            &rdata->ssigmay[ip * model->nt() * model->ny],
                            model->nt() * model->nytrue, atol, rtol, "ssigmay");
            }
            if(hdf5::locationExists(file, resultPath + "/ssigmaz") || !rdata->ssigmaz.empty()) {
                if(model->nz>0) {
                    expected = hdf5::getDoubleDataset3D(file, resultPath + "/ssigmaz", m, n, o);
                    for(int ip = 0; ip < model->nplist(); ++ip)
                        checkEqualArray(&expected[ip * model->nMaxEvent() * model->nztrue],
                                &rdata->ssigmaz[ip * model->nMaxEvent() * model->nz],
                                model->nMaxEvent() * model->nztrue, atol, rtol, "ssigmaz");
                }
            }
        /* TODO REMOVE ASAP */
        }
        /* /TODO REMOVE ASAP */

    }

    if(rdata->sensi >= SensitivityOrder::second) {
        expected = hdf5::getDoubleDataset2D(file, resultPath + "/s2llh", m, n);
        checkEqualArray(expected, rdata->s2llh, atol, rtol, "s2llh");
    } else {
        CHECK_EQUAL(0, rdata->s2llh.size());
        CHECK_EQUAL(0, rdata->s2rz.size());
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
            char *demangled = nullptr;
            int status = -1;
            if (info.dli_sname[0] == '_')
                demangled = abi::__cxa_demangle(info.dli_sname, nullptr, nullptr, &status);
            snprintf(buf, sizeof(buf), "%-3d %*p %s + %zd\n",
                     i, int(2 + sizeof(void*) * 2), callstack[i],
                     status == 0 ? demangled :
                                   info.dli_sname == nullptr ? symbols[i] : info.dli_sname,
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
