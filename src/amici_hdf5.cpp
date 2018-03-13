#include "include/amici_hdf5.h"
#include "include/amici.h"
#include "include/edata.h"
#include "include/rdata.h"
#include <amici_exception.h>

#include <hdf5_hl.h>

#include <cassert>
#include <cstring>
#include <algorithm>
#ifdef AMI_HDF5_H_DEBUG
#ifndef __APPLE__
#include <cexecinfo>
#else
#include <execinfo.h>
#endif
#endif
#include <unistd.h>
#include <cmath>


namespace amici {
namespace hdf5 {

/**
 * @brief assertMeasurementDimensionsCompatible
 * @param m
 * @param n
 * @param model
 */
void checkMeasurementDimensionsCompatible(hsize_t m, hsize_t n, Model const& model) {
    bool compatible = true;
    // if this is rank 1, n and m can be swapped
    if (n == 1) {
        compatible &= (n == (unsigned)model.nt() || n == (unsigned)model.nytrue);
        compatible &= (m == (unsigned)model.nytrue || m == (unsigned)model.nt());
        compatible &= (m * n == (unsigned)model.nytrue * model.nt());
    } else {
        compatible &= (n == (unsigned)model.nytrue);
        compatible &= (m == (unsigned)model.nt());
    }

    if(!compatible)
        throw(AmiException("HDF5 measurement data does not match model. Incompatible dimensions."));
}


/**
 * @brief assertEventDimensionsCompatible
 * @param m
 * @param n
 * @param model
 */
void checkEventDimensionsCompatible(hsize_t m, hsize_t n, Model const& model) {
    bool compatible = true;

    // if this is rank 1, n and m can be swapped
    if (n == 1) {
        compatible &= (n == (unsigned)model.nMaxEvent() ||
               n == (unsigned)model.nztrue);
        compatible &= (m == (unsigned)model.nztrue ||
               m == (unsigned)model.nMaxEvent());
        compatible &= (m * n == (unsigned)model.nytrue * model.nMaxEvent());
    } else {
        compatible &= (n == (unsigned)model.nztrue);
        compatible &= (m == (unsigned)model.nMaxEvent());
    }

    if(!compatible)
        throw(AmiException("HDF5 event data does not match model. Incompatible dimensions."));
}


void createGroup(H5::H5File &file,
                 std::string const& groupPath,
                 bool recursively) {

    hid_t groupCreationPropertyList = H5P_DEFAULT;

    if (recursively) {
        groupCreationPropertyList = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group(groupCreationPropertyList, 1);
    }

    hid_t group = H5Gcreate(file.getId(), groupPath.c_str(),
                            groupCreationPropertyList,
                            H5P_DEFAULT, H5P_DEFAULT);
    H5Pclose(groupCreationPropertyList);

    if (group < 0)
        throw(AmiException("Failed to create group in hdf5CreateGroup: %s", groupPath.c_str()));
    H5Gclose(group);
}

std::unique_ptr<ExpData> readSimulationExpData(std::string const& hdf5Filename,
                                               std::string const& hdf5Root,
                                               Model const& model) {
    H5::H5File file(hdf5Filename, H5F_ACC_RDONLY);

    hsize_t m, n;

    auto edata = std::unique_ptr<ExpData>(new ExpData(model));

    auto my = getDoubleDataset2D(file, hdf5Root + "/Y", m, n);
    if(m * n > 0) {
        checkMeasurementDimensionsCompatible(m, n, model);
        edata->my = my;
    }

    auto sigmay = getDoubleDataset2D(file, hdf5Root + "/Sigma_Y", m, n);
    if(m * n > 0) {
        checkMeasurementDimensionsCompatible(m, n, model);
        edata->sigmay = sigmay;
    }

    if (model.nz) {
        edata->mz = getDoubleDataset2D(file, hdf5Root + "/Z", m, n);
        checkEventDimensionsCompatible(m, n, model);

        edata->sigmaz = getDoubleDataset2D(file, hdf5Root + "/Sigma_Z", m, n);
        checkEventDimensionsCompatible(m, n, model);
    }

    return edata;
}

void writeReturnData(const ReturnData &rdata, H5::H5File &file, const std::string &hdf5Location)
{

    if(!locationExists(file, hdf5Location))
        createGroup(file, hdf5Location);

    if (rdata.ts.size())
        createAndWriteDouble1DDataset(file, hdf5Location + "/t",
                                      rdata.ts.data(), rdata.nt);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(), "llh", &rdata.llh, 1);
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(), "chi2", &rdata.chi2, 1);

    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(), "status", &rdata.status, 1);

    if (rdata.sllh.size())
        createAndWriteDouble1DDataset(file, hdf5Location + "/sllh",
                                      rdata.sllh.data(), rdata.nplist);

    if (rdata.x0.size())
        createAndWriteDouble1DDataset(file, hdf5Location + "/x0", rdata.x0.data(), rdata.nx);

    if (rdata.x.size())
        createAndWriteDouble2DDataset(file, hdf5Location + "/x", rdata.x.data(),
                                        rdata.nt, rdata.nx);

    if (rdata.y.size())
        createAndWriteDouble2DDataset(file, hdf5Location + "/y", rdata.y.data(),
                                        rdata.nt, rdata.ny);

    if (rdata.z.size())
        createAndWriteDouble2DDataset(file, hdf5Location + "/z", rdata.z.data(),
                                        rdata.nmaxevent, rdata.nz);

    if (rdata.rz.size())
        createAndWriteDouble2DDataset(file, hdf5Location + "/rz", rdata.rz.data(),
                                        rdata.nmaxevent, rdata.nz);

    if (rdata.sigmay.size())
        createAndWriteDouble2DDataset(file, hdf5Location + "/sigmay", rdata.sigmay.data(), rdata.nt, rdata.ny);

    if (rdata.sigmaz.size())
        createAndWriteDouble2DDataset(file, hdf5Location + "/sigmaz", rdata.sigmaz.data(), rdata.nmaxevent, rdata.nz);

    if (rdata.s2llh.size())
        createAndWriteDouble2DDataset(file, hdf5Location + "/s2llh", rdata.s2llh.data(),
                                        rdata.nJ - 1, rdata.nplist);

    if (rdata.sx0.size())
        createAndWriteDouble2DDataset(file, hdf5Location + "/sx0", rdata.sx0.data(), rdata.nx, rdata.nplist);

    if (rdata.sx.size())
        createAndWriteDouble3DDataset(file, hdf5Location + "/sx", rdata.sx.data(), rdata.nt, rdata.nplist, rdata.nx);

    if (rdata.sy.size())
        createAndWriteDouble3DDataset(file, hdf5Location + "/sy", rdata.sy.data(), rdata.nt, rdata.nplist, rdata.ny);

    if (rdata.ssigmay.size())
        createAndWriteDouble3DDataset(file, hdf5Location + "/ssigmay",
                                        rdata.ssigmay.data(), rdata.nt,
                                        rdata.nplist, rdata.ny);

    if (rdata.sz.size())
        createAndWriteDouble3DDataset(file, hdf5Location + "/sz", rdata.sz.data(),
                                        rdata.nmaxevent, rdata.nplist, rdata.nz);

    if (rdata.srz.size())
        createAndWriteDouble3DDataset(file, hdf5Location + "/srz", rdata.srz.data(),
                                        rdata.nmaxevent, rdata.nplist, rdata.nz);

    if (rdata.ssigmaz.size())
        createAndWriteDouble3DDataset(file, hdf5Location + "/ssigmaz", rdata.ssigmaz.data(),
                                      rdata.nmaxevent, rdata.nplist, rdata.nz);

    writeReturnDataDiagnosis(rdata, file, hdf5Location + "/diagnosis");
}

void writeReturnDataDiagnosis(const ReturnData &rdata,
                              H5::H5File& file,
                              const std::string& hdf5Location) {

    if(!locationExists(file, hdf5Location))
        createGroup(file, hdf5Location);

    if (rdata.xdot.size())
        createAndWriteDouble1DDataset(file, hdf5Location + "/xdot",
                                      rdata.xdot.data(), rdata.nx);

    if (rdata.numsteps.size())
        createAndWriteInt1DDataset(file, hdf5Location + "/numsteps", rdata.numsteps);

    if (rdata.numrhsevals.size())
        createAndWriteInt1DDataset(file, hdf5Location + "/numrhsevals", rdata.numrhsevals);

    if (rdata.numerrtestfails.size())
        createAndWriteInt1DDataset(file, hdf5Location + "/numerrtestfails", rdata.numerrtestfails);

    if (rdata.numnonlinsolvconvfails.size())
        createAndWriteInt1DDataset(file, hdf5Location + "/numnonlinsolvconvfails",
                                   rdata.numnonlinsolvconvfails);

    if (rdata.order.size())
        createAndWriteInt1DDataset(file, hdf5Location + "/order", rdata.order);

    if (rdata.numstepsB.size())
        createAndWriteInt1DDataset(file, hdf5Location + "/numstepsB", rdata.numstepsB);

    if (rdata.numrhsevalsB.size())
        createAndWriteInt1DDataset(file, hdf5Location + "/numrhsevalsB", rdata.numrhsevalsB);

    if (rdata.numerrtestfailsB.size())
        createAndWriteInt1DDataset(file, hdf5Location + "/numerrtestfailsB", rdata.numerrtestfailsB);

    if (rdata.numnonlinsolvconvfailsB.size())
        createAndWriteInt1DDataset(file, hdf5Location + "/numnonlinsolvconvfailsB", rdata.numnonlinsolvconvfailsB);

    H5LTset_attribute_int(file.getId(), (hdf5Location + "").c_str(), "newton_status", &rdata.newton_status, 1);

    if (rdata.newton_numsteps.size())
        createAndWriteInt1DDataset(file, hdf5Location + "/newton_numsteps", rdata.newton_numsteps);

    if (rdata.newton_numlinsteps.size())
        createAndWriteInt2DDataset(file, hdf5Location + "/newton_numlinsteps", rdata.newton_numlinsteps.data(),
                                    rdata.newton_maxsteps, 2);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(), "newton_time", &rdata.newton_time, 1);

    if (rdata.J.size())
        createAndWriteDouble2DDataset(file, hdf5Location + "/J", rdata.J.data(),
                                        rdata.nx, rdata.nx);

}


void writeReturnData(ReturnData const& rdata,
                     std::string const& hdf5Filename,
                     std::string const& hdf5Location) {
    auto file = createOrOpenForWriting(hdf5Filename);

    writeReturnData(rdata, file, hdf5Location);
}

double getDoubleScalarAttribute(H5::H5File const& file,
                                std::string const& optionsObject,
                                std::string const& attributeName) {

    double data = NAN;
    herr_t status = H5LTget_attribute_double(file.getId(), optionsObject.c_str(),
                                             attributeName.c_str(), &data);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %e\n", attributeName, data);
#endif

    if(status < 0)
        throw AmiException("Attribute %s not found for object %s.",
                           attributeName.c_str(), optionsObject.c_str());

    return data;
}

int getIntScalarAttribute(H5::H5File const& file,
                          std::string const& optionsObject,
                          std::string const& attributeName) {
    int data = 0;
    herr_t status = H5LTget_attribute_int(file.getId(), optionsObject.c_str(),
                                          attributeName.c_str(), &data);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d\n", attributeName.c_str(), data);
#endif

    if(status < 0)
        throw AmiException("Attribute %s not found for object %s.",
                           attributeName.c_str(), optionsObject.c_str());

    return data;
}


void createAndWriteInt1DDataset(H5::H5File& file,
                                     std::string const& datasetName,
                                     const int *buffer, hsize_t m) {
    H5::DataSpace dataspace(1, &m);
    auto dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_INT, dataspace);
    dataset.write(buffer, H5::PredType::NATIVE_INT);
}

void createAndWriteInt1DDataset(H5::H5File& file,
                                     std::string const& datasetName,
                                std::vector<int> const& buffer) {
    createAndWriteInt1DDataset(file, datasetName, buffer.data(), buffer.size());
}

void createAndWriteDouble1DDataset(H5::H5File& file,
                                     std::string const& datasetName,
                                     const double *buffer, hsize_t m) {
    H5::DataSpace dataspace(1, &m);
    auto dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(buffer, H5::PredType::NATIVE_DOUBLE);
}

void createAndWriteDouble2DDataset(H5::H5File& file,
                                     std::string const& datasetName,
                                     const double *buffer, hsize_t m,
                                     hsize_t n) {
    const hsize_t adims[] {m, n};
    H5::DataSpace dataspace(2, adims);
    auto dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(buffer, H5::PredType::NATIVE_DOUBLE);
}

void createAndWriteInt2DDataset(H5::H5File& file,
                                     std::string const& datasetName,
                                     const int *buffer, hsize_t m,
                                     hsize_t n) {
    const hsize_t adims[] {m, n};
    H5::DataSpace dataspace(2, adims);
    auto dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_INT, dataspace);
    dataset.write(buffer, H5::PredType::NATIVE_INT);
}

void createAndWriteDouble3DDataset(H5::H5File& file,
                                     std::string const& datasetName,
                                     const double *buffer, hsize_t m,
                                     hsize_t n, hsize_t o) {
    const hsize_t adims[] {m, n, o};
    H5::DataSpace dataspace(3, adims);
    auto dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(buffer, H5::PredType::NATIVE_DOUBLE);
}


bool attributeExists(H5::H5File const& file,
                     const std::string &optionsObject,
                     const std::string &attributeName) {
    AMICI_H5_SAVE_ERROR_HANDLER;
    int result = H5Aexists_by_name(file.getId(), optionsObject.c_str(), attributeName.c_str(), H5P_DEFAULT);
    AMICI_H5_RESTORE_ERROR_HANDLER;
    return result > 0;
}

bool attributeExists(H5::H5Object const& object,
                     const std::string &attributeName) {
    AMICI_H5_SAVE_ERROR_HANDLER;
    int result = H5Aexists(object.getId(), attributeName.c_str());
    AMICI_H5_RESTORE_ERROR_HANDLER;
    return result > 0;
}

void readSolverSettingsFromHDF5(H5::H5File const& file, Solver &solver, const std::string &datasetPath) {

    if(attributeExists(file, datasetPath, "atol")) {
        solver.setAbsoluteTolerance(getDoubleScalarAttribute(file, datasetPath, "atol"));
    }

    if(attributeExists(file, datasetPath, "rtol")) {
        solver.setRelativeTolerance(getDoubleScalarAttribute(file, datasetPath, "rtol"));
    }

    if(attributeExists(file, datasetPath, "maxsteps")) {
        solver.setMaxSteps(getIntScalarAttribute(file, datasetPath, "maxsteps"));
    }

    if(attributeExists(file, datasetPath, "lmm")) {
        solver.setLinearMultistepMethod(
                    static_cast<LinearMultistepMethod>(getIntScalarAttribute(file, datasetPath, "lmm")));
    }

    if(attributeExists(file, datasetPath, "iter")) {
        solver.setNonlinearSolverIteration(
                    static_cast<NonlinearSolverIteration>(getIntScalarAttribute(file, datasetPath, "iter")));
    }

    if(attributeExists(file, datasetPath, "stldet")) {
        solver.setStabilityLimitFlag(getIntScalarAttribute(file, datasetPath, "stldet"));
    }

    if(attributeExists(file, datasetPath, "ordering")) {
        solver.setStateOrdering(static_cast<StateOrdering>(
                                    getIntScalarAttribute(file, datasetPath, "ordering")));
    }

    if(attributeExists(file, datasetPath, "interpType")) {
        solver.setInterpolationType(static_cast<InterpolationType>(getIntScalarAttribute(file, datasetPath, "interpType")));
    }

    if(attributeExists(file, datasetPath, "sensi_meth")) {
        solver.setSensitivityMethod(static_cast<AMICI_sensi_meth>(getIntScalarAttribute(file, datasetPath, "sensi_meth")));
    }

    if(attributeExists(file, datasetPath, "sensi")) {
        solver.setSensitivityOrder(static_cast<AMICI_sensi_order>(getIntScalarAttribute(file, datasetPath, "sensi")));
    }

    if(attributeExists(file, datasetPath, "newton_maxsteps")) {
        solver.setNewtonMaxSteps(getIntScalarAttribute(file, datasetPath, "newton_maxsteps"));
    }

    if(attributeExists(file, datasetPath, "newton_preeq")) {
        solver.setNewtonPreequilibration(getIntScalarAttribute(file, datasetPath, "newton_preeq"));
    }

    if(attributeExists(file, datasetPath, "newton_precon")) {
        solver.setNewtonPreconditioner(getIntScalarAttribute(file, datasetPath, "newton_precon"));
    }

    if(attributeExists(file, datasetPath, "newton_maxlinsteps")) {
        solver.setNewtonMaxLinearSteps(getIntScalarAttribute(file, datasetPath, "newton_maxlinsteps"));
    }

    if(attributeExists(file, datasetPath, "linsol")) {
        solver.setLinearSolver(static_cast<LinearSolver>(getIntScalarAttribute(file, datasetPath, "linsol")));
    }

    if(attributeExists(file, datasetPath, "ism")) {
        solver.setInternalSensitivityMethod(static_cast<InternalSensitivityMethod>(getIntScalarAttribute(file, datasetPath, "ism")));
    }
}

void readSolverSettingsFromHDF5(const std::string &hdffile, Solver &solver, const std::string &datasetPath) {
    H5::H5File file(hdffile, H5F_ACC_RDONLY, H5P_DEFAULT);

    readSolverSettingsFromHDF5(file, solver, datasetPath);
}

void readModelDataFromHDF5(const std::string &hdffile, Model &model, const std::string &datasetPath) {
    H5::H5File file(hdffile, H5F_ACC_RDONLY, H5P_DEFAULT);

    readModelDataFromHDF5(file, model, datasetPath);
}

void readModelDataFromHDF5(const H5::H5File &file, Model &model, const std::string &datasetPath) {
    if(attributeExists(file, datasetPath, "tstart")) {
        model.setT0(getDoubleScalarAttribute(file, datasetPath, "tstart"));
    }

    if(locationExists(file, datasetPath + "/pscale")) {
        auto pscaleInt = getIntDataset1D(file, datasetPath + "/pscale");
        std::vector<AMICI_parameter_scaling> pscale(pscaleInt.size());
        for(int i = 0; (unsigned)i < pscaleInt.size(); ++i)
            pscale[i] = static_cast<AMICI_parameter_scaling>(pscaleInt[i]);
        model.setParameterScale(pscale);
    } else if (attributeExists(file, datasetPath, "pscale")) {
        // if pscale is the same for all parameters,
        // it can be set as scalar attribute for convenience
        model.setParameterScale(static_cast<AMICI_parameter_scaling>(getDoubleScalarAttribute(file,  datasetPath, "pscale")));
    }

    if(attributeExists(file, datasetPath, "nmaxevent")) {
        model.setNMaxEvent(getIntScalarAttribute(file, datasetPath, "nmaxevent"));
    }

    if(locationExists(file, datasetPath + "/qpositivex")) {
        auto qPosX = getIntDataset1D(file, datasetPath + "/qpositivex");
        model.setPositivityFlag(qPosX);
    }

    if(locationExists(file, datasetPath + "/theta")) {
        model.setParameters(getDoubleDataset1D(file, datasetPath + "/theta"));
    }

    if(locationExists(file, datasetPath + "/kappa")) {
        model.setFixedParameters(getDoubleDataset1D(file, datasetPath + "/kappa"));
    }

    if(locationExists(file, datasetPath + "/ts")) {
        model.setTimepoints(getDoubleDataset1D(file, datasetPath + "/ts"));
    }

    if(locationExists(file, datasetPath + "/sens_ind")) {
        auto sensInd = getIntDataset1D(file, datasetPath + "/sens_ind");
        model.setParameterList(sensInd);
    }

    if(locationExists(file, datasetPath + "/x0")) {
        auto x0 = getDoubleDataset1D(file, datasetPath + "/x0");
        if(x0.size())
            model.setInitialStates(x0);
    }

    if(locationExists(file, datasetPath + "/sx0")) {
        hsize_t length0 = 0;
        hsize_t length1 = 0;
        auto sx0 = getDoubleDataset2D(file, datasetPath + "/sx0", length0, length1);
        if(sx0.size()) {
            if (length0 != (unsigned) model.nx && length1 != (unsigned) model.nplist())
                throw(AmiException("Dimension mismatch when reading sx0. Expected %dx%d, got %llu, %llu.",
                                   model.nx, model.nplist(), length0, length1));
            model.setInitialStateSensitivities(sx0);
        }
    }

    if(locationExists(file, datasetPath + "/pbar")) {
        auto pbar = getDoubleDataset1D(file, datasetPath + "/pbar");
        model.setParameterScaling(pbar);
    }
}

H5::H5File createOrOpenForWriting(const std::string &hdf5filename)
{
    H5::H5File file;
    AMICI_H5_SAVE_ERROR_HANDLER;
    try {
        file = H5::H5File(hdf5filename, H5F_ACC_RDWR);
        AMICI_H5_RESTORE_ERROR_HANDLER;
    } catch(...) {
        AMICI_H5_RESTORE_ERROR_HANDLER;
        file = H5::H5File(hdf5filename, H5F_ACC_EXCL);
    }

    return file;
}

bool locationExists(const H5::H5File &file, const std::string &location)
{
    AMICI_H5_SAVE_ERROR_HANDLER;
    auto result = H5Lexists(file.getId(), location.c_str(), H5P_DEFAULT) > 0;
    AMICI_H5_RESTORE_ERROR_HANDLER;
    return result;
}

bool locationExists(const std::string &filename, const std::string &location)
{
    H5::H5File file(filename, H5F_ACC_RDONLY);
    return locationExists(file, location);
}

std::vector<int> getIntDataset1D(const H5::H5File &file,
                                 std::string const& name) {
    auto dataset = file.openDataSet(name);
    auto dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if(rank != 1)
        throw(AmiException("Expected array of rank 1 in %s", name.c_str()));

    hsize_t dim;
    dataspace.getSimpleExtentDims(&dim);
    std::vector<int> result(dim);
    if(result.size())
        dataset.read(result.data(), H5::PredType::NATIVE_INT);
    return result;
}


std::vector<double> getDoubleDataset1D(const H5::H5File &file, const std::string &name)
{
    auto dataset = file.openDataSet(name);
    auto dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if(rank != 1)
        throw(AmiException("Expected array of rank 1 in %s", name.c_str()));

    hsize_t dim;
    dataspace.getSimpleExtentDims(&dim);
    std::vector<double> result(dim);
    if(result.size())
        dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE);

    return result;

}

std::vector<double> getDoubleDataset2D(const H5::H5File &file, const std::string &name, hsize_t &m, hsize_t &n)
{
    m = n = 0;

    auto dataset = file.openDataSet(name);
    auto dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if(rank != 2)
        throw(AmiException("Expected array of rank 2 in %s", name.c_str()));

    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims);
    m = dims[0];
    n = dims[1];

    std::vector<double> result(m * n);
    if(result.size())
        dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE);

    return result;
}

std::vector<double> getDoubleDataset3D(const H5::H5File &file, const std::string &name, hsize_t &m, hsize_t &n, hsize_t &o)
{
    m = n = o = 0;

    auto dataset = file.openDataSet(name);
    auto dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if(rank != 3)
        throw(AmiException("Expected array of rank 3 in %s", name.c_str()));

    hsize_t dims[3];
    dataspace.getSimpleExtentDims(dims);
    m = dims[0];
    n = dims[1];
    o = dims[2];

    std::vector<double> result(m * n * o);
    if(result.size())
        dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE);

    return result;
}

} // namespace hdf5
} // namespace amici
