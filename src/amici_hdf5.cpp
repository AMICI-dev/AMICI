#include "include/amici_hdf5.h"
#include "include/amici.h"
#include "include/amici_model.h"
#include "include/edata.h"
#include "include/rdata.h"
#include <amici_solver.h>
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
void assertMeasurementDimensionsCompatible(hsize_t m, hsize_t n, Model const& model) {
    // if this is rank 1, n and m can be swapped
    if (n == 1) {
        assert(n == (unsigned)model.nt() || n == (unsigned)model.nytrue);
        assert(m == (unsigned)model.nytrue || m == (unsigned)model.nt());
        assert(m * n == (unsigned)model.nytrue * model.nt());
    } else {
        assert(n == (unsigned)model.nt());
        assert(m == (unsigned)model.nytrue);
    }
}


/**
 * @brief assertEventDimensionsCompatible
 * @param m
 * @param n
 * @param model
 */
void assertEventDimensionsCompatible(hsize_t m, hsize_t n, Model const& model) {
    // if this is rank 1, n and m can be swapped
    if (n == 1) {
        assert(n == (unsigned)model.nMaxEvent() ||
               n == (unsigned)model.nztrue);
        assert(m == (unsigned)model.nztrue ||
               m == (unsigned)model.nMaxEvent());
        assert(m * n == (unsigned)model.nytrue * model.nMaxEvent());
    } else {
        assert(n == (unsigned)model.nMaxEvent());
        assert(m == (unsigned)model.nztrue);
    }
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

    auto my = getDoubleArrayAttribute2D(file, hdf5Root, "Y", m, n);
    if(m * n > 0) {
        assertMeasurementDimensionsCompatible(m, n, model);
        edata->my = my;
    }

    auto sigmay = getDoubleArrayAttribute2D(file, hdf5Root, "Sigma_Y", m, n);
    if(m * n > 0) {
        assertMeasurementDimensionsCompatible(m, n, model);
        edata->sigmay = sigmay;
    }

    if (model.nz) {
        edata->mz = getDoubleArrayAttribute2D(file, hdf5Root, "Z", m, n);
        assertEventDimensionsCompatible(m, n, model);

        edata->sigmaz = getDoubleArrayAttribute2D(file, hdf5Root, "Sigma_Z", m, n);
        assertEventDimensionsCompatible(m, n, model);
    }

    return edata;
}

void writeReturnData(const ReturnData &rdata, H5::H5File &file, const std::string &hdf5Location)
{

    if(!locationExists(file, hdf5Location))
        createGroup(file, hdf5Location);

    if (rdata.ts)
        createAndWriteDouble1DDataset(file, hdf5Location + "//t",
                                      rdata.ts, rdata.nt);

    if (rdata.xdot)
        createAndWriteDouble1DDataset(file, hdf5Location + "/xdot",
                                      rdata.xdot, rdata.nx);

    if (rdata.llh)
        H5LTset_attribute_double(file.getId(), hdf5Location.c_str(), "llh", rdata.llh, 1);

    if (rdata.status)
        H5LTset_attribute_double(file.getId(), hdf5Location.c_str(), "status", rdata.status, 1);

    if (rdata.sllh)
        createAndWriteDouble1DDataset(file, hdf5Location + "/sllh",
                                      rdata.sllh, rdata.nplist);

    // are double, but should write as int:
    // TODO: change to dataset as soon as are of type int
    if (rdata.numsteps)
        setAttributeIntFromDouble(file, hdf5Location, "numsteps",
                                  rdata.numsteps, rdata.nt);

    if (rdata.numrhsevals)
        setAttributeIntFromDouble(file, hdf5Location, "numrhsevals",
                                  rdata.numrhsevals, rdata.nt);

    if (rdata.numerrtestfails)
        setAttributeIntFromDouble(file, hdf5Location,
                                  "numerrtestfails",
                                  rdata.numerrtestfails, rdata.nt);

    if (rdata.numnonlinsolvconvfails)
        setAttributeIntFromDouble(
                    file, hdf5Location, "numnonlinsolvconvfails",
                    rdata.numnonlinsolvconvfails, rdata.nt);

    if (rdata.order)
        setAttributeIntFromDouble(file, hdf5Location, "order",
                                  rdata.order, rdata.nt);

    if (rdata.numstepsB)
        setAttributeIntFromDouble(file, hdf5Location, "numstepsB",
                                  rdata.numstepsB, rdata.nt);

    if (rdata.numrhsevalsB)
        setAttributeIntFromDouble(file, hdf5Location, "numrhsevalsB",
                                  rdata.numrhsevalsB, rdata.nt);

    if (rdata.numerrtestfailsB)
        setAttributeIntFromDouble(file, hdf5Location,
                                  "numerrtestfailsB",
                                  rdata.numerrtestfailsB, rdata.nt);

    if (rdata.numnonlinsolvconvfailsB)
        setAttributeIntFromDouble(
                    file, hdf5Location, "numnonlinsolvconvfailsB",
                    rdata.numnonlinsolvconvfailsB, rdata.nt);

    auto group = file.openGroup(hdf5Location);

    if (rdata.J)
        createAndWriteDouble2DDataset(file, hdf5Location + "/J", rdata.J,
                                        rdata.nx, rdata.nx);

    if (rdata.x)
        createAndWriteDouble2DDataset(file, hdf5Location + "/x", rdata.x,
                                        rdata.nt, rdata.nx);

    if (rdata.y)
        createAndWriteDouble2DDataset(file, hdf5Location + "/y", rdata.y,
                                        rdata.nt, rdata.ny);

    if (rdata.z)
        createAndWriteDouble2DDataset(file, hdf5Location + "/z", rdata.z,
                                        rdata.nmaxevent, rdata.nz);

    if (rdata.rz)
        createAndWriteDouble2DDataset(file, hdf5Location + "/rz", rdata.rz,
                                        rdata.nmaxevent, rdata.nz);

    if (rdata.sigmay)
        createAndWriteDouble2DDataset(file, hdf5Location + "/sigmay", rdata.sigmay, rdata.nt, rdata.ny);

    if (rdata.sigmaz)
        createAndWriteDouble2DDataset(file, hdf5Location + "/sigmaz", rdata.sigmaz, rdata.nt, rdata.nz);

    if (rdata.s2llh)
        createAndWriteDouble2DDataset(file, hdf5Location + "/s2llh", rdata.s2llh,
                                        rdata.nJ, rdata.nplist);

    if (rdata.sx)
        createAndWriteDouble3DDataset(file, hdf5Location + "/sx", rdata.sx, rdata.nt, rdata.nx, rdata.nplist);

    if (rdata.sy)
        createAndWriteDouble3DDataset(file, hdf5Location + "/sy", rdata.sy, rdata.nt, rdata.ny, rdata.nplist);

    if (rdata.ssigmay)
        createAndWriteDouble3DDataset(file, hdf5Location + "/ssigmay",
                                        rdata.ssigmay, rdata.nt,
                                        rdata.ny, rdata.nplist);

    if (rdata.sz)
        createAndWriteDouble3DDataset(file, hdf5Location + "/sz", rdata.sz,
                                        rdata.nmaxevent, rdata.nz, rdata.nplist);

    if (rdata.srz)
        createAndWriteDouble3DDataset(file, hdf5Location + "/srz", rdata.srz,
                                        rdata.nmaxevent, rdata.nz, rdata.nplist);

    if (rdata.ssigmaz)
        createAndWriteDouble3DDataset(file, hdf5Location + "/ssigmaz", rdata.ssigmaz,
                                      rdata.nmaxevent, rdata.nz, rdata.nplist);
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

std::vector<double> getDoubleArrayAttribute(H5::H5File file,
                                            const std::string &optionsObject,
                                            const std::string &attributeName) {
    H5T_class_t type_class;
    size_t type_size;
    herr_t status = 0;
    hsize_t length = 0; // TODO: check why not set when reading scalar
    status = H5LTget_attribute_info(file.getId(), optionsObject.c_str(), attributeName.c_str(),
                                    &length, &type_class, &type_size);

    if (status < 0) {
        throw AmiException("Error in getDoubleArrayAttribute: Cannot read "
                           "attribute '%s' of '%s'\n",
                           attributeName.c_str(), optionsObject.c_str());
    }

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %lld: ", attributeName, length);
#endif

    std::vector<double> buffer(length);
    if(length) {
        status = H5LTget_attribute_double(file.getId(), optionsObject.c_str(), attributeName.c_str(),
                                          buffer.data());
        if (status < 0)
            throw AmiException("Error in getDoubleArrayAttribute: Cannot read "
                               "attribute '%s' of '%s'\n",
                               attributeName.c_str(), optionsObject.c_str());
    }
#ifdef AMI_HDF5_H_DEBUG
    // printfArray(*destination, length, "%e ");
    printf("\n");
#endif
    return buffer;
}

std::vector<double> getDoubleArrayAttribute2D(H5::H5File const& file,
                                              std::string const& optionsObject,
                                              std::string const& attributeName,
                                              hsize_t &m, hsize_t &n) {

    int rank = 0;
    m = n = 0;

    herr_t status = H5LTget_attribute_ndims(file.getId(), optionsObject.c_str(),
                                            attributeName.c_str(), &rank);
    if(rank == 0)
        return std::vector<double>(); // This should also throw
    if(status < 0 || rank > 2) {
        throw AmiException("Error reading attribute %s of %s. Non-existing or of wrong rank.",
                           attributeName.c_str(), optionsObject.c_str());
    }

    hsize_t dims[2];
    H5T_class_t type_class;
    size_t type_size;
    status = H5LTget_attribute_info(file.getId(), optionsObject.c_str(), attributeName.c_str(),
                                    dims, &type_class, &type_size);
    if (rank == 1) {
        m = dims[0];
        n = 1;
        return getDoubleArrayAttribute(file, optionsObject, attributeName);
    }

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %lld x %lld: ", attributeName, dims[0], dims[1]);
#endif
    m = dims[0];
    n = dims[1];

    std::vector<double> result(m * n);
    H5LTget_attribute_double(file.getId(), optionsObject.c_str(), attributeName.c_str(),
                             result.data());

#ifdef AMI_HDF5_H_DEBUG
    // printfArray(*destination, (*m) * (*n), "%e ");
    printf("\n");
#endif

    return result;
}

std::vector<double> getDoubleArrayAttribute3D(H5::H5File const&file,
                                              std::string const& optionsObject,
                                              std::string const& attributeName,
                                              hsize_t &m, hsize_t &n, hsize_t &o) {

    int rank = 0;
    m = n = o = 0;

    herr_t status = H5LTget_attribute_ndims(file.getId(), optionsObject.c_str(),
                                            attributeName.c_str(), &rank);

    if(status < 0 || rank != 3)
        throw(AmiException("Expected array of rank 3 in %s::%s", optionsObject.c_str(), attributeName.c_str()));
    hsize_t dims[3];
    H5T_class_t type_class;
    size_t type_size;
    status = H5LTget_attribute_info(file.getId(), optionsObject.c_str(), attributeName.c_str(),
                                    dims, &type_class, &type_size);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %lld x %lld x %lld", attributeName, dims[0], dims[1], dims[2]);
#endif
    m = dims[0];
    n = dims[1];
    o = dims[2];

    std::vector<double> result(m * n * o);
    H5LTget_attribute_double(file.getId(), optionsObject.c_str(), attributeName.c_str(),
                             result.data());

#ifdef AMI_HDF5_H_DEBUG
    // printfArray(*destination, (*m) * (*n), "%e ");
    printf("\n");
#endif

    return result;
}

std::vector<int> getIntArrayAttribute(H5::H5File const&file,
                                      const std::string &optionsObject,
                                      const std::string &attributeName) {
    H5T_class_t type_class;
    size_t type_size;
    herr_t status = 0;
    hsize_t length = 0; // TODO: check why not set when reading scalar
    status = H5LTget_attribute_info(file.getId(), optionsObject.c_str(), attributeName.c_str(),
                                    &length, &type_class, &type_size);

    if (status < 0) {
        throw AmiException("Error in getIntArrayAttribute: Cannot read "
                           "attribute '%s' of '%s'\n",
                           attributeName.c_str(), optionsObject.c_str());
    }

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %lld: ", attributeName, length);
#endif

    std::vector<int> buffer(length);
    if(length)
        status = H5LTget_attribute_int(file.getId(), optionsObject.c_str(), attributeName.c_str(),
                                   buffer.data());
    if (status < 0)
        throw AmiException("Error in getIntArrayAttribute: Cannot read "
                           "attribute '%s' of '%s'\n",
                           attributeName.c_str(), optionsObject.c_str());

#ifdef AMI_HDF5_H_DEBUG
    // printfArray(*destination, length, "%e ");
    printf("\n");
#endif
    return buffer;
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

void createAndWriteDouble3DDataset(H5::H5File& file,
                                     std::string const& datasetName,
                                     const double *buffer, hsize_t m,
                                     hsize_t n, hsize_t o) {
    const hsize_t adims[] {m, n, o};
    H5::DataSpace dataspace(3, adims);
    auto dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(buffer, H5::PredType::NATIVE_DOUBLE);
}


void createAndWriteDouble2DAttribute(H5::H5Object& location,
                                     std::string const& attributeName,
                                     const double *buffer, hsize_t m,
                                     hsize_t n) {
    if (attributeExists(location, attributeName)) {
        location.removeAttr(attributeName);
    }

    const hsize_t adims[] {m, n};
    H5::DataSpace space(2, adims);
    auto attr = location.createAttribute(attributeName,  H5::PredType::NATIVE_DOUBLE, space);
    attr.write(H5::PredType::NATIVE_DOUBLE, buffer);
}

void createAndWriteDouble3DAttribute(H5::H5Object& location,
                                     std::string const& attributeName,
                                     const double *buffer, hsize_t m,
                                     hsize_t n, hsize_t o) {
    if (attributeExists(location, attributeName)) {
        location.removeAttr(attributeName);
    }

    const hsize_t adims[] {m, n, o};
    H5::DataSpace space(3, adims);
    auto attr = location.createAttribute(attributeName, H5::PredType::NATIVE_DOUBLE, space);
    attr.write(H5::PredType::NATIVE_DOUBLE, buffer);
}

void setAttributeIntFromDouble(H5::H5File const& file,
                               const std::string &optionsObject,
                               const std::string &attributeName,
                               const double *bufferDouble,
                               hsize_t length) {
    int bufferInt[length];
    for (int i = 0; (unsigned)i < length; ++i)
        bufferInt[i] = bufferDouble[i];

    auto result = H5LTset_attribute_int(file.getId(), optionsObject.c_str(), attributeName.c_str(), bufferInt, length);

    if(result < 0)
        throw(AmiException("Failure writing attribute %s on %s.",
                           attributeName.c_str(), optionsObject.c_str()));
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

    if(attributeExists(file, datasetPath, "pscale")) {
        model.setParameterScale(static_cast<AMICI_parameter_scaling>(getIntScalarAttribute(file, datasetPath, "pscale")));
    }

    if(attributeExists(file, datasetPath, "nmaxevent")) {
        model.setNMaxEvent(getIntScalarAttribute(file, datasetPath, "nmaxevent"));
    }

    if(attributeExists(file, datasetPath, "qpositivex")) {
        // TODO double vs int?!
        auto dblQPosX = getDoubleArrayAttribute(file, datasetPath, "qpositivex");
        if (dblQPosX.size() == (unsigned) model.nx)
            model.setPositivityFlag(std::vector<int>(dblQPosX.begin(), dblQPosX.end()));
        else if(dblQPosX.size() != 0) { // currently not written correctly from matlab
            throw(AmiException("Failed reading qpositivex (%d).", dblQPosX.size()));
        }
    }

    if(attributeExists(file, datasetPath, "theta")) {
        model.setParameters(getDoubleArrayAttribute(file, datasetPath, "theta"));
    }

    if(attributeExists(file, datasetPath, "kappa")) {
        model.setFixedParameters(getDoubleArrayAttribute(file, datasetPath, "kappa"));
    }

    if(attributeExists(file, datasetPath, "ts")) {
        model.setTimepoints(getDoubleArrayAttribute(file, datasetPath, "ts"));
    }

    if(attributeExists(file, datasetPath, "sens_ind")) {
        auto sensInd = getIntArrayAttribute(file, datasetPath, "sens_ind");
        // currently base 1 indices are written
        for (int i = 0; (unsigned)i < sensInd.size(); ++i) {
            sensInd[i] -= 1;
        }
        model.setParameterList(sensInd);
    }

    if(attributeExists(file, datasetPath, "x0")) {
        auto x0 = getDoubleArrayAttribute(file, datasetPath, "x0");
        if(x0.size())
            model.setInitialStates(x0);
    }

    if(attributeExists(file, datasetPath, "sx0")) {
        hsize_t length0 = 0;
        hsize_t length1 = 0;
        auto sx0 = getDoubleArrayAttribute2D(file, datasetPath, "sx0", length0, length1);
        if(sx0.size()) {
            if (length0 != (unsigned) model.nx && length1 != (unsigned) model.nplist())
                throw(AmiException("Dimension mismatch when reading sx0. Expected %dx%d, got %llu, %llu.",
                                   model.nx, model.nplist(), length0, length1));
            model.setInitialStateSensitivities(sx0);
        }
    }

    if(attributeExists(file, datasetPath, "pbar")) {
        auto pbar = getDoubleArrayAttribute(file, datasetPath, "pbar");
        if(pbar.size())
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
    dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE);

    return result;
}

} // namespace hdf5
} // namespace amici
