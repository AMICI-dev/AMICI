#ifndef AMICI_HDF5_H
#define AMICI_HDF5_H

#include <string>

#include <hdf5.h>
#include <hdf5_hl.h>

#undef AMI_HDF5_H_DEBUG

namespace amici {

class ReturnData;
class ExpData;
class Model;
class Solver;

/* Functions for reading and writing AMICI data to/from HDF5 files. */
// TODO: proper type checking, exception instead of return code, c++ api; check if all fields are read and saved
/**
 * @brief Read solver options from HDF5 file
 * @param fileId hdf5 file handle to read from
 * @param solver solver to set options on
 * @param datasetPath Path inside the HDF5 file
 */
void readSolverSettingsFromHDF5(hid_t fileId, Solver& solver, std::string const& datasetPath);

/**
 * @brief Read solver options from HDF5 file
 * @param hdffile Name of HDF5 file
 * @param solver solver to set options on
 * @param datasetPath Path inside the HDF5 file
 */
void readSolverSettingsFromHDF5(std::string const& hdffile, Solver& solver, std::string const& datasetPath);

/**
 * @brief Read model data from HDF5 file
 * @param hdffile Name of HDF5 file
 * @param model model to set data on
 * @param datasetPath Path inside the HDF5 file
 */
void readModelDataFromHDF5(std::string const& hdffile, Model& model, std::string const& datasetPath);

/**
 * @brief Read model data from HDF5 file
 * @param fileId hdf5 file handle to read from
 * @param model model to set data on
 * @param datasetPath Path inside the HDF5 file
 */
void readModelDataFromHDF5(hid_t fileId, Model& model, std::string const& datasetPath);


/**
  * @brief AMI_HDF5_writeReturnData writes ReturnData struct to attributes of an
 * HDF5 dataset
  * @param rdata Data to write
  * @param hdffile Filename of HDF5 file
  * @param datasetPath Full dataset path inside the HDF5 file (will be created)
  */

void AMI_HDF5_writeReturnData(const ReturnData *rdata,
                                      const char *hdffile,
                                      const char *datasetPath);
/**
 * @brief AMI_HDF5_readSimulationExpData reads AMICI experimental data from
 * attributes in HDF5 file.
 * @param hdffile Name of HDF5 file
 * @param dataObject Path inside the HDF5 file to object having ExpData as
 * attributes
 * @return
 */

ExpData *AMI_HDF5_readSimulationExpData(const char *hdffile,
                                                const char *dataObject,
                                                Model *model);

/**
 * @brief AMI_HDF5_attributeExists Check whether an attribute with the given
 * name exists on the given dataset
 * @param fileId The HDF5 file object
 * @param datasetPath Dataset of which attributes should be checked
 * @param attributeName Name of the attribute of interest
 * @return
 */
int AMI_HDF5_attributeExists(hid_t fileId, const char *datasetPath,
                                     const char *attributeName);

// Helper functions to reading and writing HDF5 attributes:

herr_t AMI_HDF5_createAndWriteDouble2DAttribute(
    hid_t dataset, const char *attributeName, const double *buffer, hsize_t m,
    hsize_t n);

herr_t AMI_HDF5_createAndWriteDouble3DAttribute(
    hid_t dataset, const char *attributeName, const double *buffer, hsize_t m,
    hsize_t n, hsize_t o);

double AMI_HDF5_getDoubleScalarAttribute(hid_t file_id,
                                                 const char *optionsObject,
                                                 const char *attributeName, double *attributeValue);

int AMI_HDF5_getIntScalarAttribute(hid_t file_id,
                                           const char *optionsObject,
                                           const char *attributeName, int *attributeValue);

int AMI_HDF5_getDoubleArrayAttribute(hid_t file_id,
                                             const char *optionsObject,
                                             const char *attributeName,
                                             double **destination,
                                             hsize_t *length);

void AMI_HDF5_getIntArrayAttribute(hid_t file_id,
                                           const char *optionsObject,
                                           const char *attributeName,
                                           int **destination, hsize_t *length);

int AMI_HDF5_getDoubleArrayAttribute2D(hid_t file_id,
                                                const char *optionsObject,
                                                const char *attributeName,
                                                double **destination,
                                                hsize_t *m, hsize_t *n);

int AMI_HDF5_getDoubleArrayAttribute3D(hid_t file_id,
                                               const char *optionsObject,
                                               const char *attributeName,
                                               double **destination, hsize_t *m,
                                               hsize_t *n, hsize_t *o);

void AMI_HDF5_getDoubleArrayAttribute4D(
    hid_t file_id, const char *optionsObject, const char *attributeName,
    double **destination, hsize_t *m, hsize_t *n, hsize_t *o, hsize_t *pp);

void AMI_HDF5_setAttributeIntFromDouble(hid_t file_id,
                                                const char *obj_name,
                                                const char *attr_name,
                                                const double *bufferDouble,
                                                size_t size);

} // namespace amici

#endif
