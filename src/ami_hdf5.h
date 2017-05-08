#ifndef AMI_HDF5_H
#define AMI_HDF5_H

#include <hdf5.h>
#include <hdf5_hl.h>
#include <include/amici.h>

#undef AMI_HDF5_H_DEBUG

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

/* Functions for reading and writing AMICI data to/from HDF5 files. */

/**
  * @brief AMI_HDF5_writeReturnData writes ReturnData struct to attributes of an HDF5 dataset
  * @param rdata Data to write
  * @param udata UserData for model from which rdata was obtained.
  * @param hdffile Filename of HDF5 file
  * @param datasetPath Full dataset path inside the HDF5 file (will be created)
  */

EXTERNC void AMI_HDF5_writeReturnData(const ReturnData *rdata, const UserData *udata, const char* hdffile, const char* datasetPath);

/**
 * @brief AMI_HDF5_readSimulationUserDataFromFileName reads AMICI simulation options from HDF5 file.
 * @param fileName Name of HDF5 file
 * @param datasetPath Path inside the HDF5 file
 * @return
 */
EXTERNC UserData *AMI_HDF5_readSimulationUserDataFromFileName(const char* fileName, const char *datasetPath);

EXTERNC UserData *AMI_HDF5_readSimulationUserDataFromFileObject(hid_t fileId, const char *datasetPath);

EXTERNC ExpData *AMI_HDF5_readSimulationExpData(const char* hdffile, UserData *udata);

// Helper functions to reading and writing HDF5 attributes:

EXTERNC herr_t AMI_HDF5_createAndWriteDouble2DAttribute(hid_t dataset, const char *attributeName, const double *buffer, hsize_t m, hsize_t n);

EXTERNC herr_t AMI_HDF5_createAndWriteDouble3DAttribute(hid_t dataset, const char *attributeName, const double *buffer, hsize_t m, hsize_t n, hsize_t o);

EXTERNC double AMI_HDF5_getDoubleScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName);

EXTERNC int AMI_HDF5_getIntScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName);

EXTERNC int AMI_HDF5_getDoubleArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *length);

EXTERNC void AMI_HDF5_getIntArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, int **destination, hsize_t *length);

EXTERNC void AMI_HDF5_getDoubleArrayAttribute2D(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *m, hsize_t *n);

EXTERNC void AMI_HDF5_getDoubleArrayAttribute3D(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *m, hsize_t *n, hsize_t *o);

EXTERNC void AMI_HDF5_getDoubleArrayAttribute4D(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *m, hsize_t *n, hsize_t *o, hsize_t *pp);

EXTERNC void AMI_HDF5_setAttributeIntFromDouble(hid_t file_id, const char *obj_name, const char *attr_name, const double *bufferDouble, size_t size );

EXTERNC int AMI_HDF5_attributeExists(hid_t fileId, const char *datasetPath, const char *attributeName);
#endif
