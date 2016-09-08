#ifndef AMI_HDF5_H
#define AMI_HDF5_H

#include<hdf5.h>
#include<include/symbolic_functions.h>
#include<include/amici.h>

#undef AMI_HDF5_H_DEBUG

void storeSimulation(const char* fileName, ReturnData rdata);
UserData readSimulationUserData(const char* fileName);
double getDoubleScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName);

int getIntScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName);
void getDoubleArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *length);
void getIntArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, int **destination, hsize_t *length);
void getDoubleArrayAttribute2D(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *m, hsize_t *n);

#endif
