#include<hdf5_hl.h>
#include<assert.h>

#include "ami_hdf5.h"


double getDoubleScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName) {
    double doubleScalar;
    H5LTget_attribute_double(file_id, optionsObject, attributeName, &doubleScalar);
#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %e\n", attributeName, doubleScalar);
#endif
    return doubleScalar;
}

int getIntScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName) {
    int intScalar;
    H5LTget_attribute_int(file_id, optionsObject, attributeName, &intScalar);
#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d\n", attributeName, intScalar);
#endif
    return intScalar;
}

void getDoubleArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *length) {
    H5T_class_t type_class;
    size_t type_size;
    H5LTget_attribute_info(file_id, optionsObject, attributeName, length, &type_class, &type_size);
#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d: ", attributeName, *length);
#endif
    *destination = malloc((*length) * sizeof(double)); // vs. type_size
    H5LTget_attribute_double(file_id, optionsObject, attributeName, *destination);
#ifdef AMI_HDF5_H_DEBUG
    printfArray(*destination, *length, "%e ");
    printf("\n");
#endif
}

void getDoubleArrayAttribute2D(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *m, hsize_t *n) {
    int rank;
    H5LTget_attribute_ndims(file_id, optionsObject, attributeName, &rank);
    assert(rank <= 2);

    hsize_t dims[2];
    H5T_class_t type_class;
    size_t type_size;
    H5LTget_attribute_info(file_id, optionsObject, attributeName, dims, &type_class, &type_size);

    if(rank == 1) {
        *m = 1;
        *n = dims[0];
        getDoubleArrayAttribute(file_id, optionsObject, attributeName, destination, m);
    } else {
#ifdef AMI_HDF5_H_DEBUG
        printf("%s: %d x %d: ", attributeName, dims[0], dims[1]);
#endif
        *m = dims[0];
        *n = dims[1];

        *destination = malloc(type_size * (*m) * (*n));
        H5LTget_attribute_double(file_id, optionsObject, attributeName, *destination);
#ifdef AMI_HDF5_H_DEBUG
        printfArray(*destination, (*m) * (*n), "%e ");
        printf("\n");
#endif
    }
}

void getIntArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, int **destination, hsize_t *length) {
    H5T_class_t type_class;
    size_t type_size;

    H5LTget_attribute_info(file_id, optionsObject, attributeName, length, &type_class, &type_size);
#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d: ", attributeName, *length);
#endif
    *destination = malloc(sizeof(int) * (*length));
    H5LTget_attribute_int(file_id, optionsObject, attributeName, *destination);
#ifdef AMI_HDF5_H_DEBUG
    printfArray(*destination, *length, "%d ");
    printf("\n");
#endif
}

void createAndWriteDouble2DAttribute(hid_t dataset, const char *attributeName, const double *buffer, hsize_t m, hsize_t n) {
    const hsize_t adims[] = {m, n};
    hid_t space = H5Screate_simple(2, adims, NULL);
    hid_t attr = H5Acreate2(dataset, attributeName, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Awrite(attr, H5T_NATIVE_DOUBLE, buffer);
}
