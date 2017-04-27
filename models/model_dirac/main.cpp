#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "wrapfunctions.h" /* model specific functions */
#include <include/amici.h> /* AMICI API */
#include <src/ami_hdf5.h>  /* AMICI HDF5 I/O functions */
#include <include/udata_accessors.h> /* Accessor macros for UserData members */
#include <include/rdata_accessors.h> /* Accessor macros for ReturnData members */

/* This is a scaffold for a stand-alone AMICI simulation executable demonstrating
 * use of the AMICI C++ API.
 *
 * This program reads AMICI options from an HDF5 file, prints some results
 * and writes additional results to an HDF5 file. The name of the HDF5 file
 * is expected as single command line argument.
 *
 * An initial HDF5 file with the required fields can be generated using MATLAB by adding the following lines
 * in the simulate_${MODEL_NAME}.m file just before "sol = ami_${MODEL_NAME}" close to the end of the file:
 *
 *    %% Write data that is passed to AMICI to HDF5
 *    hdffile = fullfile(pwd, 'mydata.h5');
 *    structToHDF5Attribute(hdffile, '/options', options_ami);
 *    h5writeatt(hdffile, '/options', 'ts', tout);
 *    h5writeatt(hdffile, '/options', 'nt', numel(tout));
 *    h5writeatt(hdffile, '/options', 'theta', theta);
 *    h5writeatt(hdffile, '/options', 'kappa', kappa);
 *    if(~isempty(data))
 *      structToHDF5Attribute(hdffile, '/data', data);
 *    end
 *
 * ... and then running a simulation from MATLAB as usual.
 *
 */

// Function prototypes
void processReturnData(ReturnData *rdata, UserData *udata);
void printReturnData(ReturnData *rdata, UserData *udata);


int main(int argc, char **argv)
{  
    // HDF5 file to read and write data (full path)
    const char *hdffile;

    // Check command line arguments
    if(argc != 2) {
        fprintf(stderr, "Error: must provide HDF5 input file as first and only argument.\n");
        return 1;
    } else {
        hdffile = argv[1];
    }
    
    // Read UserData (AMICI settings and model parameters) from HDF5 file
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(hdffile, "/options");
    if (udata == NULL) {
        return 1;
    }

    // Read ExpData (experimental data for model) from HDF5 file
    ExpData *edata = AMI_HDF5_readSimulationExpData(hdffile, udata);
    if (edata == NULL) {
        freeUserData(udata);
        return 1;
    }

    // Run the simulation
    int status = 0;
    ReturnData *rdata = getSimulationResults(udata, edata, &status);
    if (rdata == NULL) {
        freeExpData(edata);
        freeUserData(udata);
        return 1;
    }

    // Do something with the simulation results
    processReturnData(rdata, udata);

    // Save simulation results to HDF5 file
    AMI_HDF5_writeReturnData(rdata, udata, hdffile, "/solution");

    // Free memory
    freeExpData(edata);
    freeUserData(udata);
    freeReturnData(rdata);

    return 0;
}


void processReturnData(ReturnData *rdata, UserData *udata) {
    // show some the simulation results
    printReturnData(rdata, udata);
}

void printReturnData(ReturnData *rdata, UserData *udata) {
    //Print of some the simulation results

    printf("Timepoints (tsdata): ");
    printArray(tsdata, nt);

    printf("\n\nStates (xdata):\n");
    for(int i = 0; i < nx; ++i) {
        for(int j = 0; j < nt; ++j)
            printf("%e\t", rdata->am_xdata[j +  nt * i]);
        printf("\n");
    }

    printf("\nObservables (ydata):\n");
    for(int i = 0; i < ny; ++i) {
        for(int j = 0; j < nt; ++j)
            printf("%e\t", rdata->am_ydata[j +  nt * i]);
        printf("\n");
    }

    printf("\n\ndx/dt (xdotdata):\n");
    for(int i = 0; i < nx; ++i)
        printf("%e\t", rdata->am_xdotdata[i]);

//    printf("\nJacobian (jdata)\n");
//    for(int i = 0; i < nx; ++i) {
//        for(int j = 0; j < nx; ++j)
//            printf("%e\t", rdata->am_Jdata[i + nx * j]);
//        printf("\n");
//    }

    printf("\nnumsteps: \t\t");
    printfArray(numstepsdata, nt, "%.0f ");

    printf("\nnumrhsevalsdata: \t");
    printfArray(numrhsevalsdata, nt, "%.0f ");

    printf("\norder: \t\t");
    printfArray(orderdata, nt, "%.0f ");

    printf("\n");
    printf("Loglikelihood (llhdata): %e\n", *llhdata);
}

