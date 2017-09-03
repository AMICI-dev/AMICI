#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "wrapfunctions.h"               /* model-provided functions */
#include <include/amici_hdf5.h>          /* AMICI HDF5 I/O functions */
#include <include/amici_interface_cpp.h> /* AMICI API */
#include <include/amici_model.h>

/* This is a scaffold for a stand-alone AMICI simulation executable
 * demonstrating
 * use of the AMICI C++ API.
 *
 * This program reads AMICI options from an HDF5 file, prints some results
 * and writes additional results to an HDF5 file. The name of the HDF5 file
 * is expected as single command line argument.
 *
 * An initial HDF5 file with the required fields can be generated using MATLAB
 * by adding the following lines
 * at the end of simulate_${MODEL_NAME}.m file just before the final "end":
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
 * Default UserData settings can be written to an HDF5 file with:
 *     structToHDF5Attribute('test.h5', '/options', amioption())
 */

// Function prototypes
void processReturnData(ReturnData *rdata, UserData *udata, Model *model);
void printReturnData(ReturnData *rdata, UserData *udata, Model *model);

int main(int argc, char **argv) {
    // HDF5 file to read and write data (full path)
    const char *hdffile;

    // Check command line arguments
    if (argc != 2) {
        fprintf(stderr, "Error: must provide HDF5 input file as first and only "
                        "argument.\n");
        return 1;
    } else {
        hdffile = argv[1];
    }

    Model *model = getModel();

    // Read UserData (AMICI settings and model parameters) from HDF5 file
    UserData *udata =
        AMI_HDF5_readSimulationUserDataFromFileName(hdffile, "/options", model);
    if (udata == NULL) {
        return 1;
    }

    // Read ExpData (experimental data for model) from HDF5 file
    ExpData *edata =
        AMI_HDF5_readSimulationExpData(hdffile, udata, "/data", model);

    // Run the simulation
    ReturnData *rdata = getSimulationResults(model, udata, edata);
    if (rdata == NULL) {
        if (edata)
            delete edata;
        if (udata)
            delete udata;
        return 1;
    }

    // Do something with the simulation results
    processReturnData(rdata, udata, model);

    // Save simulation results to HDF5 file
    AMI_HDF5_writeReturnData(rdata, udata, hdffile, "/solution");

    // Free memory
    delete model;
    if (edata)
        delete edata;
    if (udata)
        delete udata;
    if (rdata)
        delete rdata;

    return 0;
}

void processReturnData(ReturnData *rdata, UserData *udata, Model *model) {
    // show some the simulation results
    printReturnData(rdata, udata, model);
}

void printReturnData(ReturnData *rdata, UserData *udata, Model *model) {
    // Print of some the simulation results

    printf("Timepoints (tsdata): ");
    printArray(rdata->ts, udata->nt);

    printf("\n\nStates (xdata):\n");
    for (int i = 0; i < model->nx; ++i) {
        for (int j = 0; j < udata->nt; ++j)
            printf("%e\t", rdata->x[j + udata->nt * i]);
        printf("\n");
    }

    printf("\nObservables (ydata):\n");
    for (int i = 0; i < model->ny; ++i) {
        for (int j = 0; j < udata->nt; ++j)
            printf("%e\t", rdata->y[j + udata->nt * i]);
        printf("\n");
    }

    printf("\n\ndx/dt (xdotdata):\n");
    for (int i = 0; i < model->nx; ++i)
        printf("%e\t", rdata->xdot[i]);

    //    printf("\nJacobian (jdata)\n");
    //    for(int i = 0; i < nx; ++i) {
    //        for(int j = 0; j < nx; ++j)
    //            printf("%e\t", rdata->J[i + nx * j]);
    //        printf("\n");
    //    }

    printf("\nnumsteps: \t\t");
    printfArray(rdata->numsteps, udata->nt, "%.0f ");

    printf("\nnumrhsevalsdata: \t");
    printfArray(rdata->numrhsevals, udata->nt, "%.0f ");

    printf("\norder: \t\t");
    printfArray(rdata->order, udata->nt, "%.0f ");

    printf("\n");
    printf("Loglikelihood (llh): %e\n", *rdata->llh);
}
