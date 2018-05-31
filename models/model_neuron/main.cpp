#include <cassert>
#include <cmath>
#include <memory>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <amici/amici.h>    /* AMICI base functions */
#include <amici/hdf5.h>     /* AMICI HDF5 I/O functions */
#include "wrapfunctions.h"  /* model-provided functions */

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
void processReturnData(amici::ReturnData *rdata, amici::Model *model);
void printReturnData(amici::ReturnData *rdata, amici::Model *model);

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

    auto model = getModel();
    auto solver = model->getSolver();

    // Read AMICI settings and model parameters from HDF5 file
    amici::hdf5::readModelDataFromHDF5(hdffile, *model, "/options");
    amici::hdf5::readSolverSettingsFromHDF5(hdffile, *solver, "/options");

    // Read ExpData (experimental data for model) from HDF5 file
    auto edata = amici::hdf5::readSimulationExpData(hdffile, "/data", *model);

    // Run the simulation
    auto rdata = runAmiciSimulation(*solver, edata.get(), *model);

    // Do something with the simulation results
    processReturnData(rdata.get(), model.get());

    // Save simulation results to HDF5 file
    amici::hdf5::writeReturnData(*rdata, hdffile, "/solution");

    return 0;
}

void processReturnData(amici::ReturnData *rdata, amici::Model *model) {
    // show some the simulation results
    printReturnData(rdata, model);
}

void printReturnData(amici::ReturnData *rdata, amici::Model *model) {
    // Print of some the simulation results

    printf("\n\nStates (xdata):\n");
    for (int i = 0; i < model->nx; ++i) {
        for (int j = 0; j < model->nt(); ++j)
            printf("%e\t", rdata->x[j + model->nt() * i]);
        printf("\n");
    }

    printf("\nObservables (ydata):\n");
    for (int i = 0; i < model->ny; ++i) {
        for (int j = 0; j < model->nt(); ++j)
            printf("%e\t", rdata->y[j + model->nt() * i]);
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

    printf("\n");
    printf("Loglikelihood (llh): %e\n", rdata->llh);
}
