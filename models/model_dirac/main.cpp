#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "wrapfunctions.h" /* user functions */
#include "include/symbolic_functions.h"
#include <include/amici.h> /* amici functions */
#include <src/ami_hdf5.h>
#include <include/udata_accessors.h>
#include <include/rdata_accessors.h>

void processUserData(UserData *udata);
void processReturnData(ReturnData *rdata, UserData *udata);
void printReturnData(ReturnData *rdata, UserData *udata);

int main(int argc, char **argv)
{  
    const char *hdffile;

    if(argc != 2) {
        fprintf(stderr, "Error: must provide input file as first and only argument.");
        return 1;
    } else {
        hdffile = argv[1];
    }
    
    UserData *udata = readSimulationUserData(hdffile);
    if (udata == NULL) {
        return 1;
    }

    ExpData *edata = readSimulationExpData(hdffile, udata);
    if (edata == NULL) {
        freeUserData(udata);
        return 1;
    }

    int status = 0;
    ReturnData *rdata = getSimulationResults(udata, edata, &status);
    if (rdata == NULL) {
        freeExpData(edata);
        freeUserData(udata);
        return 1;
    }

    processReturnData(rdata, udata);
    writeReturnData(hdffile, rdata, udata);

    freeExpData(edata);
    freeUserData(udata);
    freeReturnData(rdata);

    return 0;
}


void processReturnData(ReturnData *rdata, UserData *udata) {
    printReturnData(rdata, udata);
}

void printReturnData(ReturnData *rdata, UserData *udata) {

    printf("tsdata: ");
    printArray(tsdata, nt);


    printf("\n\nxdata\n");
    for(int i = 0; i < nx; ++i) {
        for(int j = 0; j < nt; ++j)
            printf("%e\t", rdata->am_xdata[j +  nt * i]);
        printf("\n");
    }

    printf("\nydata\n");
    for(int i = 0; i < ny; ++i) {
        for(int j = 0; j < nt; ++j)
            printf("%e\t", rdata->am_ydata[j +  nt * i]);
        printf("\n");
    }

    printf("\n\nxdotdata: ");
    for(int i = 0; i < nx; ++i)
        printf("%e\t", rdata->am_xdotdata[i]);



    printf("\njdata\n");
    for(int i = 0; i < nx; ++i) {
        for(int j = 0; j < nx; ++j)
            printf("%e\t", rdata->am_Jdata[i + nx * j]);
        printf("\n");
    }

    printf("\nnumsteps: \t\t");
    printfArray(numstepsdata, nt, "%.0f ");
    printf("\nnumrhsevalsdata: \t");
    printfArray(numrhsevalsdata, nt, "%.0f ");
    printf("\norder: \t\t");
    printfArray(orderdata, nt, "%.0f ");
    printf("\n");
    printf("llh: %e\n", *llhdata);
}

