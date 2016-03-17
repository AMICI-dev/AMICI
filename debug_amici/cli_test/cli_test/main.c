//
//  main.c
//  cli_test
//
//  Created by Fabian Fröhlich on 27/01/16.
//  Copyright © nt16 Fabian Fröhlich. All rights reserved.
//

#include <stdio.h>
#include <memory.h>
#include "mex.h"
#include "matrix.h"

#define setOption(FIELD,VALUE) \
mx ## FIELD = mxCreateDoubleScalar(VALUE); \
fn = mxGetFieldNumber(mxoptions, #FIELD); \
mxSetFieldByNumber(mxoptions,0,fn,mx  ## FIELD) \


int main(int argc, const char * argv[]) {
    // insert code here...
    
    
    
    int nt = 11;
    int fn = 0;
    mwSize dimscalar[2] = {1,1};
    mxArray* prhs[9];
    mxArray* plhs[1];
    mxArray *mxsol;
    mxArray *mxstatus, *mxllh, *mxchi2, *mxt, *mxnumsteps, *mxnumrhsevals, *mxorder, *mxnumstepsS, *mxnumrhsevalsS, *mxz, *mxx, *mxy, *mxxdot, *mxJ, *mxdydx, *mxdydp, *mxdxdotdp;
    mxArray *mxoptions;
    mxArray *mxatol, *mxrtol,*mxmaxsteps, *mxsens_ind, *mxid, *mxne, *mxtstart, *mxlmm, *mxiter, *mxlinsol, *mxstldet, *mxNd, *mxinterpType, *mxlmmB, *mxiterB, *mxism, *mxsensi_meth, *mxsensi, *mxnmaxevent, *mxz2event, *mxubw, *mxlbw, *mxdata_model, *mxevent_model, *mxordering, *mxss, *mxnp, *mxnx, *mxny, *mxnz, *mxnnz;
    mxArray *mxtout, *mxtheta, *mxkappa, *mxplist, *mxpbar, *mxxscale;
    double *dtout, *dtheta, *dkappa, *dplist, *dpbar;
    mxArray *mxdata, *mxY, *mxSigma_Y, *mxZ, *mxSigma_Z;
    double *dY, *dSigma_Y, *dZ, *dSigma_Z;
    
    
    // sol
    const char *field_names_sol[] = {"status","llh","chi2","t","numsteps","numrhsevals","order","numstepsS","numrhsevalsS","z","x","y","xdot","J","dydp","dydx","dxdotdp"};
    mxsol = mxCreateStructMatrix(1,1,17,field_names_sol);
    
    mxstatus = mxCreateDoubleMatrix(1,1,mxREAL);
    mxSetField(mxsol,0,"status",mxstatus);
    mxllh = mxCreateDoubleMatrix(1,1,mxREAL);
    mxSetField(mxsol,0,"llh",mxllh);
    mxchi2 = mxCreateDoubleMatrix(1,1,mxREAL);
    mxSetField(mxsol,0,"chi2",mxchi2);
    mxt = mxCreateDoubleMatrix(1,nt,mxREAL);
    mxSetField(mxsol,0,"t",mxt);
    mxnumsteps = mxCreateDoubleMatrix(nt,1,mxREAL);
    mxSetField(mxsol,0,"numsteps",mxnumsteps);
    mxnumrhsevals = mxCreateDoubleMatrix(nt,1,mxREAL);
    mxSetField(mxsol,0,"numrhsevals",mxnumrhsevals);
    mxorder = mxCreateDoubleMatrix(nt,1,mxREAL);
    mxSetField(mxsol,0,"order",mxorder);
    mxnumstepsS = mxCreateDoubleMatrix(nt,1,mxREAL);
    mxSetField(mxsol,0,"numstepsS",mxnumstepsS);
    mxnumrhsevalsS = mxCreateDoubleMatrix(nt,1,mxREAL);
    mxSetField(mxsol,0,"numrhsevalsS",mxnumrhsevalsS);
    mxz = mxCreateDoubleMatrix(2,2,mxREAL);
    mxSetField(mxsol,0,"z",mxz);
    mxx = mxCreateDoubleMatrix(nt,3,mxREAL);
    mxSetField(mxsol,0,"x",mxx);
    mxy = mxCreateDoubleMatrix(nt,1,mxREAL);
    mxSetField(mxsol,0,"y",mxy);
    mxxdot = mxCreateDoubleMatrix(nt,3,mxREAL);
    mxSetField(mxsol,0,"xdot",mxxdot);
    mxJ = mxCreateDoubleMatrix(3,3,mxREAL);
    mxSetField(mxsol,0,"J",mxJ);
    mxdydp = mxCreateDoubleMatrix(nt,1,mxREAL);
    mxSetField(mxsol,0,"dydp",mxdydp);
    mxdydx = mxCreateDoubleMatrix(nt,1,mxREAL);
    mxSetField(mxsol,0,"dydx",mxdydx);
    mxdxdotdp = mxCreateDoubleMatrix(nt,1,mxREAL);
    mxSetField(mxsol,0,"dxdotdp",mxdxdotdp);
    
    
    
    prhs[0] = mxsol;
    
    
    // tout
    mxtout = mxCreateDoubleMatrix(0,0,mxREAL);
    double tout[] ={0,1,2,3,4,5,6,7,8,9,10};
    dtout = mxMalloc(nt*sizeof(double));
    memcpy(dtout,tout,nt*sizeof(double));
    mxSetPr(mxtout, dtout);
    mxSetM(mxtout, nt);
    mxSetN(mxtout, 1);
    
    prhs[1] = mxtout;
    
    // theta
    mxtheta = mxCreateDoubleMatrix(0,0,mxREAL);
    double theta[] ={0.5,2,0.5,0.5};
    dtheta = mxMalloc(4*sizeof(double));
    memcpy(dtheta,theta,4*sizeof(double));
    mxSetPr(mxtheta, dtheta);
    mxSetM(mxtheta, 4);
    mxSetN(mxtheta, 1);
    
    prhs[2] = mxtheta;

    // kappa
    mxkappa = mxCreateDoubleMatrix(0,0,mxREAL);
    double kappa[] ={4,8,10,4};
    dkappa = mxMalloc(4*sizeof(double));
    memcpy(dkappa,kappa,4*sizeof(double));
    mxSetPr(mxkappa, dkappa);
    mxSetM(mxkappa, 4);
    mxSetN(mxkappa, 1);
    
    prhs[3] = mxkappa;
    
    // options_ami
    const char *field_names_options[] = {"atol", "rtol","maxsteps", "sens_ind", "id", "ne", "tstart", "lmm", "iter", "linsol", "stldet", "Nd", "interpType", "lmmB", "iterB", "ism", "sensi_meth", "sensi", "nmaxevent", "z2event", "ubw", "lbw", "data_model", "event_model", "ordering", "ss", "np", "nx", "ny", "nz", "nnz"};
    mxoptions = mxCreateStructArray(2,dimscalar,31,field_names_options);
    setOption(atol, 1E-8);
    setOption(rtol, 1E-8);
    setOption(maxsteps, 10000.0);
    double sens_ind[] ={1,2,3,4};
    double* dsens_ind = mxMalloc(4*sizeof(double));
    memcpy(dsens_ind,sens_ind,4*sizeof(double));
    fn = mxGetFieldNumber(mxoptions, "sens_ind");
    mxsens_ind = mxCreateDoubleMatrix(0,0,mxREAL);
    mxSetPr(mxsens_ind, dsens_ind);
    mxSetFieldByNumber(mxoptions,0,fn,mxsens_ind);
    double id[] ={0,0,0};
    double* did = mxMalloc(3*sizeof(double));
    memcpy(did,id,3);
    fn = mxGetFieldNumber(mxoptions, "id");
    mxid = mxCreateDoubleMatrix(0,0,mxREAL);
    mxSetPr(mxid, did);
    mxSetFieldByNumber(mxoptions,0,fn,mxid);
    setOption(ne, 3.0);
    setOption(tstart, 0.0);
    setOption(lmm, 2.0);
    setOption(iter, 2.0);
    setOption(linsol, 9.0);
    setOption(stldet, 1.0);
    setOption(Nd, 1000.0);
    setOption(interpType, 1.0);
    setOption(lmmB, 2.0);
    setOption(iterB, 2.0);
    setOption(ism, 1.0);
    setOption(sensi_meth, 1.0);
    setOption(sensi, 0.0);
    setOption(nmaxevent, 2.0);
    double z2event[] ={0,0,0};
    double* dz2event = mxMalloc(3*sizeof(double));
    memcpy(dz2event,z2event,3*sizeof(double));
    fn = mxGetFieldNumber(mxoptions, "z2event");
    mxz2event = mxCreateDoubleMatrix(0,0,mxREAL);
    mxSetPr(mxz2event, dz2event);
    mxSetFieldByNumber(mxoptions,0,fn,mxz2event);
    setOption(ubw, 0.0);
    setOption(lbw, 1.0);
    setOption(data_model, 1.0);
    setOption(event_model, 1.0);
    setOption(ordering, 1.0);
    setOption(ss, 0.0);
    setOption(np, 4.0);
    setOption(nx, 3.0);
    setOption(ny, 1.0);
    setOption(nz, 2.0);
    setOption(nnz, 4.0);
    
    prhs[4] = mxoptions;
    
    // plist
    mxplist = mxCreateDoubleMatrix(0,0,mxREAL);
    double plist[] ={1,2,3,4};
    dplist = mxMalloc(4*sizeof(double));
    memcpy(dplist,plist,4*sizeof(double));
    mxSetPr(mxplist, dplist);
    mxSetM(mxplist, 4);
    mxSetN(mxplist, 1);
    
    prhs[5] = mxplist;
    
    // pbar
    mxpbar = mxCreateDoubleMatrix(0,0,mxREAL);
    double pbar[] ={1,1,1,1};
    dpbar = mxMalloc(4*sizeof(double));
    memcpy(dpbar,pbar,4*sizeof(double));
    mxSetPr(mxpbar, dpbar);
    mxSetM(mxpbar, 4);
    mxSetN(mxpbar, 1);
    
    prhs[6] = mxpbar;
    
    // xscale
    mxxscale = mxCreateDoubleMatrix(0,0,mxREAL);
    
    prhs[7] = mxxscale;
    
    // data
    
    const char *field_names_data[] = {"Y", "Sigma_Y","Z", "Sigma_Z"};
    mxdata = mxCreateStructMatrix(1,1,4,field_names_data);
    
    prhs[8] = mxdata;
    
    mxY = mxCreateDoubleMatrix(0,0,mxREAL);
    double Y[] ={1,1,1,1,1,1,1,1,1,1,1,1};
    dY = mxMalloc(nt*sizeof(double));
    memcpy(dY,Y,nt*sizeof(double));
    mxSetPr(mxY, dY);
    mxSetM(mxY, nt);
    mxSetN(mxY, 1);
    mxSetField(mxdata,0,"Y",mxY);
    
    mxSigma_Y = mxCreateDoubleMatrix(0,0,mxREAL);
    double Sigma_Y[] ={1,1,1,1,1,1,1,1,1,1,1,1};
    dSigma_Y = mxMalloc(nt*sizeof(double));
    memcpy(dSigma_Y,Sigma_Y,nt*sizeof(double));
    mxSetPr(mxSigma_Y, dSigma_Y);
    mxSetM(mxSigma_Y, nt);
    mxSetN(mxSigma_Y, 1);
    mxSetField(mxdata,0,"Sigma_Y",mxSigma_Y);
    
    mxZ = mxCreateDoubleMatrix(0,0,mxREAL);
    double Z[] ={1,1,1,1};
    dZ = mxMalloc(4*sizeof(double));
    memcpy(dZ,Z,4*sizeof(double));
    mxSetPr(mxZ, dZ);
    mxSetM(mxZ, 2);
    mxSetN(mxZ, 2);
    mxSetField(mxdata,0,"Z",mxZ);
    
    mxSigma_Z = mxCreateDoubleMatrix(0,0,mxREAL);
    double Sigma_Z[] ={1,1,1,1};
    dSigma_Z = mxMalloc(nt*sizeof(double));
    memcpy(dSigma_Z,Sigma_Z,nt*sizeof(double));
    mxSetPr(mxSigma_Z, dSigma_Z);
    mxSetM(mxSigma_Z, 2);
    mxSetN(mxSigma_Z, 2);
    mxSetField(mxdata,0,"Sigma_Z",mxSigma_Z);
    
    
    mexFunction(1,plhs, 8, prhs);
    
    
    return 0;
}
