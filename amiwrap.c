/**
 * @file   amiwrap.c
 * @brief  core routines for mex interface
 *
 * This file defines the fuction mexFunction which is executed upon calling the mex file from matlab
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES /* MS definition of PI and other constants */
#include <math.h>
#ifndef M_PI /* define PI if we still have no definition */
#define M_PI 3.14159265358979323846
#endif
#include <mex.h>
#include "wrapfunctions.h" /* user functions */
#include <include/amici.h> /* amici functions */

/*!
 * mexFunction is the main function of the mex simulation file this function carries out all numerical integration and writes results into the sol struct.
 *
 * @param[in] nlhs number of output arguments of the matlab call @type int
 * @param[out] plhs pointer to the array of output arguments @type mxArray
 * @param[in] nrhs number of input arguments of the matlab call @type int
 * @param[in] prhs pointer to the array of input arguments @type mxArray
 * @return void
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    int ip, ix, jx, iy, it, ir; /* integers for indexing in loops */
    int jp; /* integers for indexing in loops */
    int ncheck; /* the number of (internal) checkpoints stored so far */
    int nroots; /* number of found roots */
    int idisc; /* number of found roots */
    
    void *ami_mem; /* pointer to cvodes memory block */
    UserData udata; /* user data */
    ReturnData rdata; /* return data */
    ExpData edata; /* experimental data */
    TempData tdata; /* temporary data */
    int status; /* general status flag */
    double *pstatus; /* return status flag */
    int cv_status; /* status flag returned by integration method */

    long int numsteps;
    long int numrhsevals;
    long int numlinsolvsetups;
    long int numerrtestfails;
    int order;
    long int numnonlinsolviters;
    long int numjacevals;
    long int numliniters;
    long int numconvfails;
    long int numprecevals;
    long int numprecsolves;
    long int numjtimesevals;
    
    realtype tlastroot; /* storage for last found root */
    
    bool rootflag, discflag;

    udata = setupUserData(prhs);
    if (udata == NULL) goto freturn;
    
    /* solution struct */
    
    if (!prhs[0]) {
        mexErrMsgTxt("No solution struct provided!");
    }
    
    /* options */
    if (!prhs[4]) {
        mexErrMsgTxt("No options provided!");
    }
    
    tdata = (TempData) mxMalloc(sizeof *tdata);
    if (tdata == NULL) goto freturn;
    
    ami_mem = setupAMI(&status, udata, tdata);
    if (ami_mem == NULL) goto freturn;
    
    rdata = setupReturnData(prhs, udata);
    if (rdata == NULL) goto freturn;
    
    edata = setupExpData(prhs, udata);
    if (edata == NULL) goto freturn;
    
    cv_status = 0;
    
    /* FORWARD PROBLEM */
    
    ncheck = 0; /* the number of (internal) checkpoints stored so far */
    
    t = tstart;
    
    nroots = 0;
    idisc = 0;
    tlastroot = 0;
    for (it=0; it < nt; it++) {
        /* only integrate if no errors occured */
        if(cv_status == 0.0) {
            if(ts[it] > tstart) {
                if(nx>0) {
                    while (t<ts[it]) {
                        if (nr+ndisc>0) {
                            /* we have to find roots */
                            if(sensi_meth == AMI_ASA && sensi >= 1) {
                                cv_status = AMISolveF(ami_mem, RCONST(ts[it]), x, dx, &t, AMI_NORMAL, &ncheck);
                                if (cv_status==AMI_ROOT_RETURN) {
                                    if (t == tlastroot) {
                                        /* we are stuck in a root => turn off rootfinding */
                                        /* at some point we should find a more intelligent solution here, and turn on rootfinding again after some time */
                                        AMIRootInit(ami_mem, 0, NULL);
                                        cv_status = 0;
                                    }
                                    tlastroot = t;
                                    getRootDataASA(&status, &nroots, &idisc, ami_mem, udata, rdata, edata, tdata);
                                    if (t==ts[it]) {
                                        cv_status = 0;
                                    }
                                }
                            } else {
                                cv_status = AMISolve(ami_mem, RCONST(ts[it]), x, dx, &t, AMI_NORMAL);
                                if (cv_status==AMI_ROOT_RETURN) {
                                    if (t == tlastroot) {
                                        /* we are stuck in a root => turn off rootfinding */
                                        /* at some point we should find a more intelligent solution here, and turn on rootfinding again after some time */
                                        AMIRootInit(ami_mem, 0, NULL);
                                        cv_status = 0;
                                    }
                                    tlastroot = t;
                                    getRootDataFSA(&status, &nroots, ami_mem, udata, rdata, tdata);
                                    if (t==ts[it]) {
                                        cv_status = 0;
                                    }
                                }
                                if (cv_status == -22) {
                                    /* clustering of roots => turn off rootfinding */
                                    AMIRootInit(ami_mem, 0, NULL);
                                    cv_status = 0;
                                }
                            }
                        } else {
                            /* solve for forward problem and store checkpoints */
                            if(sensi_meth == AMI_ASA && sensi >= 1) {
                                cv_status = AMISolveF(ami_mem, RCONST(ts[it]), x, dx, &t, AMI_NORMAL, &ncheck);
                            } else {
                                cv_status = AMISolve(ami_mem, RCONST(ts[it]), x, dx, &t, AMI_NORMAL);
                            }
                        }
                        if (cv_status < 0) {
                            status = cv_status;
                            if (status != AMI_SUCCESS) goto freturn;
                        }
                    }
                }
            }
            
            tsdata[it] = ts[it];
            x_tmp = NV_DATA_S(x);
            for (ix=0; ix<nx; ix++) {
                xdata[it+nt*ix] = x_tmp[ix];
            }
            
            if (it == nt-1) {
                if( sensi_meth == AMI_SS) {
                    status = fxdot(t,x,dx,xdot,udata);
                    if (status != AMI_SUCCESS) goto freturn;
                    
                    xdot_tmp = NV_DATA_S(xdot);
                    
                    status = fJ(nx,ts[it],0,x,dx,xdot,Jtmp,udata,NULL,NULL,NULL);
                    if (status != AMI_SUCCESS) goto freturn;
                    
                    memcpy(xdotdata,xdot_tmp,nx*sizeof(realtype));
                    memcpy(Jdata,Jtmp->data,nx*nx*sizeof(realtype));
                    
                    status = fdxdotdp(t,x,dxdotdpdata,udata);
                    if (status != AMI_SUCCESS) goto freturn;
                    status = fdydp(ts[it],it,dydpdata,ydata,xdata,udata);
                    if (status != AMI_SUCCESS) goto freturn;
                    status = fdydx(ts[it],it,dydxdata,ydata,xdata,udata);
                    if (status != AMI_SUCCESS) goto freturn;
                }
            }
            
            if(ts[it] > tstart) {
                getDiagnosis(&status, it, ami_mem, udata, rdata);
            }
            
        } else {
            for(ix=0; ix < nx; ix++) xdata[ix*nt+it] = mxGetNaN();
        }
        
        if(cv_status == 0) {
            status = fy(ts[it],it,ydata,xdata,udata);
            if (status != AMI_SUCCESS) goto freturn;
            
            for (iy=0; iy<ny; iy++) {
                
                if(data_model != AMI_ONEOUTPUT) {
                    if (mxIsNaN(ysigma[iy*nt+it])) {
                        status =fsigma_y(t,sigma_y,udata);
                        if (status != AMI_SUCCESS) goto freturn;
                        
                    } else {
                        sigma_y[iy] = ysigma[iy*nt+it];
                    }
                }
                
                if (data_model == AMI_NORMAL) {
                    if(!mxIsNaN(my[iy*nt+it])){
                        g += 0.5*log(2*pi*pow(sigma_y[iy],2)) + 0.5*pow( ( ydata[iy*nt+it] - my[iy*nt+it] )/sigma_y[iy] , 2);
                        *chi2data += pow( ( ydata[iy*nt+it] - my[iy*nt+it] )/sigma_y[iy] , 2);
                    }
                }
                if (data_model == AMI_LOGNORMAL) {
                    if(!mxIsNaN(my[iy*nt+it])){
                        g += 0.5*log(2*pi*pow(sigma_y[iy]*ydata[iy*nt+it],2)) + 0.5*pow( ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) )/sigma_y[iy] , 2);
                        *chi2data += pow( ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) )/sigma_y[iy] , 2);
                    }
                }
                if (data_model == AMI_ONEOUTPUT) {
                    g += ydata[iy*nt+it];
                }
            }
            if (sensi >= 1) {
                if (sensi_meth == AMI_FSA) {
                    for(ip=0; ip < np; ip++) {
                        if(nx>0) {
                            if(ts[it] > tstart) {
                                status = AMIGetSens(ami_mem, &t, sx);
                                if (status != AMI_SUCCESS) goto freturn;
                            }
                            
                            sx_tmp = NV_DATA_S(sx[ip]);
                            for(ix=0; ix < nx; ix++) {
                                xSdata[(ip*nx + ix)*nt + it] = sx_tmp[ix];
                            }
                        }
                    }
                    fsy(ts[it],it,ySdata,xdata,xSdata,udata);
                }
                if (sensi_meth == AMI_ASA) {
                    status = fdydx(ts[it],it,dydx,ydata,xdata,udata);
                    if (status != AMI_SUCCESS) goto freturn;
                    status = fdydp(ts[it],it,dydp,ydata,xdata,udata);
                    if (status != AMI_SUCCESS) goto freturn;
                    for (iy=0; iy<ny; iy++) {
                        if(data_model != AMI_ONEOUTPUT) {
                            if (mxIsNaN(ysigma[iy*nt+it])) {
                                status = fsigma_y(t,sigma_y,udata);
                                if (status != AMI_SUCCESS) goto freturn;
                                status = fdsigma_ydp(t,dsigma_ydp,udata);
                                if (status != AMI_SUCCESS) goto freturn;
                            } else {
                                for (ip=0; ip<np; ip++) {
                                    dsigma_ydp[ip*ny+iy] = 0;
                                }
                                sigma_y[iy] = ysigma[iy*nt+it];
                            }
                        }
                        for (ip=0; ip<np; ip++) {
                            if(data_model == AMI_NORMAL) {
                                if(!mxIsNaN(my[iy*nt+it])){
                                    dgdp[ip] += dsigma_ydp[ip*ny+iy]/sigma_y[iy] + ( dydp[ip*ny+iy]* ( ydata[iy*nt+it] - my[iy*nt+it] ) )/pow( sigma_y[iy] , 2) - dsigma_ydp[ip*ny+iy]*pow( ( ydata[iy*nt+it] - my[iy*nt+it] ),2)/pow( sigma_y[iy] , 3);
                                }
                            }
                            if(data_model == AMI_LOGNORMAL) {
                                if(!mxIsNaN(my[iy*nt+it])){
                                    dgdp[ip] += (sigma_y[iy]*dydp[ip*ny+iy] + ydata[iy*nt+it]*dsigma_ydp[ip*ny+iy])/(sigma_y[iy]*ydata[iy*nt+it]) + ( dydp[ip*ny+iy]/ydata[iy*nt+it] * ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) ) )/pow( sigma_y[iy] , 2) -  dsigma_ydp[ip*ny+iy]*pow( ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) ),2)/pow(sigma_y[iy] , 3);
                                }
                            }
                            if(data_model == AMI_ONEOUTPUT) {
                                dgdp[ip] += dydp[ip*ny+iy];
                            }
                        }
                        for (ix=0; ix<nx; ix++) {
                            if(data_model == AMI_NORMAL) {
                                if(!mxIsNaN(my[iy*nt+it])){
                                    dgdx[it+ix*nt] += ( dydx[ix*ny+iy] * ( ydata[iy*nt+it] - my[iy*nt+it] ) )/pow( sigma_y[iy] , 2);
                                }
                            }
                            if(data_model == AMI_LOGNORMAL) {
                                if(!mxIsNaN(my[iy*nt+it])){
                                    dgdx[it+ix*nt] += 1/(2*pi)*dydx[ix*ny+iy]/ydata[iy*nt+it] + ( dydx[ix*ny+iy]/ydata[iy*nt+it] * ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) ) )/pow( sigma_y[iy] , 2);
                                }
                            }
                            if(data_model == AMI_ONEOUTPUT) {
                                dgdx[it+ix*nt] += dydx[ix*ny+iy];
                            }
                        }
                    }
                }
            }
        }
    }
    
    /* add event if we did not have one yet */
    if (nr>0) {
        if (nroots==0) {
            for (ir=0; ir<nr; ir++){
                rootdata[nroots + nmaxroot*ir] = t;
                status = froot(t, x, dx, rootvaltmp, udata);
                if (status != AMI_SUCCESS) goto freturn;
                rootvaldata[nroots + nmaxroot*ir] = rootvaltmp[ir];
                /* extract sensitivity information */
                rootidx[nroots] = ir;
                if(sensi >= 1) {
                    if(sensi_meth == AMI_FSA) {
                        status = AMIGetSens(ami_mem, &t, sx);
                        if (status != AMI_SUCCESS) goto freturn;
                        status = fsrootval(t,nroots,rootvalSdata,x,sx,udata);
                        if (status != AMI_SUCCESS) goto freturn;
                        if (sensi >= 2) {
                            status = fs2rootval(t,nroots,rootvalS2data,x,sx,udata);
                            if (status != AMI_SUCCESS) goto freturn;
                        }
                        for (ip=0; ip<np; ip++) {
                            rootSdata[nroots + nmaxroot*(ip*nr + ir)] = 0;
                            if (sensi >= 2) {
                                for (jp=0; jp<np; jp++) {
                                    rootS2data[nroots + nmaxroot*((np*ip+jp)*nr + ir)] = 0;
                                }
                            }
                        }
                    }
                }
                if(!mxIsNaN(mt[ir*nmaxroot+nroots])) {
                    if (event_model == AMI_NORMAL) {
                        r += 0.5*log(2*pi*pow(tsigma[ir*nmaxroot+nroots],2)) + 0.5*pow( ( t - mt[ir*nmaxroot+nroots] )/tsigma[ir*nmaxroot+nroots] , 2);
                        *chi2data += pow( ( t - mt[ir*nmaxroot+nroots] )/tsigma[ir*nmaxroot+nroots] , 2);
                    }
                    
                    if (event_model == AMI_NORMAL) {
                        r += 0.5*log(2*pi*pow(tsigma[ir*nmaxroot+nroots],2)) + 0.5*pow( ( rootvaltmp[ir] )/tsigma[ir*nmaxroot+nroots] , 2);
                        *chi2data += pow( ( rootvaltmp[ir] )/tsigma[ir*nmaxroot+nroots] , 2);
                    }
                    if (sensi>=1) {
                        if(sensi_meth == AMI_ASA) {
                            x_tmp = NV_DATA_S(x);
                            status = fdtdp(t,dtdp,x_tmp,udata);
                            if (status != AMI_SUCCESS) goto freturn;
                            status = fdtdx(t,dtdx,x_tmp,udata);
                            if (status != AMI_SUCCESS) goto freturn;
                            status = fdrvaldp(t,drvaldp,x_tmp,udata);
                            if (status != AMI_SUCCESS) goto freturn;
                            status = fdrvaldx(t,drvaldx,x_tmp,udata);
                            if (status != AMI_SUCCESS) goto freturn;
                            if (mxIsNaN(tsigma[ir*nmaxroot + nroots])) {
                                status = fsigma_t(t,sigma_t,udata);
                                if (status != AMI_SUCCESS) goto freturn;
                                status = fdsigma_tdp(t,dsigma_tdp,udata);
                                if (status != AMI_SUCCESS) goto freturn;
                            } else {
                                for (ip=0; ip<np; ip++) {
                                    dsigma_ydp[ip*nr+ir] = 0;
                                }
                                sigma_t[ir] = tsigma[ir*nmaxroot + nroots];
                            }
                            for (ip=0; ip<np; ip++) {
                                if(event_model == AMI_NORMAL) {
                                                                    drdp[ip] += dsigma_tdp[ip*nr+ir]/sigma_t[ir] + ( dtdp[ip*nr+ir]* ( t - mt[ir*nmaxroot+ nroots] ) )/pow( sigma_t[ir] , 2) - dsigma_tdp[ip*nr+ir]*pow( ( t - mt[ir*nmaxroot+ nroots] ),2)/pow( sigma_t[ir] , 3);                                }
/*                                if(event_model  == AMI_LOGNORMAL) {
                                      drdp[ip] += 1/(2*pi)*dtdp[ip*nr+ir]/t + ( dtdp[ip*nr+ir]/t * ( log(t) - log(mt[ir*nmaxroot+nroots]) ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                                  }*/
                            }
                            for (ix=0; ix<nx; ix++) {
                                if(event_model  == AMI_NORMAL) {
                                    drdx[nroots+ix*nmaxroot] += ( dtdx[ix*nr+ir] * ( t - mt[ir*nmaxroot+nroots] ) )/pow( sigma_t[ir] , 2);
                                }
/*                                if(event_model  == AMI_LOGNORMAL) {
                                      drdx[nroots+ix*nmaxroot] += 1/(2*pi)*dtdx[ix*nr+ir]/t + ( dtdx[ix*nr+ir]/t * ( log(t) - log(mt[ir*nmaxroot+nroots]) ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                                  }*/
                            }
                            for (ip=0; ip<np; ip++) {
                                if(event_model  == AMI_NORMAL) {
                                     drdp[ip] += dsigma_tdp[ip*nr+ir]/sigma_t[ir] + ( drvaldp[ip*nr+ir]*rootvaltmp[ir] )/pow( sigma_t[ir] , 2) - dsigma_tdp[ip*nr+ir]*pow(rootvaltmp[ir],2)/pow( sigma_t[ir] , 3);
                                }
/*                                if(event_model  == AMI_LOGNORMAL) {
                                      drdp[ip] += 1/(2*pi)*drvaldp[ip*nr+ir]/rootvaltmp[ir] + ( drvaldp[ip*nr+ir]/rootvaltmp[ir] * ( log(rootvaltmp[ir]) ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                                  }*/
                            }
                            for (ix=0; ix<nx; ix++) {
                                if(event_model  == AMI_NORMAL) {
                                    drdx[nroots+ix*nmaxroot] += ( drvaldx[ix*nr+ir] * ( rootvaltmp[ir] ) )/pow( sigma_t[ir] , 2);
                                }
/*                                if(event_model  == AMI_LOGNORMAL) {
                                      drdx[nroots+ix*nmaxroot] += 1/(2*pi)*drvaldx[ix*nr+ir]/rootvaltmp[ir] + ( drvaldx[ix*nr+ir]/rootvaltmp[ir] * ( log(rootvaltmp[ir]) ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                                  }*/
                            }
                        }
                    }
                }
            }
            nroots++;
        }
    }
    
    if (sensi >= 1) {
        /* only set output sensitivities if no errors occured */
        if(sensi_meth == AMI_ASA) {
            if(cv_status == 0) {
                
                setupAMIB(&status, ami_mem, udata, tdata);
                
                nroots--;
                idisc--;
                for (it=nt-1; it>0; it--) {
                    /* enter while */
                    rootflag = TRUE;
                    while (rootflag || discflag) {
                        rootflag = FALSE;
                        discflag = FALSE;
                        if(idisc>-1) {
                            if(discs[idisc]>=ts[it-1]) {
                                if (nroots>-1) {
                                    if (rootdata[nroots + nmaxroot*rootidx[nroots]] > discs[idisc]) {
                                        rootflag = TRUE;
                                    } else {
                                        discflag = TRUE;
                                    }
                                } else {
                                    discflag = TRUE;
                                }
                            }
                        }
                        if (discflag == FALSE && rootflag == FALSE) {
                            if (nroots>-1) {
                                if (rootdata[nroots + nmaxroot*rootidx[nroots]] >= ts[it-1]) {
                                    rootflag = TRUE;
                                }
                            }
                        }
                        
                        if (discflag) {
                            if (discs[idisc] == ts[it]) {
                                
                                status = AMIGetB(ami_mem, which, &t, xB, dxB);
                                if (status != AMI_SUCCESS) goto freturn;
                                status = AMIGetQuadB(ami_mem, which, &t, xQB);
                                if (status != AMI_SUCCESS) goto freturn;

                                status = ideltadisc(t,irdiscs[idisc],x,xB,xQB,udata);
                                if (status != AMI_SUCCESS) goto freturn;
                                status = bdeltadisc(t,irdiscs[idisc],xB,udata);
                                if (status != AMI_SUCCESS) goto freturn;
                                
                                status = AMIReInitB(ami_mem, which, t, xB, dxB);
                                if (status != AMI_SUCCESS) goto freturn;
                                status = AMIQuadReInitB(ami_mem, which, xQB);
                                if (status != AMI_SUCCESS) goto freturn;
                                
                                status = AMICalcICB(ami_mem, which, tstart, xB, dxB);
                                if (status != AMI_SUCCESS) goto freturn;

                                idisc--;
                            } else {
                                cv_status = AMISolveB(ami_mem, discs[idisc], AMI_NORMAL);

                                status = AMIGetB(ami_mem, which, &t, xB, dxB);
                                if (status != AMI_SUCCESS) goto freturn;
                                status = AMIGetQuadB(ami_mem, which, &t, xQB);
                                if (status != AMI_SUCCESS) goto freturn;
                                
                                status = ideltadisc(t,irdiscs[idisc],x,xB,xQB,udata);
                                if (status != AMI_SUCCESS) goto freturn;
                                status = bdeltadisc(t,irdiscs[idisc],xB,udata);
                                if (status != AMI_SUCCESS) goto freturn;
                                
                                if (t == ts[it-1]) {
                                    for (ix=0; ix<nx; ix++) {
                                        xB_tmp[ix] += dgdx[it-1+ix*nt];
                                    }
                                    getDiagnosisB(&status,it-1,ami_mem,udata,rdata,tdata);
                                }
                                
                                status = AMIReInitB(ami_mem, which, t, xB, dxB);
                                if (status != AMI_SUCCESS) goto freturn;
                                status = AMIQuadReInitB(ami_mem, which, xQB);
                                if (status != AMI_SUCCESS) goto freturn;
                                
                                status = AMICalcICB(ami_mem, which, tstart, xB, dxB);
                                if (status != AMI_SUCCESS) goto freturn;
                                
                                idisc--;
                            }
                        }
                        if (rootflag) {

                            cv_status = AMISolveB(ami_mem, rootdata[nroots + nmaxroot*rootidx[nroots]], AMI_NORMAL);
                            
                            status = AMIGetB(ami_mem, which, &t, xB, dxB);
                            if (status != AMI_SUCCESS) goto freturn;
                            status = AMIGetQuadB(ami_mem, which, &t, xQB);
                            if (status != AMI_SUCCESS) goto freturn;
                            
                            xB_tmp = NV_DATA_S(xB);
                            for (ix=0; ix<nx; ix++) {
                                xB_tmp[ix] += drdx[nroots + ix*nmaxroot];
                            }
                            
                            status = AMIReInitB(ami_mem, which, t, xB, dxB);
                            if (status != AMI_SUCCESS) goto freturn;
                            status = AMIQuadReInitB(ami_mem, which, xQB);
                            if (status != AMI_SUCCESS) goto freturn;
                            
                            status = AMICalcICB(ami_mem, which, tstart, xB, dxB);
                            if (status != AMI_SUCCESS) goto freturn;

                            nroots--;
                        }
                    }
                    if(cv_status == 0) {
                        if (nx>0) {
                            /* solve for backward problems */
                            if (ts[it-1] < t) {
                                cv_status = AMISolveB(ami_mem, RCONST(ts[it-1]), AMI_NORMAL);
                            
                                status = AMIGetB(ami_mem, which, &t, xB, dxB);
                                if (status != AMI_SUCCESS) goto freturn;
                                status = AMIGetQuadB(ami_mem, which, &t, xQB);
                                if (status != AMI_SUCCESS) goto freturn;

                                getDiagnosisB(&status,it-1,ami_mem,udata,rdata,tdata);
                                
                                xB_tmp = NV_DATA_S(xB);
                                for (ix=0; ix<nx; ix++) {
                                    xB_tmp[ix] += dgdx[it-1+ix*nt];
                                }
                                
                                status = AMIReInitB(ami_mem, which, t, xB, dxB);
                                if (status != AMI_SUCCESS) goto freturn;
                                status = AMIQuadReInitB(ami_mem, which, xQB);
                                if (status != AMI_SUCCESS) goto freturn;
                                
                                status = AMICalcICB(ami_mem, which, tstart, xB, dxB);
                                if (status != AMI_SUCCESS) goto freturn;
                            }
                        }
                    }
                    
                }
                
                for (ir=0; ir<nr; ir++) {
                    if (nroots>-1) {
                        if (rootdata[nroots + nmaxroot*ir]< t && rootdata[nroots + nmaxroot*ir]>ts[it-1]) {
                            cv_status = AMISolveB(ami_mem, rootdata[nroots + nmaxroot*ir], AMI_NORMAL);
                            
                            status = AMIGetQuadB(ami_mem, which, &t, xQB);
                            if (status != AMI_SUCCESS) goto freturn;
                            status = AMIGetB(ami_mem, which, &t, xB, dxB);
                            xB_tmp = NV_DATA_S(xB);
                            
                            for (ix=0; ix<nx; ix++) {
                                xB_tmp[ix] += drdx[it+ix*nt];
                            }

                            status = AMIReInitB(ami_mem, which, t, xB, dxB);
                            if (status != AMI_SUCCESS) goto freturn;
                            status = AMIQuadReInitB(ami_mem, which, xQB);
                            if (status != AMI_SUCCESS) goto freturn;
                            
                            status = AMICalcICB(ami_mem, which, tstart, xB, dxB);
                            if (status != AMI_SUCCESS) goto freturn;
                            
                            nroots--;
                        }
                    }
                }
                
                if (t>tstart) {
                    if(cv_status == 0) {
                        if (nx>0) {
                            /* solve for backward problems */
                            cv_status = AMISolveB(ami_mem, tstart, AMI_NORMAL);
                            
                            status = AMIGetQuadB(ami_mem, which, &t, xQB);
                            if (status != AMI_SUCCESS) goto freturn;
                            status = AMIGetB(ami_mem, which, &t, xB, dxB);
                            if (status != AMI_SUCCESS) goto freturn;
                        }
                    }
                }
                
                /* evaluate initial values */
                sx = N_VCloneVectorArray_Serial(np,x);
                if (sx == NULL) goto freturn;
                
                status = fx0(x,udata);
                if (status != AMI_SUCCESS) goto freturn;
                status = fdx0(x,dx,udata);
                if (status != AMI_SUCCESS) goto freturn;
                status = fsx0(sx, x, dx, udata);
                if (status != AMI_SUCCESS) goto freturn;

                if(cv_status == 0) {

                    xB_tmp = NV_DATA_S(xB);
                    
                    for (ip=0; ip<np; ip++) {
                        llhS0[ip] = 0.0;
                        sx_tmp = NV_DATA_S(sx[ip]);
                        for (ix = 0; ix < nx; ix++) {
                            llhS0[ip] = llhS0[ip] + xB_tmp[ix] * sx_tmp[ix];
                        }
                    }
                    
                    xQB_tmp = NV_DATA_S(xQB);
                    
                    for(ip=0; ip < np; ip++) {
                        llhSdata[ip] = - llhS0[ip] - dgdp[ip] - drdp[ip] - xQB_tmp[ip];
                    }
                    
                } else {
                    for(ip=0; ip < np; ip++) {
                        llhSdata[ip] = mxGetNaN();
                    }
                }
            } else {
                for(ip=0; ip < np; ip++) {
                    llhSdata[ip] = mxGetNaN();
                }
            }
        }
    }
    
    /* evaluate likelihood */
        
    *llhdata = - g - r;
    
    status = cv_status;
    
    goto freturn;
    
 freturn:
    /* Free memory */
    if(nx>0) {
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(dx);
        N_VDestroy_Serial(xdot);
        AMIFree(&ami_mem);
        DestroyMat(Jtmp);
        if (nr+ndisc>0) {
            if(rootsfound) mxFree(rootsfound);
        }
        if (nr>0) {
            if(rootvaltmp) mxFree(rootvaltmp);
            if(rootidx) mxFree(rootidx);
            if(sigma_t)    mxFree(sigma_t);
        }
        if (ndisc>0) {
            if(discs) mxFree(discs);
            if(irdiscs) mxFree(irdiscs);
        }

        if(sigma_y)    mxFree(sigma_y);
        if (sensi >= 1) {
            N_VDestroyVectorArray_Serial(sx,np);
            if (sensi_meth == AMI_FSA) {
                N_VDestroyVectorArray_Serial(sdx, np);
            }
            if (sensi_meth == AMI_ASA) {
                if(dydx)    mxFree(dydx);
                if(dydp)    mxFree(dydp);
                if(llhS0)     mxFree(llhS0);
                if(dgdp)    mxFree(dgdp);
                if(dgdx)    mxFree(dgdx);
                if (nr>0) {
                    if(dtdp)    mxFree(dtdp);
                    if(dtdx)    mxFree(dtdx);
                }
                if(drdp)    mxFree(drdp);
                if(drdx)    mxFree(drdx);
                if(drvaldp)    mxFree(drvaldp);
                if(drvaldx)    mxFree(drvaldx);
                if(dsigma_ydp)    mxFree(dsigma_ydp);
                if (nr>0) {
                    if(dsigma_tdp)    mxFree(dsigma_tdp);
                }
                if(dxB)      N_VDestroy_Serial(dxB);
                if(xB)      N_VDestroy_Serial(xB);
                if(xQB)     N_VDestroy_Serial(xQB);
            }
            
        }
    }
    
    N_VDestroy_Serial(id);
    
    if(udata)    mxFree(udata);
    if(tdata)    mxFree(tdata);
    
    if(mxGetField(prhs[0], 0 ,"status")) { pstatus = mxGetPr(mxGetField(prhs[0], 0 ,"status")); } else { mexErrMsgTxt("Parameter status not specified as field in solution struct!"); }
    *pstatus = (double) status;
    
    return;
}