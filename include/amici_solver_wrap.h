#ifndef AMICI_SOLVER_WRAP_H
#define AMICI_SOLVER_WRAP_H

int AMIAdjInit(void *mem, long int steps, int interp);
int AMIBandB(void *mem, int which, int nx, int ubw, int lbw);
int AMIBand(void *mem, int nx, int ubw, int lbw);
int AMICalcICB(void *mem, int which, realtype tout1, N_Vector xB, N_Vector dxB);
int AMICalcIC(void *mem, realtype tout1);
int AMICreateB(void *mem, int lmm, int iter, int *which);
int AMIDenseB(void *mem, int which, int nx);
int AMIDense(void *mem, int nx);
int AMIDiagB(void *mem, int which);
int AMIDiag(void *mem);
int AMIGetB(void *mem, int which, realtype *tret, N_Vector yy, N_Vector yp);
int AMIGetDky(void *mem, realtype t, int k, N_Vector dky);
int AMIGetLastOrder(void *mem, int *order);
int AMIGetNumErrTestFails(void *mem, long int *numerrtestfails);
int AMIGetNumNonlinSolvConvFails(void *mem, long int *numnonlinsolvconvfails);
int AMIGetNumRhsEvals(void *mem, long int *numrhsevals);
int AMIGetNumSteps(void *mem, long int *numsteps);
int AMIGetQuadB(void *mem, int which, realtype *tret, N_Vector qB);
int AMIGetRootInfo(void *mem, int *rootsfound);
int AMIGetSens(void *mem, realtype *tret, N_Vector *yySout);
int AMIKLUB(void *mem, int which, int nx, int nnz, int sparsetype);
int AMIKLUSetOrderingB(void *mem, int which, int ordering);
int AMIKLUSetOrdering(void *mem, int ordering);
int AMIKLU(void *mem, int nx, int nnz, int sparsetype);
int AMIQuadReInitB(void *mem, int which, N_Vector yQB0);
int AMIQuadSStolerancesB(void *mem, int which, realtype reltolQB,
                         realtype abstolQB);
int AMIReInitB(void *mem, int which, realtype tB0, N_Vector yyB0,
               N_Vector ypB0);
int AMIReInit(void *mem, realtype t0, N_Vector yy0, N_Vector yp0);
int AMIRootInit(void *mem, int nrtfn, CVRootFn ptr);
int AMISensEEtolerances(void *mem);
int AMISensReInit(void *mem, int ism, N_Vector *yS0, N_Vector *ypS0);
int AMISetErrHandlerFn(void *mem);
int AMISetId(void *mem, N_Vector id);
int AMISetMaxNumStepsB(void *mem, int which, long int mxstepsB);
int AMISetMaxNumSteps(void *mem, long int mxsteps);
int AMISetQuadErrConB(void *mem, int which, bool flag);
int AMISetSensErrCon(void *mem, bool error_corr);
int AMISetSensParams(void *mem, realtype *p, realtype *pbar, int *plist);
int AMISetStabLimDetB(void *mem, int which, int stldet);
int AMISetStabLimDet(void *mem, int stldet);
int AMISetStopTime(void *mem, realtype tstop);
int AMISetSuppressAlg(void *mem, bool flag);
int AMISetUserDataB(void *mem, int which, void *user_data);
int AMISetUserData(void *mem, void *user_data);
int AMISolveB(void *mem, realtype tBout, int itaskB);
int AMISolveF(void *mem, realtype tout, N_Vector yret, N_Vector ypret,
              realtype *tret, int itask, int *ncheckPtr);
int AMISolve(void *mem, realtype tout, N_Vector yret, N_Vector ypret,
             realtype *tret, int itask);
int AMISpbcgB(void *mem, int which, int prectype, int maxl);
int AMISpbcg(void *mem, int prectype, int maxl);
int AMISpgmrB(void *mem, int which, int prectype, int maxl);
int AMISpgmr(void *mem, int prectype, int maxl);
int AMISptfqmrB(void *mem, int which, int prectype, int maxl);
int AMISptfqmr(void *mem, int prectype, int maxl);
int AMISStolerancesB(void *mem, int which, realtype relTolB, realtype absTolB);
int AMISStolerances(void *mem, double rtol, double atol);
void *AMICreate(int lmm, int iter);
void AMIFree(void **mem);
void *AMIGetAdjBmem(void *mem, int which);

#endif // AMICI_SOLVER_WRAP_H
