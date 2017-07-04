#ifndef _am_model_jakstat_adjoint_dJzdx_h
#define _am_model_jakstat_adjoint_dJzdx_h

int dJzdx_model_jakstat_adjoint(realtype t, int ie, N_Vector x, realtype *z, realtype *mz, realtype *dzdx,  void *user_data, TempData *tdata, ReturnData *rdata);


#endif /* _am_model_jakstat_adjoint_dJzdx_h */
