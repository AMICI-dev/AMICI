#ifndef _am_model_jakstat_adjoint_dJydx_h
#define _am_model_jakstat_adjoint_dJydx_h

int dJydx_model_jakstat_adjoint(realtype t, int it, realtype *dJydx, realtype *y, N_Vector x, realtype *dydx, realtype *my, realtype *sigma_y, void *user_data);


#endif /* _am_model_jakstat_adjoint_dJydx_h */
