#ifndef _am_model_jakstat_adjoint_deltasx_h
#define _am_model_jakstat_adjoint_deltasx_h

int deltasx_model_jakstat_adjoint(realtype t, int ie, realtype *deltasx, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, void *user_data);


#endif /* _am_model_jakstat_adjoint_deltasx_h */
