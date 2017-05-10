#ifndef _am_model_jakstat_adjoint_JB_h
#define _am_model_jakstat_adjoint_JB_h

int JB_model_jakstat_adjoint(long int NeqBdot, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);


#endif /* _am_model_jakstat_adjoint_JB_h */
