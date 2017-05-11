#ifndef _am_model_jakstat_adjoint_JBandB_h
#define _am_model_jakstat_adjoint_JBandB_h

int JBandB_model_jakstat_adjoint(long int NeqBdot, long int mupper, long int mlower, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);


#endif /* _am_model_jakstat_adjoint_JBandB_h */
