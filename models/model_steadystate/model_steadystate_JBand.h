#ifndef _am_model_steadystate_JBand_h
#define _am_model_steadystate_JBand_h

int JBand_model_steadystate(long int N, long int mupper, long int mlower, realtype t, N_Vector x, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


#endif /* _am_model_steadystate_JBand_h */
