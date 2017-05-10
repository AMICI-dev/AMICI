#ifndef _am_model_jakstat_adjoint_sxdot_h
#define _am_model_jakstat_adjoint_sxdot_h

int sxdot_model_jakstat_adjoint(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip,  N_Vector sx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2);


#endif /* _am_model_jakstat_adjoint_sxdot_h */
