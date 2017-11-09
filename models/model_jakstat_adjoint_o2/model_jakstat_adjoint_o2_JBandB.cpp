
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_JB.h"
#include "model_jakstat_adjoint_o2_w.h"

using namespace amici;

void JBandB_model_jakstat_adjoint_o2(long int NeqBdot, long int mupper, long int mlower, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
return(JB_model_jakstat_adjoint_o2(NeqBdot, t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B));}


