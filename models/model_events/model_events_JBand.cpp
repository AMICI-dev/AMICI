
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void JBand_model_events(long int N, long int mupper, long int mlower, realtype t, N_Vector x, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
return(J_model_events(N, t, x, xdot, J, user_data, tmp1, tmp2, tmp3));}

