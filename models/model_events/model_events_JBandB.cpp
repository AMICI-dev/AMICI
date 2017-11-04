
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void JBandB_model_events(long int NeqBdot, long int mupper, long int mlower, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
return(JB_model_events(NeqBdot, t, x, xB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B));}

