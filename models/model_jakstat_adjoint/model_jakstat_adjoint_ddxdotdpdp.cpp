
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_dwdp.h"
#include "model_jakstat_adjoint_w.h"

int ddxdotdpdp_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
int ix;
memset(tdata->ddxdotdpdp,0,sizeof(realtype)*9*udata->nplist*udata->nplist);
status = w_model_jakstat_adjoint(t,x,NULL,tdata);
status = dwdp_model_jakstat_adjoint(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 0: {
      for (int jp=0; jp<udata->nplist; jp++) {
          switch (udata->plist[jp]) {
            case 5: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -x_tmp[0]*tdata->dwdp[0];
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = x_tmp[0]*tdata->dwdp[0];
            } break;

            case 6: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -x_tmp[0]*tdata->dwdp[1];
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = x_tmp[0]*tdata->dwdp[1];
            } break;

            case 7: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -x_tmp[0]*tdata->dwdp[2];
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = x_tmp[0]*tdata->dwdp[2];
            } break;

            case 8: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -x_tmp[0]*tdata->dwdp[3];
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = x_tmp[0]*tdata->dwdp[3];
            } break;

            case 9: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -x_tmp[0]*tdata->dwdp[4];
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = x_tmp[0]*tdata->dwdp[4];
            } break;

          }
      }

  } break;

  case 5: {
      for (int jp=0; jp<udata->nplist; jp++) {
          switch (udata->plist[jp]) {
            case 0: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -x_tmp[0]*tdata->dwdp[0];
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = x_tmp[0]*tdata->dwdp[0];
            } break;

            case 5: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(4,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(4,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 6: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(6,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(6,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 7: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(8,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(8,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 8: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 9: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

          }
      }

  } break;

  case 6: {
      for (int jp=0; jp<udata->nplist; jp++) {
          switch (udata->plist[jp]) {
            case 0: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -x_tmp[0]*tdata->dwdp[1];
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = x_tmp[0]*tdata->dwdp[1];
            } break;

            case 5: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(6,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(6,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 6: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(6,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(6,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 7: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(8,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(8,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 8: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 9: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

          }
      }

  } break;

  case 7: {
      for (int jp=0; jp<udata->nplist; jp++) {
          switch (udata->plist[jp]) {
            case 0: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -x_tmp[0]*tdata->dwdp[2];
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = x_tmp[0]*tdata->dwdp[2];
            } break;

            case 5: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(8,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(8,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 6: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(8,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(8,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 7: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(8,8,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(8,8,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 8: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,8,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,8,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 9: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,8,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,8,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

          }
      }

  } break;

  case 8: {
      for (int jp=0; jp<udata->nplist; jp++) {
          switch (udata->plist[jp]) {
            case 0: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -x_tmp[0]*tdata->dwdp[3];
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = x_tmp[0]*tdata->dwdp[3];
            } break;

            case 5: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 6: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 7: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,8,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,8,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 8: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,10,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(10,10,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 9: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,10,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,10,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

          }
      }

  } break;

  case 9: {
      for (int jp=0; jp<udata->nplist; jp++) {
          switch (udata->plist[jp]) {
            case 0: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -x_tmp[0]*tdata->dwdp[4];
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = x_tmp[0]*tdata->dwdp[4];
            } break;

            case 5: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 6: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 7: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,8,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,8,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 8: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,10,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,10,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

            case 9: {
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 0] = -tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,12,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
  tdata->ddxdotdpdp[udata->plist[ip]*model->nx*udata->nplist + udata->plist[jp]*model->nx + 1] = tdata->p[0]*x_tmp[0]*am_DDspline_pos(12,12,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);
            } break;

          }
      }

  } break;

}
}
return(status);

}


