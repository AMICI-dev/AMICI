#include "include/udata.h"
#include<include/udata_accessors.h>

void freeUserData(UserData *udata) {
    if(udata) {
#ifdef AMICI_WITHOUT_MATLAB
        if(qpositivex) delete[] qpositivex;
        if(plist) delete[] plist;
        if(p) delete[] p;
        if(k) delete[] k;
        if(ts) delete[] ts;
        if(pbar) delete[] pbar;
        if(xbar) delete[] xbar;
        if(idlist) delete[] idlist;
        if(sx0data) delete sx0data;
        if(z2event) delete[] z2event;
#endif
        if(h) delete[] h;
        if(tmp_dxdotdp) delete[] tmp_dxdotdp;
        if(w_tmp) delete[] w_tmp;
        if(dwdx_tmp) delete[] dwdx_tmp;
        if(dwdp_tmp) delete[] dwdp_tmp;
        if(M_tmp) delete[] M_tmp;
        if(dfdx_tmp) delete[] dfdx_tmp;
        if(stau_tmp) delete[] stau_tmp;
        SparseDestroyMat(tmp_J);
    }

    delete udata;
}
