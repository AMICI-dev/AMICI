#include "include/udata.h"
#include<include/udata_accessors.h>

void freeUserData(UserData udata) {
    if(udata) {
        if(qpositivex) free(qpositivex);
        if(plist) free(plist);
        if(p) free(p);
        if(k) free(k);
        if(ts) free(ts);
        if(pbar) free(pbar);
        if(xbar) free(xbar);
        if(idlist) free(idlist);
        if(sx0data) free(sx0data);
        if(z2event) free(z2event);
        if(h) free(h);
        if(tmp_dxdotdp) free(tmp_dxdotdp);
        if(w_tmp) free(w_tmp);
        if(dwdx_tmp) free(dwdx_tmp);
        if(dwdp_tmp) free(dwdp_tmp);
        if(M_tmp) free(M_tmp);
        if(dfdx_tmp) free(dfdx_tmp);
        if(stau_tmp) free(stau_tmp);
        DestroySparseMat(tmp_J);
    }

    free(udata);
}
