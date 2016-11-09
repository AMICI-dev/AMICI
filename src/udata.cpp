#include "include/udata.h"
#include<include/udata_accessors.h>
#include <stdio.h>

void freeUserData(UserData *udata) {
    if(udata) {
#ifdef AMICI_WITHOUT_MATLAB
        if(qpositivex) delete[] qpositivex;
        if(p) delete[] p;
        if(k) delete[] k;
        if(ts) delete[] ts;
        if(pbar) delete[] pbar;
        if(xbar) delete[] xbar;
        if(idlist) delete[] idlist;
        if(sx0data) delete sx0data;
        if(z2event) delete[] z2event;
#endif
        if(plist) delete[] plist;
        if(h) delete[] h;
        if(tmp_dxdotdp) delete[] tmp_dxdotdp;
        if(w_tmp) delete[] w_tmp;
        if(dwdx_tmp) delete[] dwdx_tmp;
        if(dwdp_tmp) delete[] dwdp_tmp;
        if(M_tmp) delete[] M_tmp;
        if(dfdx_tmp) delete[] dfdx_tmp;
        if(stau_tmp) delete[] stau_tmp;
        if(tmp_J) SparseDestroyMat(tmp_J);

        delete udata;
    }
}

void printUserData(UserData *udata) {
    printf("am_qpositivex: %p\n", udata->am_qpositivex);
    printf("am_plist: %p\n", udata->am_plist);
    printf("am_np: %d\n", udata->am_np);
    printf("am_ny: %d\n", udata->am_ny);
    printf("am_nytrue: %d\n", udata->am_nytrue);
    printf("am_nx: %d\n", udata->am_nx);
    printf("am_nxtrue: %d\n", udata->am_nxtrue);
    printf("am_nz: %d\n", udata->am_nz);
    printf("am_nztrue: %d\n", udata->am_nztrue);
    printf("am_ne: %d\n", udata->am_ne);
    printf("am_nt: %d\n", udata->am_nt);
    printf("am_ng: %d\n", udata->am_ng);
    printf("am_nw: %d\n", udata->am_nw);
    printf("am_ndwdx: %d\n", udata->am_ndwdx);
    printf("am_nnz: %d\n", udata->am_nnz);
    printf("am_nmaxevent: %d\n", udata->am_nmaxevent);
    printf("am_pscale: %d\n", (int)udata->am_pscale);
    printf("am_p: %p\n", udata->am_p);
    printf("am_k: %p\n", udata->am_k);
    printf("am_tstart: %e\n", udata->am_tstart);
    printf("am_ts: %p\n", udata->am_ts);
    printf("am_pbar: %p\n", udata->am_pbar);
    printf("am_xbar: %p\n", udata->am_xbar);
    printf("am_idlist: %p\n", udata->am_idlist);
    printf("am_sensi: %d\n", udata->am_sensi);
    printf("am_atol: %e\n", udata->am_atol);
    printf("am_rtol: %e\n", udata->am_rtol);
    printf("am_maxsteps: %d\n", udata->am_maxsteps);
    printf("am_ism: %d\n", udata->am_ism);
    printf("am_sensi_meth: %d\n", udata->am_sensi_meth);
    printf("am_linsol: %d\n", udata->am_linsol);
    printf("am_interpType: %d\n", udata->am_interpType);
    printf("am_lmm: %d\n", udata->am_lmm);
    printf("am_iter: %d\n", udata->am_iter);
    printf("am_stldet: %d\n", udata->am_stldet);
    printf("am_ubw: %d\n", udata->am_ubw);
    printf("am_lbw: %d\n", udata->am_lbw);
    printf("am_bx0: %d\n", udata->am_bx0);
    printf("am_bsx0: %d\n", udata->am_bsx0);
    printf("am_x0data: %p\n", udata->am_x0data);
    printf("am_sx0data: %p\n", udata->am_sx0data);
    printf("am_ordering: %d\n", udata->am_ordering);
    printf("am_z2event: %p\n", udata->am_z2event);
    printf("am_h: %p\n", udata->am_h);
}
