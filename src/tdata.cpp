#include "include/tdata.h"

#include <include/udata.h>

TempData::TempData(const UserData *udata) {
    
    xB = xB_old = dxB = xQB = xQB_old = NULL;
    x_disc = xdot_disc = xdot_old_disc = NULL;
    x = x_old = dx = dx_old = xdot = xdot_old = NULL;
    sx = sdx = NULL;
    Jtmp = NULL;
    dydx = dydp = dJydp = dJydx = dJydy = dzdp = dzdx = dJzdp = dJzdx = NULL;
    dsigmaydp = dsigmazdp = llhS0 = NULL;
    
    which = 0;
    
    nplist = udata->nplist;
    
    Jy = new realtype[udata->nJ]();
    Jz = new realtype[udata->nJ]();
    
    sigmay = new realtype[udata->ny]();
    
    if(udata->nx>0) {
        x = N_VNew_Serial(udata->nx);
        x_old = N_VNew_Serial(udata->nx);
        dx = N_VNew_Serial(udata->nx); /* only needed for idas */
        dx_old = N_VNew_Serial(udata->nx); /* only needed for idas */
        xdot = N_VNew_Serial(udata->nx);
        xdot_old = N_VNew_Serial(udata->nx);
        Jtmp = NewDenseMat(udata->nx,udata->nx);
    }
    
    /* EVENTS */
    rootsfound = new int[udata->ne]();
    rootvals = new realtype[udata->ne]();
    h = new realtype[udata->ne]();
    rootidx = new int[udata->nmaxevent*udata->ne*udata->ne]();
    nroots = new int[udata->ne]();
    discs = new realtype[udata->nmaxevent*udata->ne]();
    deltax = new realtype[udata->nx]();
    deltasx = new realtype[udata->nx*udata->nplist]();
    deltaxB = new realtype[udata->nx]();
    deltaqB = new realtype[udata->nJ*udata->nplist]();
    sigmaz = new realtype[udata->nz]();

    
    /* SENSITIVITIES */
    if(udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        dydx = new realtype[udata->ny * udata->nx]();
        dydp = new realtype[udata->ny * udata->nplist]();
        dJydp = new realtype[udata->nJ * udata->nplist * udata->nytrue]();
        dJydx = new realtype[udata->nJ * udata->nxtrue * udata->nt]();
        dJydy = new realtype[udata->nytrue * udata->nJ * udata->ny]();
        dzdp = new realtype[udata->nz*udata->nplist]();
        dzdx = new realtype[udata->nz*udata->nx]();
        dJzdp = new realtype[udata->nJ * udata->nplist * udata->nztrue * udata->nmaxevent]();
        dJzdx = new realtype[udata->nJ * udata->nx * udata->nztrue * udata->nmaxevent]();
        
        dsigmaydp = new realtype[udata->ny * udata->nplist]();
        dsigmazdp = new realtype[udata->nz * udata->nplist]();
        
        if(udata->nplist && x) {
            sx = N_VCloneVectorArray_Serial(udata->nplist, x);
            sdx = N_VCloneVectorArray_Serial(udata->nplist, x);
        }
        
        if (udata->sensi_meth == AMICI_SENSI_ASA) {
            llhS0 = new realtype[udata->nJ * udata->nplist]();
            if(udata->ne > 0 && udata->nmaxevent > 0 && x) {
                x_disc = N_VCloneVectorArray_Serial(udata->ne * udata->nmaxevent, x);
                xdot_disc = N_VCloneVectorArray_Serial(udata->ne * udata->nmaxevent, x);
                xdot_old_disc = N_VCloneVectorArray_Serial(udata->ne * udata->nmaxevent, x);
            }
            if(udata->nx > 0 && udata->nplist > 0) {
                xB = N_VNew_Serial(udata->nx);
                xB_old = N_VNew_Serial(udata->nx);
                dxB = N_VNew_Serial(udata->nx);
                xQB = N_VNew_Serial(udata->nJ * udata->nplist);
                xQB_old = N_VNew_Serial(udata->nJ * udata->nplist);
            }
        }
    }
}

TempData::~TempData() {
    
    if(x) N_VDestroy_Serial(x);
    if(dx) N_VDestroy_Serial(dx);
    if(xdot) N_VDestroy_Serial(xdot);
    if(x_old) N_VDestroy_Serial(x_old);
    if(dx_old) N_VDestroy_Serial(dx_old);
    if(xdot_old) N_VDestroy_Serial(xdot_old);
    if(Jtmp) DestroyMat(Jtmp);
    
    if(dxB) N_VDestroy_Serial(dxB);
    if(xB) N_VDestroy_Serial(xB);
    if(xB_old) N_VDestroy_Serial(xB_old);
    if(xQB) N_VDestroy_Serial(xQB);
    if(xQB_old) N_VDestroy_Serial(xQB_old);

    
    if(sx) N_VDestroyVectorArray_Serial(sx,nplist);
    if(sdx) N_VDestroyVectorArray_Serial(sdx,nplist);
    
    if(Jy) delete[] Jy;
    if(Jz) delete[] Jz;
        
    
    if(rootsfound) delete[] rootsfound;
    if(rootvals) delete[] rootvals;
    if(rootidx) delete[] rootidx;
    if(sigmaz) delete[] sigmaz;
    if(nroots) delete[] nroots;
    if(discs) delete[] discs;
        
    if(deltax) delete[] deltax;
    if(deltasx) delete[] deltasx;
    if(deltaxB) delete[] deltaxB;
    if(deltaqB) delete[] deltaqB;
    if(h) delete[] h;
    if(sigmay)    delete[] sigmay;
    if(dydx) delete[] dydx;
    if(dydp) delete[] dydp;
    if(dJydp) delete[] dJydp;
    if(dJydy) delete[] dJydy;
    if(dJydx) delete[] dJydx;
    if(dJzdp) delete[] dJzdp;
    if(dJzdx) delete[] dJzdx;
    if(dzdp) delete[] dzdp;
    if(dzdx) delete[] dzdx;
    if(dsigmaydp) delete[] dsigmaydp;
    if(dsigmazdp) delete[] dsigmazdp;
        
    if(llhS0) delete[] llhS0;
}

