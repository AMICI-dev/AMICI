#include "include/tdata.h"

#include <include/udata.h>
#include <include/amici_model.h>

TempData::TempData(const UserData *udata, Model *model) : udata(udata), model(model) {
    
    xB = xB_old = dxB = xQB = xQB_old = NULL;
    x_disc = xdot_disc = xdot_old_disc = NULL;
    x = x_old = dx = dx_old = xdot = xdot_old = NULL;
    sx = sdx = NULL;
    Jtmp = NULL;
    dydx = dydp = dJydp = dJydx = dJydy = dzdp = dzdx = drzdp = drzdx = dJzdp = dJzdx = dJzdz = dJrzdz = NULL;
    dJydsigma = dJzdsigma = dJrzdsigma = dsigmaydp = dsigmazdp = llhS0 = NULL;
    
    which = 0;
    
    nplist = udata->nplist;
    
    Jy = new realtype[udata->nJ]();
    Jz = new realtype[udata->nJ]();
    sigmay = new realtype[udata->ny]();
    sigmaz = new realtype[udata->nz]();

    if(udata->nx>0) {
        x = N_VNew_Serial(udata->nx);
        x_old = N_VNew_Serial(udata->nx);
        dx = N_VNew_Serial(udata->nx); /* only needed for idas */
        dx_old = N_VNew_Serial(udata->nx); /* only needed for idas */
        xdot = N_VNew_Serial(udata->nx);
        xdot_old = N_VNew_Serial(udata->nx);
        Jtmp = NewDenseMat(udata->nx,udata->nx);

        /* initialise temporary jacobian storage */
        J = SparseNewMat(udata->nx,udata->nx,udata->nnz,CSC_MAT);
        M = new realtype[udata->nx*udata->nx]();
        dfdx = new realtype[udata->nx*udata->nx]();
    }

    if (udata->ne>0) {
        /* initialise temporary stau storage */
        stau = new realtype[nplist]();
    }

    w = new realtype[udata->nw]();
    dwdx = new realtype[udata->ndwdx]();
    dwdp = new realtype[udata->ndwdp]();

    
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
    
    /* SENSITIVITIES */
    if(udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        /* initialise temporary dxdotdp storage */
        dxdotdp = new realtype[udata->nx*nplist]();

        dydx = new realtype[udata->ny * udata->nx]();
        dydp = new realtype[udata->ny * udata->nplist]();
        dJydp = new realtype[udata->nJ * udata->nplist]();
        dJydx = new realtype[udata->nJ * udata->nx * udata->nt]();
        dJydy = new realtype[udata->nJ * udata->nytrue * udata->ny]();
        dJydsigma = new realtype[udata->nJ * udata->nytrue * udata->ny]();
        dzdx = new realtype[udata->nz*udata->nx]();
        dzdp = new realtype[udata->nz*udata->nplist]();
        drzdx = new realtype[udata->nz*udata->nx]();
        drzdp = new realtype[udata->nz*udata->nplist]();
        dJzdp = new realtype[udata->nJ * udata->nplist]();
        dJzdx = new realtype[udata->nJ * udata->nx * udata->nmaxevent]();
        dJzdz = new realtype[udata->nJ * udata->nztrue * udata->nz]();
        dJzdsigma = new realtype[udata->nJ * udata->nztrue * udata->nz]();
        dJrzdz = new realtype[udata->nJ * udata->nztrue * udata->nz]();
        dJrzdsigma = new realtype[udata->nJ * udata->nztrue * udata->nz]();
        
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
                xB = N_VNew_Serial(udata->nxtrue * udata->nJ);
                xB_old = N_VNew_Serial(udata->nxtrue * udata->nJ);
                dxB = N_VNew_Serial(udata->nxtrue * udata->nJ);
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
    if(dJydsigma) delete[] dJydsigma;
    if(dJydx) delete[] dJydx;
    if(dJzdp) delete[] dJzdp;
    if(dJzdz) delete[] dJzdz;
    if(dJzdsigma) delete[] dJzdsigma;
    if(dJrzdz) delete[] dJrzdz;
    if(dJrzdsigma) delete[] dJrzdsigma;
    if(dJzdx) delete[] dJzdx;
    if(dzdp) delete[] dzdp;
    if(dzdx) delete[] dzdx;
    if(drzdp) delete[] drzdp;
    if(drzdx) delete[] drzdx;
    if(dsigmaydp) delete[] dsigmaydp;
    if(dsigmazdp) delete[] dsigmazdp;
        
    if(llhS0) delete[] llhS0;

    if(dxdotdp) delete[] dxdotdp;
    if(w) delete[] w;
    if(dwdx) delete[] dwdx;
    if(dwdp) delete[] dwdp;
    if(M) delete[] M;
    if(dfdx) delete[] dfdx;
    if(stau) delete[] stau;
    if(J) SparseDestroyMat(J);

}

