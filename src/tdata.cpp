#include "include/tdata.h"

#include <include/udata.h>
#include <include/rdata.h>
#include <include/amici_model.h>

TempData::TempData(const UserData *udata, Model *model, ReturnData *rdata) : udata(udata), model(model), rdata(rdata) {
    
    xB = xB_old = dxB = xQB = xQB_old = NULL;
    x_disc = xdot_disc = xdot_old_disc = NULL;
    x = x_old = dx = dx_old = xdot = xdot_old = NULL;
    sx = sdx = NULL;
    Jtmp = NULL;
    dydx = dydp = dJydp = dJydx = dJydy = dzdp = dzdx = drzdp = drzdx = dJzdp = dJzdx = dJzdz = dJrzdz = NULL;
    dJydsigma = dJzdsigma = dJrzdsigma = dsigmaydp = dsigmazdp = llhS0 = NULL;
    
    which = 0;
    
    p = new realtype[model->np];
    udata->unscaleParameters(model, p);

    nplist = udata->nplist;
    
    Jy = new realtype[model->nJ]();
    Jz = new realtype[model->nJ]();
    sigmay = new realtype[model->ny]();
    sigmaz = new realtype[model->nz]();

    if(model->nx>0) {
        x = N_VNew_Serial(model->nx);
        x_old = N_VNew_Serial(model->nx);
        dx = N_VNew_Serial(model->nx); /* only needed for idas */
        dx_old = N_VNew_Serial(model->nx); /* only needed for idas */
        xdot = N_VNew_Serial(model->nx);
        xdot_old = N_VNew_Serial(model->nx);
        Jtmp = NewDenseMat(model->nx,model->nx);

        /* initialise temporary jacobian storage */
        J = SparseNewMat(model->nx,model->nx,model->nnz,CSC_MAT);
        M = new realtype[model->nx*model->nx]();
        dfdx = new realtype[model->nx*model->nx]();
    }

    if (model->ne>0) {
        /* initialise temporary stau storage */
        stau = new realtype[nplist]();
    }

    w = new realtype[model->nw]();
    dwdx = new realtype[model->ndwdx]();
    dwdp = new realtype[model->ndwdp]();

    
    /* EVENTS */
    rootsfound = new int[model->ne]();
    rootvals = new realtype[model->ne]();
    h = new realtype[model->ne]();
    h_udata = new realtype[model->ne]();

    rootidx = new int[udata->nmaxevent*model->ne*model->ne]();
    nroots = new int[model->ne]();
    discs = new realtype[udata->nmaxevent*model->ne]();
    deltax = new realtype[model->nx]();
    deltasx = new realtype[model->nx*udata->nplist]();
    deltaxB = new realtype[model->nx]();
    deltaqB = new realtype[model->nJ*udata->nplist]();
    
    /* SENSITIVITIES */
    if(udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        /* initialise temporary dxdotdp storage */
        dxdotdp = new realtype[model->nx*nplist]();

        dydx = new realtype[model->ny * model->nx]();
        dydp = new realtype[model->ny * udata->nplist]();
        dJydp = new realtype[model->nJ * udata->nplist]();
        dJydx = new realtype[model->nJ * model->nx * udata->nt]();
        dJydy = new realtype[model->nJ * model->nytrue * model->ny]();
        dJydsigma = new realtype[model->nJ * model->nytrue * model->ny]();
        dzdx = new realtype[model->nz*model->nx]();
        dzdp = new realtype[model->nz*udata->nplist]();
        drzdx = new realtype[model->nz*model->nx]();
        drzdp = new realtype[model->nz*udata->nplist]();
        dJzdp = new realtype[model->nJ * udata->nplist]();
        dJzdx = new realtype[model->nJ * model->nx * udata->nmaxevent]();
        dJzdz = new realtype[model->nJ * model->nztrue * model->nz]();
        dJzdsigma = new realtype[model->nJ * model->nztrue * model->nz]();
        dJrzdz = new realtype[model->nJ * model->nztrue * model->nz]();
        dJrzdsigma = new realtype[model->nJ * model->nztrue * model->nz]();
        
        dsigmaydp = new realtype[model->ny * udata->nplist]();
        dsigmazdp = new realtype[model->nz * udata->nplist]();
        
        if(udata->nplist && x) {
            sx = N_VCloneVectorArray_Serial(udata->nplist, x);
            sdx = N_VCloneVectorArray_Serial(udata->nplist, x);
        }
        
        if (udata->sensi_meth == AMICI_SENSI_ASA) {
            llhS0 = new realtype[model->nJ * udata->nplist]();
            if(model->ne > 0 && udata->nmaxevent > 0 && x) {
                x_disc = N_VCloneVectorArray_Serial(model->ne * udata->nmaxevent, x);
                xdot_disc = N_VCloneVectorArray_Serial(model->ne * udata->nmaxevent, x);
                xdot_old_disc = N_VCloneVectorArray_Serial(model->ne * udata->nmaxevent, x);
            }
            if(model->nx > 0 && udata->nplist > 0) {
                xB = N_VNew_Serial(model->nxtrue * model->nJ);
                xB_old = N_VNew_Serial(model->nxtrue * model->nJ);
                dxB = N_VNew_Serial(model->nxtrue * model->nJ);
                xQB = N_VNew_Serial(model->nJ * udata->nplist);
                xQB_old = N_VNew_Serial(model->nJ * udata->nplist);
            }
        }
    }
}

TempData::~TempData() {
    if(p) delete[] p;

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
    if(h_udata) delete[] h_udata;
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

