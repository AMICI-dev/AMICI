/*
 * This file is included by both amici_interface_cpp.cpp and amici_interface_matlab.cpp
 * it depends on the initField* macros defined there.
 */

    initField2(llh,1,1);
    initField2(chi2,1,1);
    initField2(numsteps,udata->nt,1);
    initField2(numrhsevals,udata->nt,1);
    initField2(order,udata->nt,1);
    if(udata->sensi >= AMI_SENSI_ORDER_FIRST){
        initField2(numstepsS,udata->nt,1);
        initField2(numrhsevalsS,udata->nt,1);
    }
    if((udata->nz>0) & (udata->ne>0)){
        initField2(z,udata->nmaxevent,udata->nz);
        initField2(rz,udata->nmaxevent,udata->nz);
        initField2(sigmaz,udata->nmaxevent,udata->nz);
    }
    if(udata->nx>0) {
        initField2(x,udata->nt,udata->nx);
        initField2(xdot,1,udata->nx);
        initField2(J,udata->nx,udata->nx);
    }
    if(udata->ny>0) {
        initField2(y,udata->nt,udata->ny);
        initField2(sigmay,udata->nt,udata->ny);
        if (udata->sensi_meth == AMI_SENSI_SS) {
            initField2(dydp,udata->ny,udata->nplist);
            initField2(dydx,udata->ny,udata->nx);
            initField2(dxdotdp,udata->nx,udata->nplist);
        }
    }
    if(udata->sensi >= AMI_SENSI_ORDER_FIRST) {
        initField2(sllh,udata->nplist,1);
        if (udata->sensi_meth == AMI_SENSI_FSA) {
            initField3(sx,udata->nt,udata->nx,udata->nplist);
            if(udata->ny>0) {
                initField3(sy,udata->nt,udata->ny,udata->nplist);
                initField3(ssigmay,udata->nt,udata->ny,udata->nplist);
            }
            if((udata->nz>0) & (udata->ne>0)){
                initField3(srz,udata->nmaxevent,udata->nz,udata->nplist);
                if(udata->sensi >= AMI_SENSI_ORDER_SECOND){
                    initField4(s2rz,udata->nmaxevent,udata->nz,udata->nplist,udata->nplist);
                }
                initField3(sz,udata->nmaxevent,udata->nz,udata->nplist);
                initField3(ssigmaz,udata->nmaxevent,udata->nz,udata->nplist);
            }
        }
        if (udata->sensi_meth == AMI_SENSI_ASA) {
            if(udata->ny>0) {
                initField3(ssigmay,udata->nt,udata->ny,udata->nplist);
            }
            if((udata->nz>0) & (udata->ne>0)){
                initField3(ssigmaz,udata->nmaxevent,udata->nz,udata->nplist);
            }
        }
        if(udata->sensi >= AMI_SENSI_ORDER_SECOND) {
            initField2(s2llh,udata->ng-1,udata->nplist);
        }
    }
