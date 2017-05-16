/*
 * This file is included by both amici_interface_cpp.cpp and amici_interface_matlab.cpp
 * it depends on the initField* macros defined there.
 */

    initField2(llh,1,1);
    initField2(chi2,1,1);
    initField2(numsteps,nt,1);
    initField2(numrhsevals,nt,1);
    initField2(order,nt,1);
    if(sensi >= AMI_SENSI_ORDER_FIRST){
        initField2(numstepsS,nt,1);
        initField2(numrhsevalsS,nt,1);
    }
    if((nz>0) & (ne>0)){
        initField2(z,nmaxevent,nz);
        initField2(rz,nmaxevent,nz);
        initField2(sigmaz,nmaxevent,nz);
    }
    if(nx>0) {
        initField2(x,nt,nx);
        initField2(xdot,1,nx);
        initField2(J,nx,nx);
    }
    if(ny>0) {
        initField2(y,nt,ny);
        initField2(sigmay,nt,ny);
        if (sensi_meth == AMI_SENSI_SS) {
            initField2(dydp,ny,nplist);
            initField2(dydx,ny,nx);
            initField2(dxdotdp,nx,nplist);
        }
    }
    if(sensi >= AMI_SENSI_ORDER_FIRST) {
        initField2(sllh,nplist,1);
        if (sensi_meth == AMI_SENSI_FSA) {
            initField3(sx,nt,nx,nplist);
            if(ny>0) {
                initField3(sy,nt,ny,nplist);
                initField3(ssigmay,nt,ny,nplist);
            }
            if((nz>0) & (ne>0)){
                initField3(srz,nmaxevent,nz,nplist);
                if(sensi >= AMI_SENSI_ORDER_SECOND){
                    initField4(s2rz,nmaxevent,nz,nplist,nplist);
                }
                initField3(sz,nmaxevent,nz,nplist);
                initField3(ssigmaz,nmaxevent,nz,nplist);
            }
        }
        if (sensi_meth == AMI_SENSI_ASA) {
            if(ny>0) {
                initField3(ssigmay,nt,ny,nplist);
            }
            if((nz>0) & (ne>0)){
                initField3(ssigmaz,nmaxevent,nz,nplist);
            }
        }
        if(sensi >= AMI_SENSI_ORDER_SECOND) {
            initField2(s2llh,ng-1,nplist);
        }
    }
