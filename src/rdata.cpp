
#include "include/rdata.h"

ReturnData::ReturnData()
{
    ts = xdot = dxdotdp = dydx = dydp = J = z = sigmaz = sz = ssigmaz = rz = srz = s2rz = x = sx = y = sigmay = NULL;
    sy = ssigmay = numsteps = numstepsS = numrhsevals = numrhsevalsS = order = llh = chi2 = sllh = s2llh = NULL;
    status = NULL;
}

ReturnData::ReturnData(const UserData *udata)
{
    ReturnData();
    initFields(udata);
}

ReturnData::~ReturnData()
{
    if(ts) delete[] ts;
    if(xdot) delete[] xdot;
    if(dxdotdp) delete[] dxdotdp;
    if(dydx) delete[] dydx;
    if(dydp) delete[] dydp;
    if(J) delete[] J;
    if(z) delete[] z;
    if(sigmaz) delete[] sigmaz;
    if(sz) delete[] sz;
    if(ssigmaz) delete[] ssigmaz;
    if(rz) delete[] rz;
    if(srz) delete[] srz;
    if(s2rz) delete[] s2rz;
    if(x) delete[] x;
    if(sx) delete[] sx;
    if(y) delete[] y;
    if(sigmay) delete[] sigmay;
    if(sy) delete[] sy;
    if(ssigmay) delete[] ssigmay;
    if(numsteps) delete[] numsteps;
    if(numrhsevals) delete[] numrhsevals;
    if(order) delete[] order;
    if(numstepsS) delete[] numstepsS;
    if(numrhsevalsS) delete[] numrhsevalsS;
    if(llh) delete[] llh;
    if(sllh) delete[] sllh;
    if(s2llh) delete[] s2llh;
    if(chi2) delete[] chi2;
    if(status) delete[] status;
}

void ReturnData::initFields(const UserData *udata)
{
    initField1(&status, "status", 1);
    initField1(&ts, "t", udata->nt);
    initField2(&llh, "llh", 1, 1);
    initField2(&chi2, "chi", 1, 1);
    initField2(&numsteps, "numsteps",udata->nt, 1);
    initField2(&numrhsevals, "numrhsevals", udata->nt, 1);
    initField2(&order, "order", udata->nt,1);

    if((udata->nz>0) & (udata->ne>0)){
        initField2(&z, "z", udata->nmaxevent,udata->nz);
        initField2(&rz, "rz", udata->nmaxevent,udata->nz);
        initField2(&sigmaz, "sigmaz", udata->nmaxevent,udata->nz);
    }
    if(udata->nx>0) {
        initField2(&x, "x", udata->nt,udata->nx);
        initField2(&xdot, "xdot", 1,udata->nx);
        initField2(&J, "J", udata->nx,udata->nx);
    }
    if(udata->ny>0) {
        initField2(&y, "y", udata->nt,udata->ny);
        initField2(&sigmay, "sigmay", udata->nt,udata->ny);
        if (udata->sensi_meth == AMI_SENSI_SS) {
            initField2(&dydp, "dydp", udata->ny,udata->nplist);
            initField2(&dydx, "dydx", udata->ny,udata->nx);
            initField2(&dxdotdp, "dxdotdp", udata->nx,udata->nplist);
        }
    }
    if(udata->sensi >= AMI_SENSI_ORDER_FIRST) {
        initField2(&numstepsS, "numstepsS", udata->nt, 1);
        initField2(&numrhsevalsS, "numrhsevalsS", udata->nt, 1);
        initField2(&sllh, "sllh", udata->nplist,1);

        if (udata->sensi_meth == AMI_SENSI_FSA) {
            initField3(&sx, "sx", udata->nt,udata->nx,udata->nplist);
            if(udata->ny>0) {
                initField3(&sy, "sy", udata->nt,udata->ny,udata->nplist);
                initField3(&ssigmay, "ssigmay", udata->nt,udata->ny,udata->nplist);
            }
            if((udata->nz>0) & (udata->ne>0)){
                initField3(&srz, "srz", udata->nmaxevent,udata->nz,udata->nplist);
                if(udata->sensi >= AMI_SENSI_ORDER_SECOND){
                    initField4(&s2rz, "s2rz", udata->nmaxevent,udata->nz,udata->nplist,udata->nplist);
                }
                initField3(&sz, "sz", udata->nmaxevent,udata->nz,udata->nplist);
                initField3(&ssigmaz, "ssigmaz", udata->nmaxevent,udata->nz,udata->nplist);
            }
        }

        if (udata->sensi_meth == AMI_SENSI_ASA) {
            if(udata->ny>0) {
                initField3(&ssigmay, "ssigmay", udata->nt,udata->ny,udata->nplist);
            }
            if((udata->nz>0) & (udata->ne>0)){
                initField3(&ssigmaz, "ssigmaz", udata->nmaxevent,udata->nz,udata->nplist);
            }
        }

        if(udata->sensi >= AMI_SENSI_ORDER_SECOND) {
            initField2(&s2llh, "s2llh", udata->ng-1,udata->nplist);
        }
    }

}

void ReturnData::initField1(double **fieldPointer, const char *fieldName, int dim)
{
    *fieldPointer = new double[dim]();
}

void ReturnData::initField2(double **fieldPointer, const char *fieldName, int dim1, int dim2)
{
    *fieldPointer = new double[(dim1) * (dim2)]();

}

void ReturnData::initField3(double **fieldPointer, const char *fieldName, int dim1, int dim2, int dim3)
{
    *fieldPointer = new double[(dim1) * (dim2) * (dim3)]();
}

void ReturnData::initField4(double **fieldPointer, const char *fieldName, int dim1, int dim2, int dim3, int dim4)
{
    *fieldPointer = new double[(dim1) * (dim2) * (dim3) * (dim4)]();
}
