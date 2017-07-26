#include "include/rdata.h"
#include "include/udata.h"
#include "include/amici_misc.h"
#include "include/symbolic_functions.h"

ReturnData::ReturnData()
{
    setDefaults();
}

ReturnData::ReturnData(const UserData *udata)
{
    setDefaults();

    initFields(udata);
}

void ReturnData::setDefaults()
{
    ts = xdot = dxdotdp = dydx = dydp = J = z = sigmaz = sz = ssigmaz = rz = srz = s2rz = x = sx = y = sigmay = NULL;
    sy = ssigmay = numsteps = numstepsB = numrhsevals = numrhsevalsB = order = llh = chi2 = sllh = s2llh = NULL;
    numerrtestfails = numnonlinsolvconvfails = numerrtestfailsB = numnonlinsolvconvfailsB = NULL;
    newton = newton_time = newton_numsteps = newton_numlinsteps = xss = NULL;
    status = NULL;

    freeFieldsOnDestruction = true;
}

void ReturnData::invalidate(const UserData *udata)
{
    if (llh)
        *llh = amiGetNaN();

    if (sllh)
        setLikelihoodSensitivityFirstOrderNaN(udata);

    if (s2llh)
        setLikelihoodSensitivitySecondOrderNaN(udata);

}

void ReturnData::setLikelihoodSensitivityFirstOrderNaN(const UserData *udata)
{
    fillArray(sllh, udata->nplist, amiGetNaN());
}

void ReturnData::setLikelihoodSensitivitySecondOrderNaN(const UserData *udata)
{
    fillArray(s2llh, udata->nplist*(udata->nJ-1), amiGetNaN());
}

int ReturnData::applyChainRuleFactorToSimulationResults(const UserData *udata, const ExpData *edata)
{
    if (udata->pscale == AMICI_SCALING_NONE)
        return AMICI_SUCCESS;

    // chain-rule factor: multiplier for am_p
    realtype coefficient;
    realtype *pcoefficient, *augcoefficient;

    pcoefficient = new realtype[udata->nplist]();
    augcoefficient = new realtype[udata->np]();

    switch(udata->pscale) {
        case AMICI_SCALING_LOG10:
            coefficient = log(10.0);
            for(int ip = 0; ip < udata->nplist; ++ip)
                pcoefficient[ip] = udata->p[udata->plist[ip]]*log(10);
            if (udata->sensi == 2)
                if (udata->o2mode == AMICI_O2MODE_FULL)
                    for(int ip = 0; ip < udata->np; ++ip)
                        augcoefficient[ip] = udata->p[ip]*log(10);
            break;
        case AMICI_SCALING_LN:
            coefficient = 1.0;
            for(int ip = 0; ip < udata->nplist; ++ip)
                pcoefficient[ip] = udata->p[udata->plist[ip]];
            if (udata->sensi == 2)
                if (udata->o2mode == AMICI_O2MODE_FULL)
                    for(int ip = 0; ip < udata->np; ++ip)
                        augcoefficient[ip] = udata->p[ip];
            break;
        case AMICI_SCALING_NONE:
            //this should never be reached
            break;
    }

    if (udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        // recover first order sensitivies from states for adjoint sensitivity analysis
        if (udata->sensi == AMICI_SENSI_ORDER_SECOND){
            if (udata->sensi_meth == AMICI_SENSI_ASA){
                if (x)
                    if (sx)
                        for(int ip = 0; ip < udata->nplist; ++ip)
                            for(int ix = 0; ix < udata->nxtrue; ++ix)
                                for(int it = 0; it < udata->nt; ++it)
                                    sx[(ip*udata->nxtrue + ix)*udata->nt + it] = x[(udata->nxtrue + ip*udata->nxtrue + ix)*udata->nt + it];

                if (y)
                    if (sy)
                        for(int ip = 0; ip < udata->nplist; ++ip)
                            for(int iy = 0; iy < udata->nytrue; ++iy)
                                for(int it = 0; it < udata->nt; ++it)
                                    sy[(ip*udata->nytrue + iy)*udata->nt + it] = y[(udata->nytrue + ip*udata->nytrue + iy)*udata->nt + it];

                if (z)
                    if (sz)
                        for(int ip = 0; ip < udata->nplist; ++ip)
                            for(int iz = 0; iz < udata->nztrue; ++iz)
                                for(int it = 0; it < udata->nt; ++it)
                                    sz[(ip * udata->nztrue + iz)*udata->nt + it] = z[(udata->nztrue + ip*udata->nztrue + iz)*udata->nt + it];

            }
        }

        if (edata) {
            if (sllh)
                for(int ip = 0; ip < udata->nplist; ++ip)
                    sllh[ip] *= pcoefficient[ip];
        }

#define chainRule(QUANT,IND1,N1T,N1,IND2,N2) \
if (s ## QUANT ) \
for(int ip = 0; ip < udata->nplist; ++ip) \
for(int IND1 = 0; IND1 < N1T; ++IND1) \
for(int IND2 = 0; IND2 < N2; ++IND2){ \
s ## QUANT [(ip * N1 + IND1) * N2 + IND2] *= pcoefficient[ip];} \

        chainRule(x,ix,udata->nxtrue,udata->nx,it,udata->nt)
        chainRule(y,iy,udata->nytrue,udata->ny,it,udata->nt)
        chainRule(sigmay,iy,udata->nytrue,udata->ny,it,udata->nt)
        chainRule(z,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
        chainRule(sigmaz,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
        chainRule(rz,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
    }
    if (udata->sensi_meth == AMICI_SENSI_SS) {
        if (dxdotdp)
            for(int ip = 0; ip < udata->nplist; ++ip)
                for(int ix = 0; ix < udata->nx; ++ix)
                    dxdotdp[ix + ip*udata->nxtrue] *= pcoefficient[ip];

        if (dydp)
            for(int ip = 0; ip < udata->nplist; ++ip)
                for(int iy = 0; iy < udata->ny; ++iy)
                    dydp[iy + ip*udata->nytrue] *= pcoefficient[ip];
    }
    if (udata->o2mode == AMICI_O2MODE_FULL) { //full
        if (edata){
            if (s2llh) {
                if (sllh) {
                    for(int ip = 0; ip < udata->nplist; ++ip) {
                        for(int iJ = 1; iJ < udata->nJ; ++iJ) {
                            s2llh[ip*udata->nplist+(iJ-1)] *= pcoefficient[ip]*augcoefficient[iJ-1];
                            if (udata->plist[ip] == iJ-1)
                                s2llh[ip*udata->nplist+(iJ-1)] += sllh[ip]*coefficient;
                        }
                    }
                }
            }
        }

#define s2ChainRule(QUANT,IND1,N1T,N1,IND2,N2) \
if (s ## QUANT ) \
for(int ip = 0; ip < udata->nplist; ++ip) \
for(int iJ = 1; iJ < udata->nJ; ++iJ) \
for(int IND1 = 0; IND1 < N1T; ++IND1) \
for(int IND2 = 0; IND2 < N2; ++IND2){ \
s ## QUANT [(ip*N1 + iJ*N1T + IND1)*N2 + IND2] *= pcoefficient[ip]*augcoefficient[iJ-1]; \
if (udata->plist[ip]==iJ-1) \
s  ## QUANT [(ip*N1 + iJ*N1T + IND1)*N2 + IND2] += s ## QUANT [(ip*N1 + IND1)*N2 + IND2]*coefficient;}

        s2ChainRule(x,ix,udata->nxtrue,udata->nx,it,udata->nt)
        s2ChainRule(y,iy,udata->nytrue,udata->ny,it,udata->nt)
        s2ChainRule(sigmay,iy,udata->nytrue,udata->ny,it,udata->nt)
        s2ChainRule(z,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
        s2ChainRule(sigmaz,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
        s2ChainRule(rz,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
    }

    if (udata->o2mode == AMICI_O2MODE_DIR) { //directional
        if (s2llh) {
            if (sllh) {
                for(int ip = 0; ip < udata->nplist; ++ip) {
                    s2llh[ip] *= pcoefficient[ip];
                    s2llh[ip] += udata->k[udata->nk-udata->nplist+ip]*sllh[ip]/udata->p[udata->plist[ip]];
                }
            }
        }

#define s2vecChainRule(QUANT,IND1,N1T,N1,IND2,N2) \
if (s ## QUANT ) \
for(int ip = 0; ip < udata->nplist; ++ip) \
for(int IND1 = 0; IND1 < N1T; ++IND1) \
for(int IND2 = 0; IND2 < N2; ++IND2){ \
s ## QUANT [(ip*N1 + N1T + IND1)*N2 + IND2] *= pcoefficient[ip]; \
s ## QUANT [(ip*N1 + N1T + IND1)*N2 + IND2] += udata->k[udata->nk-udata->nplist+ip]*s ## QUANT [(ip*N1 + IND1)*N2 + IND2]/udata->p[udata->plist[ip]];}

        s2vecChainRule(x,ix,udata->nxtrue,udata->nx,it,udata->nt)
        s2vecChainRule(y,iy,udata->nytrue,udata->ny,it,udata->nt)
        s2vecChainRule(sigmay,iy,udata->nytrue,udata->ny,it,udata->nt)
        s2vecChainRule(z,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
        s2vecChainRule(sigmaz,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
        s2vecChainRule(rz,iz,udata->nztrue,udata->nz,ie,udata->nmaxevent)
    }

    delete[] pcoefficient;
    delete[] augcoefficient;
    return AMICI_SUCCESS;
}

ReturnData::~ReturnData()
{
    if(!freeFieldsOnDestruction)
        return;

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
    if(numerrtestfails) delete [] numerrtestfails;
    if(numnonlinsolvconvfails) delete [] numnonlinsolvconvfails;
    if(order) delete[] order;
    if(numstepsB) delete[] numstepsB;
    if(numrhsevalsB) delete[] numrhsevalsB;
    if(numerrtestfailsB) delete [] numerrtestfailsB;
    if(numnonlinsolvconvfailsB) delete [] numnonlinsolvconvfailsB;
    if(xss) delete[] xss;
    if(newton) delete[] newton;
    if(newton_numsteps) delete[] newton_numsteps;
    if(newton_numlinsteps) delete[] newton_numlinsteps;
    if(newton_time) delete[] newton_time;
    if(llh) delete[] llh;
    if(sllh) delete[] sllh;
    if(s2llh) delete[] s2llh;
    if(chi2) delete[] chi2;
    if(status) delete[] status;
}

void ReturnData::initFields(const UserData *udata)
{
    initField1(&status, "status", 1);
    if(udata) {
        initField1(&ts, "t", udata->nt);
        initField1(&llh, "llh", 1);
        initField1(&chi2, "chi2", 1);
        initField2(&numsteps, "numsteps",udata->nt, 1);
        initField2(&numrhsevals, "numrhsevals", udata->nt, 1);
        initField2(&numerrtestfails, "numerrtestfails",udata->nt, 1);
        initField2(&numnonlinsolvconvfails, "numnonlinsolvconvfails", udata->nt, 1);
        initField2(&order, "order", udata->nt,1);
        initField2(&numsteps, "numsteps",udata->nt, 1);
        
        if((udata->nz>0) & (udata->ne>0)){
            initField2(&z, "z", udata->nmaxevent,udata->nz);
            initField2(&rz, "rz", udata->nmaxevent,udata->nz);
            initField2(&sigmaz, "sigmaz", udata->nmaxevent,udata->nz);
        }
        if(udata->nx>0) {
            initField2(&x, "x", udata->nt,udata->nx);
            initField2(&xdot, "xdot", 1,udata->nx);
            initField2(&J, "J", udata->nx,udata->nx);
            initField2(&x, "xss", 1,udata->nx);
            initField2(&newton, "newton", 1,1);
            initField2(&newton_numsteps, "newton_numsteps", 1,udata->nx);
            initField2(&newton_numlinsteps, "newton_numlinsteps", udata->newton_maxsteps,2);
            initField2(&newton_time, "newton_time", 1,2);
        }
        if(udata->ny>0) {
            initField2(&y, "y", udata->nt,udata->ny);
            initField2(&sigmay, "sigmay", udata->nt,udata->ny);
            if (udata->sensi_meth == AMICI_SENSI_SS) {
                initField2(&dydp, "dydp", udata->ny,udata->nplist);
                initField2(&dydx, "dydx", udata->ny,udata->nx);
                initField2(&dxdotdp, "dxdotdp", udata->nx,udata->nplist);
            }
        }
        if(udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
            initField2(&sllh, "sllh", udata->nplist,1);
            
            if (udata->sensi_meth == AMICI_SENSI_FSA) {
                initField3(&sx, "sx", udata->nt,udata->nx,udata->nplist);
                if(udata->ny>0) {
                    initField3(&sy, "sy", udata->nt,udata->ny,udata->nplist);
                    initField3(&ssigmay, "ssigmay", udata->nt,udata->ny,udata->nplist);
                }
                if((udata->nz>0) & (udata->ne>0)){
                    initField3(&srz, "srz", udata->nmaxevent,udata->nz,udata->nplist);
                    if(udata->sensi >= AMICI_SENSI_ORDER_SECOND){
                        initField4(&s2rz, "s2rz", udata->nmaxevent,udata->nztrue,udata->nplist,udata->nplist);
                    }
                    initField3(&sz, "sz", udata->nmaxevent,udata->nz,udata->nplist);
                    initField3(&ssigmaz, "ssigmaz", udata->nmaxevent,udata->nz,udata->nplist);
                }
            }
            
            if (udata->sensi_meth == AMICI_SENSI_ASA) {
                if(udata->ny>0) {
                    initField3(&ssigmay, "ssigmay", udata->nt,udata->ny,udata->nplist);
                }
                if((udata->nz>0) & (udata->ne>0)){
                    initField3(&ssigmaz, "ssigmaz", udata->nmaxevent,udata->nz,udata->nplist);
                }
                initField2(&numstepsB, "numstepsB", udata->nt, 1);
                initField2(&numrhsevalsB, "numrhsevalsB", udata->nt, 1);
                initField2(&numerrtestfailsB, "numerrtestfailsB",udata->nt, 1);
                initField2(&numnonlinsolvconvfailsB, "numnonlinsolvconvfailsB", udata->nt, 1);
            }
            
            if(udata->sensi >= AMICI_SENSI_ORDER_SECOND) {
                initField2(&s2llh, "s2llh", udata->nJ-1,udata->nplist);
            }
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
