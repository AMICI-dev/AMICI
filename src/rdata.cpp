#include "include/rdata.h"
#include "include/udata.h"
#include "include/amici_misc.h"
#include "include/amici_model.h"
#include "include/symbolic_functions.h"
#include <cstring>

ReturnData::ReturnData() : np(0), nk(0), nx(0), nxtrue(0), ny(0), nytrue(0), nz(0), nztrue(0),
    ne(0), nJ(0), nplist(0), nmaxevent(0), nt(0), newton_maxsteps(0), pscale(AMICI_SCALING_NONE),
    o2mode(AMICI_O2MODE_NONE), sensi(AMICI_SENSI_ORDER_NONE), sensi_meth(AMICI_SENSI_NONE)
{
    setDefaults();
}

ReturnData::ReturnData(const UserData *udata, const Model *model) :
    np(model->np), nk(model->nk), nx(model->nx), nxtrue(model->nxtrue), ny(model->ny), nytrue(model->nytrue),
    nz(model->nz), nztrue(model->nztrue), ne(model->ne), nJ(model->nJ),
    nplist(udata->nplist), nmaxevent(udata->nmaxevent), nt(udata->nt), newton_maxsteps(udata->newton_maxsteps),
    pscale(model->pscale), o2mode(model->o2mode),
    sensi(udata->sensi), sensi_meth(udata->sensi_meth)
{
    setDefaults();

    initFields();

    memcpy(ts, udata->ts, nt * sizeof(realtype));
}

void ReturnData::setDefaults()
{
    ts = xdot = dxdotdp = dydx = dydp = J = z = sigmaz = sz = ssigmaz = rz = srz = s2rz = x = sx = y = sigmay = nullptr;
    sy = ssigmay = numsteps = numstepsB = numrhsevals = numrhsevalsB = order = llh = chi2 = sllh = s2llh = nullptr;
    numerrtestfails = numnonlinsolvconvfails = numerrtestfailsB = numnonlinsolvconvfailsB = nullptr;
    newton_status = newton_time = newton_numsteps = newton_numlinsteps = xss = nullptr;
    status = nullptr;

    freeFieldsOnDestruction = true;
}

void ReturnData::invalidate()
{
    if (llh)
        *llh = amiGetNaN();

    if (sllh)
        setLikelihoodSensitivityFirstOrderNaN();

    if (s2llh)
        setLikelihoodSensitivitySecondOrderNaN();

}

void ReturnData::setLikelihoodSensitivityFirstOrderNaN()
{
    fillArray(sllh, nplist, amiGetNaN());
}

void ReturnData::setLikelihoodSensitivitySecondOrderNaN()
{
    fillArray(s2llh, nplist*(nJ-1), amiGetNaN());
}

int ReturnData::applyChainRuleFactorToSimulationResults(UserData *udata)
{
    if (pscale == AMICI_SCALING_NONE)
        return AMICI_SUCCESS;

    // chain-rule factor: multiplier for am_p
    realtype coefficient;
    realtype *pcoefficient, *augcoefficient;

    pcoefficient = new realtype[nplist]();
    augcoefficient = new realtype[np]();

    switch(pscale) {
        case AMICI_SCALING_LOG10:
            coefficient = log(10.0);
            for(int ip = 0; ip < nplist; ++ip)
                pcoefficient[ip] = udata->p[udata->plist[ip]]*log(10);
            if (udata->sensi == 2)
                if (o2mode == AMICI_O2MODE_FULL)
                    for(int ip = 0; ip < np; ++ip)
                        augcoefficient[ip] = udata->p[ip]*log(10);
            break;
        case AMICI_SCALING_LN:
            coefficient = 1.0;
            for(int ip = 0; ip < nplist; ++ip)
                pcoefficient[ip] = udata->p[udata->plist[ip]];
            if (udata->sensi == 2)
                if (o2mode == AMICI_O2MODE_FULL)
                    for(int ip = 0; ip < np; ++ip)
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
                        for(int ip = 0; ip < nplist; ++ip)
                            for(int ix = 0; ix < nxtrue; ++ix)
                                for(int it = 0; it < nt; ++it)
                                    sx[(ip*nxtrue + ix)*nt + it] = x[(nxtrue + ip*nxtrue + ix)*nt + it];

                if (y)
                    if (sy)
                        for(int ip = 0; ip < nplist; ++ip)
                            for(int iy = 0; iy < nytrue; ++iy)
                                for(int it = 0; it < nt; ++it)
                                    sy[(ip*nytrue + iy)*nt + it] = y[(nytrue + ip*nytrue + iy)*nt + it];

                if (z)
                    if (sz)
                        for(int ip = 0; ip < nplist; ++ip)
                            for(int iz = 0; iz < nztrue; ++iz)
                                for(int it = 0; it < nt; ++it)
                                    sz[(ip * nztrue + iz)*nt + it] = z[(nztrue + ip*nztrue + iz)*nt + it];

            }
        }

        if (sllh)
            for(int ip = 0; ip < nplist; ++ip)
                sllh[ip] *= pcoefficient[ip];

#define chainRule(QUANT,IND1,N1T,N1,IND2,N2) \
if (s ## QUANT ) \
for(int ip = 0; ip < nplist; ++ip) \
for(int IND1 = 0; IND1 < N1T; ++IND1) \
for(int IND2 = 0; IND2 < N2; ++IND2){ \
s ## QUANT [(ip * N1 + IND1) * N2 + IND2] *= pcoefficient[ip];} \

        chainRule(x,ix,nxtrue,nx,it,nt)
        chainRule(y,iy,nytrue,ny,it,nt)
        chainRule(sigmay,iy,nytrue,ny,it,nt)
        chainRule(z,iz,nztrue,nz,ie,nmaxevent)
        chainRule(sigmaz,iz,nztrue,nz,ie,nmaxevent)
        chainRule(rz,iz,nztrue,nz,ie,nmaxevent)
    }
    if (udata->sensi_meth == AMICI_SENSI_SS) {
        if (dxdotdp)
            for(int ip = 0; ip < nplist; ++ip)
                for(int ix = 0; ix < nx; ++ix)
                    dxdotdp[ix + ip*nxtrue] *= pcoefficient[ip];

        if (dydp)
            for(int ip = 0; ip < nplist; ++ip)
                for(int iy = 0; iy < ny; ++iy)
                    dydp[iy + ip*nytrue] *= pcoefficient[ip];
    }
    if (o2mode == AMICI_O2MODE_FULL) { //full
            if (s2llh) {
                if (sllh) {
                    for(int ip = 0; ip < nplist; ++ip) {
                        for(int iJ = 1; iJ < nJ; ++iJ) {
                            s2llh[ip*nplist+(iJ-1)] *= pcoefficient[ip]*augcoefficient[iJ-1];
                            if (udata->plist[ip] == iJ-1)
                                s2llh[ip*nplist+(iJ-1)] += sllh[ip]*coefficient;
                        }
                    }
                }
            }

#define s2ChainRule(QUANT,IND1,N1T,N1,IND2,N2) \
if (s ## QUANT ) \
for(int ip = 0; ip < nplist; ++ip) \
for(int iJ = 1; iJ < nJ; ++iJ) \
for(int IND1 = 0; IND1 < N1T; ++IND1) \
for(int IND2 = 0; IND2 < N2; ++IND2){ \
s ## QUANT [(ip*N1 + iJ*N1T + IND1)*N2 + IND2] *= pcoefficient[ip]*augcoefficient[iJ-1]; \
if (udata->plist[ip]==iJ-1) \
s  ## QUANT [(ip*N1 + iJ*N1T + IND1)*N2 + IND2] += s ## QUANT [(ip*N1 + IND1)*N2 + IND2]*coefficient;}

        s2ChainRule(x,ix,nxtrue,nx,it,nt)
        s2ChainRule(y,iy,nytrue,ny,it,nt)
        s2ChainRule(sigmay,iy,nytrue,ny,it,nt)
        s2ChainRule(z,iz,nztrue,nz,ie,nmaxevent)
        s2ChainRule(sigmaz,iz,nztrue,nz,ie,nmaxevent)
        s2ChainRule(rz,iz,nztrue,nz,ie,nmaxevent)
    }

    if (o2mode == AMICI_O2MODE_DIR) { //directional
        if (s2llh) {
            if (sllh) {
                for(int ip = 0; ip < nplist; ++ip) {
                    s2llh[ip] *= pcoefficient[ip];
                    s2llh[ip] += udata->k[nk-nplist+ip]*sllh[ip]/udata->p[udata->plist[ip]];
                }
            }
        }

#define s2vecChainRule(QUANT,IND1,N1T,N1,IND2,N2) \
if (s ## QUANT ) \
for(int ip = 0; ip < nplist; ++ip) \
for(int IND1 = 0; IND1 < N1T; ++IND1) \
for(int IND2 = 0; IND2 < N2; ++IND2){ \
s ## QUANT [(ip*N1 + N1T + IND1)*N2 + IND2] *= pcoefficient[ip]; \
s ## QUANT [(ip*N1 + N1T + IND1)*N2 + IND2] += udata->k[nk-nplist+ip]*s ## QUANT [(ip*N1 + IND1)*N2 + IND2]/udata->p[udata->plist[ip]];}

        s2vecChainRule(x,ix,nxtrue,nx,it,nt)
        s2vecChainRule(y,iy,nytrue,ny,it,nt)
        s2vecChainRule(sigmay,iy,nytrue,ny,it,nt)
        s2vecChainRule(z,iz,nztrue,nz,ie,nmaxevent)
        s2vecChainRule(sigmaz,iz,nztrue,nz,ie,nmaxevent)
        s2vecChainRule(rz,iz,nztrue,nz,ie,nmaxevent)
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
    if(newton_status) delete[] newton_status;
    if(newton_numsteps) delete[] newton_numsteps;
    if(newton_numlinsteps) delete[] newton_numlinsteps;
    if(newton_time) delete[] newton_time;
    if(llh) delete[] llh;
    if(sllh) delete[] sllh;
    if(s2llh) delete[] s2llh;
    if(chi2) delete[] chi2;
    if(status) delete[] status;
}

void ReturnData::initFields()
{
    initField1(&status, "status", 1);

    initField1(&ts, "t", nt);
    initField1(&llh, "llh", 1);
    initField1(&chi2, "chi2", 1);
    initField2(&numsteps, "numsteps",nt, 1);
    initField2(&numrhsevals, "numrhsevals", nt, 1);
    initField2(&numerrtestfails, "numerrtestfails",nt, 1);
    initField2(&numnonlinsolvconvfails, "numnonlinsolvconvfails", nt, 1);
    initField2(&order, "order", nt,1);

    if((nz>0) & (ne>0)){
        initField2(&z, "z", nmaxevent,nz);
        initField2(&rz, "rz", nmaxevent,nz);
        initField2(&sigmaz, "sigmaz", nmaxevent,nz);
    }
    if(nx>0) {
        initField2(&x, "x", nt,nx);
        initField2(&xdot, "xdot", 1,nx);
        initField2(&J, "J", nx,nx);
        initField2(&xss, "xss", 1,nx);
        initField2(&newton_status, "newton_status", 1,1);
        initField2(&newton_numsteps, "newton_numsteps", 1,2);
        initField2(&newton_numlinsteps, "newton_numlinsteps", newton_maxsteps,2);
        initField2(&newton_time, "newton_time", 1,2);
    }
    if(ny>0) {
        initField2(&y, "y", nt,ny);
        initField2(&sigmay, "sigmay", nt,ny);
        if (sensi_meth == AMICI_SENSI_SS) {
            initField2(&dydp, "dydp", ny,nplist);
            initField2(&dydx, "dydx", ny,nx);
            initField2(&dxdotdp, "dxdotdp", nx,nplist);
        }
    }
    if(sensi >= AMICI_SENSI_ORDER_FIRST) {
        initField2(&sllh, "sllh", nplist,1);

        if (sensi_meth == AMICI_SENSI_FSA) {
            initField3(&sx, "sx", nt,nx,nplist);
            if(ny>0) {
                initField3(&sy, "sy", nt,ny,nplist);
                initField3(&ssigmay, "ssigmay", nt,ny,nplist);
            }
            if((nz>0) & (ne>0)){
                initField3(&srz, "srz", nmaxevent,nz,nplist);
                if(sensi >= AMICI_SENSI_ORDER_SECOND){
                    initField4(&s2rz, "s2rz", nmaxevent,nztrue,nplist,nplist);
                }
                initField3(&sz, "sz", nmaxevent,nz,nplist);
                initField3(&ssigmaz, "ssigmaz", nmaxevent,nz,nplist);
            }
        }

        if (sensi_meth == AMICI_SENSI_ASA) {
            if(ny>0) {
                initField3(&ssigmay, "ssigmay", nt,ny,nplist);
            }
            if((nz>0) & (ne>0)){
                initField3(&ssigmaz, "ssigmaz", nmaxevent,nz,nplist);
            }
            initField2(&numstepsB, "numstepsB", nt, 1);
            initField2(&numrhsevalsB, "numrhsevalsB", nt, 1);
            initField2(&numerrtestfailsB, "numerrtestfailsB",nt, 1);
            initField2(&numnonlinsolvconvfailsB, "numnonlinsolvconvfailsB", nt, 1);
        }

        if(sensi >= AMICI_SENSI_ORDER_SECOND) {
            initField2(&s2llh, "s2llh", nJ-1,nplist);
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
