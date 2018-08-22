#include "amici/model.h"
#include "amici/amici.h"

#include <cstring>
#include <cmath>
#include <typeinfo>
#include <utility>
#include <algorithm>

namespace amici {

/** Sensitivity of measurements y, total derivative sy = dydx * sx + dydp
 * @param it timepoint index
 * @param rdata pointer to return data instance
 */
void Model::fsy(const int it, ReturnData *rdata) {
    if (!ny)
        return;
    
    // copy dydp for current time to sy
    std::copy(dydp.begin(), dydp.end(), &rdata->sy[it * nplist() * ny]);

    // compute sy = 1.0*dydx*sx + 1.0*sy
    // dydx A[ny,nx] * sx B[nx,nplist] = sy C[ny,nplist]
    //        M  K          K  N              M  N
    //        lda           ldb               ldc
    amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans, ny, nplist(), nx,
                1.0, dydx.data(), ny, getsx(it,rdata), nx, 1.0,
                &rdata->sy[it*nplist()*ny], ny);
}

/** Sensitivity of z at final timepoint (ignores sensitivity of timepoint),
 * total derivative
 * @param nroots number of events for event index
 * @param ie event index
 * @param rdata pointer to return data instance
 */
void Model::fsz_tf(const int *nroots, const int ie, ReturnData *rdata) {
    
    for (int iz = 0; iz < nz; ++iz)
        if (z2event[iz] - 1 == ie)
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->sz.at((nroots[ie]*nplist()+ip)*nz + iz) = 0.0;
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy, total
 * derivative
 * @param it timepoint index
 * @param dJydx vector with values of state derivative of Jy
 * @param rdata pointer to return data instance
 */
void Model::fsJy(const int it, const std::vector<realtype>& dJydx, ReturnData *rdata) {

    // Compute dJydx*sx for current 'it'
    // dJydx        rdata->nt x nJ x nx
    // sx           rdata->nt x nx x nplist()
    std::vector<realtype> multResult(nJ * nplist(), 0);

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans, nJ,
                nplist(), nx, 1.0, &dJydx.at(it*nJ*nx), nJ, getsx(it,rdata), nx, 0.0,
                multResult.data(), nJ);

    // multResult    nJ x nplist()
    // dJydp         nJ x nplist()
    // dJydxTmp      nJ x nx
    // sxTmp         nx x nplist()

    // sJy += multResult + dJydp
    for (int iJ = 0; iJ < nJ; ++iJ) {
        if (iJ == 0)
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->sllh.at(ip) -= multResult.at(ip * nJ) + dJydp.at(ip * nJ);
        else
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->s2llh.at((iJ - 1) + ip * (nJ - 1)) -=
                        multResult.at(iJ + ip * nJ) + dJydp.at(iJ + ip * nJ);
    }
}

/** Compute sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t.
 * parameters for the given timepoint. Add result to respective fields in rdata.
 * @param it timepoint index
 * @param edata pointer to experimental data instance
 * @param rdata pointer to return data instance
 */
void Model::fdJydp(const int it, const ExpData *edata,
                   ReturnData *rdata) {

    // dJydy         nJ x nytrue x ny
    // dydp          nplist * ny
    // dJydp         nplist x nJ
    // dJydsigma

    getmy(it,edata);
    std::fill(dJydp.begin(), dJydp.end(), 0.0);

    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (isNaN(my.at(iyt)))
            continue;

        // dJydp = 1.0 * dJydp +  1.0 * dJydy * dydp
        amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
                    nJ, nplist(), ny,
                    1.0, &dJydy.at(iyt*nJ*ny), nJ,
                    dydp.data(), ny,
                    1.0, dJydp.data(), nJ);

        // dJydp = 1.0 * dJydp +  1.0 * dJydsigma * dsigmaydp
        amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
                    nJ, nplist(), ny,
                    1.0, &dJydsigma.at(iyt*nJ*ny), nJ,
                    dsigmaydp.data(), ny,
                    1.0, dJydp.data(), nJ);
    }

    if (rdata->sensi_meth != SensitivityMethod::adjoint)
        return;

    if(!ny)
        return;

    for (int iJ = 0; iJ < nJ; iJ++) {
        for (int ip = 0; ip < nplist(); ip++) {
            if (iJ == 0) {
                rdata->sllh.at(ip) -= dJydp[ip * nJ];
            } else {
                rdata->s2llh.at((iJ - 1) + ip * (nJ - 1)) -=
                        dJydp[iJ + ip * nJ];
            }
        }
    }
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t.
 * state variables
 * @param dJydx pointer to vector with values of state derivative of Jy
 * @param it timepoint index
 * @param edata pointer to experimental data instance
 * @param rdata pointer to return data instance
 */
void Model::fdJydx(std::vector<realtype> *dJydx, const int it, const ExpData *edata, const ReturnData *rdata) {

    // dJydy         nJ x ny x nytrue
    // dydx          ny x nx
    // dJydx         nJ x nx x nt
    getmy(it,edata);
    for (int iyt = 0; iyt < nytrue; ++iyt) {
        if (isNaN(my.at(iyt)))
            continue;
    // dJydy A[nyt,nJ,ny] * dydx B[ny,nx] = dJydx C[it,nJ,nx]
    //         slice                                slice
    //             M  K            K  N                M  N
    //             lda             ldb                 ldc
        amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
                    nJ, nx, ny, 1.0, &dJydy.at(iyt*ny*nJ), nJ, dydx.data(), ny, 1.0,
                    &dJydx->at(it*nx*nJ), nJ);
    }
}

/** Sensitivity of event-resolved measurement negative log-likelihood Jz, total
 * derivative
 * @param nroots event index
 * @param dJzdx vector with values of state derivative of Jz
 * @param sx pointer to state sensitivities
 * @param rdata pointer to return data instance
 */
void Model::fsJz(const int nroots, const std::vector<realtype>& dJzdx, AmiVectorArray *sx, ReturnData *rdata) {
    // sJz           nJ x nplist()
    // dJzdp         nJ x nplist()
    // dJzdx         nmaxevent x nJ x nx
    // sx            rdata->nt x nx x nplist()

    // Compute dJzdx*sx for current 'ie'
    // dJzdx        rdata->nt x nJ x nx
    // sx           rdata->nt x nx x nplist()

    std::vector<realtype> multResult(nJ * nplist(), 0);
    std::vector<realtype> sxTmp(rdata->nplist * nx, 0);
    for (int ip = 0; ip < rdata->nplist; ++ip) {
        for (int ix = 0; ix < nx; ++ix)
            sxTmp.at(ix + ip * nx) = sx->at(ix,ip);
    }

    // C := alpha*op(A)*op(B) + beta*C,
    amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans, nJ,
                nplist(), nx, 1.0, &dJzdx.at(nroots*nx*nJ), nJ, sxTmp.data(), nx, 1.0,
                multResult.data(), nJ);

    // sJy += multResult + dJydp
    for (int iJ = 0; iJ < nJ; ++iJ) {
        if (iJ == 0)
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->sllh.at(ip) -= multResult.at(ip * nJ) + dJzdp.at(ip * nJ);
        else
            for (int ip = 0; ip < nplist(); ++ip)
                rdata->s2llh.at((iJ - 1) + ip * (nJ - 1)) -=
                        multResult.at(iJ + ip * nJ) + dJzdp.at(iJ + ip * nJ);
    }
}

/** Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t.
 * parameters
 * @param nroots event index
 * @param t current timepoint
 * @param edata pointer to experimental data instance
 * @param rdata pointer to return data instance
 */
void Model::fdJzdp(const int nroots, realtype t, const ExpData *edata,
                   const ReturnData *rdata) {
    // dJzdz         nJ x nz x nztrue
    // dJzdsigma     nJ x nz x nztrue
    // dzdp          nz x nplist()
    // dJzdp         nJ x nplist()

    getmz(nroots,edata);
    std::fill(dJzdp.begin(),dJzdp.end(),0.0);
    for (int izt = 0; izt < nztrue; ++izt) {
        if (isNaN(mz.at(izt)))
            continue;

        if (t < rdata->ts.at(rdata->ts.size() - 1)) {
            // with z
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
                        nJ, nplist(), nz, 1.0, &dJzdz.at(izt*nz*nJ), nJ, dzdp.data(), nz,
                        1.0, dJzdp.data(), nJ);
        } else {
            // with rz
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nplist(), nz, 1.0,
                        &dJrzdsigma.at(izt*nz*nJ), nJ, dsigmazdp.data(), nz, 1.0,
                        dJzdp.data(), nJ);
            
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
                        nJ, nplist(), nz, 1.0, &dJrzdz.at(izt*nz*nJ), nJ, dzdp.data(), nz,
                        1.0, dJzdp.data(), nJ);
        }

        amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans, BLASTranspose::noTrans,
                    nJ, nplist(), nz, 1.0, &dJzdsigma.at(izt*nz*nJ), nJ,
                    dsigmazdp.data(), nz, 1.0, dJzdp.data(), nJ);
    }
}

/** Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t.
 * state variables
 * @param dJzdx pointer to vector with values of state derivative of Jz
 * @param nroots event index
 * @param t current timepoint
 * @param edata pointer to experimental data instance
 * @param rdata pointer to return data instance
 */
void Model::fdJzdx(std::vector<realtype> *dJzdx, const int nroots, realtype t, const ExpData *edata, const ReturnData *rdata) {
    // dJzdz         nJ x nz x nztrue
    // dzdx          nz x nx
    // dJzdx         nJ x nx x nmaxevent
    getmz(nroots,edata);
    for (int izt = 0; izt < nztrue; ++izt) {
        if (isNaN(mz.at(izt)))
            continue;
        
        if (t < rdata->ts.at(rdata->ts.size() - 1)) {
            // z
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nx, nz, 1.0, &dJzdz.at(izt*nz*nJ), nJ,
                        dzdx.data(), nz, 1.0, &dJzdx->at(nroots*nx*nJ), nJ);
        } else {
            // rz
            amici_dgemm(BLASLayout::colMajor, BLASTranspose::noTrans,
                        BLASTranspose::noTrans, nJ, nx, nz, 1.0, &dJrzdz.at(izt*nz*nJ), nJ,
                        drzdx.data(), nz, 1.0, &dJzdx->at(nroots*nx*nJ), nJ);
        }
    }
}

/** initialization of model properties
 * @param x pointer to state variables
 * @param dx pointer to time derivative of states (DAE only)
 */
void Model::initialize(AmiVector *x, AmiVector *dx) {

    initializeStates(x);
    
    fdx0(x, dx);
    
    if(ne)
        initHeaviside(x,dx);
    
}

/** initialization of initial states
 * @param x pointer to state variables
 */
void Model::initializeStates(AmiVector *x) {

    if (x0data.empty()) {
        fx0(x);
    } else {
        for (int ix = 0; ix < nx; ix++) {
            (*x)[ix] = (realtype) x0data.at(ix);
        }
    }
}

/**
 * initHeaviside initialises the heaviside variables h at the intial time t0
 * heaviside variables activate/deactivate on event occurences
 * @param x pointer to state variables
 * @param dx pointer to time derivative of states (DAE only)
 */
void Model::initHeaviside(AmiVector *x, AmiVector *dx) {
    std::vector<realtype> rootvals(ne,0.0);
    froot(tstart, x, dx, rootvals.data());
    for (int ie = 0; ie < ne; ie++) {
        if (rootvals.at(ie) < 0) {
            h.at(ie) = 0.0;
        } else if (rootvals.at(ie) == 0) {
            throw AmiException("Simulation started in an event. This could lead to "
                               "unexpected results, aborting simulation! Please "
                               "specify an earlier simulation start via "
                               "options.t0");
        } else {
            h.at(ie) = 1.0;
        }
    }
}

int Model::nplist() const {
    return plist_.size();
}

int Model::np() const {
    return originalParameters.size();
}

int Model::nk() const {
    return fixedParameters.size();
}

const double *Model::k() const{
    return fixedParameters.data();
}

int Model::nMaxEvent() const {
    return nmaxevent;
}

void Model::setNMaxEvent(int nmaxevent) {
    this->nmaxevent = nmaxevent;
}

int Model::nt() const {
    return ts.size();
}

const std::vector<ParameterScaling> &Model::getParameterScale() const {
    return pscale;
}

void Model::setParameterScale(ParameterScaling pscale) {
    this->pscale.assign(this->pscale.size(), pscale);
    unscaledParameters.resize(originalParameters.size());
    unscaleParameters(unscaledParameters.data());
}

void Model::setParameterScale(std::vector<ParameterScaling> const& pscale) {
    this->pscale = pscale;
    unscaledParameters.resize(originalParameters.size());
    unscaleParameters(unscaledParameters.data());
}

std::vector<realtype> const& Model::getParameters() const {
    return originalParameters;
}

void Model::setParameters(const std::vector<realtype> &p) {
    if(p.size() != (unsigned) this->originalParameters.size())
        throw AmiException("Dimension mismatch. Size of parameters does not match number of model parameters.");
    this->originalParameters = p;
    this->unscaledParameters.resize(originalParameters.size());
    unscaleParameters(this->unscaledParameters.data());
}

const std::vector<realtype> &Model::getUnscaledParameters() const {
    return unscaledParameters;
}

const std::vector<realtype> &Model::getFixedParameters() const {
    return fixedParameters;
}

void Model::setFixedParameters(const std::vector<realtype> &k) {
    if(k.size() != (unsigned) this->fixedParameters.size())
        throw AmiException("Dimension mismatch. Size of fixedParameters does not match number of fixed model parameters.");
    this->fixedParameters = k;
}

std::vector<realtype> const& Model::getTimepoints() const {
    return ts;
}

void Model::setTimepoints(const std::vector<realtype> &ts) {
    if (!std::is_sorted(ts.begin(), ts.end()))
        throw AmiException("Encountered non-monotonic timepoints, please order timepoints such that they are monotonically increasing!");
    this->ts = std::move(ts);
}

double Model::t(int idx) const {
    return ts.at(idx);
}

const std::vector<int> &Model::getParameterList() const {
    return plist_;
}

void Model::setParameterList(const std::vector<int> &plist) {
    for(auto const& idx: plist)
        if(idx < 0 || idx >= np())
            throw AmiException("Indices in plist must be in [0..np]");
    this->plist_ = plist;

    initializeVectors();
}

std::vector<realtype> const& Model::getInitialStates() const {
    return x0data;
}

void Model::setInitialStates(const std::vector<realtype> &x0) {
    if(x0.size() != (unsigned) nx && x0.size() != 0)
        throw AmiException("Dimension mismatch. Size of x0 does not match number of model states.");
    if (x0.size() == 0)
        this->x0data.clear();
    else
        this->x0data = x0;
}

const std::vector<realtype> &Model::getInitialStateSensitivities() const {
    return sx0data;
}

void Model::setInitialStateSensitivities(const std::vector<realtype> &sx0) {
    if(sx0.size() != (unsigned) nx * nplist() && sx0.size() != 0)
        throw AmiException("Dimension mismatch. Size of sx0 does not match number of model states * number of parameter selected for sensitivities.");
    if (sx0.size() == 0)
        this->sx0data.clear();
    else
        this->sx0data = sx0;
}

double Model::t0() const{
    return tstart;
}

void Model::setT0(double t0) {
    tstart = t0;
}

int Model::plist(int pos) const{
    return plist_.at(pos);
}


Model::Model(const int nx,
             const int nxtrue,
             const int ny,
             const int nytrue,
             const int nz,
             const int nztrue,
             const int ne,
             const int nJ,
             const int nw,
             const int ndwdx,
             const int ndwdp,
             const int nnz,
             const int ubw,
             const int lbw,
             SecondOrderMode o2mode,
             const std::vector<realtype>& p,
             std::vector<realtype> k,
             const std::vector<int>& plist,
             std::vector<realtype> idlist,
             std::vector<int> z2event)
    : nx(nx), nxtrue(nxtrue),
      ny(ny), nytrue(nytrue),
      nz(nz), nztrue(nztrue),
      ne(ne),
      nw(nw),
      ndwdx(ndwdx),
      ndwdp(ndwdp),
      nnz(nnz),
      nJ(nJ),
      ubw(ubw),
      lbw(lbw),
      o2mode(o2mode),
      z2event(std::move(z2event)),
      idlist(std::move(idlist)),
      sigmay(ny, 0.0),
      dsigmaydp(ny*plist.size(), 0.0),
      sigmaz(nz, 0.0),
      dsigmazdp(nz*plist.size(), 0.0),
      dJydp(nJ*plist.size(), 0.0),
      dJzdp(nJ*plist.size(), 0.0),
      deltax(nx, 0.0),
      deltasx(nx*plist.size(), 0.0),
      deltaxB(nx, 0.0),
      deltaqB(nJ*plist.size(), 0.0),
      dxdotdp(nx*plist.size(), 0.0),
      my(nytrue, 0.0),
      mz(nztrue, 0.0),
      dJydy(nJ*nytrue*ny, 0.0),
      dJydsigma(nJ*nytrue*ny, 0.0),
      dJzdz(nJ*nztrue*nz, 0.0),
      dJzdsigma(nJ*nztrue*nz, 0.0),
      dJrzdz(nJ*nztrue*nz, 0.0),
      dJrzdsigma(nJ*nztrue*nz, 0.0),
      dzdx(nz*nx, 0.0),
      dzdp(nz*plist.size(), 0.0),
      drzdx(nz*nx, 0.0),
      drzdp(nz*plist.size(), 0.0),
      dydp(ny*plist.size(), 0.0),
      dydx(ny*nx,0.0),
      w(nw, 0.0),
      dwdx(ndwdx, 0.0),
      dwdp(ndwdp, 0.0),
      M(nx*nx, 0.0),
      stau(plist.size(), 0.0),
      h(ne,0.0),
      unscaledParameters(p),
      originalParameters(p),
      fixedParameters(std::move(k)),
      plist_(plist),
      pscale(std::vector<ParameterScaling>(p.size(), ParameterScaling::none))
{
    J = SparseNewMat(nx, nx, nnz, CSC_MAT);
    requireSensitivitiesForAllParameters();
}

Model::Model(const Model &other)
    : nx(other.nx), nxtrue(other.nxtrue),
      ny(other.ny), nytrue(other.nytrue),
      nz(other.nz), nztrue(other.nztrue),
      ne(other.ne), nw(other.nw),
      ndwdx(other.ndwdx), ndwdp(other.ndwdp),
      nnz(other.nnz), nJ(other.nJ),
      ubw(other.ubw), lbw(other.lbw),
      o2mode(other.o2mode),
      z2event(other.z2event),
      idlist(other.idlist),

      sigmay(other.sigmay),
      dsigmaydp(other.dsigmaydp),
      sigmaz(other.sigmaz),
      dsigmazdp(other.dsigmazdp),
      dJydp(other.dJydp),
      dJzdp(other.dJzdp),
      deltax(other.deltax),
      deltasx(other.deltasx),
      deltaxB(other.deltaxB),
      deltaqB(other.deltaqB),
      dxdotdp(other.dxdotdp),
      my(other.my),
      mz(other.mz),
      dJydy(other.dJydy),
      dJydsigma(other.dJydsigma),
      dJzdz(other.dJzdz),
      dJzdsigma(other.dJzdsigma),
      dJrzdz(other.dJrzdz),
      dJrzdsigma(other.dJrzdsigma),
      dzdx(other.dzdx),
      dzdp(other.dzdp),
      drzdx(other.drzdx),
      drzdp(other.drzdp),
      dydp(other.dydp),
      dydx(other.dydx),
      w(other.w),
      dwdx(other.dwdx),
      dwdp(other.dwdp),
      M(other.M),
      stau(other.stau),
      h(other.h),
      unscaledParameters(other.unscaledParameters),
      originalParameters(other.originalParameters),
      fixedParameters(other.fixedParameters),
      plist_(other.plist_),
      x0data(other.x0data),
      sx0data(other.sx0data),
      ts(other.ts),
      nmaxevent(other.nmaxevent),
      pscale(other.pscale),
      tstart(other.tstart)
{
    J = SparseNewMat(nx, nx, nnz, CSC_MAT);
    SparseCopyMat(other.J, J);
}

Model::~Model() {
    if(J)
        SparseDestroyMat(J);
}

void Model::initializeVectors()
{
    dsigmaydp.resize(ny * nplist(), 0.0);
    dsigmazdp.resize(nz * nplist(), 0.0);
    dJydp.resize(nJ * nplist(), 0.0);
    dJzdp.resize(nJ * nplist(), 0.0);
    deltasx.resize(nx * nplist(), 0.0);
    deltaqB.resize(nJ * nplist(), 0.0);
    dxdotdp.resize(nx * nplist(), 0.0);
    dzdp.resize(nz * nplist(), 0.0);
    drzdp.resize(nz * nplist(), 0.0);
    dydp.resize(ny * nplist(), 0.0);
    stau.resize(nplist(), 0.0);
}


void Model::fx0(AmiVector *x) {
    x->reset();
    fx0(x->data(),tstart, unscaledParameters.data(),fixedParameters.data());
}

void Model::fx0_fixedParameters(AmiVector *x) {
    fx0_fixedParameters(x->data(),tstart, unscaledParameters.data(),fixedParameters.data());
}

void Model::fdx0(AmiVector *x0, AmiVector *dx0) {}

void Model::fsx0(AmiVectorArray *sx, const AmiVector *x) {
    sx->reset();
    for(int ip = 0; (unsigned)ip<plist_.size(); ip++)
        fsx0(sx->data(ip),tstart,x->data(), unscaledParameters.data(),fixedParameters.data(),plist_.at(ip));
}

void Model::fsdx0() {}

/** Sensitivity of event timepoint, total derivative
     * @param t current timepoint
     * @param ie event index
     * @param x pointer to state variables
     * @param sx pointer to state sensitivity variables
     */
void Model::fstau(const realtype t, const int ie, const AmiVector *x, const AmiVectorArray *sx) {
    std::fill(stau.begin(),stau.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++){
        fstau(&stau.at(ip),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),sx->data(ip),plist_.at(ip),ie);
    }
}

/** Observables / measurements
  * @param it timepoint index
  * @param rdata pointer to return data instance
  */
void Model::fy(int it, ReturnData *rdata) {
    if (!ny)
        return;
    
    fy(&rdata->y.at(it*ny),rdata->ts.at(it),getx(it,rdata), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/** partial derivative of observables y w.r.t. model parameters p
     * @param it timepoint index
     * @param rdata pointer to return data instance
     */
void Model::fdydp(const int it, ReturnData *rdata) {
    if (!ny)
        return;
    
    std::fill(dydp.begin(),dydp.end(),0.0);

    for(int ip = 0; (unsigned)ip < plist_.size(); ip++){
        // get dydp slice (ny) for current time and parameter
        fdydp(&dydp.at(ip*ny),
              rdata->ts.at(it),
              getx(it,rdata),
              unscaledParameters.data(),
              fixedParameters.data(),
              h.data(),
              plist_.at(ip));
    }
}

/** partial derivative of observables y w.r.t. state variables x
     * @param it timepoint index
     * @param rdata pointer to return data instance
     */
void Model::fdydx(const int it, ReturnData *rdata) {
    if (!ny)
        return;
    
    std::fill(dydx.begin(),dydx.end(),0.0);
    fdydx(dydx.data(),rdata->ts.at(it),getx(it,rdata), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/** Event-resolved output
     * @param nroots number of events for event index
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param rdata pointer to return data instance
     */
void Model::fz(const int nroots, const int ie, const realtype t, const AmiVector *x, ReturnData *rdata) {
    fz(&rdata->z.at(nroots*nz),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/** Sensitivity of z, total derivative
     * @param nroots number of events for event index
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param sx current state sensitivities
     * @param rdata pointer to return data instance
     */
void Model::fsz(const int nroots, const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata) {
    for(int ip = 0; (unsigned)ip < plist_.size();  ip++ ){
        fsz(&rdata->sz.at((nroots*nplist()+ip)*nz),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),sx->data(ip),plist_.at(ip));
    }
}

/** Event root function of events (equal to froot but does not include
     * non-output events)
     * @param nroots number of events for event index
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param rdata pointer to return data instance
     */
void Model::frz(const int nroots, const int ie, const realtype t, const AmiVector *x, ReturnData *rdata) {
    frz(&rdata->rz.at(nroots*nz),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/** Sensitivity of rz, total derivative
     * @param nroots number of events for event index
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param sx current state sensitivities
     * @param rdata pointer to return data instance
     */
void Model::fsrz(const int nroots, const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata) {
    for(int ip = 0; (unsigned)ip < plist_.size();  ip++ ){
        fsrz(&rdata->srz.at((nroots*nplist()+ip)*nz),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),sx->data(ip),plist_.at(ip));
    }
}

/** partial derivative of event-resolved output z w.r.t. to model parameters p
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
void Model::fdzdp(const realtype t, const int ie, const AmiVector *x) {
    std::fill(dzdp.begin(),dzdp.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++){
        fdzdp(dzdp.data(),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),plist_.at(ip));
    }
}

/** partial derivative of event-resolved output z w.r.t. to model states x
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
void Model::fdzdx(const realtype t, const int ie, const AmiVector *x) {
    std::fill(dzdx.begin(),dzdx.end(),0.0);
    fdzdx(dzdx.data(),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/** Sensitivity of event-resolved root output w.r.t. to model parameters p
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
void Model::fdrzdp(const realtype t, const int ie, const AmiVector *x) {
    std::fill(drzdp.begin(),drzdp.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++){
        fdrzdp(drzdp.data(),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),plist_.at(ip));
    }
}

/** Sensitivity of event-resolved measurements rz w.r.t. to model states x
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
void Model::fdrzdx(const realtype t, const int ie, const AmiVector *x) {
    std::fill(drzdx.begin(),drzdx.end(),0.0);
    fdrzdx(drzdx.data(),ie,t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/** State update functions for events
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
void Model::fdeltax(const int ie, const realtype t, const AmiVector *x,
                    const AmiVector *xdot, const AmiVector *xdot_old) {
    std::fill(deltax.begin(),deltax.end(),0.0);
    fdeltax(deltax.data(),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),ie,xdot->data(),xdot_old->data());
}

/** Sensitivity update functions for events, total derivative
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param sx current state sensitivity
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
void Model::fdeltasx(const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx,
                     const AmiVector *xdot, const AmiVector *xdot_old) {
    fw(t,x->getNVector());
    std::fill(deltasx.begin(),deltasx.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++)
        fdeltasx(&deltasx.at(nx*ip),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),w.data(),
                 plist_.at(ip),ie,xdot->data(),xdot_old->data(),sx->data(ip),&stau.at(ip));
}

/** Adjoint state update functions for events
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param xB current adjoint state
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
void Model::fdeltaxB(const int ie, const realtype t, const AmiVector *x, const AmiVector *xB,
                     const AmiVector *xdot, const AmiVector *xdot_old) {
    std::fill(deltaxB.begin(),deltaxB.end(),0.0);
    fdeltaxB(deltaxB.data(),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),ie,xdot->data(),xdot_old->data(),xB->data());
}

/** Quadrature state update functions for events
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param xB current adjoint state
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
void Model::fdeltaqB(const int ie, const realtype t, const AmiVector *x, const AmiVector *xB,
                     const AmiVector *xdot, const AmiVector *xdot_old) {
    std::fill(deltaqB.begin(),deltaqB.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++)
        fdeltaqB(deltaqB.data(),t,x->data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                 plist_.at(ip),ie,xdot->data(),xdot_old->data(),xB->data());
}

/** Standard deviation of measurements
     * @param it timepoint index
     * @param edata pointer to experimental data instance
     * @param rdata pointer to return data instance
     */
void Model::fsigmay(const int it, ReturnData *rdata, const ExpData *edata) {
    if (!ny)
        return;
    
    std::fill(sigmay.begin(),sigmay.end(),0.0);
    fsigmay(sigmay.data(),rdata->ts.at(it), unscaledParameters.data(),fixedParameters.data());
    for (int iytrue = 0; iytrue < nytrue; iytrue++) {
        /* extract the value for the standard deviation, if the data value
             is NaN, use the parameter value. Store this value in the return struct */
        if(edata){
            if (edata->isSetObservedDataStdDev(it, iytrue)) {
                auto sigmay_edata = edata->getObservedDataStdDevPtr(it);
                sigmay.at(iytrue) = sigmay_edata[iytrue];
            }
        }
        rdata->sigmay[it * rdata->ny + iytrue] = sigmay.at(iytrue);
    }
}

/** partial derivative of standard deviation of measurements w.r.t. model
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @param edata pointer to ExpData data instance holding sigma values
     */
void Model::fdsigmaydp(const int it, ReturnData *rdata, const ExpData *edata) {
    if (!ny)
        return;
    
    std::fill(dsigmaydp.begin(), dsigmaydp.end(), 0.0);

    for(int ip = 0; (unsigned)ip < plist_.size(); ip++)
        // get dsigmaydp slice (ny) for current timepoint and parameter
        fdsigmaydp(&dsigmaydp.at(ip*ny),
                    rdata->ts.at(it),
                    unscaledParameters.data(),
                    fixedParameters.data(),
                    plist_.at(ip));

    // sigmas in edata override model-sigma -> for those sigmas, set dsigmaydp to zero
    if(edata){
        for (int iy = 0; iy < nytrue; iy++) {
            if (!edata->isSetObservedDataStdDev(it, iy))
                continue;
            for (int ip = 0; ip < nplist(); ip++) {
                dsigmaydp[ip * ny + iy] = 0.0;
            }
        }
    }

    // copy dsigmaydp slice for current timepoint
    std::copy(dsigmaydp.begin(), dsigmaydp.end(), &rdata->ssigmay[it * nplist() * ny]);
}

/** Standard deviation of events
     * @param t current timepoint
     * @param ie event index
     * @param nroots array with event numbers
     * @param edata pointer to experimental data instance
     * @param rdata pointer to return data instance
     */
void Model::fsigmaz(const realtype t, const int ie, const int *nroots, ReturnData *rdata,
                    const ExpData *edata) {
    std::fill(sigmaz.begin(),sigmaz.end(),0.0);
    fsigmaz(sigmaz.data(),t, unscaledParameters.data(),fixedParameters.data());
    
    auto sigmaz_edata = edata->getObservedEventsStdDevPtr(nroots[ie]);
    for (int iztrue = 0; iztrue < nztrue; iztrue++) {
        if (z2event.at(iztrue) - 1 == ie) {
            if(edata) {
                if (edata->isSetObservedEventsStdDev(nroots[ie],iztrue)) {
                    sigmaz.at(iztrue) = sigmaz_edata[iztrue];
                }
            }
            rdata->sigmaz[nroots[ie]*rdata->nz + iztrue] = sigmaz.at(iztrue);
        }
    }
}

/** Sensitivity of standard deviation of events measurements w.r.t. model parameters p
 * @param t current timepoint
 * @param ie event index
 * @param nroots array with event numbers
 * @param rdata pointer to return data instance
 * @param edata pointer to experimental data instance
 */
void Model::fdsigmazdp(const realtype t, const int ie, const int *nroots, ReturnData *rdata, const ExpData *edata) {
    std::fill(dsigmazdp.begin(),dsigmazdp.end(),0.0);
    for(int ip = 0; (unsigned)ip < plist_.size(); ip++) {
        // get dsigmazdp slice (nz) for current event and parameter
        fdsigmazdp(&dsigmazdp.at(ip*nz),
                   t,
                   unscaledParameters.data(),
                   fixedParameters.data(),
                   plist_.at(ip));
    }
    
    // sigmas in edata override model-sigma -> for those sigmas, set dsigmazdp to zero
    if(edata) {
        for (int iz = 0; iz < nztrue; iz++) {
            if (z2event.at(iz) - 1 == ie && !edata->isSetObservedEventsStdDev(nroots[ie],iz)) {
                for(int ip = 0; (unsigned)ip < plist_.size(); ip++)
                    dsigmazdp.at(iz + nz * ip) = 0;
            }
        }
    }
    
    // copy dsigmazdp slice for current event
    std::copy(dsigmazdp.begin(), dsigmazdp.end(), &rdata->ssigmaz[nroots[ie] * nplist() * nz]);
}

/** negative log-likelihood of measurements y
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fJy(const int it, ReturnData *rdata, const ExpData *edata) {
    std::vector<realtype> nllh(nJ,0.0);
    getmy(it,edata);
    for(int iytrue = 0; iytrue < nytrue; iytrue++){
        if(!isNaN(my.at(iytrue))){
            std::fill(nllh.begin(),nllh.end(),0.0);
            fJy(nllh.data(),iytrue, unscaledParameters.data(),fixedParameters.data(),gety(it,rdata),sigmay.data(),my.data());
            rdata->llh -= nllh.at(0);
        }
    }
}

/** negative log-likelihood of event-resolved measurements z
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fJz(const int nroots, ReturnData *rdata, const ExpData *edata) {
    std::vector<realtype> nllh(nJ,0.0);
    getmz(nroots,edata);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            std::fill(nllh.begin(),nllh.end(),0.0);
            fJz(nllh.data(),iztrue, unscaledParameters.data(),fixedParameters.data(),getz(nroots,rdata),sigmaz.data(),mz.data());
            rdata->llh -= nllh.at(0);
        }
    }
}

/** regularization of negative log-likelihood with roots of event-resolved
     * measurements rz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fJrz(const int nroots, ReturnData *rdata, const ExpData *edata) {
    std::vector<realtype> nllh(nJ,0.0);
    getrz(nroots,rdata);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            std::fill(nllh.begin(),nllh.end(),0.0);
            fJrz(nllh.data(),iztrue, unscaledParameters.data(),fixedParameters.data(),getrz(nroots,rdata),sigmaz.data());
            rdata->llh -= nllh.at(0);
        }
    }
}

/** partial derivative of time-resolved measurement negative log-likelihood Jy
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fdJydy(const int it, const ReturnData *rdata,
                   const ExpData *edata) {
    // load measurements to my
    getmy(it,edata);
    std::fill(dJydy.begin(),dJydy.end(),0.0);

    for(int iytrue = 0; iytrue < nytrue; iytrue++){
        if(!isNaN(my.at(iytrue))){
            // get dJydy slice (ny) for current timepoint and observable
            fdJydy(&dJydy.at(iytrue*ny*nJ),
                   iytrue,
                   unscaledParameters.data(),
                   fixedParameters.data(),
                   gety(it,rdata),
                   sigmay.data(),
                   my.data());
        }
    }
}

/** Sensitivity of time-resolved measurement negative log-likelihood Jy
     * w.r.t. standard deviation sigma
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fdJydsigma(const int it, const ReturnData *rdata,
                       const ExpData *edata) {
    // load measurements to my
    getmy(it, edata);
    std::fill(dJydsigma.begin(),dJydsigma.end(),0.0);

    for(int iytrue = 0; iytrue < nytrue; iytrue++){
        if(!isNaN(my.at(iytrue))){
            // get dJydsigma slice (ny) for current timepoint and observable
            fdJydsigma(&dJydsigma.at(iytrue*ny*nJ),
                       iytrue,
                       unscaledParameters.data(),
                       fixedParameters.data(),
                       gety(it,rdata),
                       sigmay.data(),
                       my.data());
        }
    }
}

/** partial derivative of event measurement negative log-likelihood Jz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fdJzdz(const int nroots, const ReturnData *rdata,
                   const ExpData *edata) {
    getmz(nroots,edata);
    std::fill(dJzdz.begin(),dJzdz.end(),0.0);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            fdJzdz(&dJzdz.at(iztrue*nz*nJ),iztrue, unscaledParameters.data(),fixedParameters.data(),getz(nroots,rdata),sigmaz.data(),mz.data());
        }
    }
}

/** Sensitivity of event measurement negative log-likelihood Jz
     * w.r.t. standard deviation sigmaz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fdJzdsigma(const int nroots, const ReturnData *rdata,
                       const ExpData *edata) {
    getmz(nroots,edata);
    std::fill(dJzdsigma.begin(),dJzdsigma.end(),0.0);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            fdJzdsigma(&dJzdsigma.at(iztrue*nz*nJ),iztrue, unscaledParameters.data(),fixedParameters.data(),getz(nroots,rdata),sigmaz.data(),mz.data());
        }
    }
}

/** partial derivative of event measurement negative log-likelihood Jz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fdJrzdz(const int nroots, const ReturnData *rdata,
                    const ExpData *edata) {
    getmz(nroots,edata);
    std::fill(dJrzdz.begin(),dJrzdz.end(),0.0);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            fdJrzdz(&dJrzdz.at(iztrue*nz*nJ),iztrue, unscaledParameters.data(),fixedParameters.data(),getrz(nroots,rdata),sigmaz.data());
        }
    }
}

/** Sensitivity of event measurement negative log-likelihood Jz
     * w.r.t. standard deviation sigmaz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
void Model::fdJrzdsigma(const int nroots,const ReturnData *rdata,
                        const ExpData *edata) {
    std::fill(dJrzdsigma.begin(),dJrzdsigma.end(),0.0);
    for(int iztrue = 0; iztrue < nztrue; iztrue++){
        if(!isNaN(mz.at(iztrue))){
            fdJrzdsigma(&dJrzdsigma.at(iztrue*nz*nJ),iztrue, unscaledParameters.data(),fixedParameters.data(),getrz(nroots,rdata),sigmaz.data());
        }
    }
}

/**
 * @brief Recurring terms in xdot
 * @param t timepoint
 * @param x Vector with the states
 */
void Model::fw(const realtype t, const N_Vector x) {
    std::fill(w.begin(),w.end(),0.0);
    fw(w.data(),t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data());
}

/**
 * @brief Recurring terms in xdot, parameter derivative
 * @param t timepoint
 * @param x Vector with the states
 */
void Model::fdwdp(const realtype t, const N_Vector x) {
    fw(t,x);
    std::fill(dwdp.begin(),dwdp.end(),0.0);
    fdwdp(dwdp.data(),t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),w.data());
}

/**
 * @brief Recurring terms in xdot, state derivative
 * @param t timepoint
 * @param x Vector with the states
 */
void Model::fdwdx(const realtype t, const N_Vector x) {
    fw(t,x);
    std::fill(dwdx.begin(),dwdx.end(),0.0);
    fdwdx(dwdx.data(),t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),w.data());
}

void Model::fres(const int it, ReturnData *rdata, const ExpData *edata) {
    if (!edata || rdata->res.empty())
        return;
    
    auto observedData = edata->getObservedDataPtr(it);
    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * edata->nytrue;
        int iyt = iy + it * rdata->ny;
        if (!edata->isSetObservedData(it, iy))
            continue;
        rdata->res.at(iyt_true) = (rdata->y.at(iyt) - observedData[iy])/rdata->sigmay.at(iyt);
    }

}
    
void Model::fchi2(const int it, ReturnData *rdata) {
    if (rdata->res.empty())
        return;
    
    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * rdata->nytrue;
        rdata->chi2 += pow(rdata->res.at(iyt_true), 2);
    }
}
    
void Model::fsres(const int it, ReturnData *rdata, const ExpData *edata) {
    if (!edata || rdata->sres.empty())
        return;
    
    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * edata->nytrue;
        int iyt = iy + it * rdata->ny;
        if (!edata->isSetObservedData(it,iy))
            continue;
        for (int ip = 0; ip < nplist(); ++ip) {
            rdata->sres.at(iyt_true * nplist() + ip) =
            rdata->sy.at(iy + rdata->ny*(ip + it*nplist()))/rdata->sigmay.at(iyt);
        }
    }
}
    
void Model::fFIM(const int it, ReturnData *rdata) {
    if (rdata->sres.empty())
        return;

    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * rdata->nytrue;
        for (int ip = 0; ip < nplist(); ++ip) {
            for (int jp = 0; jp < nplist(); ++jp) {
                rdata->FIM.at(ip + nplist() * jp) +=
                rdata->sres.at(iyt_true * nplist() + ip)
                * rdata->sres.at(iyt_true * nplist() + jp);
            }
        }
    }
}
    
void Model::updateHeaviside(const std::vector<int>& rootsfound) {
    for (int ie = 0; ie < ne; ie++) {
        h.at(ie) += rootsfound.at(ie);
    }
}

void Model::updateHeavisideB(const int *rootsfound) {
    for (int ie = 0; ie < ne; ie++) {
        h.at(ie) -= rootsfound[ie];
    }
}

/** create my slice at timepoint
     * @param it timepoint index
     * @param edata pointer to experimental data instance
     */
void Model::getmy(const int it, const ExpData *edata){
    if(edata) {
        std::copy_n(edata->getObservedDataPtr(it), nytrue, my.begin());
    } else {
        std::fill(my.begin(), my.end(), getNaN());
    }
}

/** create x slice at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @return x x-slice from rdata instance
     */
const realtype *Model::getx(const int it, const ReturnData *rdata) const {
    return &rdata->x.at(it*nx);
}

/** create sx slice at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @return sx sx-slice from rdata instance
     */
const realtype *Model::getsx(const int it, const ReturnData *rdata) const {
    return &rdata->sx.at(it*nx*nplist());
}

/** create y slice at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @return y y-slice from rdata instance
     */
const realtype *Model::gety(const int it, const ReturnData *rdata) const {
    return &rdata->y.at(it*ny);
}


/** get current timepoint from index
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @return current timepoint
     */
realtype Model::gett(const int it, const ReturnData *rdata) const {
    return rdata->ts.at(it);
}

/** create mz slice at event
     * @param nroots event occurence
     * @param edata pointer to experimental data instance
     */
void Model::getmz(const int nroots, const ExpData *edata) {
    if(edata){
        std::copy_n(edata->getObservedEventsPtr(nroots), nztrue, mz.begin());
    } else {
        std::fill(mz.begin(), mz.end(), getNaN());
    }
}

/** create z slice at event
     * @param nroots event occurence
     * @param rdata pointer to return data instance
     * @return z slice
     */
const realtype *Model::getz(const int nroots, const ReturnData *rdata) const {
    return(&rdata->z.at(nroots*nz));
}

/** create rz slice at event
     * @param nroots event occurence
     * @param rdata pointer to return data instance
     * @return rz slice
     */
const realtype *Model::getrz(const int nroots, const ReturnData *rdata) const {
    return(&rdata->rz.at(nroots*nz));
}

/** create sz slice at event
     * @param nroots event occurence
     * @param ip sensitivity index
     * @param rdata pointer to return data instance
     * @return z slice
     */
const realtype *Model::getsz(const int nroots, const int ip, const ReturnData *rdata) const {
    return(&rdata->sz.at((nroots*nplist()+ip)*nz));
}

/** create srz slice at event
     * @param nroots event occurence
     * @param ip sensitivity index
     * @param rdata pointer to return data instance
     * @return rz slice
     */
const realtype *Model::getsrz(const int nroots, const int ip, const ReturnData *rdata) const {
    return(&rdata->srz.at((nroots*nplist()+ip)*nz));
}

int Model::checkFinite(const int N, const realtype *array, const char *fun) const
{
    auto result = amici::checkFinite(N, array, fun);

    if(result != AMICI_SUCCESS) {
        amici::checkFinite(ts.size(), ts.data(), "ts");
        amici::checkFinite(fixedParameters.size(), fixedParameters.data(), "k");
        amici::checkFinite(unscaledParameters.size(), unscaledParameters.data(), "p");
        amici::checkFinite(w.size(), w.data(), "w");
    }

    return result;
}

void Model::unscaleParameters(double *bufferUnscaled) const
{
    /**
         * unscaleParameters removes parameter scaling according to the parameter
         * scaling in pscale
         *
         * @param[out] bufferUnscaled unscaled parameters are written to the array
         * @type double
         *
         * @return status flag indicating success of execution @type int
         */
    for (int ip = 0; ip < np(); ++ip) {
        switch (pscale[ip]) {
        case ParameterScaling::log10:
            bufferUnscaled[ip] = pow(10, originalParameters[ip]);
            break;
        case ParameterScaling::ln:
            bufferUnscaled[ip] = exp(originalParameters[ip]);
            break;
        case ParameterScaling::none:
            bufferUnscaled[ip] = originalParameters[ip];
            break;
        }
    }
}

void Model::requireSensitivitiesForAllParameters() {
    plist_.resize(np());
    std::iota(plist_.begin(), plist_.end(), 0);
    initializeVectors();
}

bool operator ==(const Model &a, const Model &b)
{
    if (typeid(a) != typeid(b))
        return false;

    return (a.nx == b.nx)
            && (a.nxtrue == b.nxtrue)
            && (a.ny == b.ny)
            && (a.nytrue == b.nytrue)
            && (a.nz == b.nz)
            && (a.nztrue == b.nztrue)
            && (a.ne == b.ne)
            && (a.nw == b.nw)
            && (a.ndwdx == b.ndwdx)
            && (a.ndwdp == b.ndwdp)
            && (a.nnz == b.nnz)
            && (a.nJ == b.nJ)
            && (a.ubw == b.ubw)
            && (a.lbw == b.lbw)
            && (a.o2mode == b.o2mode)
            && (a.z2event == b.z2event)
            && (a.idlist == b.idlist)
            && (a.h == b.h)
            && (a.unscaledParameters == b.unscaledParameters)
            && (a.originalParameters == b.originalParameters)
            && (a.fixedParameters == b.fixedParameters)
            && (a.plist_ == b.plist_)
            && (a.x0data == b.x0data)
            && (a.sx0data == b.sx0data)
            && (a.ts == b.ts)
            && (a.nmaxevent == b.nmaxevent)
            && (a.pscale == b.pscale)
            && (a.tstart == b.tstart);
}



} // namespace amici
