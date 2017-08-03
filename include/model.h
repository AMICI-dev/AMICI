#ifndef MODEL_H
#define MODEL_H

#include <include/cvodewrap.h>

class UserData;
class ExpData;

/**
 * @brief The Model class represents an AMICI ODE model
 */
class Model
{
public:
    Model() {}
    // TODO model dimensions constructors

    virtual int fdx0(N_Vector x0, N_Vector dx0, void *user_data);
    virtual int fdx0(const realtype *k, realtype *x) {}

    virtual int fsx0(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data) {}

    virtual int fsdx0(N_Vector *sdx0, N_Vector x, N_Vector dx, void *user_data) {}

    virtual int fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {}

    virtual int fJB(long int NeqBdot, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {}

    virtual int froot(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data) {}

    virtual int frz(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, ReturnData *rdata) {}

    virtual int fsrz(realtype t, int ie, N_Vector x, N_Vector *sx, void *user_data, TempData *tdata, ReturnData *rdata) {}

    virtual int fstau(realtype t, int ie, N_Vector x, N_Vector *sx, void *user_data, TempData *tdata) {}

    virtual int fy(realtype t, int it, N_Vector x, void *user_data, ReturnData *rdata) {}

    virtual int fdydp(realtype t, int it, N_Vector x, void *user_data, TempData *tdata) {}

    virtual int fdydx(realtype t, int it, N_Vector x, void *user_data, TempData *tdata) {}

    virtual int fz(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, ReturnData *rdata) {}

    virtual int fsz(realtype t, int ie, N_Vector x, N_Vector *sx, void *user_data, TempData *tdata, ReturnData *rdata) {}

    virtual int fdzdp(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata) {}

    virtual int fdzdx(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata) {}

    virtual int fdrzdp(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata) {}

    virtual int fdrzdx(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata) {}

    virtual int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data) {}

    virtual int fxBdot(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data) {}

    virtual int fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data) {}

    virtual int fdxdotdp(realtype t, N_Vector x, N_Vector dx, void *user_data) {}

    virtual int fdeltax(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, void *user_data, TempData *tdata) {}

    virtual int fdeltasx(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, void *user_data, TempData *tdata) {}

    virtual int fdeltaxB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, void *user_data, TempData *tdata) {}

    virtual int fdeltaqB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, void *user_data, TempData *tdata) {}

    virtual int fsigma_y(realtype t, void *user_data, TempData *tdata) {}

    virtual int fdsigma_ydp(realtype t, void *user_data, TempData *tdata) {}

    virtual int fsigma_z(realtype t, int ie, void *user_data, TempData *tdata) {}

    virtual int fdsigma_zdp(realtype t, int ie, void *user_data, TempData *tdata) {}

    virtual int fJy(realtype t, int it, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {}

    virtual int fJz(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {}

    virtual int fJrz(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {}

    virtual int fdJydy(realtype t, int it, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {}

    virtual int fdJydsigma(realtype t, int it, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {}

    virtual int fdJzdz(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {}

    virtual int fdJzdsigma(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {}

    virtual int fdJrzdz(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {}

    virtual int fdJrzdsigma(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {}


    // Generic implementations
    static int fsy(int it, UserData *udata, TempData *tdata, ReturnData *rdata);

    static int fsz_tf(int ie, UserData *udata, TempData *tdata, ReturnData *rdata);

    static  int fsJy(int it, UserData *udata, TempData *tdata, const ExpData *edata, ReturnData *rdata);

    static  int fdJydp(int it, UserData *udata, TempData *tdata, const ExpData *edata, ReturnData *rdata);

    static  int fdJydx(int it, UserData *udata, TempData *tdata, const ExpData *edata);

    static  int fsJz(int ie, UserData *udata, TempData *tdata, const ExpData *edata, ReturnData *rdata);

    static  int fdJzdp(int ie, UserData *udata, TempData *tdata, const ExpData *edata, ReturnData *rdata);

    static  int fdJzdx(int ie, UserData *udata, TempData *tdata, const ExpData *edata);


    UserData *udata;
    ExpData *edata;
};

#endif // MODEL_H
