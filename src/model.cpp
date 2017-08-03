#include "../include/model.h"

Model::Model()
{

}

int Model::fdx0(N_Vector x0, N_Vector dx0, void *user_data)
{
    UserData *udata = (UserData*) user_data;
    realtype *x0_tmp = N_VGetArrayPointer(x0);

    return fdx0(udata->k, x0_tmp);
}

int Model::fsy(int it, UserData *udata, TempData *tdata, ReturnData *rdata) {}

int Model::fsz_tf(int ie, UserData *udata, TempData *tdata, ReturnData *rdata) {}

int Model::fsJy(int it, UserData *udata, TempData *tdata, const ExpData *edata, ReturnData *rdata) {}

int Model::fdJydp(int it, UserData *udata, TempData *tdata, const ExpData *edata, ReturnData *rdata) {}

int Model::fdJydx(int it, UserData *udata, TempData *tdata, const ExpData *edata) {}

int Model::fsJz(int ie, UserData *udata, TempData *tdata, const ExpData *edata, ReturnData *rdata) {}

int Model::fdJzdp(int ie, UserData *udata, TempData *tdata, const ExpData *edata, ReturnData *rdata) {}

int Model::fdJzdx(int ie, UserData *udata, TempData *tdata, const ExpData *edata) {}
