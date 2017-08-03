#ifndef AMICI_MODEL_FUNCTIONS_H
#define AMICI_MODEL_FUNCTIONS_H

class UserData;
class ReturnData;
class TempData;
class ExpData;
class Solver;
class Model;

UserData getUserData();
Solver *getSolver();
Model *getModel(UserData*, const ExpData *);

#endif // AMICI_MODEL_FUNCTIONS_H
