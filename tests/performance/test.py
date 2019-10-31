#!/usr/bin/env python3

import amici
import sys

import CS_Signalling_ERBB_RAS_AKT_petab as model_module


def main():
    arg = sys.argv[1]

    model = model_module.getModel()
    solver = model.getSolver()
    # TODO
    edata = amici.ExpData(model)
    edata.setTimepoints([1e8])
    edata.setObservedData([1.0])

    if arg == 'forward_simulation':
        solver.setSensitivityMethod(amici.SensitivityMethod_none)
        solver.setSensitivityOrder(amici.SensitivityOrder_none)
    elif arg == 'forward_sensitivities':
        solver.setSensitivityMethod(amici.SensitivityMethod_forward)
        solver.setSensitivityOrder(amici.SensitivityOrder_first)
    elif arg == 'adjoint_sensitivities':
        solver.setSensitivityMethod(amici.SensitivityMethod_adjoint)
        solver.setSensitivityOrder(amici.SensitivityOrder_first)
    else:
        print("Unknown argument:", arg)
        sys.exit(1)
    rdata = amici.runAmiciSimulation(model, solver, edata)
    assert rdata['status'] == amici.AMICI_SUCCESS


if __name__ == '__main__':
    main()
