#!/usr/bin/env python3

import amici

import CS_Signalling_ERBB_RAS_AKT_petab as model_module


def main():
    model = model_module.getModel()
    solver = model.getSolver()
    rdata = amici.runAmiciSimulation(model, solver)
    assert rdata['status'] == amici.AMICI_SUCCESS


if __name__ == '__main__':
    main()
