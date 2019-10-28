#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import os

from example import AmiciExample

class ExampleJakStatAdjoint(AmiciExample):

    def __init__(self):
        AmiciExample.__init__( self )

        self.numX = 9
        self.numP = 17
        self.numK = 2

        curPath = os.path.dirname(os.path.realpath(__file__))
        dataPath = curPath + "/../../matlab/examples/example_jakstat_adjoint/pnas_data_original.xls"
        xls = pd.ExcelFile(dataPath).parse()
        self.modelOptions['ts'] = xls.time
        self.modelOptions['theta'] = np.array([0.60, 3, -0.95, -0.0075, 0,
                                               -2.8, -0.26, -0.075, -0.41, -5,
                                               -0.74, -0.64, -0.11, 0.027, -0.5,
                                               0, -0.5])
        self.modelOptions['kappa'] = [1.4, 0.45]
        self.modelOptions['pscale'] = 2

        self.solverOptions['atol'] = 1e-16
        self.solverOptions['maxsteps'] = 1e4
        self.solverOptions['nmaxevent'] = 10
        self.solverOptions['rtol'] = 1e-12
        self.solverOptions['sensi'] = 0
        self.solverOptions['sensi_meth'] = 1

        self.data['Y'] = np.array(xls.loc[:,['pSTAT_au', 'tSTAT_au', 'pEpoR_au']])
        self.data['Sigma_Y'] = np.full(self.data['Y'].shape, np.nan)
        self.data['Sigma_Z'] = []
        self.data['Z'] = []
        self.data['condition'] = self.modelOptions['kappa']
        self.data['t'] = self.modelOptions['ts']


def writeNoSensi(filename):
    ex = ExampleJakStatAdjoint()

    ex.writeToFile(filename, '/model_jakstat_adjoint/nosensi/')


def writeSensiForward(filename):
    ex = ExampleJakStatAdjoint()

    ex.modelOptions['theta'] = 0.1 + ex.modelOptions['theta']
    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 1

    ex.writeToFile(filename, '/model_jakstat_adjoint/sensiforward/')

def writeSensiAdjoint(filename):
    ex = ExampleJakStatAdjoint()

    ex.modelOptions['theta'] = 0.1 + ex.modelOptions['theta']

    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 1
    ex.solverOptions['sensi_meth'] = 2

    ex.writeToFile(filename, '/model_jakstat_adjoint/sensiadjoint/')


def writeSensiForwardEmptySensInd(filename):
    ex = ExampleJakStatAdjoint()

    ex.solverOptions['sens_ind'] = np.array([])
    ex.solverOptions['sensi'] = 1
    ex.solverOptions['sensi_meth'] = 1

    ex.writeToFile(filename, '/model_jakstat_adjoint/sensiforwardemptysensind/')


def writeSensiAdjointEmptySensInd(filename):
    ex = ExampleJakStatAdjoint()

    ex.solverOptions['sens_ind'] = np.array([])
    ex.solverOptions['sensi'] = 1
    ex.solverOptions['sensi_meth'] = 2

    ex.writeToFile(filename, '/model_jakstat_adjoint/sensiadjointemptysensind/')


def writeSensi2Forward(filename):
    ex = ExampleJakStatAdjoint()

    ex.modelOptions['theta'] = 0.1 + ex.modelOptions['theta']

    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 2
    ex.solverOptions['sensi_meth'] = 1

    ex.writeToFile(filename, '/model_jakstat_adjoint/sensi2forward/')


def writeSensi2Adjoint(filename):
    ex = ExampleJakStatAdjoint()

    ex.modelOptions['theta'] = 0.1 + ex.modelOptions['theta']

    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 2
    ex.solverOptions['sensi_meth'] = 2

    ex.writeToFile(filename, '/model_jakstat_adjoint/sensi2adjoint/')

def writeSensiForwardLogParam(filename):
    ex = ExampleJakStatAdjoint()

    ex.modelOptions['theta'] = np.log(np.power(10.0, ex.modelOptions['theta'] + 0.1))
    ex.modelOptions['pscale'] = 1

    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 1
    ex.solverOptions['sensi_meth'] = 1

    ex.writeToFile(filename, '/model_jakstat_adjoint/sensiforwardlogparam/')

def writeSensi2ForwardLogParam(filename):
    ex = ExampleJakStatAdjoint()

    ex.modelOptions['theta'] = 0.1 + ex.modelOptions['theta']
    ex.modelOptions['pscale'] = 1

    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 2
    ex.solverOptions['sensi_meth'] = 1

    ex.writeToFile(filename, '/model_jakstat_adjoint/sensi2forwardlogparam/')

def main():
    if len(sys.argv) < 2:
        print("Error: Must provide output file as first and only argument.")
        sys.exit(1)
    filename = sys.argv[1]

    writeNoSensi(filename)
    writeSensiForward(filename)
    writeSensiAdjoint(filename)
    writeSensi2Forward(filename)
    writeSensi2Adjoint(filename)
    writeSensiForwardLogParam(filename)
    writeSensi2ForwardLogParam(filename)
    writeSensiForwardEmptySensInd(filename)
    writeSensiAdjointEmptySensInd(filename)


if __name__ == "__main__":
    main()

