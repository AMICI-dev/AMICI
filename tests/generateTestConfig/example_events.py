#!/usr/bin/env python3

import sys
import numpy as np

from example import AmiciExample

class ExampleEvents(AmiciExample):

    def __init__(self):
        AmiciExample.__init__( self )

        self.numZ = 2
        self.numX = 3
        self.numP = 4
        self.numK = 4

        self.modelOptions['theta'] = np.log10([0.5, 2, 0.5, 0.5])
        self.modelOptions['kappa'] = [4.0, 8.0, 10.0, 4.0]
        self.modelOptions['ts'] = np.linspace(0.0, 10.0, 20)
        self.modelOptions['pscale'] = 2

        self.solverOptions['atol'] = 1e-16
        self.solverOptions['maxsteps'] = 1e4
        self.solverOptions['nmaxevent'] = 2
        self.solverOptions['rtol'] = 1e-8
        self.solverOptions['sens_ind'] = []
        self.solverOptions['sensi'] = 0
        self.solverOptions['sensi_meth'] = 1

        self.data['Y'] = np.full((len(self.modelOptions['ts']), 1), np.nan)
        self.data['Sigma_Y'] = np.full((len(self.modelOptions['ts']), 1), np.nan)

        self.data['Z'] = np.full((self.solverOptions['nmaxevent'], self.numZ ), np.nan)
        self.data['Sigma_Z'] = np.full((self.solverOptions['nmaxevent'], self.numZ ), np.nan)

        self.data['condition'] = self.modelOptions['kappa']
        self.data['t'] = self.modelOptions['ts']


def writeNoSensi(filename):
    ex = ExampleEvents()

    ex.writeToFile(filename, '/model_events/nosensi/')


def writeSensiForward(filename):
    ex = ExampleEvents()

    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 1
    ex.solverOptions['sensi_meth'] = 1

    ex.writeToFile(filename, '/model_events/sensiforward/')



def main():
    if len(sys.argv) < 2:
        print("Error: Must provide output file as first and only argument.")
        sys.exit(1)
    filename = sys.argv[1]
    writeNoSensi(filename)
    writeSensiForward(filename)

if __name__ == "__main__":
    main()
