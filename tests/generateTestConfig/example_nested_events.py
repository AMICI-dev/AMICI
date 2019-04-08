#!/usr/bin/env python3

import sys
import numpy as np

from example import AmiciExample

class ExampleNestedEvents(AmiciExample):

    def __init__(self):
        AmiciExample.__init__( self )

        self.numX = 1
        self.numP = 5
        self.numK = 0
        self.numZ = 0

        self.modelOptions['theta'] = np.log10([0.1, 1000, 2, 8e-1, 1.6])
        self.modelOptions['ts'] = np.linspace(0, 20, 100)
        self.modelOptions['pscale'] = 2

        self.solverOptions['atol'] = 1e-12
        self.solverOptions['maxsteps'] = 1e4
        self.solverOptions['nmaxevent'] = 2
        self.solverOptions['rtol'] = 1e-14
        self.solverOptions['sens_ind'] = []
        self.solverOptions['sensi'] = 0
        self.solverOptions['sensi_meth'] = 1

        self.data['Y'] = np.full((len(self.modelOptions['ts']), 1), np.nan)
        self.data['Sigma_Y'] = np.full((len(self.modelOptions['ts']), 1), np.nan)

        self.data['Z'] = np.full((self.numZ, self.solverOptions['nmaxevent']), np.nan)
        self.data['Sigma_Z'] = np.full((self.numZ, self.solverOptions['nmaxevent']), 0.5)

        self.data['condition'] = self.modelOptions['kappa']
        self.data['t'] = self.modelOptions['ts']


def writeNoSensi(filename):
    ex = ExampleNestedEvents()

    ex.writeToFile(filename, '/model_nested_events/nosensi/')


def writeSensiForward(filename):
    ex = ExampleNestedEvents()

    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 1

    ex.writeToFile(filename, '/model_nested_events/sensiforward/')


def main():
    if len(sys.argv) < 2:
        print("Error: Must provide output file as first and only argument.")
        sys.exit(1)
    filename = sys.argv[1]

    writeNoSensi(filename)
    writeSensiForward(filename)

if __name__ == "__main__":
    main()
