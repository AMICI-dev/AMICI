#!/usr/bin/env python3

import sys
import numpy as np
from example import AmiciExample

class ExampleDirac(AmiciExample):

    def __init__(self):
        AmiciExample.__init__( self )

        self.numX = 2
        self.numP = 4
        self.numK = 0

        self.modelOptions['theta'] = np.log10([1, 0.5, 2, 3])
        self.modelOptions['ts'] = np.linspace(0, 3, 1001)
        self.modelOptions['pscale'] = 2

        self.solverOptions['atol'] = 1e-16
        self.solverOptions['maxsteps'] = 1e4
        self.solverOptions['nmaxevent'] = 10
        self.solverOptions['rtol'] = 1e-8
        self.solverOptions['sens_ind'] = []
        self.solverOptions['sensi'] = 0
        self.solverOptions['sensi_meth'] = 1


def writeNoSensi(filename):
    ex = ExampleDirac()

    ex.writeToFile(filename, '/model_dirac/nosensi/')


def writeSensiForward(filename):
    ex = ExampleDirac()

    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 1

    ex.writeToFile(filename, '/model_dirac/sensiforward/')

def main():
    if len(sys.argv) < 2:
        print("Error: Must provide output file as first and only argument.")
        sys.exit(1)
    filename = sys.argv[1]

    writeNoSensi(filename)
    writeSensiForward(filename)


if __name__ == "__main__":
    main()
