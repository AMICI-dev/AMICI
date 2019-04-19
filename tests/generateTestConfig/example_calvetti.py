#!/usr/bin/env python3
import sys
import numpy as np
from example import AmiciExample


class ExampleCalvetti(AmiciExample):

    def __init__(self):
        AmiciExample.__init__( self )

        self.numX = 6
        self.numP = 0
        self.numK = 6

        self.modelOptions['theta'] = []
        self.modelOptions['kappa'] = [0.29, 0.74, 0.44, 0.08, 0.27, 0.18]
        self.modelOptions['ts'] = np.linspace(0, 20, 201)
        self.modelOptions['pscale'] = 0

        self.solverOptions['atol'] = 1e-6
        self.solverOptions['rtol'] = 1e-4
        self.solverOptions['sens_ind'] = []
        self.solverOptions['sensi'] = 0
        self.solverOptions['sensi_meth'] = 1


def writeNoSensi(filename):
    ex = ExampleCalvetti()

    ex.writeToFile(filename, '/model_calvetti/nosensi/')


def main():
    if len(sys.argv) < 2:
        print("Error: Must provide output file as first and only argument.")
        sys.exit(1)
    filename = sys.argv[1]

    writeNoSensi(filename)


if __name__ == "__main__":
    main()
