#!/usr/bin/env python3
import sys
import numpy as np
from example import AmiciExample


class ExampleRobertson(AmiciExample):

    def __init__(self):
        AmiciExample.__init__( self )

        self.numX = 3
        self.numP = 3
        self.numK = 1


        self.modelOptions['theta'] = np.log10([0.04, 1e4, 3e7])
        self.modelOptions['kappa'] = [0.9]
        # self.modelOptions['ts'] = np.append(0, 4 * np.logspace(-6, 6))
        # logspace output from matlab slightly different:
        self.modelOptions['ts'] = [0.0, 4.0E-6, 7.030042499419172E-6, 1.2355374385909914E-5, 2.1714701757295438E-5, 3.8163819053999774E-5, 6.70733174724404E-5, 1.1788206810207238E-4, 2.071789871692485E-4, 3.641192711966091E-4, 6.39943487842423E-4, 0.0011247074791896924, 0.001976685344529533, 0.003474045495005412, 0.006105671868700933, 0.010730783181118898, 0.018859465453829577, 0.03314571091418737, 0.05825393910004978, 0.10238191690798133, 0.17993730675877775, 0.31624172843630804, 0.5557981977492555, 0.97682123781946, 1.7167737040515112, 3.017248025341849, 5.302845462360432, 9.319807242061488, 16.37966024952171, 28.787426920046055, 50.59420867421183, 88.91985930104782, 156.27759748218483, 274.6595380017199, 482.71705625573054, 848.380355168077, 1491.0374881259752, 2620.5142274381983, 4605.5815973057925, 8094.358590900622, 14225.921224892543, 25002.207701095904, 43941.64567950229, 77227.90915533014, 135728.87087581318, 238544.93266378547, 419245.253661875, 736827.9877306866, 1294983.017127056, 2275946.411607322, 4000000.0]
        self.modelOptions['pscale'] = 2

        self.solverOptions['atol'] = 1e-12
        self.solverOptions['rtol'] = 1e-8
        self.solverOptions['sens_ind'] = []
        self.solverOptions['sensi'] = 0
        self.solverOptions['sensi_meth'] = 1


def writeNoSensi(filename):
    ex = ExampleRobertson()

    ex.writeToFile(filename, '/model_robertson/nosensi/')


def writeSensiForward(filename):
    ex = ExampleRobertson()

    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 1

    ex.writeToFile(filename, '/model_robertson/sensiforward/')


def writeSensiForwardDense(filename):
    ex = ExampleRobertson()

    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 1
    ex.solverOptions['linsol'] = 1

    ex.writeToFile(filename, '/model_robertson/sensiforwarddense/')


def writeSensiForwardSPBCG(filename):
    ex = ExampleRobertson()

    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 1
    ex.solverOptions['linsol'] = 7
    ex.solverOptions['atol'] = 1e-14
    ex.solverOptions['rtol'] = 1e-12

    ex.writeToFile(filename, '/model_robertson/sensiforwardSPBCG/')


def main():
    if len(sys.argv) < 2:
        print("Error: Must provide output file as first and only argument.")
        sys.exit(1)
    filename = sys.argv[1]

    writeNoSensi(filename)
    writeSensiForward(filename)
    writeSensiForwardDense(filename)
    writeSensiForwardSPBCG(filename)


if __name__ == "__main__":
    main()
