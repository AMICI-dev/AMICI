#!/usr/bin/env python3
import sys
import h5py
import numpy as np

from example import AmiciExample, dict2attrs

class ExampleSteadystate(AmiciExample):
  
    def __init__(self):
        AmiciExample.__init__( self )

        self.numX = 3
        self.numP = 5
        self.numK = 4
        
        self.modelOptions['theta'] = np.log10([1, 0.5, 0.4, 2, 0.1])
        self.modelOptions['kappa'] = [0.1, 0.4, 0.7, 1]
        self.modelOptions['ts'] = np.linspace(0, 100, 50)
        self.modelOptions['pscale'] = 2
        self.modelOptions['qpositivex'] = [0] * self.numX 
            
def writeNoSensi(filename):
    ex = ExampleSteadystate()

    ex.modelOptions['ts'] = np.append(np.linspace(0, 100, 50), np.inf)
    ex.writeToFile(filename, '/model_steadystate/nosensi/')


def writeSensiForward(filename):
    ex = ExampleSteadystate()

    ex.modelOptions['ts'] = np.append(np.linspace(0, 100, 50), np.inf)
    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP) + 1
    ex.solverOptions['pbar'] = [1.0] * len(ex.solverOptions['sens_ind'])
    ex.solverOptions['sensi'] = 1

    ex.writeToFile(filename, '/model_steadystate/sensiforward/')


def writeSensiForwardPlist(filename):
    ex = ExampleSteadystate()

    ex.modelOptions['ts'] = np.append(np.linspace(0, 100, 50), np.inf)
    ex.solverOptions['sens_ind'] = [3, 1, 2, 4]
    ex.solverOptions['pbar'] = [1.0] * len(ex.solverOptions['sens_ind'])
    ex.solverOptions['sensi'] = 1

    ex.writeToFile(filename, '/model_steadystate/sensiforwardplist/')

def writeSensiForwardDense(filename):
    ex = ExampleSteadystate()
    
    ex.modelOptions['ts'] = np.append(np.linspace(0, 100, 50), np.inf)
    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP) + 1
    ex.solverOptions['pbar'] = [1.0] * len(ex.solverOptions['sens_ind'])
    ex.solverOptions['sensi'] = 1
    ex.solverOptions['linsol'] = 1
    ex.writeToFile(filename, '/model_steadystate/sensiforwarddense/')

    
def writeNosensiSPBCG(filename):
    ex = ExampleSteadystate()

    ex.modelOptions['ts'] = np.append(np.linspace(0, 100, 50), np.inf)
    ex.solverOptions['sensi'] = 0
    ex.solverOptions['linsol'] = 7

    ex.writeToFile(filename, '/model_steadystate/nosensiSPBCG/')

def writeSensiForwardErrorInt(filename):
    ex = ExampleSteadystate()

    ex.modelOptions['ts'] = np.append(np.linspace(0, 100, 50), np.inf)
    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP) + 1
    ex.solverOptions['pbar'] = [1.0] * len(ex.solverOptions['sens_ind'])
    ex.solverOptions['sensi'] = 1
    ex.solverOptions['linsol'] = 9
    ex.solverOptions['maxsteps'] = 100

    ex.writeToFile(filename, '/model_steadystate/sensiforwarderrorint/')


def writeSensiForwardErrorNewt(filename):
    ex = ExampleSteadystate()
   
    ex.modelOptions['ts'] = [0, np.inf]
    ex.solverOptions['sensi'] = 1
    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP) + 1
    ex.solverOptions['pbar'] = [1.0] * len(ex.solverOptions['sens_ind'])
    ex.solverOptions['linsol'] = 9
    ex.solverOptions['maxsteps'] = 100
    ex.solverOptions['newton_maxsteps'] = 2
    
    ex.writeToFile(filename, '/model_steadystate/sensiforwarderrornewt/')


def main():
    if len(sys.argv) < 2:
        print("Error: Must provide output file as first and only argument.")
        sys.exit(1)
    filename = sys.argv[1]

    writeNoSensi(filename)
    writeSensiForward(filename)
    writeSensiForwardPlist(filename)
    writeSensiForwardDense(filename)
    writeNosensiSPBCG(filename)
    writeSensiForwardErrorInt(filename)
    writeSensiForwardErrorNewt(filename)
    
if __name__ == "__main__":
    main()
    