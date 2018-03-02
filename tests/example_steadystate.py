import h5py
import numpy as np
from doctest import Example

def dict2attrs(object, dictionary):
    for key, value in dictionary.items():
        object.attrs[key] = value
        
class ExampleSteadystate:
  
    def __init__(self):
        self.numX = 3
        self.numP = 5
        self.numK = 4
        
        self.modelOptions = {
            #'pscale' : []
            'theta' : np.log10([1, 0.5, 0.4, 2, 0.1]),
            'kappa' : [0.1, 0.4, 0.7, 1],
            'ts' : np.linspace(0, 100, 50),
            'tstart' : 0.0
        }
        self.solverOptions = {
            'atol' : 1e-16,
            'interpType' : 1,
            'ism': 1,
            'iter' : 2,
            'linsol': 9,
            'lmm': 2,
            'maxsteps' : 1e4,
            'newton_maxlinsteps': 100,
            'newton_maxsteps' : 40,
            # newton_precon ununsed
            'newton_preeq' : 0,
            'nmaxevent' : 10,
            'ordering' : 0,
            # pbar
            # 'qpositivex'
            'rtol' : 1e-8,
            'sens_ind' : [],
            'sensi': 0,
            'sensi_meth' : 1,
            'ss' : 0, # ?
            'stldet' : 1,
            # sx0
            # x0
            'z2event' : 0.0
        }
        self.modelOptions['pscale'] = 2
        self.modelOptions['qpositivex'] = [0] * self.numX 
        #self.modelOptions['sx0'] = [1.0] * np 
        #self.modelOptions['x0'] = [1.0] * np 


    def writeToFile(self, filename, root = "/"):
        with h5py.File(filename, "a") as f:
            g = f.require_group(root + 'options')
            
            dict2attrs(g, self.modelOptions)
            dict2attrs(g, self.solverOptions)

            
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
    filename = "testOptions.h5"
    writeNoSensi(filename)
    writeSensiForward(filename)
    writeSensiForwardPlist(filename)
    writeSensiForwardDense(filename)
    writeNosensiSPBCG(filename)
    writeSensiForwardErrorInt(filename)
    writeSensiForwardErrorNewt(filename)
    
if __name__ == "__main__":
    main()
    