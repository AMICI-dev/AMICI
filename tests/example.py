import h5py
import numpy as np

def dict2attrs(object, dictionary):
    for key, value in dictionary.items():
        object.attrs[key] = value
        
class AmiciExample:
  
    def __init__(self):
        
        self.modelOptions = {
            #'pscale' : []
            'theta' : [],
            'kappa' : [],
            'ts' : [],
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
        self.modelOptions['pscale'] = 0
        #self.modelOptions['qpositivex'] = [0] * self.numX 
        #self.modelOptions['sx0'] = [1.0] * np 
        #self.modelOptions['x0'] = [1.0] * np 
        self.data = {}

    def writeToFile(self, filename, root = "/"):
        with h5py.File(filename, "a") as f:
            g = f.require_group(root + '/options')
            dict2attrs(g, self.modelOptions)
            dict2attrs(g, self.solverOptions)

            if 'data' in self.__dict__ and len(self.data):
                g = f.require_group(root + '/data')
                dict2attrs(g, self.data)
    