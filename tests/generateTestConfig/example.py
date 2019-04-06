import h5py
import numpy as np
import pandas as pd

def dict2hdf5(object, dictionary):
    for key, value in dictionary.items():
        if isArray(value):
            a = np.array(value)
            if not len(value):
                dtype = 'f8'
            elif isArray(value[0]):
                if isinstance(value[0][0], (np.float64, float)):
                    dtype = 'f8'
                else:
                    dtype = '<i4'
            elif isinstance(value[0], (np.float64, float)):
                dtype = 'f8'
            else:
                dtype = '<i4'
            object.require_dataset(name=key,
                                   data=a,
                                   shape=a.shape,
                                   dtype=dtype)
        else:
            object.attrs[key] = value

def dict2attrs(object, dictionary):
    for key, value in dictionary.items():
        object.attrs[key] = value

def isArray(var):
    return isinstance(var, (list,
                            tuple,
                            np.ndarray,
                            pd.core.frame.DataFrame,
                            pd.core.series.Series))

class AmiciExample:

    def __init__(self):

        self.modelOptions = {
            #'pscale' : []
            'theta' : [],
            'kappa' : [],
            'ts' : [],
            'tstart' : 0.0,
            'pscale' : 0
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
            'newton_preeq' : 0,
            'nmaxevent' : 10,
            'ordering' : 0,
            'rtol' : 1e-8,
            'sens_ind' : [],
            'sensi': 0,
            'sensi_meth' : 1,
            'stldet' : 1,
            # sx0
            # x0
        }
        #self.modelOptions['sx0'] = [1.0] * np
        #self.modelOptions['x0'] = [1.0] * np
        self.data = {}

    def writeToFile(self, filename, root = "/"):
        with h5py.File(filename, "a") as f:
            g = f.require_group(root + '/options')
            dict2hdf5(g, self.modelOptions)
            dict2hdf5(g, self.solverOptions)

            if 'data' in self.__dict__ and len(self.data):
                g = f.require_group(root + '/data')
                dict2hdf5(g, self.data)

