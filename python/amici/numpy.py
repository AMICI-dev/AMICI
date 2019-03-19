import numpy as np
import copy
import collections


class SwigClassView(collections.abc.Mapping):
    _swigptr = None
    _field_names = []
    _field_dimensions = []
    _cache = dict()

    def __getitem__(self, item):
        if len(self._field_names) == 0 or len(self._field_dimensions) == 0 or \
                self._swigptr is None:
            raise NotImplementedError('Cannot get items from abstract class.')

        if item == 'ptr':
            return self._swigptr

        if item in self._cache:
            return self._cache[item]

        if item not in self._field_names or item not in self._field_dimensions:
            self.__missing__(item)

        value = fieldAsNumpy(
            self._field_dimensions[item], self._field_names[item], self._swigptr
        )
        self._cache[item] = value
        return value

    def __missing__(self, key):
        raise KeyError(f'Unknown field name {key}.')

    def __getattr__(self, item):
        return self.__getitem__(item)

    def __init__(self, swigptr):
        self._swigptr = swigptr
        super(SwigClassView, self).__init__()

    def __len__(self):
        return len(self._field_names)

    def __iter__(self):
        return iter(self._field_names)

    def __copy__(self):
        other = SwigClassView(self._swigptr)
        other._field_names = self._field_names
        other._field_dimensions = self._field_dimensions
        other._cache = self._cache
        return other

    def __deepcopy__(self):
        other = SwigClassView(self._swigptr)
        other._field_names = copy.deepcopy(self._field_names)
        other._field_dimensions = copy.deepcopy(self._field_dimensions)
        other._cache = copy.deepcopy(self._cache)
        return other


class ReturnDataView(SwigClassView):
    """ Interface class for C++ Return Data objects that avoids possibly costly
    copies of member data.
    """

    _field_names = ['t', 'x', 'x0', 'x_ss', 'sx', 'sx0', 'sx_ss', 'y', 'sigmay',
                   'sy', 'ssigmay', 'z', 'rz', 'sigmaz', 'sz', 'srz',
                   'ssigmaz', 'sllh', 's2llh', 'J', 'xdot', 'status', 'llh',
                   'chi2', 'res', 'sres', 'FIM', 'wrms_steadystate',
                   't_steadystate', 'newton_numlinsteps', 'newton_numsteps',
                   'numsteps', 'numrhsevals', 'numerrtestfails',
                   'numnonlinsolvconvfails', 'order', 'numstepsB',
                    'numrhsevalsB', 'numerrtestfailsB',
                    'numnonlinsolvconvfailsB']

    def __init__(self, rdata):
        self.nt = rdata.nt
        self.nx = rdata.nx
        self.nxtrue = rdata.nxtrue
        self.nxsolver = rdata.nx_solver
        self.ny = rdata.ny
        self.nytrue = rdata.nytrue
        self.nz = rdata.nz
        self.nztrue = rdata.nztrue
        self.np = rdata.np
        self.nplist = rdata.nplist
        self.nJ = rdata.nJ
        self.nmaxevent = rdata.nmaxevent
        self.newton_maxsteps = self.newton_maxsteps

        self.field_dimensions = {
            'ts': [rdata.nt],
            'x': [rdata.nt, rdata.nx],
            'x0': [rdata.nx],
            'x_ss': [rdata.nx],
            'sx': [rdata.nt, rdata.nplist, rdata.nx],
            'sx0': [rdata.nplist, rdata.nx],
            'sx_ss': [rdata.nplist, rdata.nx],

            # observables
            'y': [rdata.nt, rdata.ny],
            'sigmay': [rdata.nt, rdata.ny],
            'sy': [rdata.nt, rdata.nplist, rdata.ny],
            'ssigmay': [rdata.nt, rdata.nplist, rdata.ny],

            # event observables
            'z': [rdata.nmaxevent, rdata.nz],
            'rz': [rdata.nmaxevent, rdata.nz],
            'sigmaz': [rdata.nmaxevent, rdata.nz],
            'sz': [rdata.nmaxevent, rdata.nplist, rdata.nz],
            'srz': [rdata.nmaxevent, rdata.nplist, rdata.nz],
            'ssigmaz': [rdata.nmaxevent, rdata.nplist, rdata.nz],

            # objective function
            'sllh': [rdata.nplist],
            's2llh': [rdata.np, rdata.nplist],

            'res': [rdata.nt * rdata.nytrue],
            'sres': [rdata.nt * rdata.nytrue, rdata.nplist],
            'FIM': [rdata.nplist, rdata.nplist],

            # diagnosis
            'J': [rdata.nx_solver, rdata.nx_solver],
            'xdot': [rdata.nx_solver],
            'newton_numlinsteps': [rdata.newton_maxsteps, 2],
            'newton_numsteps': [1, 3],
            'numsteps': [rdata.nt],
            'numrhsevals': [rdata.nt],
            'numerrtestfails': [rdata.nt],
            'numnonlinsolvconvfails': [rdata.nt],
            'order': [rdata.nt],
            'numstepsB': [rdata.nt],
            'numrhsevalsB': [rdata.nt],
            'numerrtestfailsB': [rdata.nt],
            'numnonlinsolvconvfailsB': [rdata.nt],
        }
        super(ReturnDataView, self).__init__(rdata)

    def __getitem__(self, item):
        if item == 't':
            item = 'ts'
        super(ReturnDataView, self).__init__(item)


class ExpDataView(SwigClassView):
    """ Interface class for C++ Exp Data objects that avoids possibly costly
    copies of member data.
    """

    _field_names = ['observedData', 'observedDataStdDev', 'observedEvents',
                   'observedEventsStdDev', 'fixedParameters',
                   'fixedParametersPreequilibration',
                   'fixedParametersPresimulation']

    def __init__(self, edata):
        self.field_dimensions = {  # observables
            'observedData': [edata.nt(), edata.nytrue()],
            'observedDataStdDev': [edata.nt(), edata.nytrue()],

            # event observables
            'observedEvents': [edata.nmaxevent(), edata.nztrue()],
            'observedEventsStdDev': [edata.nmaxevent(), edata.nztrue()],

            # fixed parameters
            'fixedParameters': [len(edata.fixedParameters)],
            'fixedParametersPreequilibration': [
                len(edata.fixedParametersPreequilibration)],
            'fixedParametersPresimulation': [
                len(edata.fixedParametersPreequilibration)],
        }
        edata.observedData = edata.getObservedData()
        edata.observedDataStdDev = edata.getObservedDataStdDev()
        edata.observedEvents = edata.getObservedEvents()
        edata.observedEventsStdDev = edata.getObservedEventsStdDev()
        super(ExpDataView, self).__init__(edata)

    def __getitem__(self, item):
        if item == 't':
            item = 'ts'
        super(ExpDataView, self).__init__(item)


def fieldAsNumpy(field_dimensions, field, data):
    """ Convert data object field to numpy array with dimensions according to
    specified field dimensions

    Arguments:
        field_dimensions: dimension specifications
            dict({field: list([dim1, dim2, ...])})
        data: object with fields
        field: Name of field

    Returns:
        Field Data as numpy array with dimensions according to specified field
        dimensions

    Raises:

    """
    attr = getattr(data, field)
    if field in field_dimensions.keys():
        if len(attr) == 0:
            return None
        else:
            return np.array(attr).reshape(field_dimensions[field])
    else:
        return float(attr)
