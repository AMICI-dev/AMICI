"""
C++ object views
----------------
This module provides views on C++ objects for efficient access.
"""

import numpy as np
import copy
import collections

from . import ExpDataPtr, ReturnDataPtr, ExpData, ReturnData
from typing import Union, List, Dict, Iterator


class SwigPtrView(collections.abc.Mapping):
    """
    Interface class to expose std::vector<double> and scalar members of
    swig wrapped C++ objects as numpy array attributes and fields. This
    class is memory efficient as copies of the underlying C++ objects is
    only created when respective fields are accessed for the first time.
    Cached copies are used for all subsequent calls.

    :ivar _swigptr: pointer to the c++ object
    :ivar _field_names: names of members that will be exposed as numpy arrays
    :ivar _field_dimensions: dimensions of numpy arrays
    :ivar _cache: dictionary with cached values
    """

    _swigptr = None
    _field_names: List[str] = []
    _field_dimensions: Dict[str, List[int]] = dict()

    def __getitem__(self, item: str) -> Union[np.ndarray, float]:
        """
        Access to field names, copies data from C++ object into numpy
        array, reshapes according to field dimensions and stores values in
        cache.

        :param item: field name
        :return: value
        """
        if self._swigptr is None:
            raise NotImplementedError('Cannot get items from abstract class.')

        if item == 'ptr':
            return self._swigptr

        if item in self._cache:
            return self._cache[item]

        if item == 'id':
            return getattr(self._swigptr, item)

        if item not in self._field_names:
            self.__missing__(item)

        value = field_as_numpy(
            self._field_dimensions, item, self._swigptr
        )
        self._cache[item] = value
        return value

    def __missing__(self, key: str) -> None:
        """
        Default behaviour for missing keys

        :param key: field name
        """
        raise KeyError(f'Unknown field name {key}.')

    def __getattr__(self, item) -> Union[np.ndarray, float]:
        """
        Attribute accessor for field names

        :param item: field name

        :returns: value
        """
        return self.__getitem__(item)

    def __init__(self, swigptr):
        """
        Constructor

        :param swigptr: pointer to the C++ object
        """
        self._swigptr = swigptr
        self._cache = dict()
        super(SwigPtrView, self).__init__()

    def __len__(self) -> int:
        """
        Returns the number of available keys/fields

        :returns: length of _field_names
        """
        return len(self._field_names)

    def __iter__(self) -> Iterator:
        """
        Create an iterator of the keys/fields

        :returns: iterator over _field_names
        """
        return iter(self._field_names)

    def __copy__(self):
        """
        Create a shallow copy

        :return: SwigPtrView shallow copy
        """
        other = SwigPtrView(self._swigptr)
        other._field_names = self._field_names
        other._field_dimensions = self._field_dimensions
        other._cache = self._cache
        return other

    def __contains__(self, item) -> bool:
        """
        Faster implementation of __contains__ that avoids copy of the field

        :param item: item to check for

        :returns: whether item is available as key
        """
        return item in self._field_names

    def __deepcopy__(self, memo):
        """
        Create a deep copy

        :param memo: dict with id-to-object mapping

        :returns: SwigPtrView deep copy
        """
        other = SwigPtrView(self._swigptr)
        other._field_names = copy.deepcopy(self._field_names)
        other._field_dimensions = copy.deepcopy(self._field_dimensions)
        other._cache = copy.deepcopy(self._cache)
        return other


class ReturnDataView(SwigPtrView):
    """
    Interface class for C++ Return Data objects that avoids possibly costly
    copies of member data.
    """

    _field_names = [
        'ts', 'x', 'x0', 'x_ss', 'sx', 'sx0', 'sx_ss', 'y', 'sigmay',
        'sy', 'ssigmay', 'z', 'rz', 'sigmaz', 'sz', 'srz',
        'ssigmaz', 'sllh', 's2llh', 'J', 'xdot', 'status', 'llh',
        'chi2', 'res', 'sres', 'FIM', 'w', 'preeq_wrms', 'preeq_t',
        'preeq_numsteps', 'preeq_numstepsB', 'preeq_status', 'preeq_cpu_time',
        'preeq_cpu_timeB', 'posteq_wrms', 'posteq_t', 'posteq_numsteps',
        'posteq_numstepsB', 'posteq_status', 'posteq_cpu_time',
        'posteq_cpu_timeB', 'numsteps', 'numrhsevals',
        'numerrtestfails', 'numnonlinsolvconvfails', 'order', 'cpu_time',
        'numstepsB', 'numrhsevalsB', 'numerrtestfailsB',
        'numnonlinsolvconvfailsB', 'cpu_timeB', 'cpu_time_total'
    ]

    def __init__(self, rdata: Union[ReturnDataPtr, ReturnData]):
        """
        Constructor

        :param rdata: pointer to the ReturnData instance
        """
        if not isinstance(rdata, (ReturnDataPtr, ReturnData)):
            raise TypeError(f'Unsupported pointer {type(rdata)}, must be'
                            f'amici.ExpDataPtr!')
        self._field_dimensions = {
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

            'res': [rdata.nt * rdata.nytrue *
                    (2 if rdata.sigma_res else 1)],
            'sres': [rdata.nt * rdata.nytrue *
                     (2 if rdata.sigma_res else 1), rdata.nplist],
            'FIM': [rdata.nplist, rdata.nplist],

            # diagnosis
            'J': [rdata.nx_solver, rdata.nx_solver],
            'w': [rdata.nt, rdata.nw],
            'xdot': [rdata.nx_solver],
            'preeq_numlinsteps': [rdata.newton_maxsteps, 2],
            'preeq_numsteps': [1, 3],
            'preeq_status': [1, 3],
            'posteq_numlinsteps': [rdata.newton_maxsteps, 2],
            'posteq_numsteps': [1, 3],
            'posteq_status': [1, 3],
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

    def __getitem__(self, item: str) -> Union[np.ndarray, ReturnDataPtr,
                                              ReturnData, float]:
        """
        Custom getitem implementation shim to map `t` to `ts`

        :param item: field/attribute key

        :returns: self[item]
        """
        if item == 't':
            item = 'ts'
        return super(ReturnDataView, self).__getitem__(item)


class ExpDataView(SwigPtrView):
    """
    Interface class for C++ Exp Data objects that avoids possibly costly
    copies of member data.
    """

    _field_names = [
        'observedData', 'observedDataStdDev', 'observedEvents',
        'observedEventsStdDev', 'fixedParameters',
        'fixedParametersPreequilibration',
        'fixedParametersPresimulation'
    ]

    def __init__(self, edata: Union[ExpDataPtr, ExpData]):
        """
        Constructor

        :param edata: pointer to the ExpData instance
        """
        if not isinstance(edata, (ExpDataPtr, ExpData)):
            raise TypeError(f'Unsupported pointer {type(edata)}, must be'
                            f'amici.ExpDataPtr!')
        self._field_dimensions = {  # observables
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


def field_as_numpy(field_dimensions: Dict[str, List[int]],
                   field: str, data: SwigPtrView) -> Union[np.ndarray,
                                                           float,
                                                           None]:
    """
    Convert data object field to numpy array with dimensions according to
    specified field dimensions

    :param field_dimensions: dimension specifications
            dict({field: list([dim1, dim2, ...])})
    :param data: object with fields
    :param field: Name of field

    :returns: Field Data as numpy array with dimensions according to
    specified field dimensions
    """
    attr = getattr(data, field)
    if field in field_dimensions:
        if len(attr) == 0:
            return None
        else:
            return np.array(attr).reshape(field_dimensions[field])
    else:
        return float(attr)
