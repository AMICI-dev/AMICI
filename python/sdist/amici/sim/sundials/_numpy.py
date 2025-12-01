"""
C++ object views
----------------
This module provides views on C++ objects for efficient access.
"""

from __future__ import annotations

import collections
import copy
import itertools
from collections.abc import Iterator
from numbers import Number
from typing import Literal

import numpy as np
import sympy as sp
import xarray as xr
from sympy.abc import _clash

from . import (
    ExpData,
    ExpDataPtr,
    Model,
    ReturnData,
    ReturnDataPtr,
    SteadyStateStatus,
    simulation_status_to_str,
)

__all__ = [
    "ReturnDataView",
    "ExpDataView",
    "evaluate",
]

StrOrExpr = str | sp.Expr


class XArrayFactory:
    """
    Factory class to create xarray DataArrays for fields of a
    SwigPtrView instance.

    Currently, only ReturnDataView is supported.
    """

    def __init__(self, svp: SwigPtrView):
        """
        Constructor

        :param svp: SwigPtrView instance to create DataArrays from.
        """
        self._svp = svp

    def __getattr__(self, name: str) -> xr.DataArray:
        """
        Create xarray DataArray for field name

        :param name: field name

        :returns: xarray DataArray
        """
        data = getattr(self._svp, name)
        if data is None:
            return xr.DataArray(name=name)

        if isinstance(data, float):
            return xr.DataArray(data, name=name)

        if not isinstance(data, np.ndarray):
            raise TypeError(
                f"Cannot create xarray DataArray for field {name} of type"
                f" {type(data)}"
            )

        dims = None

        match name:
            case "x":
                coords = {
                    "time": self._svp.ts,
                    "state": list(self._svp.state_ids),
                }
            case "x0" | "x_ss":
                coords = {
                    "state": list(self._svp.state_ids),
                }
            case "xdot":
                coords = {
                    "state": list(self._svp.state_ids_solver),
                }
            case "y" | "sigmay":
                coords = {
                    "time": self._svp.ts,
                    "observable": list(self._svp.observable_ids),
                }
            case "sy" | "ssigmay":
                coords = {
                    "time": self._svp.ts,
                    "free_parameter": [
                        self._svp.free_parameter_ids[i]
                        for i in self._svp.plist
                    ],
                    "observable": list(self._svp.observable_ids),
                }
            case "w":
                coords = {
                    "time": self._svp.ts,
                    "expression": list(self._svp.expression_ids),
                }
            case "sx0":
                coords = {
                    "free_parameter": [
                        self._svp.free_parameter_ids[i]
                        for i in self._svp.plist
                    ],
                    "state": list(self._svp.state_ids),
                }
            case "sx":
                coords = {
                    "time": self._svp.ts,
                    "free_parameter": [
                        self._svp.free_parameter_ids[i]
                        for i in self._svp.plist
                    ],
                    "state": list(self._svp.state_ids),
                }
                dims = ("time", "free_parameter", "state")
            case "sllh":
                coords = {
                    "free_parameter": [
                        self._svp.free_parameter_ids[i]
                        for i in self._svp.plist
                    ]
                }
            case "FIM":
                coords = {
                    "free_parameter_1": [
                        self._svp.free_parameter_ids[i]
                        for i in self._svp.plist
                    ],
                    "free_parameter_2": [
                        self._svp.free_parameter_ids[i]
                        for i in self._svp.plist
                    ],
                }
            case "J":
                coords = {
                    "state_variable_1": list(self._svp.state_ids_solver),
                    "state_variable_2": list(self._svp.state_ids_solver),
                }
            case _:
                dims = tuple(f"dim_{i}" for i in range(data.ndim))
                coords = {
                    f"dim_{i}": np.arange(dim)
                    for i, dim in enumerate(data.shape)
                }
        arr = xr.DataArray(
            data,
            dims=dims,
            coords=coords,
            name=name,
        )
        return arr

    def __dir__(self):
        return sorted(
            set(
                itertools.chain(
                    dir(super()), self.__dict__, self._svp._field_names
                )
            )
        )


class SwigPtrView(collections.abc.Mapping):
    """
    Interface class to expose ``std::vector<double>`` and scalar members of
    swig wrapped C++ objects as numpy array attributes and fields. This
    class is memory efficient as copies of the underlying C++ objects is
    only created when respective fields are accessed for the first time.
    Cached copies are used for all subsequent calls.

    :ivar _swigptr: pointer to the C++ object
    :ivar _field_names: names of members that will be exposed as numpy arrays
    :ivar _field_dimensions: dimensions of numpy arrays
    :ivar _cache: dictionary with cached values
    """

    _swigptr = None
    _field_names: list[str] = []
    _field_dimensions: dict[str, list[int]] = dict()

    def __getitem__(self, item: str) -> np.ndarray | float:
        """
        Access to field names, copies data from C++ object into numpy
        array, reshapes according to field dimensions and stores values in
        cache.

        :param item: field name
        :return: value
        """
        if self._swigptr is None:
            raise NotImplementedError("Cannot get items from abstract class.")

        if item == "ptr":
            return self._swigptr

        if item in self._cache:
            return self._cache[item]

        if item in self._field_names:
            value = _field_as_numpy(
                self._field_dimensions, item, self._swigptr
            )
            self._cache[item] = value

            return value

        if not item.startswith("_") and hasattr(self._swigptr, item):
            return getattr(self._swigptr, item)

        self.__missing__(item)

    def __missing__(self, key: str) -> None:
        """
        Default behaviour for missing keys

        :param key: field name
        """
        raise KeyError(f"Unknown field name {key}.")

    def __getattr__(self, item) -> np.ndarray | float:
        """
        Attribute accessor for field names

        :param item: field name

        :returns: value
        """
        try:
            return self.__getitem__(item)
        except KeyError as e:
            raise AttributeError(item) from e

    def __init__(self, swigptr):
        """
        Constructor

        :param swigptr: pointer to the C++ object
        """
        self._swigptr = swigptr
        self._cache = {}

        super().__init__()

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
        Faster implementation of ``__contains__`` that avoids copy of the field

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
        # We assume we have a copy-ctor for the swigptr object
        other = self.__class__(copy.deepcopy(self._swigptr))
        other._field_names = copy.deepcopy(self._field_names)
        other._field_dimensions = copy.deepcopy(self._field_dimensions)
        other._cache = copy.deepcopy(self._cache)
        return other

    def __repr__(self):
        """
        String representation of the object

        :returns: string representation
        """
        return f"<{self.__class__.__name__}({self._swigptr})>"

    def __eq__(self, other):
        """
        Equality check

        :param other: other object

        :returns: whether other object is equal to this object
        """
        if not isinstance(other, self.__class__):
            return False
        return self._swigptr == other._swigptr

    def __dir__(self):
        return sorted(
            set(
                itertools.chain(dir(super()), self.__dict__, self._field_names)
            )
        )


class ReturnDataView(SwigPtrView):
    """
    Interface class for C++ :class:`ReturnData` objects that avoids
    possibly costly copies of member data.
    """

    _field_names = [
        "ts",
        "x",
        "x0",
        "x_ss",
        "sx",
        "sx0",
        "sx_ss",
        "y",
        "sigmay",
        "sy",
        "ssigmay",
        "z",
        "rz",
        "sigmaz",
        "sz",
        "srz",
        "ssigmaz",
        "sllh",
        "s2llh",
        "J",
        "xdot",
        "status",
        "llh",
        "chi2",
        "res",
        "sres",
        "FIM",
        "w",
        "preeq_wrms",
        "preeq_t",
        "preeq_numsteps",
        "preeq_numsteps_b",
        "preeq_status",
        "preeq_cpu_time",
        "preeq_cpu_time_b",
        "posteq_wrms",
        "posteq_t",
        "posteq_numsteps",
        "posteq_numsteps_b",
        "posteq_status",
        "posteq_cpu_time",
        "posteq_cpu_time_b",
        "numsteps",
        "num_rhs_evals",
        "num_err_test_fails",
        "num_non_lin_solv_conv_fails",
        "order",
        "cpu_time",
        "numsteps_b",
        "num_rhs_evals_b",
        "num_err_test_fails_b",
        "num_non_lin_solv_conv_fails_b",
        "cpu_time_b",
        "cpu_time_total",
        "messages",
        "t_last",
    ]

    def __init__(self, rdata: ReturnDataPtr | ReturnData):
        """
        Constructor

        :param rdata: pointer to the ``ReturnData`` instance
        """
        if not isinstance(rdata, (ReturnDataPtr | ReturnData)):
            raise TypeError(
                f"Unsupported pointer {type(rdata)}, must be"
                f"ReturnDataPtr or ReturnData!"
            )
        self._field_dimensions = {
            "ts": [rdata.nt],
            "x": [rdata.nt, rdata.nx_rdata],
            "x0": [rdata.nx_rdata],
            "x_ss": [rdata.nx_rdata],
            "sx": [rdata.nt, rdata.nplist, rdata.nx_rdata],
            "sx0": [rdata.nplist, rdata.nx_rdata],
            "sx_ss": [rdata.nplist, rdata.nx_rdata],
            # observables
            "y": [rdata.nt, rdata.ny],
            "sigmay": [rdata.nt, rdata.ny],
            "sy": [rdata.nt, rdata.nplist, rdata.ny],
            "ssigmay": [rdata.nt, rdata.nplist, rdata.ny],
            # event observables
            "z": [rdata.nmaxevent, rdata.nz],
            "rz": [rdata.nmaxevent, rdata.nz],
            "sigmaz": [rdata.nmaxevent, rdata.nz],
            "sz": [rdata.nmaxevent, rdata.nplist, rdata.nz],
            "srz": [rdata.nmaxevent, rdata.nplist, rdata.nz],
            "ssigmaz": [rdata.nmaxevent, rdata.nplist, rdata.nz],
            # objective function
            "sllh": [rdata.nplist],
            "s2llh": [rdata.np, rdata.nplist],
            "res": [rdata.nt * rdata.nytrue * (2 if rdata.sigma_res else 1)],
            "sres": [
                rdata.nt * rdata.nytrue * (2 if rdata.sigma_res else 1),
                rdata.nplist,
            ],
            "FIM": [rdata.nplist, rdata.nplist],
            # diagnosis
            "J": [rdata.nx_solver, rdata.nx_solver],
            "w": [rdata.nt, rdata.nw],
            "xdot": [rdata.nx_solver],
            "preeq_numlinsteps": [rdata.newton_maxsteps, 2],
            "preeq_numsteps": [3],
            "posteq_numlinsteps": [rdata.newton_maxsteps, 2],
            "posteq_numsteps": [3],
            "numsteps": [rdata.nt],
            "num_rhs_evals": [rdata.nt],
            "num_err_test_fails": [rdata.nt],
            "num_non_lin_solv_conv_fails": [rdata.nt],
            "order": [rdata.nt],
            "numsteps_b": [rdata.nt],
            "num_rhs_evals_b": [rdata.nt],
            "num_err_test_fails_b": [rdata.nt],
            "num_non_lin_solv_conv_fails_b": [rdata.nt],
        }
        self.xr = XArrayFactory(self)
        super().__init__(rdata)

    def __getitem__(
        self, item: str
    ) -> np.ndarray | ReturnDataPtr | ReturnData | float | int | list:
        """
        Access fields by name.

        Custom ``__getitem__`` implementation shim to map ``t`` to ``ts``.

        :param item: field/attribute key

        :returns: self[item]
        """
        if item == "status":
            return int(super().__getitem__(item))

        if item in ("preeq_status", "posteq_status"):
            return list(map(SteadyStateStatus, super().__getitem__(item)))

        if item == "t":
            item = "ts"

        return super().__getitem__(item)

    def __repr__(self):
        status = simulation_status_to_str(self._swigptr.status)
        return f"<{self.__class__.__name__}(id={self._swigptr.id!r}, status={status})>"

    def by_id(
        self, entity_id: str, field: str = None, model: Model = None
    ) -> np.array:
        """
        Get the value of a given field for a named entity.

        :param entity_id: The ID of the model entity that is to be extracted
            from ``field`` (e.g. a state ID).
        :param field: The requested field, e.g. 'x' for model states. This is
            optional if field would be one of ``{'x', 'y', 'w'}``
        :param model: The model from which this ReturnDataView was generated.
            This is optional if this ReturnData was generated with
            ``solver.getReturnDataReportingMode() == RDataReporting.full``.
        """
        if field is None:
            field = _entity_type_from_id(entity_id, self, model)

        if field in {"x", "x0", "x_ss", "sx", "sx0", "sx_ss"}:
            ids = (model and model.get_state_ids()) or self._swigptr.state_ids
        elif field in {"w"}:
            ids = (
                model and model.get_expression_ids()
            ) or self._swigptr.expression_ids
        elif field in {"y", "sy", "sigmay"}:
            ids = (
                model and model.get_observable_ids()
            ) or self._swigptr.observable_ids
        elif field in {"sllh"}:
            ids = (
                model and model.get_free_parameter_ids()
            ) or self._swigptr.free_parameter_ids
        else:
            raise NotImplementedError(
                f"Subsetting `{field}` by ID (`{entity_id}`) "
                "is not implemented or not possible."
            )
        col_index = ids.index(entity_id)
        return getattr(self, field)[:, ..., col_index]


class ExpDataView(SwigPtrView):
    """
    Interface class for C++ Exp Data objects that avoids possibly costly
    copies of member data.

    NOTE: This currently assumes that the underlying :class:`ExpData`
    does not change after instantiating an :class:`ExpDataView`.
    """

    _field_names = [
        "ts",
        "measurements",
        "noise_scale",
        "event_measurements",
        "event_noise_scale",
        "fixed_parameters",
        "fixed_parameters_pre_equilibration",
        "fixed_parameters_presimulation",
    ]

    def __init__(self, edata: ExpDataPtr | ExpData):
        """
        Constructor

        :param edata: pointer to the ExpData instance
        """
        if not isinstance(edata, (ExpDataPtr | ExpData)):
            raise TypeError(
                f"Unsupported pointer {type(edata)}, must be ExpDataPtr!"
            )
        self._field_dimensions = {
            "ts": [edata.nt()],
            # observables
            "measurements": [edata.nt(), edata.nytrue()],
            "noise_scale": [edata.nt(), edata.nytrue()],
            # event observables
            "event_measurements": [edata.nmaxevent(), edata.nztrue()],
            "event_noise_scale": [edata.nmaxevent(), edata.nztrue()],
            # fixed parameters
            "fixed_parameters": [len(edata.fixed_parameters)],
            "fixed_parameters_pre_equilibration": [
                len(edata.fixed_parameters_pre_equilibration)
            ],
            "fixed_parameters_presimulation": [
                len(edata.fixed_parameters_presimulation)
            ],
        }
        edata.ts = edata.timepoints
        edata.measurements = edata.get_measurements()
        edata.noise_scale = edata.get_noise_scale()
        edata.event_measurements = edata.get_event_measurements()
        edata.event_noise_scale = edata.get_event_noise_scale()
        super().__init__(edata)


def _field_as_numpy(
    field_dimensions: dict[str, list[int]], field: str, data: SwigPtrView
) -> np.ndarray | float | None:
    """
    Convert data object field to numpy array with dimensions according to
    specified field dimensions

    :param field_dimensions: dimension specifications
                ``dict({field: list([dim1, dim2, ...])})``
    :param data: object with fields
    :param field: Name of field

    :returns: Field Data as numpy array with dimensions according to
    specified field dimensions
    """
    attr = getattr(data, field)
    if field_dim := field_dimensions.get(field, None):
        return None if len(attr) == 0 else np.array(attr).reshape(field_dim)

    if isinstance(attr, Number):
        return float(attr)

    return attr


def _entity_type_from_id(
    entity_id: str,
    rdata: ReturnData | ReturnDataView = None,
    model: Model = None,
) -> Literal["x", "y", "w", "p", "k"]:
    """Guess the type of some entity by its ID."""
    for entity_type, symbol in (
        ("state", "x"),
        ("observable", "y"),
        ("expression", "w"),
        ("free_parameter", "p"),
        ("fixed_parameter", "k"),
    ):
        if model:
            if entity_id in getattr(model, f"get_{entity_type}_ids")():
                return symbol
        else:
            if entity_id in getattr(
                rdata if isinstance(rdata, ReturnData) else rdata._swigptr,
                f"{entity_type.lower()}_ids",
            ):
                return symbol

    raise KeyError(f"Unknown symbol {entity_id}.")


def evaluate(expr: StrOrExpr, rdata: ReturnDataView) -> np.array:
    """Evaluate a symbolic expression based on the given simulation outputs.

    :param expr:
        A symbolic expression, e.g. a sympy expression or a string that can be sympified.
        Can include state variable, expression, and observable IDs, depending on whether
        the respective data is available in the simulation results.
        Parameters are not yet supported.
    :param rdata:
        The simulation results.

    :return:
        The evaluated expression for the simulation output timepoints.
    """
    from sympy.utilities.lambdify import lambdify

    if isinstance(expr, str):
        expr = sp.sympify(expr, locals=_clash)

    arg_names = list(sorted(expr.free_symbols, key=lambda x: x.name))
    func = lambdify(arg_names, expr, "numpy")
    args = [rdata.by_id(arg.name) for arg in arg_names]
    return func(*args)
