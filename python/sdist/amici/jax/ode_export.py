"""
JAX Export
----------
This module provides all necessary functionality to specify an ordinary
differential equation model and generate executable jax simulation code.
The user generally won't have to directly call any function from this module
as this will be done by
:py:func:`amici.importers.pysb.pysb2jax`,
:py:func:`amici.importers.sbml.SbmlImporter.` and
:py:func:`amici.petab_import.import_model`.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import sympy as sp
import numpy as np

from amici import (
    amiciModulePath,
)
from amici._symbolic.de_model import DEModel
from amici._symbolic.sympy_utils import (
    _monkeypatch_sympy,
)
from amici.exporters.sundials.de_export import is_valid_identifier
from amici.exporters.template import apply_template
from amici.jax.jaxcodeprinter import AmiciJaxCodePrinter, _jnp_array_str
from amici.jax.model import JAXModel
from amici.jax.nn import generate_equinox
from amici.logging import get_logger, log_execution_time, set_log_level

#: python log manager
logger = get_logger(__name__, logging.ERROR)


def _jax_variable_assignments(
    model: DEModel, sym_names: tuple[str, ...]
) -> dict:
    return {
        f"{sym_name.upper()}_SYMS": "".join(
            f"{s.name}, " for s in model.sym(sym_name)
        )
        if model.sym(sym_name)
        else "_"
        for sym_name in sym_names
    }


def _jax_variable_equations(
    model: DEModel,
    code_printer: AmiciJaxCodePrinter,
    eq_names: tuple[str, ...],
    subs: dict,
    indent: int,
) -> dict:
    return {
        f"{eq_name.upper()}_EQ": "\n".join(
            code_printer._get_sym_lines(
                (s.name for s in model.sym(eq_name)),
                # sp.Matrix to support event assignments which are lists
                sp.Matrix(model.eq(eq_name)).subs(subs), 
                indent,
            )
        )[indent:]  # remove indent for first line
        for eq_name in eq_names
    }


def _jax_return_variables(
    model: DEModel,
    eq_names: tuple[str, ...],
) -> dict:
    return {
        f"{eq_name.upper()}_RET": _jnp_array_str(
            s.name for s in model.sym(eq_name)
        )
        if model.sym(eq_name) and sp.Matrix(model.eq(eq_name)).shape[0]
        else "jnp.array([])"
        for eq_name in eq_names
    }


def _jax_variable_ids(model: DEModel, sym_names: tuple[str, ...]) -> dict:
    return {
        f"{sym_name.upper()}_IDS": "".join(
            f'"{s.name}", ' for s in model.sym(sym_name)
        )
        if model.sym(sym_name)
        else "tuple()"
        for sym_name in sym_names
    }


class ODEExporter:
    """
    The ODEExporter class generates AMICI jax files for a model as
    defined in _symbolic expressions.

    :ivar model:
        DE definition

    :ivar verbose:
        more verbose output if True

    :ivar model_name:
        name of the model that will be used for compilation

    :ivar model_path:
        path to the generated model specific files

    :ivar _code_printer:
        Code printer to generate JAX code
    """

    def __init__(
        self,
        ode_model: DEModel,
        output_dir: Path | str | None = None,
        verbose: bool | int | None = False,
        model_name: str | None = "model",
        hybridization: dict[str, dict] = None,
    ):
        """
        Generate AMICI jax files for the ODE provided to the constructor.

        :param ode_model:
            DE model definition

        :param output_dir:
            see :meth:`amici.de_export.DEExporter.set_paths`

        :param verbose:
            verbosity level for logging, ``True``/``False`` default to
            :data:`logging.Error`/:data:`logging.DEBUG`

        :param model_name:
            name of the model to be used during code generation

        :param hybridization:
            dict representation of the hybridization information in the PEtab YAML file, see
            https://petab-sciml.readthedocs.io/latest/format.html#problem-yaml-file
        """
        set_log_level(logger, verbose)

        if ode_model.has_algebraic_states():
            raise NotImplementedError(
                "The JAX backend does not support models with algebraic states."
            )

        if ode_model.has_priority_events():
            raise NotImplementedError(
                "The JAX backend does not support event priorities."
            )
        
        if ode_model.has_implicit_event_assignments():
            raise NotImplementedError(
                "The JAX backend does not support event assignments with implicit triggers."
            )

        self.verbose: bool = logger.getEffectiveLevel() <= logging.DEBUG

        self.model_path: Path = Path()

        self.set_name(model_name)
        self.set_paths(output_dir)

        self.model: DEModel = ode_model

        self.hybridization = hybridization if hybridization is not None else {}

        self._code_printer = AmiciJaxCodePrinter()

    @_monkeypatch_sympy
    @log_execution_time("generating jax code", logger)
    def generate_model_code(self) -> None:
        """
        Generates the jax code for the loaded model
        """
        self._prepare_model_folder()
        self._generate_jax_code()
        self._generate_nn_code()

    def _prepare_model_folder(self) -> None:
        """
        Create model directory or remove all files if the output directory
        already exists.
        """
        self.model_path.mkdir(parents=True, exist_ok=True)

        for file in self.model_path.glob("*"):
            if file.is_file():
                file.unlink()

    @log_execution_time("generating jax code", logger)
    def _generate_jax_code(self) -> None:
        eq_names = (
            "xdot",
            "w",
            "x0",
            "y",
            "sigmay",
            "Jy",
            "x_solver",
            "x_rdata",
            "total_cl",
            "eroot",
            "iroot",
            "deltax",
        )
        sym_names = (
            "p",
            "np",
            "op",
            "x",
            "tcl",
            "ih",
            "allh",
            "w",
            "my",
            "y",
            "sigmay",
            "x_rdata",
            "iroot",
            "x_old",
        )

        indent = 8

        # replaces Heaviside variables with corresponding functions
        subs_heaviside = dict(
            zip(
                self.model.sym("eh"),
                [sp.Heaviside(x) for x in self.model.eq("eroot")],
                strict=True,
            )
        )

        # replaces observables with a generic my variable
        subs_observables = dict(
            zip(
                self.model.sym("my"),
                [sp.Symbol("my")] * len(self.model.sym("my")),
                strict=True,
            )
        )
        subs = subs_heaviside | subs_observables

        tpl_data = {
            # assign named variable using corresponding algebraic formula (function body)
            **_jax_variable_equations(
                self.model, self._code_printer, eq_names, subs, indent
            ),
            # create jax array from concatenation of named variables
            **_jax_return_variables(self.model, eq_names),
            # assign named variables from a jax array
            **_jax_variable_assignments(self.model, sym_names),
            # tuple of variable names (ids as they are unique)
            **_jax_variable_ids(self.model, ("p", "k", "y", "w", "x_rdata")),
            "P_VALUES": _jnp_array_str(self.model.val("p")),
            "ROOTS": _jnp_array_str(
                {
                    _print_trigger_root(root)
                    for e in self.model._events
                    for root in e.get_trigger_times()
                }
            ),
            "N_IEVENTS": str(len(self.model.get_implicit_roots())),
            "N_EEVENTS": str(len(self.model.get_explicit_roots())),
            "EVENT_INITIAL_VALUES": _jnp_array_str(
                [
                    e.get_initial_value() for e in self.model._events
                ]
            ),
            **{
                "MODEL_NAME": self.model_name,
                # keep track of the API version that the model was generated with so we
                # can flag conflicts in the future
                "MODEL_API_VERSION": f"'{JAXModel.MODEL_API_VERSION}'",
            },
            "NET_IMPORTS": "\n".join(
                f"{net} = _module_from_path('{net}', Path(__file__).parent / '{net}.py')"
                for net in self.hybridization.keys()
            ),
            "NETS": ",\n".join(
                f'"{net}": {net}.net(jr.PRNGKey(0))'
                for net in self.hybridization.keys()
            ),
        }

        apply_template(
            Path(amiciModulePath) / "jax" / "jax.template.py",
            self.model_path / "__init__.py",
            tpl_data,
        )

    def _generate_nn_code(self) -> None:
        for net_name, net in self.hybridization.items():
            generate_equinox(
                net["model"],
                self.model_path / f"{net_name}.py",
                net["frozen_layers"],
            )

    def _implicit_roots(self) -> list[sp.Expr]:
        """Return root functions that require rootfinding."""
        roots = []
        for root in self.model.get_implicit_roots():
            if any(
                sp.simplify(root + r) == 0 or sp.simplify(root - r) == 0
                for r in roots
            ):
                continue
            roots.append(root)
        return roots

    def set_paths(self, output_dir: str | Path | None = None) -> None:
        """
        Set output paths for the model and create if necessary

        :param output_dir:
            relative or absolute path where the generated model
            code is to be placed. If ``None``, this will default to
            ``amici-{self.model_name}`` in the current working directory.
            will be created if it does not exist.

        """
        if output_dir is None:
            output_dir = Path(os.getcwd()) / f"amici-{self.model_name}"

        self.model_path = Path(output_dir).resolve()
        self.model_path.mkdir(parents=True, exist_ok=True)

    def set_name(self, model_name: str) -> None:
        """
        Sets the model name

        :param model_name:
            name of the model (may only contain upper and lower case letters,
            digits and underscores, and must not start with a digit)
        """
        if not is_valid_identifier(model_name):
            raise ValueError(
                f"'{model_name}' is not a valid model name. "
                "Model name may only contain upper and lower case letters, "
                "digits and underscores, and must not start with a digit."
            )

        self.model_name = model_name

def _print_trigger_root(root: sp.Expr) -> str:
    """Convert a trigger root expression into a string representation.

    :param root: The trigger root expression.
    :return: A string representation of the trigger root.
    """
    if root.is_number:
        return float(root)
    return str(root).replace(" ", "")
