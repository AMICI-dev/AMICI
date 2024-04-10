"""
BNGL Import
------------
This module provides all necessary functionality to import a model specified
in the :term:`BNGL` format.
"""


from pysb.importers.bngl import model_from_bngl

from .pysb_import import pysb2amici


def bngl2amici(bngl_model: str, *args, **kwargs) -> None:
    r"""
    Generate AMICI C++ files for the provided model.

    :param bngl_model:
        bngl model file, model name will determine the name of the generated
        module

    :param args:
        see :func:`amici.pysb_import.pysb2amici` for additional arguments

    :param kwargs:
        see :func:`amici.pysb_import.pysb2amici` for additional arguments

    """
    if "model" in kwargs:
        raise ValueError("model argument not allowed")
    pysb_model = model_from_bngl(bngl_model)
    pysb2amici(pysb_model, *args, **kwargs)
