from pysb.importers.bngl import model_from_bngl

from .pysb_import import pysb2amici


def bngl2amici(bngl_model:  str, **kwargs) -> None:
    r"""
    Generate AMICI C++ files for the provided model.

    :param bngl_model:
        bngl model, :attr:`pysb.Model.name` will determine the name of the
        generated module

    :param kwargs:
        see :meth:`amici.ode_export.ODEExporter.set_paths` for additional
        arguments

    """
    if 'model' in kwargs:
        raise ValueError('model argument not allowed')
    pysb_model = model_from_bngl(bngl_model)
    pysb2amici(model=pysb_model, **kwargs)
