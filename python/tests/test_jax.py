import pytest
import amici
import amici.jax

import jax.numpy as jnp
import numpy as np

from amici.pysb_import import pysb2amici
from numpy.testing import assert_allclose

pysb = pytest.importorskip("pysb")


def test_simulation():
    pysb.SelfExporter.cleanup()  # reset pysb
    pysb.SelfExporter.do_export = True

    model = pysb.Model('conversion')
    a = pysb.Monomer('A')
    b = pysb.Monomer('B')
    pysb.Initial(a(), pysb.Parameter('a0', 1.2))
    pysb.Rule(
        'conv',
        a() >> b(), pysb.Parameter('kcat', 0.05)
    )
    pysb.Observable('b', b())

    outdir = model.name
    pysb2amici(model, outdir, verbose=True,
               observables=['b'])

    model_module = amici.import_model_module(module_name=model.name,
                                             module_path=outdir)

    amici_model = model_module.getModel()

    ts = jnp.linspace(0, 1, 10)
    amici_model.setTimepoints(np.asarray(ts, dtype=np.float64))
    sol_amici_ref = amici.runAmiciSimulation(amici_model, amici_model.getSolver())

    jax_model = model_module.get_jax_model()
    jax_solver = jax_model.get_solver()

    p = jnp.stack((1.0, 0.1), axis=-1)
    k = jnp.empty((0,))

    amici_model.setParameters(np.asarray(p, dtype=np.float64))
    amici_model.setFixedParameters(np.asarray(k, dtype=np.float64))
    edata = amici.ExpData(sol_amici_ref, 1.0, 1.0)
    edata.parameters = amici_model.getParameters()
    edata.fixedParameters = amici_model.getFixedParameters()
    edata.pscale = amici_model.getParameterScale()
    r_amici = amici.runAmiciSimulation(
        amici_model,
        amici_model.getSolver(),
        edata
    )

    r_jax = amici.jax.runAmiciSimulationJAX(
        jax_model,
        jax_solver,
        edata
    )
    for field in ['x', 'y', 'sigmay', 'llh']:
        assert_allclose(
            actual=r_amici[field],
            desired=r_jax[field],
            atol=1e-6,
            rtol=1e-6
        )
