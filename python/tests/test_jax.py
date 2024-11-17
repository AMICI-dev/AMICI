import pytest
import amici

pytest.importorskip("jax")
import amici.jax

import jax.numpy as jnp
import jax
import diffrax
import numpy as np
from beartype import beartype

from amici.pysb_import import pysb2amici
from numpy.testing import assert_allclose

pysb = pytest.importorskip("pysb")

jax.config.update("jax_enable_x64", True)


def test_conversion():
    pysb.SelfExporter.cleanup()  # reset pysb
    pysb.SelfExporter.do_export = True

    model = pysb.Model("conversion")
    a = pysb.Monomer("A", sites=["s"], site_states={"s": ["a", "b"]})
    pysb.Initial(a(s="a"), pysb.Parameter("aa0", 1.2))
    pysb.Rule("conv", a(s="a") >> a(s="b"), pysb.Parameter("kcat", 0.05))
    pysb.Observable("ab", a(s="b"))

    outdir = model.name
    pysb2amici(model, outdir, verbose=True, observables=["ab"])

    model_module = amici.import_model_module(
        module_name=model.name, module_path=outdir
    )

    ts = tuple(np.linspace(0, 1, 10))
    p = jnp.stack((1.0, 0.1), axis=-1)
    k = tuple()
    _test_model(model_module, ts, p, k)


@pytest.mark.filterwarnings(
    "ignore:Model does not contain any initial conditions"
)
def test_dimerization():
    pysb.SelfExporter.cleanup()  # reset pysb
    pysb.SelfExporter.do_export = True

    model = pysb.Model("dimerization")
    a = pysb.Monomer("A", sites=["b"])
    b = pysb.Monomer("B", sites=["a"])

    pysb.Rule(
        "turnover_a",
        a(b=None) | None,
        pysb.Parameter("kdeg_a", 10),
        pysb.Parameter("ksyn_a", 0.1),
    )
    pysb.Rule(
        "turnover_b",
        b(a=None) | None,
        pysb.Parameter("kdeg_b", 0.1),
        pysb.Parameter("ksyn_b", 10),
    )
    pysb.Rule(
        "dimer",
        a(b=None) + b(a=None) | a(b=1) % b(a=1),
        pysb.Parameter("kon", 1.0),
        pysb.Parameter("koff", 0.1),
    )

    pysb.Observable("a_obs", a())
    pysb.Observable("b_obs", b())

    outdir = model.name
    pysb2amici(
        model,
        outdir,
        verbose=True,
        observables=["a_obs", "b_obs"],
        constant_parameters=["ksyn_a", "ksyn_b"],
    )

    model_module = amici.import_model_module(
        module_name=model.name, module_path=outdir
    )

    ts = tuple(np.linspace(0, 1, 10))
    p = jnp.stack((5, 0.5, 0.5, 0.5), axis=-1)
    k = (0.5, 5)
    _test_model(model_module, ts, p, k)


def _test_model(model_module, ts, p, k):
    amici_model = model_module.getModel()

    amici_model.setTimepoints(np.asarray(ts, dtype=np.float64))
    sol_amici_ref = amici.runAmiciSimulation(
        amici_model, amici_model.getSolver()
    )

    jax_model = model_module.get_jax_model()

    amici_model.setParameters(np.asarray(p, dtype=np.float64))
    amici_model.setFixedParameters(np.asarray(k, dtype=np.float64))
    edata = amici.ExpData(sol_amici_ref, 1.0, 1.0)
    edata.parameters = amici_model.getParameters()
    edata.fixedParameters = amici_model.getFixedParameters()
    edata.pscale = amici_model.getParameterScale()
    amici_solver = amici_model.getSolver()
    amici_solver.setSensitivityMethod(amici.SensitivityMethod.forward)
    amici_solver.setSensitivityOrder(amici.SensitivityOrder.first)
    rs_amici = amici.runAmiciSimulations(amici_model, amici_solver, [edata])

    check_fields_jax(
        rs_amici, jax_model, edata, ["x", "y", "llh", "res", "x0"]
    )

    check_fields_jax(
        rs_amici,
        jax_model,
        edata,
        ["sllh", "sx0", "sx", "sres", "sy"],
        sensi_order=amici.SensitivityOrder.first,
    )


def check_fields_jax(
    rs_amici,
    jax_model,
    edata,
    fields,
    sensi_order=amici.SensitivityOrder.none,
):
    r_jax = dict()
    ts = np.array(edata.getTimepoints())
    my = np.array(edata.getObservedData()).reshape(len(ts), -1)
    ts = np.repeat(ts.reshape(-1, 1), my.shape[1], axis=1)
    iys = np.repeat(np.arange(my.shape[1]).reshape(1, -1), len(ts), axis=0)
    my = my.flatten()
    ts = ts.flatten()
    iys = iys.flatten()

    ts_preeq = ts[ts == 0]
    ts_dyn = ts[ts > 0]
    ts_posteq = np.array([])
    p = jnp.array(list(edata.parameters) + list(edata.fixedParameters))
    args = (
        jnp.array([]),  # p_preeq
        jnp.array(ts_preeq),  # ts_preeq
        jnp.array(ts_dyn),  # ts_dyn
        jnp.array(ts_posteq),  # ts_posteq
        jnp.array(my),  # my
        jnp.array(iys),  # iys
        diffrax.Kvaerno5(),  # solver
        diffrax.PIDController(atol=1e-8, rtol=1e-8),  # controller
        diffrax.RecursiveCheckpointAdjoint(),  # adjoint
        2**8,  # max_steps
    )
    fun = beartype(jax_model.simulate_condition)

    for output in ["llh", "x0", "x", "y", "res"]:
        oargs = (*args[:-2], diffrax.DirectAdjoint(), 2**8, output)
        if sensi_order == amici.SensitivityOrder.none:
            r_jax[output] = fun(p, *oargs)[0]
        if sensi_order == amici.SensitivityOrder.first:
            if output == "llh":
                r_jax[f"s{output}"] = jax.grad(fun, has_aux=True)(p, *args)[0]
            else:
                r_jax[f"s{output}"] = jax.jacfwd(fun, has_aux=True)(p, *oargs)[
                    0
                ]

    for field in fields:
        for r_amici, r_jax in zip(rs_amici, [r_jax]):
            actual = r_jax[field]
            desired = r_amici[field]
            if field == "x":
                actual = actual[iys == 0, :]
            if field == "y":
                actual = np.stack(
                    [actual[iys == iy] for iy in sorted(np.unique(iys))],
                    axis=1,
                )
            elif field == "sllh":
                actual = actual[: len(edata.parameters)]
            elif field == "sx":
                actual = np.permute_dims(
                    actual[iys == 0, :, : len(edata.parameters)], (0, 2, 1)
                )
            elif field == "sy":
                actual = np.permute_dims(
                    np.stack(
                        [
                            actual[iys == iy, : len(edata.parameters)]
                            for iy in sorted(np.unique(iys))
                        ],
                        axis=1,
                    ),
                    (0, 2, 1),
                )
            elif field == "sx0":
                actual = actual[:, : len(edata.parameters)].T
            elif field == "sres":
                actual = actual[:, : len(edata.parameters)]

            assert_allclose(
                actual=actual,
                desired=desired,
                atol=1e-5,
                rtol=1e-5,
                err_msg=f"field {field} does not match",
            )
