import os

import numpy as np
import pytest

pysb = pytest.importorskip("pysb")

from contextlib import suppress

from amici import import_model_module
from amici.importers.bngl import bngl2amici
from amici.sim.sundials import run_simulation
from amici.testing import TemporaryDirectoryWinSafe, skip_on_valgrind
from numpy.testing import assert_allclose
from pysb.importers.bngl import model_from_bngl
from pysb.simulator import ScipyOdeSimulator

tests = [
    "CaOscillate_Func",
    "deleteMolecules",
    "empty_compartments_block",
    "gene_expr",
    "gene_expr_func",
    "gene_expr_simple",
    "isomerization",
    "Motivating_example_cBNGL",
    "motor",
    "simple_system",
    "test_compartment_XML",
    "test_setconc",
    "test_synthesis_cBNGL_simple",
    "test_synthesis_complex",
    "test_synthesis_complex_0_cBNGL",
    "test_synthesis_complex_source_cBNGL",
    "test_synthesis_simple",
    "univ_synth",
    "Repressilator",
    "test_paramname",
    "tlmr",
]


@skip_on_valgrind
@pytest.mark.parametrize("example", tests)
def test_compare_to_pysb_simulation(example):
    from amici.importers.utils import RESERVED_SYMBOLS, MeasurementChannel

    # allow "NULL" as model symbol
    # (used in CaOscillate_Func and Repressilator examples)
    with suppress(ValueError):
        RESERVED_SYMBOLS.remove("NULL")

    atol = 1e-12
    rtol = 1e-10

    model_file = os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "ThirdParty",
        "BioNetGen-2.7.0",
        "Validate",
        f"{example}.bngl",
    )

    pysb_model = model_from_bngl(model_file)

    # pysb part
    tspan = np.linspace(0, 100, 101)
    sim = ScipyOdeSimulator(
        pysb_model,
        tspan=tspan,
        integrator_options={"rtol": rtol, "atol": atol},
        compiler="python",
    )
    pysb_simres = sim.run()

    # amici part
    cl = example not in ["Motivating_example_cBNGL", "univ_synth"]

    kwargs = {
        "compute_conservation_laws": cl,
        "observation_model": list(
            map(MeasurementChannel, pysb_model.observables.keys())
        ),
    }

    with TemporaryDirectoryWinSafe(prefix=pysb_model.name) as outdir:
        if not cl:
            with pytest.raises(ValueError, match="Conservation laws"):
                bngl2amici(
                    model_file,
                    output_dir=outdir,
                    compute_conservation_laws=True,
                )

        if example in ["empty_compartments_block", "motor"]:
            with pytest.raises(ValueError, match="Cannot add"):
                bngl2amici(model_file, output_dir=outdir, **kwargs)
            return
        else:
            bngl2amici(model_file, output_dir=outdir, **kwargs)

        amici_model_module = import_model_module(pysb_model.name, outdir)

        model_amici = amici_model_module.get_model()

        model_amici.set_timepoints(tspan)

        solver = model_amici.create_solver()
        solver.set_max_steps(10**6)
        solver.set_absolute_tolerance(atol)
        solver.set_relative_tolerance(rtol)
        rdata = run_simulation(model_amici, solver)

        # check agreement of species simulation
        assert_allclose(rdata.x, pysb_simres.species, 1e-4, 1e-4)
