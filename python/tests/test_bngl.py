import pytest
import os
import amici
import shutil
import numpy as np

pysb = pytest.importorskip("pysb")

from amici.bngl_import import bngl2amici
from pysb.simulator import ScipyOdeSimulator
from pysb.importers.bngl import model_from_bngl

tests = [
    'CaOscillate_Func', 'deleteMolecules', 'empty_compartments_block',
    'gene_expr', 'gene_expr_func', 'gene_expr_simple', 'isomerization',
    'Motivating_example_cBNGL', 'motor', 'simple_system',
    'test_compartment_XML', 'test_setconc', 'test_synthesis_cBNGL_simple',
    'test_synthesis_complex', 'test_synthesis_complex_0_cBNGL',
    'test_synthesis_complex_source_cBNGL', 'test_synthesis_simple',
    'univ_synth', 'Repressilator', 'test_paramname', 'tlmr'
]


@pytest.mark.parametrize('example', tests)
def test_compare_to_pysb_simulation(example):

    atol = 1e-8
    rtol = 1e-8

    model_file = os.path.join(os.path.dirname(__file__), '..', '..',
                              'ThirdParty', 'BioNetGen-2.3.2', 'Validate',
                              f'{example}.bngl')

    pysb_model = model_from_bngl(model_file)

    # pysb part
    tspan = np.linspace(0, 100, 101)
    sim = ScipyOdeSimulator(
        pysb_model,
        tspan=tspan,
        integrator_options={'rtol': rtol, 'atol': atol},
        compiler='python'
    )
    pysb_simres = sim.run()

    # amici part

    outdir = pysb_model.name

    cl = example not in ['Motivating_example_cBNGL', 'univ_synth']

    kwargs = {
        'compute_conservation_laws': cl,
        'observables': list(pysb_model.observables.keys())
    }
    if not cl:
        with pytest.raises(ValueError, match="Conservation laws"):
            bngl2amici(model_file, outdir, compute_conservation_laws=True)

    if example in ['empty_compartments_block', 'motor']:
        with pytest.raises(ValueError, match="Cannot add"):
            bngl2amici(model_file, outdir, **kwargs)
        return
    else:
        bngl2amici(model_file, outdir, **kwargs)

    amici_model_module = amici.import_model_module(pysb_model.name,
                                                   outdir)

    model_amici = amici_model_module.getModel()

    model_amici.setTimepoints(tspan)

    solver = model_amici.getSolver()
    solver.setMaxSteps(10**6)
    solver.setAbsoluteTolerance(atol)
    solver.setRelativeTolerance(rtol)
    rdata = amici.runAmiciSimulation(model_amici, solver)

    # check agreement of species simulation
    assert np.isclose(rdata.x, pysb_simres.species, 1e-4, 1e-4).all()

    shutil.rmtree(outdir, ignore_errors=True)
