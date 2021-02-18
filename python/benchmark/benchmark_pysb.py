import os
import pysb
import amici
from amici.pysb_import import pysb2amici
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import importlib
import timeit
import sys

from pysb.simulator import ScipyOdeSimulator

sys.path.insert(0, os.path.join('..', 'tests'))
from test_pysb import pysb_models

simulation_times = dict()

N_REPEATS = 100

atol = 1e-8
rtol = 1e-8

for example in pysb_models:
    simulation_times[example] = dict()

    with amici.add_path(os.path.dirname(pysb.examples.__file__)):
        with amici.add_path(os.path.join(os.path.dirname(__file__), '..',
                                         'tests', 'pysb_test_models')):

            pysb.SelfExporter.cleanup()  # reset pysb
            pysb.SelfExporter.do_export = True

            module = importlib.import_module(example)
            pysb_model = module.model
            pysb_model.name = pysb_model.name.replace('pysb.examples.', '')
            # avoid naming clash for custom pysb models
            pysb_model.name += '_amici'

            # pysb part
            tspan = np.linspace(0, 100, 101)
            sim = ScipyOdeSimulator(
                pysb_model,
                tspan=tspan,
                integrator_options={'rtol': rtol, 'atol': atol},
            )
            time_pysb = timeit.Timer(
                'pysb_simres = sim.run()',
                globals={'sim': sim}
            ).timeit(number=N_REPEATS)/N_REPEATS

            simulation_times[example]['pysb'] = time_pysb
            print(f'PySB average simulation time {example}: {time_pysb}')

            # amici part
            outdir = pysb_model.name

            if pysb_model.name in ['move_connected_amici']:
                compute_conservation_laws = False
            else:
                compute_conservation_laws = True

            pysb2amici(
                pysb_model,
                outdir,
                compute_conservation_laws=compute_conservation_laws,
                observables=list(pysb_model.observables.keys())
            )

            amici_model_module = amici.import_model_module(pysb_model.name,
                                                           outdir)

            model_pysb = amici_model_module.getModel()

            model_pysb.setTimepoints(tspan)

            solver = model_pysb.getSolver()
            solver.setMaxSteps(int(1e6))
            solver.setAbsoluteTolerance(atol)
            solver.setRelativeTolerance(rtol)
            time_amici = timeit.Timer(
                'rdata = amici.runAmiciSimulation(model, solver)',
                globals={'model': model_pysb, 'solver': solver,
                         'amici': amici}
            ).timeit(number=N_REPEATS)/N_REPEATS
            simulation_times[example]['amici'] = time_amici
            print(f'AMICI average simulation time {example}: {time_amici}')

times = pd.DataFrame(simulation_times)

ax = times.T.plot(kind='scatter', x='pysb', y='amici')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_aspect('equal')
xy_min = np.min([ax.get_xlim()[0], ax.get_ylim()[0]])
xy_max = np.max([ax.get_xlim()[1], ax.get_ylim()[1]])

ax.set_xlim([xy_min, xy_max])
ax.set_ylim([xy_min, xy_max])
ax.set_ylabel('simulation time AMICI [s]')
ax.set_xlabel('simulation time PySB [s]')
ax.plot([xy_min, xy_max], [xy_min, xy_max], 'k:')
plt.tight_layout()
plt.savefig('benchmark_pysb.eps')
