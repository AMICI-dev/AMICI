#!/usr/bin/env python3

import amici
import sys
import petab
import subprocess
import os
import re
import shutil

from amici.petab_import import import_model


def main():
    arg = sys.argv[1]

    if '_' in arg and re.match(r'O[0-2]', arg.split("_")[-1]):
        optim = arg.split("_")[-1]
        os.environ['AMICI_CXXFLAGS'] = f'-{optim}'
        suffix = f'_{optim}'
        arg = '_'.join(arg.split("_")[:-1])
    else:
        suffix = ''

    model_dir = os.path.join(os.curdir, 'CS_Signalling_ERBB_RAS_AKT',
                             'CS_Signalling_ERBB_RAS_AKT_petab' + suffix)
    model_name = 'CS_Signalling_ERBB_RAS_AKT_petab'
    if arg == 'import':
        git_dir = os.path.join(os.curdir, 'CS_Signalling_ERBB_RAS_AKT')
        if not os.path.exists(git_dir):
            subprocess.run([
                'git', 'clone', '--depth', '1',
                'https://github.com/ICB-DCM/CS_Signalling_ERBB_RAS_AKT']
            )
        os.chdir(os.path.join(os.curdir, 'CS_Signalling_ERBB_RAS_AKT'))

        pp = petab.Problem.from_yaml(
            'FroehlichKes2018/PEtab/FroehlichKes2018.yaml'
        )
        petab.lint_problem(pp)
        import_model(model_name=model_name,
                     sbml_model=pp.sbml_model,
                     condition_table=pp.condition_df,
                     observable_table=pp.observable_df,
                     measurement_table=pp.measurement_df,
                     compile=False,
                     verbose=True)

        return
    elif arg == 'compile':
        if model_name != os.path.basename(model_dir):
            shutil.copytree(
                os.path.join(os.curdir, 'CS_Signalling_ERBB_RAS_AKT',
                             model_name),
                model_dir
            )

        subprocess.run(['python', 'setup.py',
                        'build_ext', f'--build-lib=.', '--force'],
                       cwd=model_dir)
        return
    else:
        model_module = amici.import_model_module(model_name, model_dir)
        model = model_module.getModel()
        solver = model.getSolver()
        # TODO
        edata = amici.ExpData(model)
        edata.setTimepoints([1e8])
        edata.setObservedData([1.0])
        edata.setObservedDataStdDev([1.0])

    if arg == 'forward_simulation':
        solver.setSensitivityMethod(amici.SensitivityMethod.none)
        solver.setSensitivityOrder(amici.SensitivityOrder.none)
    elif arg == 'forward_sensitivities':
        model.setParameterList(list(range(100)))
        solver.setSensitivityMethod(amici.SensitivityMethod.forward)
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
    elif arg == 'adjoint_sensitivities':
        solver.setSensitivityMethod(amici.SensitivityMethod.adjoint)
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
    elif arg == 'forward_simulation_non_optimal_parameters':
        tmpPar = model.getParameters()
        model.setParameters([0.1 for _ in tmpPar])
        solver.setSensitivityMethod(amici.SensitivityMethod.none)
        solver.setSensitivityOrder(amici.SensitivityOrder.none)
    elif arg == 'adjoint_sensitivities_non_optimal_parameters':
        tmpPar = model.getParameters()
        model.setParameters([0.1 for _ in tmpPar])
        solver.setSensitivityMethod(amici.SensitivityMethod.adjoint)
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
    elif arg == 'forward_steadystate_sensitivities_non_optimal_parameters':
        tmpPar = model.getParameters()
        model.setParameters([0.1 for _ in tmpPar])
        solver.setSensitivityMethod(amici.SensitivityMethod.forward)
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
        edata.setTimepoints([float('inf')])
    elif arg == 'adjoint_steadystate_sensitivities_non_optimal_parameters':
        tmpPar = model.getParameters()
        model.setParameters([0.1 for _ in tmpPar])
        solver.setSensitivityMethod(amici.SensitivityMethod.adjoint)
        solver.setSensitivityOrder(amici.SensitivityOrder.first)
        edata.setTimepoints([float('inf')])
    else:
        print("Unknown argument:", arg)
        sys.exit(1)
    rdata = amici.runAmiciSimulation(model, solver, edata)

    diagnostics = ['numsteps', 'numstepsB', 'numrhsevals', 'numrhsevalsB',
                   'numerrtestfails', 'numerrtestfailsB',
                   'numnonlinsolvconvfails', 'numnonlinsolvconvfailsB',
                   'preeq_cpu_time', 'preeq_cpu_timeB',
                   'cpu_time', 'cpu_timeB',
                   'posteq_cpu_time', 'posteq_cpu_timeB']
    for d in diagnostics:
        print(d, rdata[d])
    assert rdata['status'] == amici.AMICI_SUCCESS


if __name__ == '__main__':
    main()
