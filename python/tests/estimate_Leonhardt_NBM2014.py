import pypesto
import pypesto.petab
import pypesto.optimize as optimize
import pypesto.visualize as visualize
import amici
import fides
import petab

import os
import numpy as np
import matplotlib.pyplot as plt

folder_base = "/home/paulstapor/Dokumente/Projekte/adjoints_and_events/petab_models"

# a collection of models that can be simulated
model_name = "Leonhardt_NBM2014"

# the yaml configuration file links to all needed files
yaml_config = os.path.join(folder_base, model_name, model_name + '.yaml')

# create a petab and the pypesto problem
petab_problem = petab.Problem.from_yaml(yaml_config)
importer = pypesto.petab.PetabImporter.from_yaml(yaml_config)
problem = importer.create_problem()

objective = problem.objective
x = petab_problem.x_nominal_free_scaled

# forward, do the optimization
optimizer = optimize.FidesOptimizer(hessian_update=fides.BFGS())
engine = pypesto.engine.MultiThreadEngine()
objective.amici_solver.setSensitivityMethod(1)
result_fsa = optimize.minimize(problem=problem, optimizer=optimizer,
                               n_starts=20, engine=engine)
fname = f'opimization_Leonhardt_NBM2014_Fides_BFGS_Forward.h5'
pypesto.store.write_result(result_fsa,
                           filename=fname,
                           problem=True, optimize=True,
                           sample=False, profile=False)

fvals_fsa = result_fsa.optimize_result.get_for_key('fval')
xs_fsa = result_fsa.optimize_result.get_for_key('x')
ret_fsa = objective(xs_fsa[0], sensi_orders=(0,1))

visualize.waterfall(result_fsa, scale_y='lin')
visualize.parameters(result_fsa)

# some sanity checking
objective = problem.objective
x = petab_problem.x_nominal_free_scaled
print('\n\n\n\n===========AMICI objective gradient check===========')
df = objective.check_grad_multi_eps(x)
print('\n___________RESULT____________')
print(df)

print('\n\n\n\n===========AMICI objective call with sensi (0,1)===========')
call = objective(x, sensi_orders=(0, 1), return_dict=True)
print('forward:', call)

print('\n\n\n\n===========AMICI objective call with sensi (0,1)===========')
call = objective(x, sensi_orders=(0, 1), return_dict=True)
objective.amici_solver.setSensitivityMethod(2)
print('adjoint:', call)

print('\n\n\n\n===========AMICI objective call with sensi (0,1)===========')
call = objective(x, sensi_orders=(0, 1), return_dict=True)
objective.amici_solver.setSensitivityMethod(2)
objective.amici_solver.setInterpolationType(2)
print('adjoint, polynomial interpolation:', call)

print('\n\n\n\n===========AMICI objective call with sensi (0,1)===========')
call = objective(x, sensi_orders=(0, 1), return_dict=True)
objective.amici_solver.setSensitivityMethod(2)
objective.amici_solver.setInterpolationType(2)
objective.amici_solver.setSensitivityMethod(2)
objective.amici_solver.setInterpolationType(2)
objective.amici_solver.setRelativeTolerance(1e-12)
objective.amici_solver.setRelativeToleranceB(1e-12)
objective.amici_solver.setAbsoluteToleranceQuadratures(1e-12)
objective.amici_solver.setRelativeToleranceQuadratures(1e-10)
print('adjoint, polynomial interpolation, high acc.:', call)

optimizer = optimize.FidesOptimizer(hessian_update=fides.BFGS())
engine = pypesto.engine.MultiThreadEngine()

# adjoint, do the optimization
objective.amici_solver.setSensitivityMethod(2)
objective.amici_solver.setInterpolationType(2)
objective.amici_solver.setRelativeTolerance(1e-12)
objective.amici_solver.setRelativeToleranceB(1e-12)
objective.amici_solver.setAbsoluteToleranceQuadratures(1e-12)
objective.amici_solver.setRelativeToleranceQuadratures(1e-10)
result_asa = optimize.minimize(problem=problem, optimizer=optimizer,
                               n_starts=100, engine=engine)
fname = f'opimization_Leonhardt_NBM2014_Fides_BFGS_Adjoint_strict.h5'
pypesto.store.write_result(result_asa,
                           filename=fname,
                           problem=True, optimize=True,
                           sample=False, profile=False)

fvals_asa = result_asa.optimize_result.get_for_key('fval')
xs_asa = result_asa.optimize_result.get_for_key('x')
ret_asa = objective(xs_asa[0], sensi_orders=(0,1))


visualize.waterfall(result_asa, scale_y='log')
visualize.parameters(result_asa)
plt.show()

# forward, do the optimization
optimizer = optimize.FidesOptimizer(hessian_update=fides.BFGS())
engine = pypesto.engine.MultiThreadEngine()
objective.amici_solver.setSensitivityMethod(1)
result_fsa = optimize.minimize(problem=problem, optimizer=optimizer,
                               n_starts=100, engine=engine)
fname = f'opimization_Leonhardt_NBM2014_Fides_BFGS_Forward.h5'
pypesto.store.write_result(result_fsa,
                           filename=fname,
                           problem=True, optimize=True,
                           sample=False, profile=False)

fvals_fsa = result_fsa.optimize_result.get_for_key('fval')
xs_fsa = result_fsa.optimize_result.get_for_key('x')
ret_fsa = objective(xs_fsa[0], sensi_orders=(0,1))

visualize.waterfall(result_fsa, scale_y='lin')
visualize.parameters(result_fsa)

