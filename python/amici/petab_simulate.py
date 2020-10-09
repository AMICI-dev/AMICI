from amici import AmiciModel
from amici.petab_import import import_petab_problem
from amici.petab_objective import simulate_petab, rdatas_to_measurement_df
from petab.simulate import Simulator


class PetabSimulator(Simulator):
    """
    amici_model:
        Amici model instance.
    """
    def _simulate_without_noise(
            self,
            amici_model: AmiciModel = None,
            **kwargs
    ):
        """
        See :py:func:`petab.simulate.Simulator.simulate()` docstring.

        :param amici_model:
            AMICI Model assumed to be compatible with ``petab_problem``.
        """
        if amici_model is None and 'amici_model' not in dir(self):
            if 'model_output_dir' not in kwargs:
                import os
                kwargs['model_output_dir'] = os.path.join(self.working_dir,
                                                          'amici_models')
            relevant_kwargs_import = [
                'model_output_dir',
                'model_name',
                'force_compile',
            ]
            # TODO don't compute sensitivities
            # TODO allow amici_model to be passed as argument
            self.amici_model = import_petab_problem(
                self.petab_problem,
                **{k: v
                   for k, v in kwargs.items()
                   if k in relevant_kwargs_import},
            )
        else:
            # TODO don't save kwarg amici_model in `self` state?
            self.amici_model = amici_model

        # due to the `simulate_petab` decorator (I think), the inspect module
        # cannot be used to identify this automatically with e.g.:
        # inspect.getfullargspec(simulate_petab).args
        # `amici_model` is skipped here, as it is handled separately
        relevant_kwargs_simulate = [
            'solver',
            'problem_parameters',
            'simulation_conditions',
            'edatas',
            'parameter_mapping',
            'scaled_parameters',
            'log_level',
        ]

        # TODO allow specification of solver?
        result = simulate_petab(
            self.petab_problem,
            self.amici_model,
            **{k: v
               for k, v in kwargs.items()
               if k in relevant_kwargs_simulate},
        )

        # TODO use `rdatas_to_simulation_df` instead?
        simulation_df = rdatas_to_measurement_df(
            result['rdatas'],
            self.amici_model,
            self.petab_problem.measurement_df,
        )

        return simulation_df
