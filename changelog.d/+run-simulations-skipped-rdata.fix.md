Fixed `run_simulations` with `failfast=True` returning results for the
skipped conditions that had the model's timepoints and parameter list instead
of the `ExpData`'s, and no `id`. Results of skipped conditions and of
conditions that raised an error now carry the `ExpData` `id`, and their
states, observables, residuals, likelihood and respective sensitivities are
NaN instead of 0.0, as for failed simulations (#3456).
