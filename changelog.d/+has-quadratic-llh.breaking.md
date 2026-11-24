When subclassing `amici::Model` in C++, override
`Model::get_noise_distribution` (and `get_event_noise_distribution`)
instead of `has_quadratic_llh` or `get_observable_scaling`. Both are no
longer virtual: overrides marked `override` no longer compile, others are
silently ignored. The undocumented
`amici.importers.utils.ObservableTransformation` and
`amici.importers.utils.noise_distribution_to_observable_transformation`
were removed
([#3315](https://github.com/AMICI-dev/AMICI/issues/3315)).
