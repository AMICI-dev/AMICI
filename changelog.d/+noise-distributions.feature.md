Models now store the noise distribution of each observable. They are
available via `Model.get_noise_distributions()` and
`Model.get_event_noise_distributions()` as `NoiseDistribution` values
(`NoiseDistribution.custom` for user-provided cost functions), and per
observable via `Model.get_noise_distribution(iy)` and
`Model.get_event_noise_distribution(iz)` (as integers).
`Model.is_quadratic_llh(edata)` tells whether residuals and the
Gauss-Newton FIM are computed
([#3315](https://github.com/AMICI-dev/AMICI/issues/3315)).
