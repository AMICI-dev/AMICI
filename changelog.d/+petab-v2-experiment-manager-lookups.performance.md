`amici.sim.sundials.petab.ExperimentManager` and `PetabSimulator` now look
up experiments, conditions, measurements, and parameters in constant time.
Before, creating the `ExpData` objects for a PEtab v2 problem was quadratic
in the number of experiments and in the number of parameters, which made it
impractically slow for large problems
([#3519](https://github.com/AMICI-dev/AMICI/issues/3519)).
