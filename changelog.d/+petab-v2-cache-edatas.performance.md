`amici.sim.sundials.petab.PetabSimulator.simulate` now creates the `ExpData`
objects for a PEtab v2 problem only on its first call, and afterwards only
updates their parameter-dependent values. Before, all `ExpData` objects were
recreated on every call, which took about half of an objective function
evaluation for problems with many experiments
([#3521](https://github.com/AMICI-dev/AMICI/issues/3521)).
`ExperimentManager.create_edatas` now accepts the problem parameters to apply.
