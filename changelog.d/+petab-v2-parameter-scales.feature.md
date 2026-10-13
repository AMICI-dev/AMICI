The PEtab v2 simulators now support parameter scales. Since PEtab v2 no
longer specifies the scale on which parameters are estimated, the
`PetabSimulator` and the JAX `JAXProblem` so far only worked on linear
scale. Now, `"lin"`, `"log"`, or `"log10"` scales can be chosen per problem
parameter via the new `parameter_scales` argument of
`PetabImporter.create_simulator`, `ExperimentManager`, or `JAXProblem`.
Parameter values passed to the `PetabSimulator` are then expected on these
scales, and the aggregated sensitivities `sllh`, `sres`, and `s2llh` are
with respect to the scaled parameters.
`ExperimentManager.scale_parameters` and
`ExperimentManager.unscale_parameters` convert between scales. As a
consequence, `ExperimentManager.apply_parameters` now always sets
`ExpData.pscale`, overriding any parameter scale set on the `Model`.
For JAX models, `JAXProblem.parameters` (initialised to the nominal values)
and `JAXProblem.update_parameters` are on these scales, and so are
gradients computed with respect to `JAXProblem.parameters`. Scales are
preserved by `JAXProblem.save` and `JAXProblem.load` (#3311).
