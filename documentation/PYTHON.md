# Python Interface {#python_interface}

In the following we will give a detailed overview how to specify models in
Python and how to call the generated simulation files.

## Model Definition

This guide will guide the user on how to specify models to import and simulate
them using the Python interface.
For example implementations see the examples in the `python/examples`
directory.

### SBML import

[SBML](http://sbml.org/) is a commonly used format for specifying systems
biology models. To import an SBML model into AMICI, first, load an SBML file
using the `amici.sbml_import.SbmlImporter` class:

    import amici
    sbml_importer = amici.SbmlImporter('model_steadystate_scaled.sbml')

the SBML model as imported by [libSBML](http://sbml.org/Software/libSBML)
is available as

    sbml_model = sbml_importer.sbml

#### Constants

Model parameters that should be considered constants can be specified in a list
of strings specifying the SBML ID of the respective parameter, e.g.:

    constant_parameters=['k4']

#### Observables

Observables are specified as a dictionary with observable ID as key and
observable formula as value.

A convenient way for specifying observables for an SBML model is storing them
as `AssignmentRule`s. Assignment rules that should be considered as observables
can then be extracted using the `amici.assignmentRules2observables` function,
e.g.:

    observables = amici.assignmentRules2observables(sbml, filter_function=lambda variable: 
                                                    variable.getId().startswith('observable_') and not variable.getId().endswith('_sigma'))

#### Standard Deviations

Standard deviations can be specified as dictionaries ...

    sigmas = {'observable_x1withsigma': 'observable_x1withsigma_sigma'}


#### Model Compilation

To generate a Python module from the SBML model, call the method
`amici.sbml_import.SbmlImporter.sbml2amici`, passing all the previously defined
model specifications:

    sbml_importer.sbml2amici('test', 'test',
                             observables=observables,
                             constant_parameters=constant_parameters,
                             sigmas=sigmas)

### PySB import

[PySB](http://pysb.org/) is a tool for specifying rule-based systems biology
models as Python code. AMICI can import PySB models via
[`amici.pysb_import.pysb2amici`](https://amici.readthedocs.io/en/latest/generated/amici.pysb_import.html#amici.pysb_import.pysb2amici).

[BioNetGen](https://www.csb.pitt.edu/Faculty/Faeder/?page_id=409) and
[Kappa](https://kappalanguage.org/) models can be imported into AMICI using
PySB.

### PEtab import

[PEtab](https://github.com/PEtab-dev/PEtab) is a format for specifying
parameter estimation problems. It is based on an SBML model and tab-separated
value files specifying the observation mdoel and experimental conditions.

AMICI can import PEtab-based model definitions and run simulations for the
specified simulations conditions. For usage, see
[python/examples/example_petab/petab.ipynb](https://github.com/AMICI-dev/AMICI/blob/develop/python/examples/example_petab/petab.ipynb).

### Importing plain ODEs

The AMICI Python interface does not currently support direct import of ODEs.
However, it is straightforward to encode them as RateRules in an SBML model.
The [yaml2sbml](https://github.com/martamatos/yaml2sbml) package may come in
handy, as it facilitates generating SBML models from a YAML-based specification
of an ODE model. Besides the SBML model it can also create
[PEtab](https://github.com/PEtab-dev/PEtab) files.

## Model Simulation 

AMICI model import creates a Python module for simulation of the respective
model. To use the model module, the model directory has to be manually added to
the python path:

    import sys
    sys.path.insert(0, 'test')

the compiled model can then be imported as

    import test as model_module

It is usually more convenient to use `amici.import_model_module()` for that
purpose.

To obtain a model instance call the `getModel()` method. This model instance
will be instantiated using the default parameter values specified in the
imported model:

    model = model_module.getModel()

Specify the simulation timepoints via `amici.Model.setTimepoints`:

    model.setTimepoints(np.linspace(0, 60, 60)) 

For simulation, we need to generate a solver instance:

    solver = model.getSolver()

The model simulation can now be carried out using `amici.runAmiciSimulation`:

    rdata = amici.runAmiciSimulation(model, solver)

## Building models with OpenMP support

AMICI can be built with OpenMP support, which allows to parallelize model
simulations for multiple experimental conditions.

On Linux and OSX this is enabled by default. This can be verified using

    import amici
    amici.compiledWithOpenMP()

If not already enabled by default, you can enable OpenMP support by setting
the environment variables `AMICI_CXXFLAGS` and `AMICI_LDFLAGS` to the
correct OpenMP flags of your compiler and linker, respectively. This has to be
done for both AMICI package installation *and* model compilation. When using
`gcc` on Linux, this would be:

    # on your shell:
    AMICI_CXXFLAGS=-fopenmp AMICI_LDFLAGS=-fopenmp pip3 install amici

    # in python, before model compilation:
    import os
    os.environ['AMICI_CXXFLAGS'] = '-fopenmp'
    os.environ['AMICI_LDFLAGS'] = '-fopenmp'
