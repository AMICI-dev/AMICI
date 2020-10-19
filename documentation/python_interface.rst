.. _python_interface:

******************************
Using AMICI's Python interface
******************************

In the following we will give a detailed overview how to specify models in
Python and how to call the generated simulation files.

Model definition
================

This guide will guide the user on how to specify models to import and simulate
them using the Python interface. Further examples are available in the AMICI
repository in the
`python/examples <https://github.com/AMICI-dev/AMICI/tree/master/python/examples>`_
directory.

SBML import
-----------

AMICI can import :term:`SBML` models via the
:py:func:`amici.sbml_import.SbmlImporter` class.

Status of SBML support in Python-AMICI
++++++++++++++++++++++++++++++++++++++

Python-AMICI currently **passes 756 out of the 1781 (~42%) test cases** from
the semantic
`SBML Test Suite <https://github.com/sbmlteam/sbml-test-suite/>`_
(`current status <https://github.com/AMICI-dev/AMICI/actions>`_).

In addition, we currently plan to add support for the following features
(see corresponding issues for details and progress):

- Events (currently Matlab-only)
- Algebraic rules
- Models without species

contributions are welcome.

However, the following features are unlikely to be supported:

- SBML extensions
- `factorial()`, `ceil()`, `floor()`, due to incompatibility with
  symbolic sensitivity computations
- initial assignments for parameters
- `delay()` due to missing SUNDIALS solver support


How to import an SBML model into AMICI
++++++++++++++++++++++++++++++++++++++

To import an :term:`SBML` model into AMICI, first, load an SBML file
using the :py:func:`amici.sbml_import.SbmlImporter` class::

    import amici
    sbml_importer = amici.SbmlImporter('model_steadystate_scaled.sbml')

the SBML model as imported by `libSBML <http://sbml.org/Software/libSBML>`_
is available as::

    sbml_model = sbml_importer.sbml

Constants
^^^^^^^^^

Model parameters that should be considered :term:`constants <fixed parameters>`
can be specified in a list
of strings specifying the SBML ID of the respective parameter, e.g.::

    constant_parameters=['k4']

Observables
^^^^^^^^^^^

Observables are specified as a dictionary with observable ID as key and
observable formula as value.

A convenient way for specifying observables for an SBML model is storing them
as ``AssignmentRule``\ s. Assignment rules that should be considered as observables
can then be extracted using the :py:func:`amici.sbml_import.assignmentRules2observables`
function, e.g.::

    observables = amici.assignmentRules2observables(sbml, filter_function=lambda variable:
                                                    variable.getId().startswith('observable_') and not variable.getId().endswith('_sigma'))

Standard deviations
^^^^^^^^^^^^^^^^^^^

Standard deviations can be specified as dictionaries, such as::

    sigmas = {'observable_x1withsigma': 'observable_x1withsigma_sigma'}

Noise distributions
^^^^^^^^^^^^^^^^^^^

Various noise distributions including normal and Laplace and discrete
distributions, and scale transformations including linear, log and log10
are supported::

    noise_distributions = {'observable_x1withsigma': 'log-normal'}

Find details in :py:func:`amici.sbml_import.noise_distribution_to_cost_function`.

Model compilation
^^^^^^^^^^^^^^^^^

To generate a Python module from the SBML model, call the method
:py:func:`amici.sbml_import.SbmlImporter.sbml2amici`, passing all the
previously defined model specifications::

    sbml_importer.sbml2amici('test', 'test',
                             observables=observables,
                             constant_parameters=constant_parameters,
                             sigmas=sigmas)

Full example
^^^^^^^^^^^^

See `here <ExampleSteadystate.ipynb>`_ for a full example.

PySB import
-----------

AMICI can import :term:`PySB` models via
:py:func:`amici.pysb_import.pysb2amici`.

`BioNetGen <https://www.csb.pitt.edu/Faculty/Faeder/?page_id=409>`_ and
`Kappa <https://kappalanguage.org/>`_ models can be imported into AMICI using
PySB.

PEtab import
------------

AMICI can import :term:`PEtab`-based model definitions and run simulations for
the specified simulations conditions. For usage, see
`python/examples/example_petab/petab.ipynb <petab.ipynb>`_.

Importing plain ODEs
--------------------

The AMICI Python interface does not currently support direct import of ODEs.
However, it is straightforward to encode them as RateRules in an SBML model.
The `yaml2sbml <https://github.com/martamatos/yaml2sbml>`_ package may come in
handy, as it facilitates generating SBML models from a YAML-based specification
of an ODE model. Besides the SBML model it can also create
`PEtab <https://github.com/PEtab-dev/PEtab>`_ files.

SED-ML import
-------------

We also plan to implement support for the `Simulation Experiment Description Markup Language (SED-ML) <https://sed-ml.org/>`_.

Model simulation
================

AMICI model import creates a Python module for simulation of the respective
model. To use the model module, the model directory has to be manually added to
the python path::

    import sys
    sys.path.insert(0, 'test')

the compiled model can then be imported as::

    import test as model_module

It is usually more convenient to use :py:func:`amici.import_model_module` for
that purpose.

To obtain a model instance call the `getModel()` method. This model instance
will be instantiated using the default parameter values specified in the
imported model::

    model = model_module.getModel()

Specify the simulation timepoints via :py:func:`amici.Model.setTimepoints`::

    model.setTimepoints(np.linspace(0, 60, 60))

For simulation, we need to generate a solver instance::

    solver = model.getSolver()

The model simulation can now be carried out using
:py:func:`amici.runAmiciSimulation`::

    rdata = amici.runAmiciSimulation(model, solver)


Examples
========

.. toctree::
   :maxdepth: 1

   ExampleSteadystate.ipynb
   petab.ipynb
   model_presimulation.ipynb
   ExampleEquilibrationLogic.ipynb


Miscellaneous
=============

.. _amici_python_openmp:

OpenMP support for parallelized simulation for multiple experimental conditions
-------------------------------------------------------------------------------

AMICI can be built with OpenMP support, which allows to parallelize model
simulations for multiple experimental conditions.

On Linux and OSX this is enabled by default. This can be verified using::

    import amici
    amici.compiledWithOpenMP()

If not already enabled by default, you can enable OpenMP support by setting
the environment variables ``AMICI_CXXFLAGS`` and ``AMICI_LDFLAGS`` to the
correct OpenMP flags of your compiler and linker, respectively. This has to be
done for both AMICI package installation *and* model compilation. When using
``gcc`` on Linux, this would be::

    # on your shell:
    AMICI_CXXFLAGS=-fopenmp AMICI_LDFLAGS=-fopenmp pip3 install amici

    # in python, before model compilation:
    import os
    os.environ['AMICI_CXXFLAGS'] = '-fopenmp'
    os.environ['AMICI_LDFLAGS'] = '-fopenmp'
