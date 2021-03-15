.. _python_interface:

******************************
Using AMICI's Python interface
******************************

In the following we will give a detailed overview how to specify models in
Python and how to call the generated simulation files.

Model definition
================

This document provides an overview of different interfaces to import models
in AMICI. Further examples are available in the AMICI repository in the
`python/examples <https://github.com/AMICI-dev/AMICI/tree/master/python/examples>`_
directory.

SBML import
-----------

AMICI can import :term:`SBML` models via the
:py:func:`amici.sbml_import.SbmlImporter` class.

.. _amici_python_sbml_support:

Status of SBML support in Python-AMICI
++++++++++++++++++++++++++++++++++++++

Python-AMICI currently **passes 862 out of the 1780 (~48%) test cases** from
the semantic
`SBML Test Suite <https://github.com/sbmlteam/sbml-test-suite/>`_
(`current status <https://github.com/AMICI-dev/AMICI/actions>`_).

The following SBML test suite tags are currently supported
(i.e., at least one test case with the respective test passes;
`tag descriptions <https://github.com/sbmlteam/sbml-test-suite/blob/master/docs/tags-documentation/all-tags.txt>`_):

**Component tags:**

* AssignmentRule
* Compartment
* CSymbolAvogadro
* CSymbolTime
* FunctionDefinition
* InitialAssignment
* Parameter
* RateRule
* Reaction
* Species

**Test tags:**

* 0D-Compartment
* Amount
* AssignedConstantStoichiometry
* AssignedVariableStoichiometry
* BoolNumericSwap
* BoundaryCondition
* Concentration
* ConstantSpecies
* ConversionFactors
* HasOnlySubstanceUnits
* InitialValueReassigned
* L3v2MathML
* LocalParameters
* MultiCompartment
* NoMathML
* NonConstantCompartment
* NonConstantParameter
* NonUnityCompartment
* NonUnityStoichiometry
* ReversibleReaction
* SpeciesReferenceInMath
* UncommonMathML
* VolumeConcentrationRates

In addition, we currently plan to add support for the following features
(see corresponding `issues <https://github.com/AMICI-dev/AMICI/milestone/14>`_
for details and progress):

- Events (currently Matlab-only) (`#757 <https://github.com/AMICI-dev/AMICI/issues/757>`_)
- Algebraic rules (`#760 <https://github.com/AMICI-dev/AMICI/issues/760>`_)

However, the following features are unlikely to be supported:

- any SBML extensions
- `factorial()`, `ceil()`, `floor()`, due to incompatibility with
  symbolic sensitivity computations
- `delay()` due to missing :term:`SUNDIALS` solver support

Tutorials
+++++++++

A basic tutorial on how to import and simulate SBML models is available in the
`Getting Started notebook <GettingStarted.ipynb>`_, while a more detailed example
including customized import and sensitivity computation is available in the
`Example Steadystate notebook <ExampleSteadystate.ipynb>`_.

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
The `yaml2sbml <https://github.com/yaml2sbml-dev/yaml2sbml>`_ package may come in
handy, as it facilitates generating SBML models from a YAML-based specification
of an ODE model. Besides the SBML model it can also create
`PEtab <https://github.com/PEtab-dev/PEtab>`_ files.

SED-ML import
-------------

We also plan to implement support for the
`Simulation Experiment Description Markup Language (SED-ML) <https://sed-ml.org/>`_.

Examples
========

.. toctree::
   :maxdepth: 1

   GettingStarted.ipynb
   ExampleSteadystate.ipynb
   petab.ipynb
   ExampleExperimentalConditions.ipynb
   ExampleEquilibrationLogic.ipynb


Miscellaneous
=============

.. _amici_python_openmp:

OpenMP support for parallelized simulation for multiple experimental conditions
-------------------------------------------------------------------------------

AMICI can be built with OpenMP support, which allows to parallelize model
simulations for multiple experimental conditions.

On Linux and OSX this is enabled by default. This can be verified using:

.. code-block:: python

   import amici
   amici.compiledWithOpenMP()

If not already enabled by default, you can enable OpenMP support by setting
the environment variables ``AMICI_CXXFLAGS`` and ``AMICI_LDFLAGS`` to the
correct OpenMP flags of your compiler and linker, respectively. This has to be
done for both AMICI package installation *and* model compilation. When using
``gcc`` on Linux, this would be:

.. code-block:: bash

   # on your shell:
   AMICI_CXXFLAGS=-fopenmp AMICI_LDFLAGS=-fopenmp pip3 install amici

.. code-block:: python

   # in python, before model compilation:
   import os
   os.environ['AMICI_CXXFLAGS'] = '-fopenmp'
   os.environ['AMICI_LDFLAGS'] = '-fopenmp'
