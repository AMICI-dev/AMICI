# Python Interface {#python_interface}

In the following we will give a detailed overview how to specify models in  Python and how to call the generated simulation files.

## Model Definition

This guide will guide the user on how to specify models in Python using SBML. For example implementations see the examples in the python/examples directory.

### SBML input

First, import an sbml file using the `amici.sbml_import.SbmlImporter` class:

    import amici
    sbmlImporter = amici.SbmlImporter('model_steadystate_scaled.sbml')
    
the sbml document as imported by [libSBML](http://sbml.org/Software/libSBML) is available as 

    sbml = sbmlImporter.sbml

### Constants 

parameters that should be considered constants can be specified in a list of strings specifying the respective SbmlId of a parameter.

    constantParameters=['k4']

### Observables

assignment rules that should be considered as observables can extracted using the `amici.assignmentRules2observables` function

    observables = amici.assignmentRules2observables(sbml, filter=lambda variableId: 
                                                    variableId.startswith('observable_') and not variableId.endswith('_sigma'))

### Standard Deviations

standard deviations can be specified as dictionaries ...

    sigmas = {'observable_x1withsigma': 'observable_x1withsigma_sigma'}


## Model Compilation

to compile the sbml as python module, the user has to call the method `amici.sbml_import.SbmlImporter.sbml2amici`, passing all the previously defined model specifications

    sbmlImporter.sbml2amici('test', 'test', 
                            observables=observables,
                            constantParameters=constantParameters,
                            sigmas=sigma)

## Model Simulation 

currently the model folder has to be manually added to the python path
    
    import sys
    sys.path.insert(0, 'test')
    
the compiled model can now be imported as python module
    
    import test as modelModule

to obtain a model instance call the `getModel()` method. This model instance will be instantiated using the defautl parameter values specified in the sbml.

    model = modelModule.getModel()

then pass the simulation timepoints as `amici.DoubleVector` to `amici.Model.setTimepoints`

    model.setTimepoints(amici.DoubleVector(np.linspace(0, 60, 60))) 
    
for simulation we need to generate a solver instance 

    solver = model.getSolver()
    
the model simulation can now be carried out using `amici.runAmiciSimulation`
    
    rdata = amici.runAmiciSimulation(model, solver)
