from .ode_export import (
    ODEExporter, ODEModel, State, Constant, Parameter, Observable, SigmaY,
    Expression, LogLikelihood
)

import symengine as sp
import sympy

try:
    import pysb.bng
    ## bool indicating whether pysb is available
    pysb_available = True
except ImportError:
    pysb_available = False


def pysb2amici(model,
               output_dir=None,
               observables=None,
               constant_parameters=None,
               sigmas=None,
               verbose=False,
               assume_pow_positivity=False,
               compiler=None
               ):
    """Generate AMICI C++ files for the model provided to the constructor.

    Arguments:
        model: pysb model, model.name will determine the name of the
        generated module @type pysb.Model

        output_dir: see sbml_import.setPaths()  @type str

        observables: list of pysb.Expressions names that should be
        interpreted as observables @type list

        sigmas: list of pysb.Expressions names that should be
        interpreted as sigmas @type list

        constantParameters: list of pysb.Parameter names that should be
        interpreted as constantParameters

        verbose: more verbose output if True @type bool

        assume_pow_positivity: if set to true, a special pow function is
        used to avoid problems with state variables that may become
        negative due to numerical errors @type bool

        compiler: distutils/setuptools compiler selection to build the
        python extension @type str

    Returns:

    Raises:

    """
    if observables is None:
        observables = []

    if constant_parameters is None:
        constant_parameters = []

    if sigmas is None:
        sigmas = {}

    ode_model = ODEModel_from_pysb_importer(
        model, constants=constant_parameters, observables=observables,
        sigmas=sigmas,
    )
    exporter = ODEExporter(
        ode_model,
        outdir=output_dir,
        verbose=verbose,
        assume_pow_positivity=assume_pow_positivity,
        compiler=compiler,
    )
    exporter.setName(model.name)
    exporter.setPaths(output_dir)
    exporter.generateModelCode()
    exporter.compileModel()


def ODEModel_from_pysb_importer(model, constants=None,
                                observables=None, sigmas=None):
    """Creates an ODEModel instance from a pysb.Model instance.


    Arguments:
        model: pysb model @type pysb.Model

        constants: list of Parameters that should be constants @type list

        observables: list of Expressions that should be observables
        @type list

        sigmas: dict with observable Expression name as key and sigma
        Expression name as expression @type list


    Returns:
    New ODEModel instance according to pysbModel

    Raises:

    """

    ODE = ODEModel()

    if not pysb_available:
        raise ImportError(
            "This function requires an installation of pysb."
        )

    if constants is None:
        constants = []

    if observables is None:
        observables = []

    if sigmas is None:
        sigmas = {}

    pysb.bng.generate_equations(model)

    process_pysb_species(model, ODE)
    process_pysb_parameters(model, ODE, constants)
    process_pysb_expressions(model, ODE, observables, sigmas)
    process_pysb_observables(model, ODE)

    ODE.generateBasicVariables()

    return ODE


def process_pysb_species(model, ODE):
    """Converts pysb Species into States and adds them to the ODEModel
    instance


    Arguments:
        model: pysb model @type pysb.Model

        ODE: ODEModel instance @type ODEModel


    Returns:

    Raises:

    """
    xdot = sp.DenseMatrix(sympy.Matrix(model.odes))

    for ix, specie in enumerate(model.species):
        init = sp.sympify('0.0')
        for ic in model.odes.model.initial_conditions:
            if str(ic[0]) == str(specie):
                # we don't want to allow expressions in initial conditions
                if ic[1] in model.expressions:
                    init = sp.sympify(
                        model.expressions[ic[1].name].expand_expr()
                    )
                else:
                    init = sp.Symbol(ic[1].name)

        ODE.add_component(
            State(
                sp.Symbol(f'__s{ix}'),
                f'{specie}',
                init,
                xdot[ix])
        )


def process_pysb_parameters(model, ODE, constants):
    """Converts pysb parameters into Parameters or Constants and adds them to
    the ODEModel instance


    Arguments:
        model: pysb model @type pysb.Model

        constants: list of Parameters that should be constants @type list

        ODE: ODEModel instance @type ODEModel


    Returns:

    Raises:

    """
    for par in model.parameters:
        if par.name in constants:
            comp = Constant
        else:
            comp = Parameter

        ODE.add_component(
            comp(sp.Symbol(f'{par.name}'), f'{par.name}', par.value)
        )


def process_pysb_expressions(model, ODE, observables, sigmas):
    """Converts pysb expressions into Observables (with corresponding
    standard deviation SigmaY and LogLikelihood) or Expressions and adds them
    to the ODEModel instance

    Arguments:
        model: pysb model @type pysb.Model

        observables: list of Expressions that should be observables
        @type list

        sigmas: dict with observable Expression name as key and sigma
        Expression name as expression @type list

        ODE: ODEModel instance @type ODEModel


    Returns:

    Raises:

    """
    for exp in model.expressions:
        if exp.name in observables:
            # here we do not define a new Expression from the
            # pysb.Expression but define an observable, so we do not need
            # to expand observables as these can be defined as Expressions
            ODE.add_component(
                Observable(
                    sp.Symbol(f'{exp.name}'),
                    f'{exp.name}',
                    sp.sympify(exp.expand_expr(expand_observables=False)))
            )

            sigma_name, sigma_value = get_sigma_name_and_value(
                model, exp.name, sigmas
            )


            ODE.add_component(
                SigmaY(
                    sp.Symbol(sigma_name),
                    f'{sigma_name}',
                    sigma_value
                )
            )

            ODE.add_component(
                LogLikelihood(
                    sp.Symbol(f'llh_{exp.name}'),
                    f'llh_{exp.name}',
                    sp.sympify(
                        f'0.5*log(2*pi*{sigma_name}**2)'
                        f' + 0.5*(({exp.name} - m{exp.name})'
                        f' / {sigma_name})**2'
                    )
                )
            )

        elif exp.name in sigmas.values():
            # do nothing
            pass
        else:
            # here we do define a new Expression from the pysb.Expression
            # so we do need to expand observables as these would otherwise
            # lead to dependencies between different Expressions
            ODE.add_component(
                Expression(
                    sp.Symbol(f'{exp.name}'),
                    f'{exp.name}',
                    sp.sympify(exp.expand_expr(expand_observables=True)))
            )


def get_sigma_name_and_value(model, name, sigmas):
    """Tries to extract standard deviation symbolic identifier and formula
    for a given observable name from the pysb model and if no specification is
    available sets default values


    Arguments:
        model: pysb model @type pysb.Model

        name: name of the observable

        sigmas: dict with observable Expression name as key and sigma
        Expression name as expression @type list


    Returns:
    tuple containing symbolic identifier and formula for the specified
    observable

    Raises:

    """
    if name in sigmas:
        if sigmas[name] not in model.expressions:
            raise Exception(f'value of sigma {name} is not a '
                            f'valid expression.')
        sigma_name = model.expressions[sigmas[name]].name
        sigma_value = sp.sympify(
            model.expressions[sigmas[name]].expand_expr()
        )
    else:
        sigma_name = f'sigma_{name}'
        sigma_value = sp.sympify(1.0)

    return sigma_name, sigma_value


def process_pysb_observables(model, ODE):
    """Converts pysb observables into Expressions and adds them to the ODEModel
    instance


    Arguments:
        model: pysb model @type pysb.Model

        ODE: ODEModel instance @type ODEModel


    Returns:

    Raises:

    """
    # only add those pysb observables that occur in the added
    # Observables as expressions
    for obs in model.observables:
        if sp.Symbol(obs.name) in ODE.eq('y').free_symbols:
            ODE.add_component(
                Expression(
                    sp.Symbol(f'{obs.name}'),
                    f'{obs.name}',
                    sp.sympify(obs.expand_obs()))
            )