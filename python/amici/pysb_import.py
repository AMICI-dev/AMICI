from .ode_export import (
    ODEExporter, ODEModel, State, Constant, Parameter, Observable, SigmaY,
    Expression, LogLikelihood, sanitize_basic_sympy
)

import sympy as sp

try:
    import pysb.bng
    ## bool indicating whether pysb is available
    pysb_available = True
    from pysb.pattern import SpeciesPatternMatcher
except ImportError:
    pysb_available = False


def pysb2amici(model,
               output_dir=None,
               observables=None,
               constant_parameters=None,
               sigmas=None,
               verbose=False,
               assume_pow_positivity=False,
               compiler=None,
               compute_conservation_laws=True,
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

        compute_conservation_laws: if set to true, conservation laws are
        automatically computed and applied such that the state-jacobian of
        the ODE right-hand-side has full rank. This option should be set to
        True when using the newton algorithm to compute steadystates @type bool

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
        sigmas=sigmas, compute_conservation_laws=compute_conservation_laws,
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
                                observables=None, sigmas=None,
                                compute_conservation_laws=True):
    """Creates an ODEModel instance from a pysb.Model instance.


    Arguments:
        model: pysb model @type pysb.Model

        constants: list of Parameters that should be constants @type list

        observables: list of Expressions that should be observables
        @type list

        sigmas: dict with observable Expression name as key and sigma
        Expression name as expression @type list

        compute_conservation_laws: see pysb2amici


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
    if compute_conservation_laws:
        process_pysb_conservation_laws(model, ODE)
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
    xdot = sp.Matrix(model.odes)

    for ix, specie in enumerate(model.species):
        init = sp.sympify('0.0')
        for ic in model.odes.model.initial_conditions:
            if str(ic[0]) == str(specie):
                # we don't want to allow expressions in initial conditions
                if ic[1] in model.expressions:
                    init = model.expressions[ic[1].name].expand_expr()
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
            comp(par, f'{par.name}', par.value)
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
            y = exp
            ODE.add_component(
                Observable(
                    y,
                    f'{exp.name}',
                    exp.expand_expr(expand_observables=False))
            )

            sigma_name, sigma_value = get_sigma_name_and_value(
                model, exp.name, sigmas
            )

            sy = sp.Symbol(sigma_name)
            ODE.add_component(
                SigmaY(
                    sy,
                    f'{sigma_name}',
                    sigma_value
                )
            )

            my = sp.Symbol(f'm{exp.name}')
            pi = sp.sympify('pi')
            ODE.add_component(
                LogLikelihood(
                    sp.Symbol(f'llh_{exp.name}'),
                    f'llh_{exp.name}',
                    0.5*sp.log(2*pi*sy**2) + 0.5*((y - my)/sy)**2
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
                    exp,
                    f'{exp.name}',
                    exp.expand_expr(expand_observables=True))
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
        sigma_value = model.expressions[sigmas[name]].expand_expr()
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
                    obs,
                    f'{obs.name}',
                    obs.expand_obs()
                )
            )


def process_pysb_conservation_laws(model, ODE):
    """Removes species according to conservation laws to ensure that the
    jacobian has full rank


    Arguments:
        model: pysb model @type pysb.core.Model

        ODE: ODEModel instance @type ODEModel


    Returns:

    Raises:

    """

    monomers_without_conservation_law = set()
    for rule in model.rules:
        monomers_without_conservation_law |= get_unconserved_monomers(rule,
                                                                      model)

    if hasattr(model, 'initial_conditions_fixed'):
        for monomer in model.monomers:
            # check if monomer has an initial condition that is fixed (means
            #  that corresponding state is constant and all conservation
            # laws are broken)
            if any([
                model.initial_conditions_fixed[ix] # true or false
                for ix, cp in enumerate(model.initial_conditions)
                if monomer.name in extract_monomers(cp)
            ]):
                monomers_without_conservation_law |= {monomer.name}

    for monomer in model.monomers:
        monomer_species = [
            specie
            for specie in model.species
            if monomer.name in extract_monomers(specie)
        ]
        # we cannot reduce species according to conservation laws if there
        # only is (less than) a single species in the first place
        if len(monomer_species) <= 1:
            monomers_without_conservation_law |= {monomer.name}

    conservation_laws = []
    for monomer in model.monomers:
        if monomer.name not in monomers_without_conservation_law:

            compartments = [
                str(mp.compartment) # string based comparison as
                # compartments are not hashable
                for cp in model.species
                for mp in cp.monomer_patterns
                if mp.monomer.name == monomer.name
            ]

            if len(set(compartments)) > 1:
                raise Exception('Conservation laws involving species in '
                                'multiple compartments are currently not '
                                'supported! Please run pysb2amici with '
                                'compute_conservation_laws=False')
                # TODO: implement this, multiply species by the volume of
                # their respective compartment and allow total_cl to depend
                # on parameters + constants and update the respective symbolic
                # derivative accordingly


            target_index = next((
                ix
                for ix, specie in enumerate(model.species)
                if extract_monomers(specie)[0] == monomer.name
                and not ODE.state_has_conservation_law(ix)),
                None
            )
            if target_index is None:
                raise Exception(f'Cannot compute suitable conservation laws '
                                f'for this model as there are cyclic '
                                f'dependencies between states involved in '
                                f'conservation laws')
            target_expression = sanitize_basic_sympy(sum([
                sp.Symbol(f'__s{ix}')
                * extract_monomers(specie).count(monomer.name)
                for ix, specie in enumerate(model.species)
                if ix != target_index
            ]))
            target_state = sp.Symbol(f'__s{target_index}')
            conservation_laws.append({
                'state': target_state,
                'law': target_expression,
            })

    # flatten conservation laws
    unflattened_conservation_laws = \
        get_unflattened_conservation_laws(conservation_laws)
    # TODO: are circular dependencies possible? if yes, how do we
    # automatically check/prevent them?
    while len(unflattened_conservation_laws):
        for cl in conservation_laws:
            cl['law'] = cl['law'].subs(unflattened_conservation_laws)
        unflattened_conservation_laws = \
            get_unflattened_conservation_laws(conservation_laws)

    for cl in conservation_laws:
        ODE.add_conservation_law(cl['state'], cl['law'])


def get_unflattened_conservation_laws(conservation_laws):
    """Removes species according to conservation laws to ensure that the
    jacobian has full rank


    Arguments:
        model: pysb model @type pysb.core.Model

        ODE: ODEModel instance @type ODEModel


    Returns:

    Raises:

    """
    free_symbols_cl = conservation_law_variables(conservation_laws)
    return [
        (cl['state'], cl['law']) for cl in conservation_laws
        if cl['state'] in free_symbols_cl
    ]


def conservation_law_variables(conservation_laws):
    """Construct the set of all free variables from a list of conservation laws


    Arguments:
        conservation_laws: list of conservation laws (sympy.Basic)  @type list

    Returns:
    Set union of all free_symbols

    Raises:

    """
    variables = set()
    for cl in conservation_laws:
        variables |= cl['law'].free_symbols
    return variables


def has_fixed_parameter_ic(specie, model, ODE):
    """Wrapper to interface ODE.state_has_fixed_parameter_initial_condition
    from a pysb specie/model arguments

    Arguments:
        specie: pysb species @type pysb.core.ComplexPattern
        model: pysb model @type pysb.core.Model
        ODE: ODE model @type amici.ODE

    Returns:
    False if the species does not have an initial condition at all.
    Otherwise the return value of
    ODE.state_has_fixed_parameter_initial_condition

    Raises:

    """
    # ComplexPatterns are not hashable, so we have to compare by string
    ic_strs = [str(ic[0]) for ic in model.initial_conditions]
    if str(specie) not in ic_strs:
        # no initial condition at all
        return False
    else:
        return ODE.state_has_fixed_parameter_initial_condition(
            ic_strs.index(str(specie))
        )


def extract_monomers(complex_patterns):
    """Constructs a list of monomer names contained in complex patterns.
    Multiplicity of names corresponds to the stoichiometry in the complex.

    Arguments:
        specie: (list of) complex pattern(s) @type pysb.core.ComplexPattern

    Returns:
    list of monomer names

    Raises:

    """
    if not isinstance(complex_patterns, list):
        complex_patterns = [complex_patterns]
    return [
        mp.monomer.name
        for cp in complex_patterns
        if cp is not None
        for mp in cp.monomer_patterns
    ]


def get_unconserved_monomers(rule, model):
    """Constructs the set of monomer names for which the rule changes the
    stoichiometry of the monomer in the specified model.

    Arguments:
        rule: rule @type pysb.core.Rule
        model: model @type pysb.core.Model

    Returns:
    set of monomer names for which the stoichiometry is not conserved

    Raises:

    """
    unconserved_monomers = set()

    if not rule.delete_molecules \
            and len(rule.product_pattern.complex_patterns) == 0:
        # if delete_molecules is not True but we have a degradation rule,
        # we have to actually go through the reactions that are created by
        # the rule
        for reaction in [r for r in model.reactions if rule.name in r['rule']]:
            unconserved_monomers |= get_changed_stoichiometries(
                [model.species[ix] for ix in reaction['reactants']],
                [model.species[ix] for ix in reaction['products']]
            )
    else:
        # otherwise we can simply extract all information for the rule
        # itself, which is computationally much more efficient
        unconserved_monomers |= get_changed_stoichiometries(
            rule.reactant_pattern.complex_patterns,
            rule.product_pattern.complex_patterns
        )

    return unconserved_monomers


def get_changed_stoichiometries(reactants, products):
    """Constructs the set of monomer names which have different
    stoichiometries in reactants and products.

    Arguments:
        reactants: (list of) complex pattern(s) @type pysb.core.ComplexPattern
        products: (list of) complex pattern(s) @type pysb.core.ComplexPattern

    Returns:
    set of monomer name for which the stoichiometry changed

    Raises:

    """

    changed_stoichiometries = set()

    reactant_monomers = extract_monomers(
        reactants
    )

    product_monomers = extract_monomers(
        products
    )

    for monomer in set(reactant_monomers + product_monomers):
        if reactant_monomers.count(monomer) != product_monomers.count(monomer):
            changed_stoichiometries |= {monomer}

    return changed_stoichiometries

