from .ode_export import (
    ODEExporter, ODEModel, State, Constant, Parameter, Observable, SigmaY,
    Expression, LogLikelihood, sanitize_basic_sympy
)

import sympy as sp
import numpy as np
import itertools


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

    ode_model = _ODEModel_from_pysb_importer(
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


def _ODEModel_from_pysb_importer(model, constants=None,
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

    _process_pysb_species(model, ODE)
    _process_pysb_parameters(model, ODE, constants)
    if compute_conservation_laws:
        _process_pysb_conservation_laws(model, ODE)
    _process_pysb_expressions(model, ODE, observables, sigmas)
    _process_pysb_observables(model, ODE)

    ODE.generateBasicVariables()

    return ODE


def _process_pysb_species(model, ODE):
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


def _process_pysb_parameters(model, ODE, constants):
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


def _process_pysb_expressions(model, ODE, observables, sigmas):
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

            sigma_name, sigma_value = _get_sigma_name_and_value(
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


def _get_sigma_name_and_value(model, name, sigmas):
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


def _process_pysb_observables(model, ODE):
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


def _process_pysb_conservation_laws(model, ODE):
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
        monomers_without_conservation_law |= \
            _get_unconserved_monomers(rule, model)

    monomers_without_conservation_law |= \
        _compute_monomers_with_fixed_initial_conditions(model)

    cl_prototypes = _generate_cl_prototypes(
        monomers_without_conservation_law, model, ODE
    )
    conservation_laws = _construct_conservation_from_prototypes(
        cl_prototypes, model
    )

    _flatten_conservation_laws(conservation_laws)

    for cl in conservation_laws:
        ODE.add_conservation_law(cl['state'], cl['law'])


def _compute_monomers_with_fixed_initial_conditions(model):
    """Computes the set of monomers in a model with species that have fixed
    initial conditions

    Arguments:
        model: pysb model @type pysb.core.Model

    Returns:
        list of monomer names

    Raises:

    """
    monomers_with_fixed_initial_conditions = set()

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
                monomers_with_fixed_initial_conditions |= {monomer.name}

    return monomers_with_fixed_initial_conditions


def _generate_cl_prototypes(excluded_monomers, model, ODE):
    """Constructs a dict that contains preprocessed information for the
    construction of conservation laws

    Arguments:
        model: pysb model @type pysb.core.Model

        ODE: ODEModel instance @type ODEModel

    Returns:
        dict('monomer.name':{'possible_indices': ..., 'target_indices': ...}

    Raises:

    """
    cl_prototypes = dict()

    _compute_possible_indices(cl_prototypes, model, excluded_monomers)
    _compute_target_index(cl_prototypes, ODE)

    return cl_prototypes


def _compute_possible_indices(cl_prototypes, model, excluded_monomers):
    """Computes viable choices for target_index, ie species that could be
    removed and replaced by an algebraic expression according to the
    conservation law

    Arguments:
        cl_prototypes: dict in which possible indices will be written @type
        dict

        model: pysb model @type pysb.core.Model

        excluded_monomers: monomers for which no conservation laws will be
        computed @type list


    Returns:

    Raises:

    """
    for monomer in model.monomers:
        if monomer.name not in excluded_monomers:

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

            cl_prototypes[monomer.name] = dict()
            cl_prototypes[monomer.name]['possible_indices'] = [
                ix
                for ix, specie in enumerate(model.species)
                if extract_monomers(specie)[0] == monomer.name
            ]
            cl_prototypes[monomer.name]['species_count'] = len(
                cl_prototypes[monomer.name]['possible_indices']
            )


def _compute_target_index(cl_prototypes, ODE):
    """Computes the target index for every monomer


    Arguments:
        cl_prototypes: dict that contains possible indices for every monomer
        @type dict

        ODE: ODEModel instance @type ODEModel

    Returns:

    Raises:

    """
    possible_indices = list(set(list(itertools.chain(*[
        cl_prototypes[monomer]['possible_indices']
        for monomer in cl_prototypes
    ]))))

    appearance_counts = ODE.get_appearance_counts(possible_indices)

    for monomer in cl_prototypes:
        prototype = cl_prototypes[monomer]
        # extract monomer specific appearance counts
        prototype['appearance_counts'] = \
            [
                appearance_counts[possible_indices.index(idx)]
                for idx in prototype['possible_indices']
            ]
        # select target index as possible index with minimal appearance count
        idx = np.argmin(prototype['appearance_counts'])

        # remove entries from possible indices and appearance counts so we
        # do not consider them again in later iterations
        prototype['target_index'] = prototype['possible_indices'].pop(idx)
        prototype['appearance_count'] = prototype['appearance_counts'].pop(idx)

        prototype['fillin'] = \
            prototype['appearance_count'] * prototype['species_count']


    # we might end up with the same index for multiple monomers, so loop until
    # we have a set of unique target indices
    while len(cl_prototypes) > len(set(get_target_indices(cl_prototypes))):
        _greedy_target_index_update(cl_prototypes)


def _greedy_target_index_update(cl_prototypes):

    target_indices = get_target_indices(cl_prototypes)

    for monomer in cl_prototypes:
        prototype = cl_prototypes[monomer]
        if target_indices.count(prototype['target_index']) > 1:
            # compute how much fillin the next best target_index would yield

            # we exclude already existing target indices to avoid that
            # updating the target index removes uniqueness from already unique
            # target indices, this may slightly reduce chances of finding a
            # solution but prevents infinite loops
            for taken_index in list(set(target_indices)):
                local_idx = prototype['possible_indices'].index(
                    taken_index, None
                )
                if local_idx:
                    del prototype['possible_indices'][local_idx]
                    del prototype['appearance_counts'][local_idx]

            if len(prototype['possible_indices']) == 0:
                # we have exhausted all possible target indices, nothing left
                # we can do
                raise Exception('Could not compute a valid set of '
                                'conservation laws for this model')

            idx = np.argmin(prototype['appearance_counts'])

            prototype['local_index'] = idx
            prototype['alternate_target_index'] = \
                prototype['possible_indices'][idx]
            prototype['alternate_appearance_count'] = \
                prototype['appearance_counts'][idx]

            prototype['alternate_fillin'] = \
                prototype['alternate_appearance_count'] \
                * prototype['species_count']
            prototype['diff_fillin'] = \
                prototype['alternate_fillin'] - prototype['fillin']
        else:
            prototype['diff_fillin'] = 0

    # this puts prototypes with high diff_fillin last
    cl_prototypes = sorted(
        cl_prototypes.items(), key=lambda kv: kv[1]['diff_fillin']
    )

    for monomer in cl_prototypes:
        prototype = cl_prototypes[monomer]
        # we check that we
        # A) have an alternative index computed, i.e. that
        # that monomer originally had a non-unique target_index
        # B) that the target_index still is not unique. due to the sorting,
        # this will always be the monomer with the highest diff_fillin (note
        #  that the target indice counts are recomputed on the fly

        if prototype['diff_fillin'] > 0 \
                and get_target_indices(cl_prototypes).count(
                    prototype['target_index']
                ) > 1:

            prototype['fillin'] = prototype['alternate_fillin']
            prototype['target_index'] = prototype['alternate_target_index']
            prototype['appearance_count'] = \
                prototype['alternate_appearance_count']

            del prototype['possible_indices'][prototype['local_index']]
            del prototype['appearance_counts'][prototype['local_index']]


def get_target_indices(cl_prototypes):
    """Computes the list target indices for the current
    conservation law prototype

    Arguments:
        cl_prototypes: dict that contains target indices for every monomer
        @type dict

    Returns:

    Raises:

    """
    return [
        cl_prototypes[monomer]['target_index'] for monomer in cl_prototypes
    ]



def _construct_conservation_from_prototypes(cl_prototypes, model):
    """Computes the algebraic expression for the total amount of a given
    monomer

    Arguments:
        cl_prototype: pysb model @type pysb.core.Model

        ODE: ODEModel instance @type ODEModel


    Returns:

    Raises:

    """
    conservation_laws = []
    for monomer_name in cl_prototypes:
        target_index = cl_prototypes[monomer_name]['target_index']

        target_expression = sanitize_basic_sympy(
            sum([
                sp.Symbol(f'__s{ix}')
                * extract_monomers(specie).count(monomer_name)
                for ix, specie in enumerate(model.species)
                if ix != target_index
            ])
            /
            extract_monomers(model.species[target_index]).count(monomer_name)
        ) # normalize by the stoichiometry of the target species
        target_state = sp.Symbol(f'__s{target_index}')
        conservation_laws.append({
            'state': target_state,
            'law': target_expression,
        })

    return conservation_laws


def _flatten_conservation_laws(conservation_laws):
    unflattened_conservation_laws = \
        _get_unflattened_conservation_laws(conservation_laws)
    while len(unflattened_conservation_laws):
        for cl in conservation_laws:
            cl['law'] = cl['law'].subs(unflattened_conservation_laws)
        unflattened_conservation_laws = \
            _get_unflattened_conservation_laws(conservation_laws)


def _get_unflattened_conservation_laws(conservation_laws):
    """ returns a list of (state, law) tuples for conservation laws that still
    appear in other conservation laws

    Arguments:
        conservation_laws: dict({'state':state, 'law':law}) @type dict

    Returns:
        list((law,tuple))

    Raises:

    """
    free_symbols_cl = _conservation_law_variables(conservation_laws)
    return [
        (cl['state'], cl['law']) for cl in conservation_laws
        if cl['state'] in free_symbols_cl
    ]


def _conservation_law_variables(conservation_laws):
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


def _get_unconserved_monomers(rule, model):
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
            unconserved_monomers |= _get_changed_stoichiometries(
                [model.species[ix] for ix in reaction['reactants']],
                [model.species[ix] for ix in reaction['products']]
            )
    else:
        # otherwise we can simply extract all information for the rule
        # itself, which is computationally much more efficient
        unconserved_monomers |= _get_changed_stoichiometries(
            rule.reactant_pattern.complex_patterns,
            rule.product_pattern.complex_patterns
        )

    return unconserved_monomers


def _get_changed_stoichiometries(reactants, products):
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

