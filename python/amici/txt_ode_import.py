import sys
from libsbml import *
import re
import warnings


def parse_vars_and_values(line):

    try:
        id, value, unit = re.findall('(\S+)\s*=\s*(\d+\.?\d+|\w+)\s*(\S*)', line)[0]
    except IndexError:
        raise Exception('Error in parsing: line \n' + line + '\n is not correctly formated! \n'
                        'Make sure it has the format <variable> = <value> [<unit>] (where unit is optional)')
    return id, value, unit


def parse_equations(line):

    try:
        id, arguments, formula = re.findall('(\w+)\((\S+)\)=(\S+)', line.replace(' ', ''))[0]
    except IndexError:
        raise Exception('Error in parsing: line \n' + line + '\n is not correctly formated! \n'
                        'Make sure it has the format <function name> ( <arguments> ) = <formula>')
    return id, arguments, formula


def create_species(model, id, initial_amount,
                   substance_units):

    s = model.createSpecies()
    s.setId(id)
    s.setInitialAmount(float(initial_amount))

    substance_units = substance_units if substance_units else 'mole'  # if not specified set units to mole

    s.setSubstanceUnits(substance_units)
    return s


def create_parameter(model, id, constant, value, units):
    k = model.createParameter()
    k.setId(id)
    k.setName(id)
    k.setConstant(constant)
    k.setValue(float(value))

    units = units if units else 'dimensionless'  # if not specified set units to dimensionless

    k.setUnits(units)
    return k


def create_functions(model, id, arguments, formula):
    lambda_sting = 'lambda('+arguments+', '+formula
    f = model.createFunctionDefinition()
    f.setId(id)
    math_ast = parseL3Formula(formula)
    f.setMath(math_ast)
    return f


def create_rate_rule(model, id, species, formula):
    r = model.createRateRule()
    r.setMetaId('d/' + id + '_'+species)
    math_ast = parseL3Formula(formula)
    r.setMath(math_ast)
    return r


def read_time_block(model, line):

    if line.strip() != '':
        try:
            time_var, time_unit = re.findall('(\S+)\s*(\S*)', line)[0]
            return time_var, time_unit
        except IndexError:
            raise('Error in parsing: line \n' + line + '\n is not correctly formated! \n'
                        'Make sure it has the format <time_variable> [<unit>] (where unit is optional)')
        time_unit = time_unit if time_unit else "seconds"
        model.setTimeUnits(time_unit)


def read_constants_block(model, line):
    if line.strip() != '':

        id, value, unit = parse_vars_and_values(line)

        create_parameter(model, id, True, value, unit)


def read_parameters_block(model, line):
    if line.strip() != '':

        id, value, unit = parse_vars_and_values(line)

        create_parameter(model, id, False, value, unit)


def read_species_block(model, line):
    if line.strip() != '':

        id, inital_amount, unit = parse_vars_and_values(line)
        create_species(model, id, inital_amount, unit)


def read_functions_block(model, line):
    if line.strip() != '':

        id, arguments, formula = parse_equations(line)
        create_functions(model, id, arguments, formula)


def read_odes_block(model, line):
    if line.strip() != '':

        id, arguments, formula = parse_equations(line)
        create_rate_rule(model, id, arguments, formula)

def read_observables_block(model, line):
    warnings.warn('Observables nor supported yet')

def read_noise_block(model, line):
    warnings.warn('Noise nor supported yet')

def read_events_block(model, line):
    warnings.warn('Events nor supported yet')


def parse_txt(text_file):

    # create an SBML

    document = SBMLDocument(3, 1)
    model = document.createModel()

    function_dict = {'time': read_time_block ,
                     'constants': read_constants_block,
                     'parameters': read_parameters_block,
                     'species': read_species_block,
                     'functions': read_functions_block,
                     'odes': read_odes_block,
                     'observables': read_observables_block,
                     'noise': read_noise_block,
                     'events': read_events_block}

    with open(text_file, 'r') as f_in:
        line = f_in.readline()

        current_block = ''

        while line:

            if line.strip().replace(':', '') in function_dict.keys():

                current_block = line.strip().replace(':', '')

            else:
                function_dict[current_block](model, line)

            line = f_in.readline()

    return writeSBMLToString(document)



def reimport_from_SBML(sbmlfile):
    pass


def import_from_txt(text_file, sbml_file):

    print(parse_txt(text_file))
    sbml_as_string = parse_txt(text_file)
    with open(sbml_file, 'w') as f_out:
        f_out.write(sbml_as_string)
    # reimport_from_SBML()


import_from_txt('/Users/jvanhoefer/Dokumente/GitHub/DCM Software/Hackathon Frankfurt /ode_input.txt',
                '/Users/jvanhoefer/Dokumente/GitHub/DCM Software/Hackathon Frankfurt /sbml_out.xml')




"""
###################################################

def create_unit_definitions(model):
    unitDefinition = model.createUnitDefinition()
    unitDefinition.setId('litre_per_second')
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_LITRE)
    unit.setExponent(1)
    unit.setScale(0)
    unit.setMultiplier(1)
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_SECOND)
    unit.setExponent(-1)
    unit.setScale(0)
    unit.setMultiplier(1)

    unitDefinition = model.createUnitDefinition()
    unitDefinition.setId('mole_per_litre')
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_MOLE)
    unit.setExponent(1)
    unit.setScale(0)
    unit.setMultiplier(1)
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_LITRE)
    unit.setExponent(-1)
    unit.setScale(0)
    unit.setMultiplier(1)


    unitDefinition = model.createUnitDefinition()
    unitDefinition.setId('mole_per_second')
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_MOLE)
    unit.setExponent(1)
    unit.setScale(0)
    unit.setMultiplier(1)
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_SECOND)
    unit.setExponent(-1)
    unit.setScale(0)
    unit.setMultiplier(1)

    unitDefinition = model.createUnitDefinition()
    unitDefinition.setId('litre2_per_mole_per_second')
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_LITRE)
    unit.setExponent(2)
    unit.setScale(0)
    unit.setMultiplier(1)
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_MOLE)
    unit.setExponent(-1)
    unit.setScale(0)
    unit.setMultiplier(1)
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_SECOND)
    unit.setExponent(-1)
    unit.setScale(0)
    unit.setMultiplier(1)


def create_parameter(model, id, constant, value, units):
    k = model.createParameter()
    k.setId(id)
    k.setName(id)
    k.setConstant(constant)
    k.setValue(value)
    k.setUnits(units)
    return k

def create_reaction(model, id, reactants, products, formula, reversible = False, fast = False):
    r = model.createReaction()
    r.setId(id)
    r.setReversible(reversible)
    r.setFast(fast)

    for (coeff, name) in reactants:
        species_ref = r.createReactant()
        species_ref.setSpecies(name)
        species_ref.setConstant(True) # TODO ?
        species_ref.setStoichiometry(coeff)
    for (coeff, name) in products:
        species_ref = r.createProduct()
        species_ref.setSpecies(name)
        species_ref.setConstant(True) # TODO ?
        species_ref.setStoichiometry(coeff)

    math_ast = parseL3Formula(formula)
    kinetic_law = r.createKineticLaw()
    kinetic_law.setMath(math_ast)
    return r


def create_assigment_rule(model, name, formula):
    rule = model.createAssignmentRule()
    rule.setId(name)
    rule.setName(name)
    rule.setVariable(name)
    rule.setFormula(formula)
    return rule


def add_sigma(model, observable):
    # corresponding sigma
    p = create_parameter(model, f'sigma_{observable}', False, 1.0,
                         'mole_per_litre')
    p = create_parameter(model, f'noiseParameter1_{observable}', True, 0.2,
                         'mole_per_litre')
    rule = create_assigment_rule(model, f'sigma_{observable}',
                                 f'noiseParameter1_{observable}')


def create_model():
    document = SBMLDocument(3, 1)
    model = document.createModel()
    model.setTimeUnits("second")
    model.setExtentUnits("mole")
    model.setSubstanceUnits('mole')

    create_unit_definitions(model)

    c1 = model.createCompartment()
    c1.setId('c1')
    c1.setConstant(True)
    c1.setSize(1)
    c1.setSpatialDimensions(3)
    c1.setUnits('litre')

    s1 = create_species(model, 'x1', 'c1', False, 0.1)
    s2 = create_species(model, 'x2', 'c1', False, 0.4)
    s3 = create_species(model, 'x3', 'c1', False, 0.7)

    #TODO: initial amounts should be parameters
    p1 = create_parameter(model, 'p1', True, 1.0, 'litre2_per_mole_per_second')
    p2 = create_parameter(model, 'p2', True, 0.5, 'litre2_per_mole_per_second')
    p3 = create_parameter(model, 'p3', True, 0.4, 'litre_per_second')
    p4 = create_parameter(model, 'p4', True, 2.0, 'litre_per_second')
    p5 = create_parameter(model, 'p5', True, 0.1, 'mole_per_second')
    k0 = create_parameter(model, 'k0', True, 1.0, 'litre_per_second')


    create_reaction(model, 'r1', [(2, 'x1')], [(1, 'x2')], 'p1 * x1^2')
    create_reaction(model, 'r2', [(1, 'x1'), (1, 'x2')], [(1, 'x3')], 'p2 * x1 * x2')
    create_reaction(model, 'r3', [(1, 'x2')], [(2, 'x1')], 'p3 * x2')
    create_reaction(model, 'r4', [(1, 'x3')], [(1, 'x1'), (1, 'x2')], 'p4 * x3')
    create_reaction(model, 'r5', [(1, 'x3')], [], 'k0 * x3')
    create_reaction(model, 'r6', [], [(1, 'x1')], 'p5')

    #    S -2 -1  2  1  0 1  v1: p1 x1^2
    #       1 -1 -1  1  0 0  v2: p2 x1 x2
    #       0  1  0 -1 -1 0  v3: p3 x2
    #                        v4: p4 x3
    #                        v5: k0 x3
    #                        v6: p5
    # R1: 2X1         ->       X2
    # R2:  X1 + X2    ->           X3
    # R3:       X2    -> 2X1
    # R4:          X3 ->  X1 + X2
    # R5:          X3 ->
    # R6:             ->  X1

    # write observables
    for i in range(1, 4):
        # all states fully observed
        observable = f'observable_x{i}'
        p = create_parameter(model, observable, False, 1.0, 'mole_per_litre')
        rule = create_assigment_rule(model, observable, 'x%d' % i)
        add_sigma(model, f'x{i}')

    # Scaled x1
    p = create_parameter(model, 'observableParameter1_x1_scaled', True, 2.0, 'dimensionless')
    p = create_parameter(model, 'observable_x1_scaled', False, 1.0, 'mole_per_litre')
    rule = create_assigment_rule(model, 'observable_x1_scaled', 'observableParameter1_x1_scaled * x1')
    add_sigma(model, 'x1_scaled')

    # x2 with offset
    p = create_parameter(model, 'observableParameter1_x2_offsetted', True, 3.0, 'mole_per_litre')
    p = create_parameter(model, 'observable_x2_offsetted', False, 1.0, 'mole_per_litre')
    rule = create_assigment_rule(model, 'observable_x2_offsetted', 'observableParameter1_x2_offsetted + x2')
    add_sigma(model, 'x2_offsetted')

    # Fully observed with sigma parameter
    # observable
    p = create_parameter(model, 'observable_x1withsigma', False, 1.0,
                         'mole_per_litre')
    rule = create_assigment_rule(model, 'observable_x1withsigma', 'x1')
    add_sigma(model, 'x1withsigma')

    return writeSBMLToString(document)


if __name__ == '__main__':
    print(create_model())
"""
