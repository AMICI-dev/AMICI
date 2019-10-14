import sys
from libsbml import parseL3Formula, SBMLDocument, writeSBMLToString
import re
import warnings
# from . import SbmlImporter


def parse_time(line:str) -> (str, str):
    """
    Parses the time block with the form <time_variable> [<time_unit>].

    Args:
        line: a line of text from the file containing the ODEs.

    Returns:
        time_var: str, the time variable.
        time_unit: str, time unit. Optional
    """
    try:
        time_var, time_unit = re.findall('(\S+)\s*(\S*)', line)[0]
    except IndexError:
        raise Exception('Error in parsing: line \n' + line + '\n is not correctly formated! \n'
                        'Make sure it has the format <time_variable> [<time_unit>]')

    return time_var, time_unit


def parse_vars_and_values(line: str) -> (str, str, str):
    """
    Parses lines of the type: <variable> = <value> [<unit>], where unit is optional and returns the <variable> value
    as id, the <value> as value, and <unit> as unit.
    In particular, this is used to parse species, constants, and parameters.

    Args:
        line: a line of text from the file containing the ODEs.

    Returns:
        id: str, the variable id.
        value: str, the value of the variable.
        unit: str, the units for the variable's value.

    """

    try:
        id, value, unit = re.findall('(\S+)\s*=\s*(\d+\.?\d+|\w+)\s*(\S*)', line)[0]
    except IndexError:
        raise Exception('Error in parsing: line \n' + line + '\n is not correctly formated! \n'
                        'Make sure it has the format <variable> = <value> [<unit>] (where unit is optional)')

    return id, value, unit


def parse_equations(line: str) -> (str, str, str):
    """
    Parses equations with the form <left_hand_side>(<arguments>) = <equation>.
    In particular, it parses functions and odes defined in the ODE text file.

    Args:
        line: a line of text from the file containing the ODEs.

    Returns:
        id: str, the variable id.
        arguments: str, the value of the variable.
        formula: str, the formula for the right hand side.
    """

    try:
        id, arguments, formula = re.findall('(\w+)\((\S+)\)=(\S+)', line.replace(' ', ''))[0]
    except IndexError:
        raise Exception('Error in parsing: line \n' + line + '\n is not correctly formated! \n'
                        'Make sure it has the format <function name> ( <arguments> ) = <formula>')

    return id, arguments, formula


def parse_assignments_and_observable(line: str) -> (str, str):
    """
    Parses an observable with the form <observable name> = <equation> and returns the observable name as id
    and the formula as name.

    Args:
        line: a line of text from the file containing the ODEs.

    Returns:
        id: str, the observable id.
        formula: str, the formula for the observable.
    """

    try:
        id, formula = re.findall('(\S+)=(\S+)', line.replace(' ', ''))[0]
    except IndexError:
        raise Exception('Error in parsing line \n' + line + '\n is not correctly formated! \n'
                        'Make sure it has the format <observable name> = <formula>')

    return id, formula


def create_time(model, time_var: str, time_unit: str = None):
    """
    Creates the time variable,  add assignment to 'time'

    Args:
        model: the SBML model to which the species will be added.
        time_var: str, the time variable
        time_unit: str, the time unit
    """
    if time_unit is not None:

        time_unitdef = create_unit_definition(model, time_unit)
        model.setTimeUnits(time_unitdef)

    if time_var is not 'time':

        time_parameter = model.createParameter()
        time_parameter.setId(time_var)
        time_parameter.setName(time_var)
        time_parameter.setConstant(False)

        time_assignment = model.createAssignmentRule()
        time_assignment.setVariable(time_var)
        time_assignment.setMath(parseL3Formula('time'))


def create_species(model, id: str, initial_amount: str, unit_name: str):
    """
    Creates a species and adds it to the given SBML model.

    Args:
        model: the SBML model to which the species will be added.
        id: the species ID
        initial_amount: the species initial amount
        unit_name: the units of the species amount

    Returns:
        s: the SBML species

    """

    s = model.createSpecies()
    s.setId(id)
    s.setInitialAmount(float(initial_amount))
    s.setConstant(False)
    s.setBoundaryCondition(False)
    s.setHasOnlySubstanceUnits(False)
    s.setCompartment('Compartment')

    if unit_name:
        unitId = create_unit_definition(model, unit_name)
        s.setUnits(unitId)

    else:
        s.setSubstanceUnits('dimensionless')

    return s


def create_parameter(model, id: str, constant: bool, value: str, unit_name: str):
    """
    Creates a parameter or constant and adds it to the given SBML model.
    The difference between parameter and constant is only whether the parameter is set as constant or not. If it is
    set as constant then it is a constant, otherwise it is a parameter.
    If the parameter/constant units are not specified, units are set as dimensionless.

    Args:
        model: the SBML model to which the parameter/constant will be added.
        id: the parameter/constant ID
        constant: whether the parameter is actually a constant or not.
        value: the parameter or constant value
        unit_name: the units of the parameter or constant

    Returns:
        k: the SBML parameter or constant

    """

    k = model.createParameter()
    k.setId(id)
    k.setName(id)
    k.setConstant(constant)
    k.setValue(float(value))

    if unit_name:
        unitId = create_unit_definition(model, unit_name)
        k.setUnits(unitId)
    else:
        k.setUnits('dimensionless')

    return k


def create_functions(model, id: str, arguments: str, formula: str):
    """
    Creates a functionDefinition and adds it to the given SBML model.

    Args:
        model: SBML model to which the function will be added.
        id: the function id/name
        arguments: the arguments of the function (species AND parameters)
        formula: the formula of the function

    Returns:
        f: the SBML function definition

    """

    f = model.createFunctionDefinition()
    f.setId(id)
    math = parseL3Formula('lambda(' + arguments + ', ' + formula + ')')
    f.setMath(math)

    return f


def create_rate_rule(model, id: str, species: str, formula: str):
    """
    Creates a SBML rateRule for a species and adds it to the given model.
    This is where the ODEs from the text file are encoded.

    Args:
        model: SBML model to which the rate rule will be added.
        id: the function id/name
        species: the species name of the ODE
        formula: the right-hand-side of the ODE

    Returns:
        r: the SBML rateRule
    """

    r = model.createRateRule()
    r.setId('d/' + id + '_'+species)
    r.setVariable(species)
    math_ast = parseL3Formula(formula)
    r.setMath(math_ast)

    return r


def create_assignment(model, id: str, formula: str):
    """
    Creates an  assignment rule, that assigns id to formula.

    Args:
        model: SBML model to which the assignment rule will be added.
        id: str, the id of the assignment rule
        formula: str: contains the equation for the assignment rule
    """

    assignment_parameter = model.createParameter()
    assignment_parameter.setId(id)
    assignment_parameter.setName(id)
    assignment_parameter.setConstant(False)
    assignment_parameter.setUnits('dimensionless')

    assignment_rule = model.createAssignmentRule()
    assignment_rule.setVariable(id)
    assignment_rule.setMath(parseL3Formula(formula))


def create_observable(model, id: str, formula: str):
    """
    Creates an patameter with the name observalble_[Id] and an assignment rule, that assigns the parameter to
    the equation given in formula.

    Args:
        model: SBML model to which the rate rule will be added.
        id: str, the id of the observable
        formula: str: contains the equation for the observable
    """

    obs_param = model.createParameter()
    obs_param.setId('observable_' + id)
    obs_param.setName(id)
    obs_param.setConstant(False)
    obs_param.setUnits('dimensionless')

    obs_assignment_rule = model.createAssignmentRule()
    obs_assignment_rule.setVariable('observable_' + id)
    obs_assignment_rule.setMath(parseL3Formula(formula))


def read_time_block(model, line: str):
    """
    Reads and processes lines in the time block.
    If time units are not given they're set to seconds.

    Args:
        model: SBML model to which the rate rule will be added.
        line: a line in the time block in the ODE text file.

    Returns:
        None
    """

    if line.strip() != '':
        time_var, time_unit = parse_time(line)
        create_time(model, time_var, time_unit)


def read_constants_block(model, line: str):
    """
    Reads and processes lines in the constants block in the ODE text file.
    In particular, it reads the constants and adds them to the given SBML model.
    The expected format for constant definition is: <constant_name> = <value> [<units>], where units are optional.

    Args:
        model: the SBML model
        line: a line containing a constant definition

    Returns:
        None
    """

    if line.strip() != '':
        id, value, unit = parse_vars_and_values(line)
        create_parameter(model, id, True, value, unit)


def read_parameters_block(model, line: str):
    """
    Reads and processes lines in the parameters block in the ODE text file.
    In particular, it reads the parameters and adds them to the given SBML file.
    The expected format for parameter definition is: <parameter_name> = <value> [<units>], where units are optional.

    Args:
        model: the SBML model
        line: a line containing a parameter definition

    Returns:
        None
    """

    if line.strip() != '':
        id, value, unit = parse_vars_and_values(line)
        create_parameter(model, id, True, value, unit)


def read_species_block(model, line: str):
    """
    Reads and processes lines in the species block in the ODE text file.
    In particular, it reads the species and adds them to the given SBML file.
    The expected format of a species definition is: <species_name> = <initial_amount> [<units>], where units are
     optional.

    Args:
        model: the SBML model
        line: a line containing a species definition

    Returns:
        None
    """

    if line.strip() != '':
        id, inital_amount, unit = parse_vars_and_values(line)
        create_species(model, id, inital_amount, unit)


def read_assignments_block(model, line: str):

    if line.strip() != '':
        id, formula = parse_assignments_and_observable(line)
        create_assignment(model, id, formula)


def read_functions_block(model, line: str):
    """
    Reads and processes lines in the functions block in the ODE text file.
    In particular, it reads the functions and adds them to the given SBML file as functionDefinitions.
    The expected format of a function definition is: <function_name>(<arguments>) = <formula>.

    Args:
        model: a SBML model
        line: a line containing a function definition

    Returns:
        None
    """

    if line.strip() != '':
        id, arguments, formula = parse_equations(line)
        create_functions(model, id, arguments, formula)


def read_odes_block(model, line: str):
    """
    Reads and processes lines in the odes block in the ODE text file.
    In particular, it reads the odes and adds them to the given SBML file as rateRules.
    The expected format of an ode definition is: <d/d<time_variable>>(<species_name>) = <formula>

    Args:
        model: a SBML model
        line: a line containing an ODE definition

    Returns:
        None
    """

    if line.strip() != '':
        id, arguments, formula = parse_equations(line)
        create_rate_rule(model, id, arguments, formula)


def read_observables_block(model, line):
    """
    Reads an processes the lines in the observables-block in the ODE text file.
    In particular it generates the Observables in the SBML file.

    Args:
        model: SBML model (libsbml)
        line: a line containing an observable definition.

    Returns:
        None
    """
    if line.strip() != '':
        id, formula = parse_assignments_and_observable(line)
        create_observable(model, id, formula)


def read_noise_block(model, line):
    warnings.warn('Noise not supported yet')


# TODO read_events_block
def read_events_block(model, line):
    warnings.warn('Events not supported yet')


def create_unit_definition(model, unit_name: str)->str:
    """
    Checks, if unit with name unit_name exists, creates one if necessary and returns it`s Id

    Args:
        model: libsbml model
        unit_name: Unit Name(string)
        obj: libsbml object, which gets the unit assigned (parameter, species, ...)

    Retruns:
        unitId: str
    """

    predefined_units = {'ampere',
                        'avogadro',
                        'gram',
                        'katal',
                        'metre',
                        'second',
                        'watt',
                        'becquerel',
                        'gray',
                        'kelvin',
                        'mole',
                        'siemens',
                        'weber',
                        'candela',
                        'henry',
                        'kilogram',
                        'newton',
                        'sievert',
                        'coulomb',
                        'hertz',
                        'litre',
                        'ohm',
                        'steradian',
                        'dimensionless',
                        'item',
                        'lumen',
                        'pascal',
                        'tesla',
                        'farad',
                        'joule',
                        'lux',
                        'radian',
                        'volt'}

    if unit_name in predefined_units:
        return unit_name

    unit = None

    for u in model.getListOfUnitDefinitions():
        if u.getName() == unit_name:
            unit = u
            break

    if not unit:

        unit = model.createUnitDefinition()
        unit.setName(unit_name)

        unit_name = unit_name.replace('*', '_mult_').replace('/', '_div_').replace('^', '_pow_')
        unit.setId(unit_name)

    return unit.getId()


def parse_txt(text_file: str) -> str:
    """
    Takes in a text file with the specification of ODEs, parses it, and converts it into the SBML format.

    Args:
        text_file: path to the text file with the ODEs specification

    Returns:
        sbml_string: a string containing the ODEs conversion from text to SBML

    """

    document = SBMLDocument(3, 1)
    model = document.createModel()
    c = model.createCompartment()
    c.setId('Compartment')
    c.setConstant(True)

    function_dict = {'time': read_time_block,
                     'constants': read_constants_block,
                     'parameters': read_parameters_block,
                     'species': read_species_block,
                     'assignments': read_assignments_block,
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

    sbml_string = writeSBMLToString(document)

    return sbml_string


def txt2amici(model_name: str, text_file_dir: str, amici_output_dir: str = None, sbml_output_dir: str = None):
    """
    Takes in a text file with the specification of ODEs, parses it, converts it into the SBML format, and reimports
    the SBML file, so that it can be used with AMICI.

    Args:
        model_name: model_name
        text_file_dir : path to the text file with the ODEs specification
        amici_output_dir: path where the C++ files are written to
        sbml_output_dir: path to the SBML file to be written out

    Returns:

    """

    if not sbml_output_dir:
        sbml_output_dir = model_name + '.xml'

    sbml_as_string = parse_txt(text_file_dir)

    # write sbml file
    with open(sbml_output_dir, 'w') as f_out:
        f_out.write(sbml_as_string)

    # reimport files
    sbml_importer = SbmlImporter(sbml_output_dir)
    sbml_importer.sbml2amici(model_name,
                             amici_output_dir)


