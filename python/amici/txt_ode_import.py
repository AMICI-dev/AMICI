import sys
from libsbml import parseL3Formula, SBMLDocument, writeSBMLToString
import re
import warnings


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
        unit: str, the units for the variable's value.
    """

    try:
        id, arguments, formula = re.findall('(\w+)\((\S+)\)=(\S+)', line.replace(' ', ''))[0]
    except IndexError:
        raise Exception('Error in parsing: line \n' + line + '\n is not correctly formated! \n'
                        'Make sure it has the format <function name> ( <arguments> ) = <formula>')

    return id, arguments, formula


def create_species(model, id: str, initial_amount: str, substance_units: str):
    """
    Creates a species and adds it to the given SBML model.

    Args:
        model: the SBML model to which the species will be added.
        id: the species ID
        initial_amount: the species initial amount
        substance_units: the units of the species amount

    Returns:
        s: the SBML species

    """

    s = model.createSpecies()
    s.setId(id)
    s.setInitialAmount(float(initial_amount))

    substance_units = substance_units if substance_units else 'mole'  # if not specified set units to mole
    s.setSubstanceUnits(substance_units)

    return s


def create_parameter(model, id: str, constant: bool, value: str, units: str):
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
        units: the units of the parameter or constant

    Returns:
        k: the SBML parameter or constant

    """

    k = model.createParameter()
    k.setId(id)
    k.setName(id)
    k.setConstant(constant)
    k.setValue(float(value))

    units = units if units else 'dimensionless'  # if not specified set units to dimensionless
    k.setUnits(units)

    return k


def create_functions(model, id: str, arguments: str, formula: str):
    """
    Creates a functionDefinition and adds it to the given SBML model.

    Args:
        model: SBML model to which the function will be added.
        id: the function id/name
        arguments: the arguments of the function
        formula: the formula of the function

    Returns:
        f: the SBML function definition

    """

    f = model.createFunctionDefinition()
    f.setId(id)
    math_ast = parseL3Formula(formula)
    f.setMath(math_ast)

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
    r.setMetaId('d/' + id + '_'+species)
    math_ast = parseL3Formula(formula)
    r.setMath(math_ast)

    return r


def read_time_block(model, line: str):
    """
    Reads and processes lines in the time block.
    If time units are not given they're set to seconds.

    Args:
        model: the SBML model
        line: a line in the time block in the ODE text file.

    Returns:
        None
    """

    if line.strip() != '':
        try:
            time_var, time_unit = re.findall('(\S+)\s*(\S*)', line)[0]
            return time_var, time_unit
        except IndexError:
            raise('Error in parsing: line \n' + line + '\n is not correctly formated! \n'
                        'Make sure it has the format <time_variable> [<unit>] (where unit is optional)')
        time_unit = time_unit if time_unit else "seconds"
        model.setTimeUnits(time_unit)


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


def read_parameters_block(model, line:str ):
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
        create_parameter(model, id, False, value, unit)


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


# TODO read_observables_block
def read_observables_block(model, line):
    warnings.warn('Observables not supported yet')


# TODO read_noise_block
def read_noise_block(model, line):
    warnings.warn('Noise not supported yet')


# TODO read_events_block
def read_events_block(model, line):
    warnings.warn('Events not supported yet')


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

    function_dict = {'time': read_time_block,
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

    sbml_string = writeSBMLToString(document)

    return sbml_string


# TODO reimport_from_SBML
def reimport_from_SBML(sbml_file):
    pass


def import_from_txt(text_file: str, sbml_file: str):
    """
    Takes in a text file with the specification of ODEs, parses it, converts it into the SBML format, and reimports
    the SBML file, so that it can be used with AMICI.

    Args:
        text_file : path to the text file with the ODEs specification
        sbml_file: path to the SBML file to be written out

    Returns:

    """

    sbml_as_string = parse_txt(text_file)

    with open(sbml_file, 'w') as f_out:
        f_out.write(sbml_as_string)

    reimport_from_SBML(sbml_file)


