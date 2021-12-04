"""
SBML Utilities
------------
This module provides helper functions for SBML files.
"""

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .sbml_import import SbmlImporter
    from typing import Optional, Union

import itertools as itt
import xml.dom.minidom
import libsbml
import sympy as sp

from sympy.printing.mathml import MathMLContentPrinter
from sympy.core.parameters import evaluate

from sympy.logic.boolalg import BooleanTrue as spTrue
from sympy.logic.boolalg import BooleanFalse as spFalse


################################################################################


sbml_time_symbol = sp.Symbol('time', real=True)
amici_time_symbol = sp.Symbol('t', real=True)


annotation_namespace = 'http://github.com/AMICI-dev/AMICI'


class SbmlException(Exception):
    "Exception class for SBML-related errors."
    pass


################################################################################


def createSbmlModel(modelId: str, level: int = 2, version: int = 5):
    """
    Helper for creating an empty SBML model.

        :param modelId:
        SBML ID of the new model.

        :param level:
        Level of the new SBML document.

        :param size:
        Version of the new SBML document.

        :return:
        a tuple containing the newly created `libsbml.SBMLDocument`
        and `libsbml.Model`.
    """
    doc = libsbml.SBMLDocument(level, version)
    model = doc.createModel()
    model.setId(modelId)
    return doc, model


################################################################################

def hasSpecies(model: libsbml.Model, speciesId) -> bool:
    """
    Check whether a species with SBML ID `speciesId` is present in the
    SBML model `model`.
    """
    speciesId = str(speciesId)
    for sp in model.getListOfSpecies():
        if sp.getIdAttribute() == speciesId:
            return True
    return False


def hasParameter(model: libsbml.Model, parameterId) -> bool:
    """
    Check whether a parameter with SBML ID `variableId` is present in the
    SBML model `model`.
    """
    parameterId = str(parameterId)
    for par in model.getListOfParameters():
        if par.getIdAttribute() == parameterId:
            return True
    return False


def hasRule(model: libsbml.Model, variableId) -> bool:
    """
    Check whether a SBML rule for quantity `variableId` is present in the
    SBML model `model`.
    """
    variableId = str(variableId)
    for rule in model.getListOfRules():
        if rule.getVariable() == variableId:
            return True
    return False


def hasRuleWithId(model: libsbml.Model, ruleId) -> bool:
    """
    Check whether a rule with SBML ID `ruleId` is present in the
    SBML model `model`.
    """
    ruleId = str(ruleId)
    for rule in model.getListOfRules():
        if rule.getIdAttribute() == ruleId:
            return True
    return False

def hasReaction(model: libsbml.Model, reactionId) -> bool:
    """
    Check whether a reaction with SBML ID `reactionId` is present in the
    SBML model `model`.
    """
    reactionId = str(reactionId)
    for r in model.getListOfReactions():
        if r.getIdAttribute() == reactionId:
            return True
    return False


################################################################################


def addCompartment(
        model,
        compartmentId,
        *,
        size: float = 1.0
    ) -> libsbml.Species:
    """
    Helper for adding a compartment to a SBML model.

        :param model:
        SBML model to which the compartment is to be added.

        :param compartmentId:
        SBML ID of the new compartment.

        :param size:
        Size of the new compartment. Defaults to `1.0`.

        :return:
        the new compartment as a `libsbml.Compartment` object.
    """
    compartmentId = str(compartmentId)

    # Check whether a compartment with the same ID already exists
    # TODO the resulting SBML may still be invalid
    #      if other types of objects (e.g., parameter) have the same ID
    if model.getCompartment(compartmentId):
        raise SbmlException(
            f'A compartment with ID {compartmentId} has already been defined'
        )

    cmp = model.createCompartment()
    if cmp.setId(compartmentId) != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise SbmlException(f'{compartmentId} is not a valid SBML ID')
    cmp.setSize(size)

    return cmp


def addSpecies(
        model,
        speciesId,
        *,
        compartmentId: Optional[str] = None,
        name: Union[bool, str, None] = None,
        initial_amount: float = 0.0,
        units: Optional[str] = None
    ) -> libsbml.Species:
    """
    Helper for adding a species to a SBML model.

        :param model:
        SBML model to which the species is to be added.

        :param speciesId:
        SBML ID of the new species.

        :param compartmentId:
        Compartment ID for the new species.
        If there is only one compartment it can be auto-selected.

        :param initial_amount:
        Initial amount of the new species.

        :param units:
        Units attribute for the new species.

        :return:
        the new species as a `libsbml.Species` object.
    """
    speciesId = str(speciesId)
    if name is True:
        name = speciesId

    # Check whether a species with the same ID already exists
    # TODO the resulting SBML may still be invalid
    #      if other types of objects (e.g., parameter) have the same ID
    if hasSpecies(model, speciesId):
        raise SbmlException(
            f'A species with ID {speciesId} has already been defined'
        )

    if compartmentId is None:
        compartments = model.getListOfCompartments()
        if len(compartments) != 1:
            raise SbmlException(
                'Compartment auto-selection is possible '
                'only if there is one and only one compartment.'
            )
        compartmentId = compartments[0].getId()
    elif not model.getCompartment(compartmentId):
        raise SbmlException(f'No compartment with ID {compartmentId}')

    sp = model.createSpecies()
    if sp.setIdAttribute(speciesId) != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise SbmlException(f'{speciesId} is not a valid SBML ID')
    sp.setCompartment(compartmentId)
    sp.setInitialAmount(float(initial_amount))
    if units is not None:
        sp.setUnits(str(units))
    if isinstance(name, str):
        sp.setName(name)

    return sp


def addParameter(
        model,
        parameterId,
        *,
        name: Union[bool, str, None] = None,
        value: Optional[float] = None,
        units: Optional[str] = None,
        constant: Optional[bool] = None
    ) -> libsbml.Parameter:
    """
    Helper for adding a parameter to a SBML model.

        :param model:
        SBML model to which the parameter is to be added.

        :param parameterId:
        SBML ID of the new parameter.

        :param value:
        Value attribute for the new parameter.

        :param units:
        Units attribute for the new parameter.

        :param constant:
        Constant attribute for the new parameter.

        :return:
        the new parameter as a `libsbml.Parameter` object.
    """
    parameterId = str(parameterId)
    if name is True:
        name = parameterId

    # Check whether a parameter with the same ID already exists
    # TODO the resulting SBML may still be invalid
    #      if other types of objects (e.g., species) have the same ID
    if hasParameter(model, parameterId):
        raise SbmlException(
            f'A parameter with ID {parameterId} has already been defined'
        )

    par = model.createParameter()
    if par.setIdAttribute(parameterId) != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise SbmlException(f'{parameterId} is not a valid SBML ID')
    if units is not None:
        par.setUnits(str(units))
    if constant is not None:
        par.setConstant(bool(constant))
    if value is not None:
        par.setValue(float(value))
    if isinstance(name, str):
        par.setName(name)

    return par


def addAssignmentRule(
        model: libsbml.Model,
        variableId,
        formula,
        ruleId: Optional[str] = None
    ) -> libsbml.AssignmentRule:
    """
    Helper for adding an assignment rule to a SBML model.

        :param model:
        SBML model to which the assignment rule is to be added.

        :param variableId:
        SBML ID of the quantity for which the assignment rule is to be added.

        :param formula:
        Formula for the assignment rule (it will be sympified).

        :param ruleId:
        SBML ID of the new assignment rule.
        Defaults to `'assignment_' + variableId`.

        :return:
        the assignment rule as a `libsbml.AssignmentRule` object.
    """
    variableId = str(variableId)
    if ruleId is None:
        ruleId = 'assignment_' + variableId

    # Check whether rules exists for this parameter or with the same name
    # TODO the resulting SBML may still be invalid
    #      if other types of objects (e.g., species) have the same ID
    if hasRule(model, variableId):
        raise SbmlException(
            f'A rule for parameter {variableId} has already been defined.'
        )
    if hasRuleWithId(model, ruleId):
        raise SbmlException(
            f'A rule with SBML ID {ruleId} has already been defined.'
        )

    rule = model.createAssignmentRule()
    if rule.setVariable(variableId) != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise SbmlException(f'{variableId} is not a valid SBML ID')
    if rule.setIdAttribute(ruleId) != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise SbmlException(f'{ruleId} is not a valid SBML ID')
    setSbmlMath(rule, formula)

    return rule


def addRateRule(
        model: libsbml.Model,
        variableId,
        formula,
        ruleId: Optional[str] = None
    ) -> libsbml.RateRule:
    """
    Helper for adding a rate rule to a SBML model.

        :param model:
        SBML model to which the rate rule is to be added.

        :param variableId:
        SBML ID of the quantity for which the rate rule is to be added.

        :param formula:
        Formula for the rate rule (it will be sympified).

        :param ruleId:
        SBML ID of the new rate rule.
        Defaults to `'rate_' + variableId`.

        :return:
        the new rate rule as a `libsbml.RateRule` object.
    """
    variableId = str(variableId)
    if ruleId is None:
        ruleId = 'rate_' + variableId

    # Check whether rules exists for this parameter or with the same name
    # TODO the resulting SBML may still be invalid
    #      if other types of objects (e.g., species) have the same ID
    if hasRule(model, variableId):
        raise SbmlException(
            f'A rule for parameter {variableId} has already been defined.'
        )
    if hasRuleWithId(model, ruleId):
        raise SbmlException(
            f'A rule with SBML ID {ruleId} has already been defined.'
        )

    rule = model.createRateRule()
    if rule.setVariable(variableId) != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise SbmlException(f'{variableId} is not a valid SBML ID')
    if rule.setIdAttribute(ruleId) != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise SbmlException(f'{ruleId} is not a valid SBML ID')
    setSbmlMath(rule, formula)

    return rule


def addInflow(model, speciesId, rate, *, reactionId: Optional[str] = None, reversible: bool = False):
    speciesId = str(speciesId)
    if reactionId is None:
        reactionId = f'inflow_of_{speciesId}'

    if hasReaction(model, reactionId):
        raise SbmlException(
            f'A reaction with SBML ID {reactionId} has already been defined.'
        )

    reaction = model.createReaction()
    if reaction.setId(reactionId) != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise SbmlException(f'{reactionId} is not a valid SBML ID')
    reaction.setReversible(reversible)

    spr = reaction.createProduct()
    spr.setSpecies(speciesId)

    kl = reaction.createKineticLaw()
    compartmentId = model.getSpecies(speciesId).getCompartment()
    setSbmlMath(kl, sp.Symbol(compartmentId) * rate)

    return reaction


################################################################################


def getSbmlUnits(model: libsbml.Model, x) -> Union[None, str]:
    """
    Try to get the units for expression `x`.
        :param model:
        SBML model.
        :param x:
        Expression to get the units of.
        :return:
        A string if the units could be determined, otherwise `None`.
    """
    # TODO can the SBML unit inference machinery be used?
    x = sp.sympify(x)
    if x.is_Symbol:
        if x.name == 'time':
            if model.isSetTimeUnits():
                return model.getTimeUnits()
            else:
                return None
        else:
            par = model.getParameter(x.name)
            if par is None:
                return None
            else:
                units = par.getUnits()
                if units == '':
                    return None
                else:
                    return units
    else:
        return None


################################################################################
# SymPy to SBML MathML/AST conversion


def pretty_xml(ugly_xml: str):
    dom = xml.dom.minidom.parseString(ugly_xml)
    pretty_xml = dom.toprettyxml()
    # We must delete the first line (xml header)
    return pretty_xml[pretty_xml.index('\n')+1:]


class MathMLSbmlPrinter(MathMLContentPrinter):
    """Prints a SymPy expression to a MathML expression parsable by libSBML.
    Differences from `sympy.MathMLContentPrinter`:
    1. underscores in symbol names are not converted to subscripts
    2. symbols with name 'time' are converted to the SBML time symbol
    """
    def _print_Symbol(self, sym):
        ci = self.dom.createElement(self.mathml_tag(sym))
        ci.appendChild(self.dom.createTextNode(sym.name))
        return ci
    def doprint(self, expr, *, pretty: bool = False):
        mathml = '<math xmlns="http://www.w3.org/1998/Math/MathML">'
        mathml += super().doprint(expr)
        mathml += '</math>'
        mathml = mathml.replace(
            f'<ci>time</ci>',
            '<csymbol encoding="text" definitionURL="http://www.sbml.org/sbml/symbols/time"> time </csymbol>'
        )
        return pretty_xml(mathml) if pretty else mathml


def sbmlMathML(expr, *, replace_time: bool = False, pretty: bool = False, **settings) -> str:
    """
    Prints a SymPy expression to a MathML expression parsable by libSBML.

        :param expr:
            expression to be converted to MathML (will be sympified).

        :param replace_time:
            replace the AMICI time symbol with the SBML time symbol.

        :param pretty:
            prettify the resulting MathML.
    """
    with evaluate(False):
        expr = sp.sympify(expr)
        if replace_time:
            expr = expr.subs(amici_time_symbol, sbml_time_symbol)
    return MathMLSbmlPrinter(settings).doprint(expr, pretty=pretty)


def sbmlMathAST(expr, **kwargs) -> libsbml.ASTNode:
    """
    Convert a SymPy expression to SBML math AST.

        :param expr:
            expression to be converted (will be sympified).

        :param kwargs:
            extra options for MathML conversion.
    """
    mathml = sbmlMathML(expr, **kwargs)
    ast = libsbml.readMathMLFromString(mathml)
    if ast is None:
        raise SbmlException(
            f'error while converting the following expression to SBML AST.\n'
            f'expression:\n{expr}\n'
            f'MathML:\n{pretty_xml(mathml)}'
        )
    return ast


def setSbmlMath(obj: libsbml.SBase, expr, **kwargs) -> None:
    """
    Set the math attribute of a SBML node using a SymPy expression.

        :param obj:
            SBML node supporting `setMath` method.

        :param expr:
            expression to which the math attribute of `obj` should be se to
            (will be sympified).

        :param kwargs:
            extra options for MathML conversion.
    """
    mathml = sbmlMathAST(expr, **kwargs)
    if obj.setMath(mathml) != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise SbmlException(
            f'Could not set math attribute of SBML object {obj}\n'
            f'expression:\n{expr}\n'
            f'MathML:\n{pretty_xml(mathml)}'
        )


################################################################################
# MathML to Sympy conversion


def mathml2sympy(mathml: str, *, evaluate: bool = False, locals=None, expression_type=None):
    ast = libsbml.readMathMLFromString(mathml)
    if ast is None:
        raise ValueError(
            f'libSBML could not parse MathML string:\n{pretty_xml(mathml)}'
        )

    formula = _parse_logical_operators(libsbml.formulaToL3String(ast))

    if evaluate:
        expr = sp.sympify(formula, locals=locals)
    else:
        with sp.core.parameters.evaluate(False):
            expr = sp.sympify(formula, locals=locals)

    expr = _parse_special_functions(expr)

    if expression_type is not None:
        _check_unsupported_functions(expr, expression_type)

    return expr


def _check_unsupported_functions(sym: sp.Basic,
                                 expression_type: str,
                                 full_sym: sp.Basic = None):
    """
    Recursively checks the symbolic expression for unsupported symbolic
    functions

    :param sym:
        symbolic expressions

    :param expression_type:
        type of expression, only used when throwing errors
    """
    if full_sym is None:
        full_sym = sym

    unsupported_functions = (
        sp.functions.factorial, sp.functions.ceiling, sp.functions.floor,
        sp.functions.sec, sp.functions.csc, sp.functions.cot,
        sp.functions.asec, sp.functions.acsc, sp.functions.acot,
        sp.functions.acsch, sp.functions.acoth,
        sp.Mod, sp.core.function.UndefinedFunction
    )

    unsupp_fun_type = next(
        (
            fun_type
            for fun_type in unsupported_functions
            if isinstance(sym.func, fun_type)
        ),
        None
    )
    if unsupp_fun_type:
        raise SBMLException(f'Encountered unsupported expression '
                            f'"{sym.func}" of type '
                            f'"{unsupp_fun_type}" as part of a '
                            f'{expression_type}: "{full_sym}"!')
    for fun in list(sym._args) + [sym]:
        unsupp_fun_type = next(
            (
                fun_type
                for fun_type in unsupported_functions
                if isinstance(fun, fun_type)
            ),
            None
        )
        if unsupp_fun_type:
            raise SBMLException(f'Encountered unsupported expression '
                                f'"{fun}" of type '
                                f'"{unsupp_fun_type}" as part of a '
                                f'{expression_type}: "{full_sym}"!')
        if fun is not sym:
            _check_unsupported_functions(fun, expression_type)


def _parse_special_functions(sym: sp.Basic, toplevel: bool = True) -> sp.Basic:
    """
    Recursively checks the symbolic expression for functions which have be
    to parsed in a special way, such as piecewise functions

    :param sym:
        symbolic expressions

    :param toplevel:
        as this is called recursively, are we in the top level expression?
    """
    args = tuple(_parse_special_functions(arg, False) for arg in sym._args)

    if sym.__class__.__name__ == 'abs':
        return sp.Abs(sym._args[0])
    elif sym.__class__.__name__ == 'xor':
        return sp.Xor(*sym.args)
    elif sym.__class__.__name__ == 'piecewise':
        # how many condition-expression pairs will we have?
        return sp.Piecewise(*grouper(args, 2, True))
    elif isinstance(sym, (sp.Function, sp.Mul, sp.Add)):
        sym._args = args
    elif toplevel:
        # Replace boolean constants by numbers so they can be differentiated
        #  must not replace in Piecewise function. Therefore, we only replace
        #  it the complete expression consists only of a Boolean value.
        if isinstance(sym, spTrue):
            sym = sp.Float(1.0)
        elif isinstance(sym, spFalse):
            sym = sp.Float(0.0)

    return sym


def _parse_logical_operators(math_str: Union[str, float, None]
                             ) -> Union[str, float, None]:
    """
    Parses a math string in order to replace logical operators by a form
    parsable for sympy

    :param math_str:
        str with mathematical expression
    :param math_str:
        parsed math_str
    """
    if not isinstance(math_str, str):
        return math_str

    if ' xor(' in math_str or ' Xor(' in math_str:
        raise SBMLException('Xor is currently not supported as logical '
                            'operation.')

    return (math_str.replace('&&', '&')).replace('||', '|')


def grouper(iterable: Iterable, n: int,
            fillvalue: Any = None) -> Iterable[Iterable]:
    """
    Collect data into fixed-length chunks or blocks

    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"

    :param iterable:
        any iterable

    :param n:
        chunk length

    :param fillvalue:
        padding for last chunk if length < n

    :return: itertools.zip_longest of requested chunks
    """
    args = [iter(iterable)] * n
    return itt.zip_longest(*args, fillvalue=fillvalue)
