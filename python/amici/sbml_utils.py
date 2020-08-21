"""
SBML Utilities
------------
This module provides helper functions for SBML files.
"""


import libsbml
import sympy as sp

from typing import Optional, Union
from sympy.printing.mathml import MathMLContentPrinter
from sympy.core.parameters import evaluate


################################################################################


sbml_time_symbol = sp.Symbol('time', real=True)
amici_time_symbol = sp.Symbol('t', real=True)


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


def hasCompartment(model: libsbml.Model, compartmentId) -> bool:
    """
    Check whether a compartment with SBML ID `compartmentId` is present in the
    SBML model `model`.
    """
    compartmentId = str(compartmentId)
    for cmp in model.getListOfCompartments():
        if cmp.getId() == compartmentId:
            return True
    return False


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
    Check whether a rule with SBML ID `variableId` is present in the
    SBML model `model`.
    """
    ruleId = str(ruleId)
    for rule in model.getListOfRules():
        if rule.getIdAttribute() == ruleId:
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
    if hasCompartment(model, compartmentId):
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
    elif not hasCompartment(model, compartmentId):
        raise SbmlException(
            f'No compartment with ID {compartmentId}'
        )

    sp = model.createSpecies()
    if sp.setIdAttribute(speciesId) != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise SbmlException(f'{speciesId} is not a valid SBML ID')
    sp.setCompartment(compartmentId)
    sp.setInitialAmount(float(initial_amount))
    if units is not None:
        sp.setUnits(str(units))

    return sp


def addParameter(
        model,
        parameterId,
        *,
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

    # Check whether a parameter with the same ID already exists
    # TODO the resulting SBML may still be invalid
    #      if other types of objects (e.g., species) have the same ID
    if hasParameter(model, parameterId):
        raise SbmlException(
            f'A parameter with ID {parameterID} has already been defined'
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
    if hasRuleWithId(mode, ruleId):
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
    if hasRuleWithId(mode, ruleId):
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


class SbmlException(Exception):
    "Exception class for SBML-related errors."
    pass


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
    def doprint(self, expr):
        mathml = '<math xmlns="http://www.w3.org/1998/Math/MathML">'
        mathml += super().doprint(expr)
        mathml += '</math>'
        mathml = mathml.replace(
            f'<ci>time</ci>',
            '<csymbol encoding="text" definitionURL="http://www.sbml.org/sbml/symbols/time"> time </csymbol>'
        )
        return mathml


def sbmlMathML(expr, *, replace_time: bool = False, **settings) -> str:
    """
    Prints a SymPy expression to a MathML expression parsable by libSBML.

        :param expr:
            expression to be converted to MathML (will be sympified).

        :param replace_time:
            replace the AMICI time symbol with the SBML time symbol.
    """
    with evaluate(False):
        expr = sp.sympify(expr)
        if replace_time:
            expr = expr.subs(amici_time_symbol, sbml_time_symbol)
    return MathMLSbmlPrinter(settings).doprint(expr)


def sbmlMathAST(expr, **kwargs) -> libsbml.ASTNode:
    """
    Convert a SymPy expression to SBML math AST.

        :param expr:
            expression to be converted (will be sympified).

        :param kwargs:
            extra options for MathML conversion.
    """
    _mathml = sbmlMathML(expr, **kwargs)
    ast = libsbml.readMathMLFromString(_mathml)
    if ast is None:
        raise SbmlException(
            f'error while converting the following expression to SBML AST.\n'
            'expression:\n{expr}\n'
            'MathML:\n{_mathml}'
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
            'expression:\n{expr}\n'
            'MathML:\n{mathml}'
        )
