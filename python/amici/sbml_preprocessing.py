import libsbml


def flatten_sbml(sbml_doc: libsbml.SBMLDocument):
    """Reformulate parts of the SBML, in order to allow import of some
    SBML features.

    Arguments:
        sbml_doc: Instance of libsbml.Model`.

    Returns:

    Raises:
    """

    process_rate_rules(sbml_doc)
    process_species_assignment_rules(sbml_doc)


def process_rate_rules(sbml_doc: libsbml.SBMLDocument):
    """ Reformulate rate rules by an reaction of the form
     R: -> X with rate f(x).

    Arguments:
        sbml_doc: Instance of libsbml.Model`.

    Returns:

    Raises:
    """

    # loop backwards over rules, so that one can delete rules inside the loop!
    for i in range(sbml_doc.getModel().num_rules - 1, -1, -1):

        rule = sbml_doc.getModel().getRule(i)

        rule_variable = sbml_doc.getModel().getElementBySId(rule.getVariable())

        if rule.isRate() and type(rule_variable) == libsbml.Species:

            # create dummy reaction R: -> X with rate f(X)
            reaction = sbml_doc.getModel().createReaction()
            reaction.setId('d_dt_' + rule_variable.getId())
            reaction.setName('d_dt_' + rule_variable.getName())
            reaction.setReversible(False)
            reaction.setFast(False)

            product = reaction.createProduct()
            product.setSpecies(rule.getVariable())
            product.setConstant(False)

            k = reaction.createKineticLaw()
            k.setMath(rule.getMath())

            # add modifiers
            for species in sbml_doc.getModel().getListOfSpecies():

                if species.getId() in k.getFormula():
                    reaction.addModifier(
                        sbml_doc.getModel().getElementBySId(species.getId())
                    )
            # remove rule
            rule.removeFromParentAndDelete()


def process_species_assignment_rules(sbml_doc: libsbml.SBMLDocument):
    """ Reformulate species assignment rules by replacing the species by
    a parameter.

    Arguments:
        sbml_doc: Instance of libsbml.Model`.

    Returns:

    Raises:
    """

    # loop backwards over rules, so that one can delete rules inside the loop!
    for i in range(sbml_doc.getModel().num_rules - 1, -1, -1):

        rule = sbml_doc.getModel().getListOfRules()[i]

        rule_variable = sbml_doc.getModel().getElementBySId(rule.getVariable())

        if rule.isAssignment() and type(rule_variable) == libsbml.Species:

            p = sbml_doc.getModel().createParameter()
            p.setId(rule_variable.getId())
            p.setName(rule_variable.getName())
            p.setUnits(rule_variable.getUnits())
            p.setConstant(False)

            print(rule_variable.getId())
            # remove modifier, if species is modifier of a reaction

            for reaction in sbml_doc.getModel().getListOfReactions():
                modifier_list = reaction.getListOfModifiers()
                for modifier in modifier_list:
                    if rule_variable.getId() == modifier.getSpecies():
                        modifier.removeFromParentAndDelete()

            # delete species and rule
            rule_variable.removeFromParentAndDelete()
            # rule.removeFromParentAndDelete()
