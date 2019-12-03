import libsbml


def flatten_sbml(sbml: libsbml.Model):
    """Reformulate parts of the SBML, in order to allow
        import of some SBML features.

    Arguments:
        sbml: Instance of libsbml.Model`.

    Returns:

    Raises:
    """

    process_rate_rules(sbml)
    process_species_assignment_rules(sbml)


def process_rate_rules(sbml: libsbml.SBMLDocument):
    """ Reformulate rate rules by an reaction of the form
     R: -> X with rate f(x).

    Arguments:
        sbml: Instance of libsbml.Model`.

    Returns:

    Raises:
    """

    # loop backwards over rules, so that one can delete rules inside the loop!
    for i in range(sbml.num_rules - 1, -1, -1):

        rule = sbml.getRule(i)
        rule_variable = sbml.getElementBySId(rule.getVariable())

        if rule.isRate() and isinstance(rule_variable, libsbml.Species):

            # create dummy reaction R: -> X with rate f(X)
            reaction = sbml.createReaction()
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
            for species in sbml.getListOfSpecies():

                if species.getId() in k.getFormula():
                    reaction.addModifier(
                        sbml.getElementBySId(species.getId())
                    )
            # remove rule
            rule.removeFromParentAndDelete()


def process_species_assignment_rules(sbml: libsbml.Model):
    """ Reformulate species assignment rules by replacing the species by
    a parameter.

    Arguments:
        sbml: Instance of libsbml.Model`.

    Returns:

    Raises:
    """

    # store the species ids to remove them as modifiers later...
    assignment_rule_species_ids = set()

    # loop backwards over rules, so that one can delete rules inside the loop!
    for i in range(sbml.num_rules - 1, -1, -1):

        rule = sbml.getRule(i)

        rule_variable = sbml.getElementBySId(rule.getVariable())

        if rule.isAssignment() and isinstance(rule_variable, libsbml.Species):
            
            assignment_rule_species_ids.add(rule_variable.getId())
            
            p = sbml.createParameter()
            p.setId(rule_variable.getId())
            p.setName(rule_variable.getName())
            p.setUnits(rule_variable.getUnits())
            p.setConstant(False)

            # delete species and rule
            rule_variable.removeFromParentAndDelete()
            rule.removeFromParentAndDelete()

    # remove modifier, if a species with an assignment rule was modifier
    for reaction in sbml.getListOfReactions():
        
        for modifier in reaction.getListOfModifiers():
            if modifier.getSpecies() in assignment_rule_species_ids:
                modifier.removeFromParentAndDelete()
