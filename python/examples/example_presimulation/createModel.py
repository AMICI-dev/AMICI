from pysb import *
from pysb.export import sbml
from pysb.util import alias_model_components
import os

import libsbml

Model()

Monomer('PROT', ['kin','drug','phospho'], {'phospho': ['u', 'p']})
Parameter('PROT_0', 10)
alias_model_components()
Initial(PROT(phospho='u', drug=None, kin=None),
        Expression('initProt',PROT_0))

Monomer('DRUG', ['bound'])
Parameter('DRUG_0', 9)
alias_model_components()
Initial(DRUG(bound=None),
        Expression('initDrug', DRUG_0))

Parameter('kon_prot_drug', 0.1)
Parameter('koff_prot_drug', 0.1)
alias_model_components()

Rule('PROT_DRUG_bind',
     DRUG(bound=None) + PROT(phospho='u', drug=None, kin=None) |
     DRUG(bound=1) % PROT(phospho='u', drug=1, kin=None),
     kon_prot_drug,
     koff_prot_drug
     )

Monomer('KIN', ['bound'])
Parameter('KIN_0', 1)
alias_model_components()
Initial(KIN(bound=None),
        Expression('initKin', KIN_0))

Parameter('kon_prot_kin', 0.1)
alias_model_components()

Rule('PROT_KIN_bind',
     KIN(bound=None) + PROT(phospho='u', drug=None, kin=None) >>
     KIN(bound=1) % PROT(phospho='u', drug=None, kin=1),
     kon_prot_kin,
     )

Parameter('kphospho_prot_kin', 0.1)
alias_model_components()

Rule('PROT_KIN_phospho',
     KIN(bound=1) % PROT(phospho='u', drug=None, kin=1) >>
     KIN(bound=None) + PROT(phospho='p', drug=None, kin=None),
     kphospho_prot_kin,
     )

Parameter('kdephospho_prot', 0.1)
alias_model_components()

Rule('PROT_dephospho',
     PROT(phospho='p', drug=None, kin=None) >>
     PROT(phospho='u', drug=None, kin=None),
     kdephospho_prot,
     )

Observable('pPROT', PROT(phospho='p'))
Observable('tPROT', PROT())
alias_model_components()

sbml_output = sbml.SbmlExporter(model).export()

outfile = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                       'model_presimulation.xml')
with open(outfile, 'w') as f:
    f.write(sbml_output)


#sbml postProcessing
SBMLreader = libsbml.SBMLReader()
sbml_doc = SBMLreader.readSBML(outfile)

sbml_model = sbml_doc.getModel()

for par in sbml_model.getListOfParameters():
    if par.getName() == 'pPROT':
        phosphoId = par.getId()

        phosphoRule = \
            sbml_model.getAssignmentRuleByVariable(
                phosphoId
            ).getFormula()

        totalId = [par.getId() for par in sbml_model.getListOfParameters()
                   if par.getName() == 'tPROT'][0]

        totalRule = \
            sbml_model.getAssignmentRuleByVariable(
                totalId
            ).getFormula()

        new_phosphoRule = '(' + phosphoRule + ')/(' + totalRule + ')'
        sbml_model.getAssignmentRuleByVariable(
            phosphoId
        ).setFormula(new_phosphoRule)

sbml_model.getAssignmentRuleByVariable(
            totalId
        ).removeFromParentAndDelete()
sbml_model.getParameter(totalId).removeFromParentAndDelete()



libsbml.writeSBML(sbml_doc, outfile)