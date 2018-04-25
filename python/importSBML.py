import libsbml as sbml
import numpy as np
import symengine as sp

def getLocalParameters(reaction):
	kinLaw = reaction.getKineticLaw()
	return(kinLaw.getListOfLocalParameters())


reader = sbml.SBMLReader()
doc = reader.readSBML('/Users/F.Froehlich/Downloads/Speedy_v3_r403445_v1.sbml')
model = doc.getModel()

reactions = model.getListOfReactions()
n_r = len(reactions)

species = model.getListOfSpecies()
n_x = len(species)
speciesIndex = {species_id.getId(): species_index for species_index, species_id in enumerate(list(species))}
sym_species = sp.zeros(n_x,1)
for index, species in enumerate(list(species)):
    sym_species[index] = sp.symbols(species.getId())

parameters = model.getListOfParameters()
n_p = len(parameters)
sym_pars = sp.zeros(n_p,1)
for index, par in enumerate(list(parameters)):
    sym_pars[index] = sp.symbols(par.getId())

localParameters = []
[localParameters.extend(getLocalParameters(reaction)) for reaction in reactions]
n_ploc = len(localParameters)
sym_locpars = sp.zeros(n_ploc,1)
for index, par in enumerate(list(localParameters)):
    sym_locpars[index] = sp.symbols(par.getId())




stoichiometry_matrix = sp.zeros(n_x,n_r)
flux = sp.zeros(n_r,1)
w = sp.zeros(n_r,1)

for reaction_index, reaction in enumerate(model.getListOfReactions()):
        reactants = {r.getSpecies(): r.getStoichiometry() for r in reaction.getListOfReactants()}
        products = {p.getSpecies(): p.getStoichiometry() for p in reaction.getListOfProducts()}

        for reactant in reactants.keys():
            try:
                stoichiometry_matrix[speciesIndex[product], reaction_index] -= sp.sympify(reactants[reactant])
            except:
                stoichiometry_matrix[speciesIndex[product], reaction_index] = sp.sympify(reactants[reactant])

        for product in products.keys():
            try:
                stoichiometry_matrix[speciesIndex[product], reaction_index] += sp.sympify(products[product])
            except:
                stoichiometry_matrix[speciesIndex[product], reaction_index] = sp.sympify(products[product])

        w[reaction_index] = sp.sympify(reaction.getKineticLaw().getFormula())
        flux[reaction_index] = sp.sympify('w' + str(reaction_index))

rhs = stoichiometry_matrix*flux

jac = rhs.jacobian(sym_species)
dxdotdp = rhs.jacobian(sym_pars)

dwdx = w.jacobian(sym_species)
dwdp = w.jacobian(sym_pars)

xdot_lines = [
    '#include "amici/symbolic_functions.h"',
    '#include "amici/defines.h" //realtype definition',
    'typedef amici::realtype realtype;',
    ''#include <cmath> ',
    '#include species.cpp',
    '#include parameters.cpp',
    '#include fluxes.cpp',
    '',
    'void xdot_model_steadystate(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w)']

xdot_lines.append('')

[xdot_lines.append('    xdot[' + str(index) + '] = ' + str(math) + ';') if not math == 0 else None for index, math in enumerate(rhs)]

xdot_lines.append('}')

file_path = 'xdot.cpp'
open(file_path, 'w').write('\n'.join(xdot_lines))

parameters_lines = []
[parameters_lines.append('#define ' + str(symbol) + ' p[' + str(index) + '];') for index, symbol in enumerate(sym_pars)]
file_path = 'parameters.cpp'
open(file_path, 'w').write('\n'.join(parameters_lines))

species_lines = []
[species_lines.append('#define ' + str(symbol) + ' x[' + str(index) + '];') for index, symbol in enumerate(sym_species)]
file_path = 'species.cpp'
open(file_path, 'w').write('\n'.join(species_lines))

reactions_lines = []
[reactions_lines.append('#define w' + str(index) + ' w[' + str(index) + '];') for index in range(0,len(reactions))]
file_path = 'reactions.cpp'
open(file_path, 'w').write('\n'.join(reactions_lines))