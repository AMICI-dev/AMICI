SBML initial assignments of parameters, compartments, conversion factors, and
species references that reference a parameter or compartment changed by an
event assignment now use its initial value. Previously, such quantities
followed the event-assigned value after the event (#3506). This also affects
PEtab v2 problems where such a parameter or compartment is changed in a later
experiment period, including the main period after pre-equilibration.
Also fixed an import error for initial assignments referencing a rate-rule
compartment without initial assignment (#3620).
