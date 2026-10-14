Fixed the sensitivities of the likelihood with respect to parameters of
event noise formulas: the derivative of the noise parameter was discarded
if the noise scale was taken from the model, and kept if it was taken from
the measurements.
