A Lotka-Volterra model, based on the model provided as an example in the `yaml2sbml` package: https://github.com/yaml2sbml-dev/yaml2sbml

The PEtab problem can be created by running `bash model/convert_to_petab.sh` (test on Ubuntu 20.04) inside the `model` directory, in a Python 3 environment with `yaml2sbml`: https://pypi.org/project/yaml2sbml/

The model is augmented with new parameters `departure_prey` and `arrival_predator`, to allow for a steady-state under certain conditions. This results in a model that can pre-equilibrate, then switch to the expected oscillations of a Lotka-Volterra model.
