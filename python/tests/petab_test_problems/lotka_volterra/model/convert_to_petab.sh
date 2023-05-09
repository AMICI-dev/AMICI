#!/usr/bin/env bash
set -eou pipefail
script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
petab_path=$script_dir/../petab

mkdir -p "${petab_path}"

# Validate input yaml2sbml model
yaml2sbml_validate lotka_volterra.yaml

# Copy measurements to PEtab directory
cp "${script_dir}/measurements.tsv" "$petab_path/measurements.tsv"

# Write the PEtab problem
python "${script_dir}/writer.py"

# Remove condition parameters from PEtab parameters table
for condition_parameter in beta delta departure_prey arrival_predator
do
	sed -i "/^${condition_parameter}/d" "${petab_path}/parameter*"
done

# Validate the PEtab problem
petablint -vy "${petab_path}/problem.yaml"
