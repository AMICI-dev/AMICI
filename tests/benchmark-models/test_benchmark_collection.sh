#!/bin/bash
# Import and run selected benchmark models with nominal parameters and check
# agreement with reference values
#
# Expects environment variable BENCHMARK_COLLECTION to provide path to
# benchmark collection model directory

# Confirmed to be working
# FIXME only test one model until sensitivities appear correct, to save compilation/simulation time
#Beer_MolBioSystems2014
#Boehm_JProteomeRes2014
#Borghans_BiophysChem1997
#Brannmark_JBC2010
#Crauste_CellSystems2017
#Elowitz_Nature2000
#Fiedler_BMC2016
#Fujita_SciSignal2010
#Schwen_PONE2014
#Sneyd_PNAS2002
#Weber_BMC2015
models="
Zheng_PNAS2012"

#
# Not matching reference for unclear reasons
# Lucarelli_CellSystems2018
#
# PEtab needs fixing: Bachmann_MSB2011
#
# Unsupported:
#
# Becker_Science2010: multiple models
#
# no reference value:
# Alkan_SciSignal2018
# Blasi_CellSystems2016
# Hass_PONE2017
# Korkut_eLIFE2015
# Perelson_Science1996
# Bruno_JExpBio2016
#
# yaml missing:
# Casaletto_PNAS2019
#
# Model missing:
# Merkle_PCB2016
#
# SBML extensions:
# Parmar_PCB2019
#
# Events:
# Swameye_PNAS2003
#
# state-dependent sigmas:
# Raia_CancerResearch2011
#
# Evaluation is known to be inconsistent:
# Chen_MSB2009
#
# Integration Failure:
# Isensee_JCB2018

set -e

[[ -n "${BENCHMARK_COLLECTION}" ]] && model_dir="${BENCHMARK_COLLECTION}"

function show_help() {
  echo "-h: this help; -n: dry run, print commands; -b path_to_models_dir"
}

OPTIND=1
while getopts "h?nb:" opt; do
  case "$opt" in
  h | \?)
    show_help
    exit 0
    ;;
  n)
    dry_run=1
    ;;
  b)
    model_dir=$OPTARG
    ;;
  esac
done

script_path=$(dirname "$BASH_SOURCE")
script_path=$(cd "$script_path" && pwd)

for model in $models; do
  yaml="${model_dir}"/"${model}"/"${model}".yaml
  amici_model_dir=test_bmc/"${model}"
  mkdir -p "$amici_model_dir"
  cmd_import="amici_import_petab --verbose -y ${yaml} -o ${amici_model_dir} -n ${model} --flatten"
  cmd_run="$script_path/test_petab_model.py --verbose -y ${yaml} -d ${amici_model_dir} -m ${model} -c -s"

  printf '=%.0s' {1..40}
  printf "   %s   " "${model}"
  printf '=%.0s' {1..40}
  echo

  if [[ -z "$dry_run" ]]; then
    $cmd_import
    $cmd_run
  else
    echo "$cmd_import"
    echo "$cmd_run"
  fi

  printf '=%.0s' {1..100}
  echo
  echo
done
