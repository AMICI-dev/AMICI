#!/bin/bash
# Import and run selected benchmark models with nominal parameters and check
# agreement with reference values
#
# Expects environment variable BENCHMARK_COLLECTION to provide path to
# benchmark collection model directory

# Confirmed to be working
models="
Bachmann_MSB2011
Beer_MolBioSystems2014
Boehm_JProteomeRes2014
Borghans_BiophysChem1997
Brannmark_JBC2010
Bruno_JExpBot2016
Crauste_CellSystems2017
Elowitz_Nature2000
Fiedler_BMCSystBiol2016
Fujita_SciSignal2010
Isensee_JCB2018
Lucarelli_CellSystems2018
Schwen_PONE2014
Smith_BMCSystBiol2013
Sneyd_PNAS2002
Weber_BMC2015
Zheng_PNAS2012"

#
# not merged:
# Becker_Science2010 (multiple models)             https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Becker_Science2010
# Hass_PONE2017 (???)                              https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Hass_PONE2017
# Korkut_eLIFE2015 (???)                           https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Korkut_eLIFE2015
# Casaletto_PNAS2019 (yaml missing)                https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Casaletto_PNAS2019
# Merkle_PCB2016 (model missing)                   https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Merkle_PCB2016
# Parmar_PCB2019 (SBML extensions)                 https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Parmar_PCB2019
# Swameye_PNAS2003 (splines)                       https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Swameye_PNAS2003
# Sobotta_Frontiers2017 (???)                      https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Sobotta_Frontiers2017
# Raia_CancerResearch2011 (state dependent sigmas) https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/model_Raia_CancerResearch2011
#
# no reference value:
# Alkan_SciSignal2018 (d2d: Alkan_DDRP_SciSignal2018)
# Bertozzi_PNAS2020 (gh: vanako, doi: https://doi.org/10.1073/pnas.2006520117, code/data: https://github.com/gomohler/pnas2020 (R))
# Blasi_CellSystems2016 (gh: Leonard Schmiester, doi: https://doi.org/10.1016/j.cels.2016.01.002, code/data: not available)
# Giordano_Nature2020 (gh: Paul Jonas Jost, doi: https://doi.org/10.1038/s41591-020-0883-7, code/data: http://users.dimi.uniud.it/~giulia.giordano/docs/papers/SIDARTHEcode.zip (MATLAB))
# Laske_PLOSComputBiol2019 (gh: Clemens Peiter, doi: https://doi.org/10.1128/JVI.00080-12 (?), code/data: ???)
# Okuonghae_ChaosSolitonsFractals2020 (gh: Paul Jonas Jost, doi: https://doi.org/10.1016/j.chaos.2020.110032, code/data: ???)
# Oliveira_NatCommun2021 (gh: lorenapohl, doi: https://doi.org/10.1038/s41467-020-19798-3, code: https://github.com/cidacslab/Mathematical-and-Statistical-Modeling-of-COVID19-in-Brazil (python) data: https://infovis.sei.ba.gov.br/covid19/ )
# Perelson_Science1996 (gh: Philipp Staedter, doi: https://doi.org/10.1126/science.271.5255.1582, code/data: ???)
# Rahman_MBS2016 (gh: Yannik Schaelte, doi: https://doi.org/10.1016/j.mbs.2016.07.009, code: not available, data: table in paper ...)
# Raimundez_PCB2020 (gh: Elba Raimundez, doi: https://doi.org/10.1371/journal.pcbi.1007147, code/data: https://zenodo.org/record/2908234#.Y5hUUS8w3yw (d2d))
# SalazarCavazos_MBoC2020 (gh: Dilan Pathirana, doi: https://doi.org/10.1091/mbc.E19-09-0548, code/data: supplement (BNGL))
# Zhao_QuantBiol2020 (gh: Iva Ewert, doi: https://doi.org/10.1007/s40484-020-0199-0, code: not available, data: table in supp)
#
# covered by performance test:
# Froehlich_CellSystems2018
#
# Unknown reasons:
# Chen_MSB2009
#

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

  # different naming scheme
  if [[ "$model" == "Bertozzi_PNAS2020" ]]; then
    yaml="${model_dir}"/"${model}"/problem.yaml
  fi

  # problems we need to flatten
  to_flatten=(
    "Bruno_JExpBot2016" "Chen_MSB2009" "Crauste_CellSystems2017"
    "Fiedler_BMCSystBiol2016" "Fujita_SciSignal2010" "SalazarCavazos_MBoC2020"
  )
  flatten=""
  for item in "${to_flatten[@]}"; do
    if [[ "$item" == "$model" ]]; then
      flatten="--flatten"
      break
    fi
  done

  amici_model_dir=test_bmc/"${model}"
  mkdir -p "$amici_model_dir"
  cmd_import="amici_import_petab ${yaml} -o ${amici_model_dir} -n ${model} ${flatten}"
  cmd_run="$script_path/test_petab_model.py -y ${yaml} -d ${amici_model_dir} -m ${model} -c"

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

cd "$script_path" && python evaluate_benchmark.py

# Test deprecated import from individual PEtab files
model="Zheng_PNAS2012"
problem_dir="${model_dir}/${model}"
amici_model_dir=test_bmc/"${model}-deprecated"
cmd_import="amici_import_petab -s "${problem_dir}/model_${model}.xml" \
  -m "${problem_dir}/measurementData_${model}.tsv" \
  -c "${problem_dir}/experimentalCondition_${model}.tsv" \
  -p "${problem_dir}/parameters_${model}.tsv" \
  -b "${problem_dir}/observables_${model}.tsv" \
  -o ${amici_model_dir} -n ${model} --no-compile"

if [[ -z "$dry_run" ]]; then
  $cmd_import
else
  echo "$cmd_import"
fi
