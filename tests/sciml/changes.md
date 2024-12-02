rename `ioId` in mapping table to `petabEntityId`
rename `mapping_table` in problem.yaml to `mapping_files` and turn into list
change `format_version` in problem.yaml to `2.0.0`
rename `model_sbml` in problem.yaml to `model_files` and turn into dict with fields location (model_sbml) and language (sbml)
change `net_files` to absolute paths
append `_o` to observable ids in observable table and measurements table to ensure uniqueness
rename `ioId` in `mapping_table` to `modelEntityId`
rename `ioValue` in `mapping_table` to `petabEntityId`
change files in `mapping_table` to absolute paths
turned `parameter_file` in problem.yaml into a list and added nn parameters
renamed `value` column in nn parameters to `nominalValue`
parameter ids in nn parameters table need to be mapped?
inputs to neural networks should have names
