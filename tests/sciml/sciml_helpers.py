import pandas as pd
from amici.sim.jax.petab import _try_float
from petab import v1, v2


def _process_prior_params(prior_params):
    if isinstance(prior_params, float):
        return prior_params
    else:
        return [float(param) for param in prior_params.split(";")]


def _v2_sciml_problem_helper(yaml_config, base_path):
    config = v2.ProblemConfig(**yaml_config)

    parameter_tables = []
    for f in config.parameter_files:
        df = pd.read_csv(f, sep="\t")
        df.nominalValue = df.nominalValue.apply(_try_float)
        if "priorParameters" in df.columns:
            df.priorParameters = df.priorParameters.apply(
                _process_prior_params
            )
        parameters = [
            v2.Parameter.model_construct(**row.to_dict())
            for _, row in df.reset_index().iterrows()
        ]
        parameter_tables.append(v2.ParameterTable(elements=parameters))

    models = [
        v1.models.model.model_factory(
            model_info.location,
            base_path=base_path,
            model_language=model_info.language,
            model_id=model_id,
        )
        for model_id, model_info in (config.model_files or {}).items()
    ]

    measurement_tables = (
        [
            v2.MeasurementTable.from_tsv(f, base_path)
            for f in config.measurement_files
        ]
        if config.measurement_files
        else None
    )

    experiment_tables = (
        [
            v2.ExperimentTable.from_tsv(f, base_path)
            for f in config.experiment_files
        ]
        if config.experiment_files
        else None
    )

    condition_tables = (
        [
            v2.ConditionTable.from_tsv(f, base_path)
            for f in config.condition_files
        ]
        if config.condition_files
        else None
    )

    if condition_tables is None:
        cond_ids = [
            cid
            for exp_table in experiment_tables
            for exp in exp_table.elements
            for p in exp.periods
            for cid in p.condition_ids
        ]
        condition_tables = [
            v2.ConditionTable(elements=[v2.Condition(id=cid, changes=[])])
            for cid in set(cond_ids)
        ]

    observable_tables = (
        [
            v2.ObservableTable.from_tsv(f, base_path)
            for f in config.observable_files
        ]
        if config.observable_files
        else None
    )

    mapping_tables = (
        [v2.MappingTable.from_tsv(f, base_path) for f in config.mapping_files]
        if config.mapping_files
        else None
    )

    return v2.Problem(
        config=config,
        models=models,
        condition_tables=condition_tables,
        experiment_tables=experiment_tables,
        observable_tables=observable_tables,
        measurement_tables=measurement_tables,
        parameter_tables=parameter_tables,
        mapping_tables=mapping_tables,
    )
