"""
PEtab Import
------------
Import a model in the :mod:`petab` (https://github.com/PEtab-dev/PEtab) format
into AMICI.
"""

import logging
import os
import re
import shutil
from pathlib import Path

import pandas as pd
import petab.v1 as petab
import petab.v2 as petabv2
from petab.v1.models import MODEL_TYPE_PYSB, MODEL_TYPE_SBML

import amici
from amici.logging import get_logger
from amici.sim.jax import reformat_petab_v2_to_v1, add_events_to_sbml

from .import_helpers import (
    _can_import_model,
    _create_model_name,
    _get_package_name_and_path,
    check_model,
)
from .sbml_import import import_model_sbml

try:
    from .pysb_import import import_model_pysb
except ModuleNotFoundError:
    # pysb not available
    import_model_pysb = None


__all__ = ["import_petab_problem"]

logger = get_logger(__name__, logging.WARNING)


def import_petab_problem(
    petab_problem: petab.Problem,
    output_dir: str | Path | None = None,
    *,
    model_name: str = None,
    compile_: bool = None,
    non_estimated_parameters_as_constants=True,
    jax=False,
    **kwargs,
) -> "amici.sim.sundials.Model | amici.jax.JAXProblem":
    """
    Create an AMICI model for a PEtab problem.

    :param petab_problem:
        A petab problem containing all relevant information on the model.

    :param output_dir:
        Directory to write the model code to. It will be created if it doesn't
        exist. Defaults to :func:`amici.get_model_dir`.

    :param model_name:
        Name of the generated model module. Defaults to the ID of the model
        or the model file name without the extension.

    :param compile_:
        If ``True``, the model will be compiled. If ``False``, the model will
        not be compiled. If ``None``, the model will be compiled if it cannot
        be imported.

    :param non_estimated_parameters_as_constants:
        Whether parameters marked as non-estimated in PEtab should be
        considered constant in AMICI. Setting this to ``True`` will reduce
        model size and simulation times. If sensitivities with respect to those
        parameters are required, this should be set to ``False``.

    :param jax:
        Whether to create a JAX-based problem. If ``True``, returns a
        :class:`amici.jax.JAXProblem` instance. If ``False``, returns a
        standard AMICI model.

    :param kwargs:
        Additional keyword arguments to be passed to
        :meth:`amici.importers.sbml.SbmlImporter.sbml2amici` or
        :func:`amici.importers.pysb.pysb2amici`, depending on the model type.

    :return:
        The imported model (if ``jax=False``) or JAX problem (if ``jax=True``).
    """
    if isinstance(petab_problem, petabv2.Problem):
        petab_problem = add_events_to_sbml(petab_problem)
        petab_problem = reformat_petab_v2_to_v1(petab_problem)

    if petab_problem.model.type_id not in (MODEL_TYPE_SBML, MODEL_TYPE_PYSB):
        raise NotImplementedError(
            "Unsupported model type " + petab_problem.model.type_id
        )

    model_name = model_name or petab_problem.model.model_id

    if petab_problem.model.type_id == MODEL_TYPE_PYSB and model_name is None:
        model_name = petab_problem.pysb_model.name
    elif model_name is None and output_dir:
        model_name = _create_model_name(output_dir)

    # generate folder and model name if necessary
    if output_dir is None:
        output_dir = amici.get_model_dir(model_name, jax=jax).absolute()
    else:
        output_dir = Path(output_dir).absolute()

    output_dir.mkdir(parents=True, exist_ok=True)

    # check if compilation necessary
    if compile_ or (
        compile_ is None and not _can_import_model(model_name, output_dir, jax)
    ):
        # check if folder exists
        if os.listdir(output_dir) and not compile_:
            raise ValueError(
                f"Cannot compile to {output_dir}: not empty. "
                "Please assign a different target or set `compile_` to `True`."
            )

        # remove folder if exists
        if not jax and os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        logger.info(f"Compiling model {model_name} to {output_dir}.")

        if "sciml" in petab_problem.extensions_config:
            from petab_sciml.standard import NNModelStandard

            config = petab_problem.extensions_config["sciml"]
            # TODO: only accept YAML format for now
            hybridizations = [
                pd.read_csv(hf, sep="\t")
                for hf in config["hybridization_files"]
            ]
            hybridization_table = pd.concat(hybridizations)

            input_mapping = dict(
                zip(
                    hybridization_table["targetId"],
                    hybridization_table["targetValue"],
                )
            )
            output_mapping = dict(
                zip(
                    hybridization_table["targetValue"],
                    hybridization_table["targetId"],
                )
            )
            observable_mapping = dict(
                zip(
                    petab_problem.observable_df["observableFormula"],
                    petab_problem.observable_df.index,
                )
            )
            hybridization = {
                net_id: {
                    "model": NNModelStandard.load_data(
                        Path(net_config["location"])
                    ),
                    "input_vars": [
                        input_mapping[petab_id]
                        for petab_id, model_id in petab_problem.mapping_df.loc[
                            petab_problem.mapping_df[petab.MODEL_ENTITY_ID]
                            .str.split(".")
                            .str[0]
                            == net_id,
                            petab.MODEL_ENTITY_ID,
                        ]
                        .to_dict()
                        .items()
                        if model_id.split(".")[1].startswith("input")
                        and petab_id in input_mapping.keys()
                    ],
                    "output_vars": {
                        output_mapping[petab_id]: _get_net_index(model_id)
                        for petab_id, model_id in petab_problem.mapping_df.loc[
                            petab_problem.mapping_df[petab.MODEL_ENTITY_ID]
                            .str.split(".")
                            .str[0]
                            == net_id,
                            petab.MODEL_ENTITY_ID,
                        ]
                        .to_dict()
                        .items()
                        if model_id.split(".")[1].startswith("output")
                        and petab_id in output_mapping.keys()
                    },
                    "observable_vars": {
                        observable_mapping[petab_id]: _get_net_index(model_id)
                        for petab_id, model_id in petab_problem.mapping_df.loc[
                            petab_problem.mapping_df[petab.MODEL_ENTITY_ID]
                            .str.split(".")
                            .str[0]
                            == net_id,
                            petab.MODEL_ENTITY_ID,
                        ]
                        .to_dict()
                        .items()
                        if model_id.split(".")[1].startswith("output")
                        and petab_id in observable_mapping.keys()
                    },
                    "frozen_layers": dict(
                        [
                            _get_frozen_layers(model_id)
                            for petab_id, model_id in petab_problem.mapping_df.loc[
                                petab_problem.mapping_df[petab.MODEL_ENTITY_ID]
                                .str.split(".")
                                .str[0]
                                == net_id,
                                petab.MODEL_ENTITY_ID,
                            ]
                            .to_dict()
                            .items()
                            if petab_id in petab_problem.parameter_df.index
                            and petab_problem.parameter_df.loc[
                                petab_id, petab.ESTIMATE
                            ]
                            == 0
                        ]
                    ),
                    **net_config,
                }
                for net_id, net_config in config["neural_nets"].items()
            }
            if not jax or petab_problem.model.type_id != MODEL_TYPE_SBML:
                raise NotImplementedError(
                    "petab_sciml extension is currently only supported for sbml models"
                )
        else:
            hybridization = None

        # compile the model
        if petab_problem.model.type_id == MODEL_TYPE_PYSB:
            import_model_pysb(
                petab_problem,
                model_name=model_name,
                output_dir=output_dir,
                jax=jax,
                **kwargs,
            )
        else:
            import_model_sbml(
                petab_problem=petab_problem,
                model_name=model_name,
                output_dir=output_dir,
                non_estimated_parameters_as_constants=non_estimated_parameters_as_constants,
                hybridization=hybridization,
                jax=jax,
                **kwargs,
            )

    # import model
    model_module = amici.import_model_module(
        *_get_package_name_and_path(model_name, output_dir, jax=jax)
    )

    if jax:
        from amici.jax import JAXProblem

        model = model_module.Model()

        logger.info(
            f"Successfully loaded jax model {model_name} from {output_dir}."
        )

        # Create and return JAXProblem
        logger.info(f"Successfully created JAXProblem for {model_name}.")
        return JAXProblem(model, petab_problem)

    model = model_module.get_model()
    check_model(amici_model=model, petab_problem=petab_problem)

    logger.info(f"Successfully loaded model {model_name} from {output_dir}.")

    return model


def _get_net_index(model_id: str):
    matches = re.findall(r"\[(\d+)\]", model_id)
    if matches:
        return int(matches[-1])


def _get_frozen_layers(model_id):
    layers = re.findall(r"\[(.*?)\]", model_id)
    array_attr = model_id.split(".")[-1]
    layer_id = layers[0] if len(layers) else None
    array_attr = array_attr if array_attr in ("weight", "bias") else None
    return layer_id, array_attr
