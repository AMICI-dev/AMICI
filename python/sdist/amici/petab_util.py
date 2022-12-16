import petab
import re
import libsbml
import pandas as pd
import numpy as np

from typing import Optional, Union, Dict

from petab.models import MODEL_TYPE_SBML, MODEL_TYPE_PYSB
from petab.C import SIMULATION_CONDITION_ID

# ID of model parameter that is to be added to SBML model to indicate
#  preequilibration
PREEQ_INDICATOR_ID = 'preequilibration_indicator'


def resolve_mapping(element: str, mapping_df: Optional[pd.DataFrame]) -> str:
    """Resolve mapping for a given element.

    :param element:
        Element to resolve.

    :param mapping_df:
        Mapping table.

    :return:
        Resolved element.
    """
    if mapping_df is None:
        return element
    if element in mapping_df.index:
        return mapping_df.loc[element, petab.MODEL_ENTITY_ID]
    return element


def get_states_in_condition_table(
        petab_problem: petab.Problem,
        condition: Union[Dict, pd.Series] = None,
) -> Dict[str, Union[float, str]]:
    """Get list of states in the condition table"""
    if petab_problem.model.type_id not in (MODEL_TYPE_SBML, MODEL_TYPE_PYSB):
        raise NotImplementedError()
    species_check_funs = {
        MODEL_TYPE_SBML:
            lambda x: _element_is_sbml_state(petab_problem.sbml_model, x),
        MODEL_TYPE_PYSB:
            lambda x: _element_is_pysb_pattern(petab_problem.model.model, x),
    }
    states = {
        resolve_mapping(col, petab_problem.mapping_df):
            np.NaN if condition is None else
            petab_problem.condition_df.loc[
                condition[SIMULATION_CONDITION_ID], col
            ]
        for col in petab_problem.condition_df
        if species_check_funs[petab_problem.model.type_id](
            resolve_mapping(col, petab_problem.mapping_df)
        ) and (condition is None or not pd.isna(
            petab_problem.condition_df.loc[
                condition[SIMULATION_CONDITION_ID], col
            ]
        ))
    }
    if petab_problem.model.type_id == MODEL_TYPE_PYSB:
        import pysb.pattern
        spm = pysb.pattern.SpeciesPatternMatcher(
            model=petab_problem.model.model
        )

        # expose model components as variables so we can evaluate patterns
        for c in petab_problem.model.model.components:
            globals()[c.name] = c

        states = {
            f'__s{ix}': value
            for pattern, value in states.items()
            for ix in spm.match(eval(pattern), index=True)
        }
    return states


def _element_is_pysb_pattern(model: 'pysb.Model', element: str) -> bool:
    """Check if element is a pysb pattern"""
    if match := re.match(r'[a-zA-Z_][\w_]*\(', element):
        return match.group(0)[:-1] in [
            m.name for m in model.monomers
        ]
    return False


def _element_is_sbml_state(sbml_model: libsbml.Model, sbml_id: str) -> bool:
    """Does the element with ID `sbml_id` correspond to a state variable?"""
    if sbml_model.getCompartment(sbml_id) is not None:
        return True
    if sbml_model.getSpecies(sbml_id) is not None:
        return True
    if (rule := sbml_model.getRuleByVariable(sbml_id)) is not None \
            and rule.getTypeCode() == libsbml.SBML_RATE_RULE:
        return True

    return False
