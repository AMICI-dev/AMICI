"""Tests for model entities whose id collides with an AMICI-internally
-derived symbol name (a name AMICI itself constructs by string-formatting
a pattern around an existing entity's id, e.g. a reaction's flux symbol
`flux_{reaction_id}`, a conservation-law total `tcl_{state_id}`, or an
observable's measurement/sigma symbols `m{obs_id}`/`sigma_{obs_id}`),
rather than a fixed reserved name (see test_reserved_symbols.py).

See https://github.com/AMICI-dev/AMICI/issues/3240.
"""

from pathlib import Path

import numpy as np
from amici.importers.antimony import antimony2amici
from amici.importers.utils import MeasurementChannel
from amici.sim.sundials import AMICI_SUCCESS
from amici.testing import skip_on_valgrind


def _generated_source(tempdir, filename: str) -> str:
    return (Path(tempdir) / filename).read_text()


@skip_on_valgrind
def test_flux_vs_species_collision(tempdir):
    """A species literally named `flux_r1` collides with reaction `r1`'s
    auto-generated flux symbol (both would otherwise be minted as the
    exact same, `real=True` sympy symbol) -- the issue's original
    reproducer. Before the fix, this silently aliased the species' own
    state to the reaction's flux value.

    `flux_r1`'s own rate law (`flux_r1' = 0.01*flux_r1`) doesn't depend on
    `A`/`B` at all, so its trajectory has a simple closed form
    (`flux_r1(0) * exp(0.01*t)`) that only holds if it truly reads its own
    state rather than the reaction rate `0.1*A`."""
    model_name = "flux_vs_species_test"
    amici_model = antimony2amici(
        """
        model test
            compartment comp = 1;
            species A in comp = 10;
            species B in comp = 0;
            species flux_r1 in comp = 5;
            r1: A -> B; 0.1*A;
            flux_r1' = 0.01*flux_r1;
        end
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=None,
        compute_conservation_laws=False,
    )

    assert amici_model.get_state_ids() == ("A", "B", "flux_r1")
    # the flux symbol itself was disambiguated, not the user's species
    assert amici_model.get_expression_ids() == ("flux_r1_2",)

    timepoints = [0.0, 1.0, 2.0]
    amici_model.set_timepoints(timepoints)
    rdata = amici_model.simulate()
    assert rdata.status == AMICI_SUCCESS
    assert np.all(np.isfinite(rdata.x))

    expected_flux_r1 = 5.0 * np.exp(0.01 * np.array(timepoints))
    np.testing.assert_allclose(rdata.x[:, 2], expected_flux_r1, rtol=1e-6)


@skip_on_valgrind
def test_conservation_law_and_observable_symbol_collisions(tempdir):
    """Model entities literally named after AMICI's derived
    conservation-law-total (`tcl_{state}`) and observable-related
    (`m{obs}`, `sigma_{obs}`) symbols must not collide with them -- the
    derived symbol is disambiguated instead.

    Checked from generated source rather than by compiling: the
    duplicate-declaration failure this guards against is already visible
    in the generated C++ text (see test_underscore_collision_disambiguation
    in test_reserved_symbols.py for the same rationale)."""
    model_name = "cl_and_obs_collision_test"
    antimony2amici(
        """
        model test
            compartment comp = 1;
            species A in comp = 10;
            species B in comp = 0;
            const species C in comp = 3;
            r1: A -> B; 0.1*A;
            tcl_C = 1;
            mobs1 = 1;
            sigma_obs1 = 1;
        end
        """,
        model_name=model_name,
        output_dir=tempdir,
        observation_model=[MeasurementChannel(id_="obs1", formula="A")],
        compute_conservation_laws=True,
        compile=False,
    )

    # original ids of the colliding parameters are preserved
    model_cpp = _generated_source(tempdir, f"{model_name}.cpp")
    assert '"tcl_C"' in model_cpp
    assert '"mobs1"' in model_cpp
    assert '"sigma_obs1"' in model_cpp

    # the conservation-law total for species C was disambiguated, and
    # consistently so everywhere it's declared
    for filename in ("w.cpp", "x_rdata.cpp"):
        source = _generated_source(tempdir, filename)
        assert "tcl_C_2_ = tcl[0]" in source
    assert "tcl_C_ " not in _generated_source(tempdir, "w.cpp")

    # the observable's measurement/sigma symbols were disambiguated too,
    # and every use is backed by an actual declared local (not just the
    # cost-function formula referencing a name nothing declares, which
    # would be a compile error) -- regression check for the disambiguated
    # symbol also being threaded into `Observable.__init__` rather than
    # silently re-derived (undisambiguated) later by `DEModel`'s own "my"
    # array construction
    jy = _generated_source(tempdir, "Jy.cpp")
    assert "mobs1_2_ = my[0]" in jy
    assert "sigma_obs1_2_" in jy
    assert "mobs1_ " not in jy
