"""Tests for model entities whose id collides with a C++ keyword, a
standard-library macro, or one of amici's own reserved argument names,
and for the time-symbol."""

from pathlib import Path

import numpy as np
import pytest
from amici import import_model_module
from amici.importers.antimony import antimony2amici, antimony2sbml
from amici.importers.utils import MeasurementChannel
from amici.sim.sundials import AMICI_SUCCESS
from amici.testing import skip_on_valgrind


def _antimony_model_with_species(species_ids: list[str]) -> str:
    """An Antimony model with one compartment and one species per given id,
    each following ``d[id]/dt = -0.1 * [id]``, starting at 1."""
    lines = ["model test", "    compartment comp = 1;"]
    lines += [f"    species {sid} in comp = 1;" for sid in species_ids]
    lines += [f"    {sid}' = -0.1*{sid};" for sid in species_ids]
    lines.append("end")
    return "\n".join(lines)


def _generated_source(tempdir, filename: str) -> str:
    return (Path(tempdir) / filename).read_text()


def _without_comments(source: str) -> str:
    """Strip trailing ``//`` comments so `__`-in-comment (original id,
    kept human-readable) isn't mistaken for `__` in an actual identifier."""
    return "\n".join(line.split("//")[0] for line in source.splitlines())


# Names that are unsafe in the generated C++ -- amici's own fixed
# array-parameter names (x, p, k, h, w, y), real C++ keywords and known
# standard-library macros (int, class, EOF, INFINITY, NULL) -- plus the
# underscore edge cases from #2226's design. None of these is reserved at
# the model level: they're all handled by the code printer mangling every
# identifier it prints.
# (Not `NAN`/`NaN`/`nan`: Antimony reserves those spellings as its own
# pre-defined floating-point-NaN constant, rejected at parse time regardless
# of context -- `NULL`/`EOF`/`INFINITY` already cover this class of case.)
RESERVED_SPECIES_IDS = [
    "x",
    "p",
    "k",
    "h",
    "w",
    "y",
    "NULL",
    "int",
    "class",
    "EOF",
    "INFINITY",
    "k_",
    "my__species",
]

# The same for the generated JAX module: amici's fixed array-parameter
# names, every other argument and local the generated methods use, the
# module-level imports they call into, and Python keywords. Plus `t`, the
# one name that *is* still renamed at the model level (see RESERVED_SYMBOLS)
# and must be reported back under its original id.
JAX_UNSAFE_SPECIES_IDS = [
    # amici's fixed array-parameter names
    "x",
    "p",
    "k",
    "h",
    "w",
    "y",
    # further arguments/locals of the generated methods
    "t",
    "tcl",
    "op",
    "np",
    "my",
    "iy",
    "args",
    "self",
    # module-level imports the generated code calls into
    "jnp",
    "jr",
    "eqx",
    "oo",
    "safe_log",
    "safe_div",
    "Path",
    "JAXModel",
    # Python keywords
    "class",
    "lambda",
    "None",
    "def",
    "return",
]


@skip_on_valgrind
def test_reserved_species_ids(tempdir):
    """Species using names that used to be rejected, silently renamed, or
    unsafe as C++ identifiers import, compile, and simulate correctly,
    reporting their original ids.

    All ids share a single model/compilation -- compiling one model per id
    would be needlessly expensive, and each id's rate rule is independent of
    every other's, so nothing is lost by combining them."""
    amici_model = antimony2amici(
        _antimony_model_with_species(RESERVED_SPECIES_IDS),
        model_name="reserved_species_test",
        output_dir=tempdir,
        observation_model=None,
        compute_conservation_laws=False,
    )

    assert amici_model.get_state_ids() == tuple(RESERVED_SPECIES_IDS)

    amici_model.set_timepoints([0, 1, 2])
    rdata = amici_model.simulate()
    assert rdata.status == AMICI_SUCCESS
    assert np.all(np.isfinite(rdata.x))


@skip_on_valgrind
def test_underscore_collision_disambiguation(tempdir):
    """Two species whose ids only differ in underscore-run length collapse
    to the same C++ identifier base and must be disambiguated, not aliased
    to the same state.

    This is a property of the generated source (are there two distinct,
    `__`-free declared locals, and does each retain its own public id), not
    of runtime behavior -- checking the generated code directly is just as
    conclusive as compiling and avoids the extra build."""
    model_name = "underscore_collision_test"
    antimony2amici(
        _antimony_model_with_species(["a__b", "a_b"]),
        model_name=model_name,
        output_dir=tempdir,
        observation_model=None,
        compute_conservation_laws=False,
        compile=False,
    )

    xdot = _generated_source(tempdir, "xdot.cpp")
    assert "__" not in _without_comments(xdot)
    # two distinct locals were actually declared, not one silently reused
    # for both (which a naming collision, rather than disambiguation,
    # would produce)
    assert "a_b_ =" in xdot
    assert "a_b_2 =" in xdot

    model_cpp = _generated_source(tempdir, f"{model_name}.cpp")
    assert '"a__b"' in model_cpp
    assert '"a_b"' in model_cpp


@skip_on_valgrind
def test_species_named_t(tempdir):
    """A species literally named `t` must not be confused with simulation
    time (#2461): it keeps its own id everywhere it's reported, and its
    generated rate law reads from its own state slot rather than from the
    time argument.

    Checked from generated source rather than by compiling and simulating:
    the bug this guards against (`_process_time()` running before the
    pre-import rename, merging the species into the time symbol) would mean
    there's no `x[0]`-bound state for `t` at all, or import/codegen raising
    outright -- both are visible without a build."""
    model_name = "species_named_t_test"
    antimony2amici(
        _antimony_model_with_species(["t"]),
        model_name=model_name,
        output_dir=tempdir,
        observation_model=None,
        compute_conservation_laws=False,
        compile=False,
    )

    # the public id must be exactly "t", never "amici_t" -- the #2461 fix
    model_cpp = _generated_source(tempdir, f"{model_name}.cpp")
    assert '"t", // x_rdata[0]' in model_cpp

    # the state is bound to its own array slot (x[0]), i.e. it's a real,
    # distinct state -- not merged away into the time argument `t`
    xdot = _generated_source(tempdir, "xdot.cpp")
    assert "amici_t_ = x[0];" in xdot
    assert "-1.0/10.0*amici_t_" in xdot


@skip_on_valgrind
def test_reserved_species_ids_jax(tempdir):
    """JAX counterpart of ``test_reserved_species_ids``: species whose ids
    are unsafe as identifiers in the generated JAX module import and
    evaluate correctly, and report their original id.

    The generated methods destructure their array arguments into per-entry
    locals, so an entity id that isn't mangled would shadow the argument it
    was unpacked from (or a module-level import, or be a syntax error) --
    hence checking that the right-hand sides still evaluate, not just that
    the id lists look right."""
    pytest.importorskip("jax")
    import jax.numpy as jnp
    from amici.importers.sbml import SbmlImporter

    species_ids = JAX_UNSAFE_SPECIES_IDS
    model_name = "reserved_species_test_jax"
    sbml_str = antimony2sbml(_antimony_model_with_species(species_ids))
    SbmlImporter(sbml_str, from_file=False).sbml2jax(
        model_name,
        output_dir=Path(tempdir) / model_name,
        observation_model=[],
        compute_conservation_laws=False,
    )
    model = import_model_module(model_name, tempdir).Model()

    assert tuple(model.state_ids) == tuple(species_ids)

    # each species follows d[id]/dt = -0.1 * [id], starting at 1
    x0 = model._x0(0.0, model.parameters)
    empty = jnp.array([])
    xdot = model._xdot(0.0, x0, (model.parameters, empty, empty))
    assert np.allclose(x0, 1.0)
    assert np.allclose(xdot, -0.1)


@skip_on_valgrind
def test_observable_named_my_jax(tempdir):
    """An observable named `my` must not shadow the generated `_nllh`'s
    measurement argument, which is called exactly that.

    Unlike a shadowed state, this one produces no error at all: the
    negative log-likelihood would just silently be computed against the
    simulated observable instead of the measurement, and come out
    independent of the data."""
    pytest.importorskip("jax")
    import jax.numpy as jnp
    from amici.importers.sbml import SbmlImporter

    model_name = "observable_named_my_test_jax"
    sbml_str = antimony2sbml(_antimony_model_with_species(["A"]))
    SbmlImporter(sbml_str, from_file=False).sbml2jax(
        model_name,
        output_dir=Path(tempdir) / model_name,
        observation_model=[
            MeasurementChannel(
                id_="my",
                formula="A",
                noise_distribution="normal",
                sigma=1.0,
            )
        ],
        compute_conservation_laws=False,
    )
    model = import_model_module(model_name, tempdir).Model()

    assert tuple(model.observable_ids) == ("my",)

    x0 = model._x0(0.0, model.parameters)  # A(0) = 1, so `my` = 1
    empty = jnp.array([])

    def nllh(measurement):
        return model._nllh(
            0.0,
            x0,
            model.parameters,
            empty,
            empty,
            jnp.array([measurement]),
            0,
            empty,
            empty,
        )

    # normal noise with sigma=1: 0.5*(y - my)**2 + 0.5*log(2*pi)
    assert np.isclose(nllh(1.0), 0.5 * np.log(2 * np.pi))
    assert np.isclose(nllh(5.0), 0.5 * 4**2 + 0.5 * np.log(2 * np.pi))
