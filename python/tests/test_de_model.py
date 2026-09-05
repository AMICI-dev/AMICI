import libsbml
import pytest
import sympy as sp
from amici import MeasurementChannel, SbmlImporter
from amici._symbolic.de_model_components import Event, FreeParameter
from amici.importers.antimony import antimony2sbml
from amici.importers.sbml.splines import CubicHermiteSpline
from amici.importers.utils import amici_time_symbol
from amici.testing import skip_on_valgrind


def _build_spline_model(*, with_spline_dependent_event: bool):
    """Build a minimal `DEModel` with a single state and a single spline,
    optionally with an event whose trigger depends on the spline's value.

    Kept intentionally tiny and built via `_build_ode_model` directly (no
    codegen/compilation) so these tests run fast.
    """
    event_line = (
        "e1: at u > 1: x = x + 1;\n" if with_spline_dependent_event else ""
    )
    ant_str = rf"""
    model spline_model
        p1 = 1;
        sp1 = 0;
        sp2 = 2;
        species x = 1;
        x' = -p1*x + u;
        var u = 0;
        {event_line}
    end
    """
    sbml_str = antimony2sbml(ant_str)
    sbml_model = libsbml.SBMLReader().readSBMLFromString(sbml_str).getModel()

    spline = CubicHermiteSpline(
        sbml_id="u",
        evaluate_at=amici_time_symbol,
        nodes=[0, 10],
        values_at_nodes=[sp.Symbol("sp1"), sp.Symbol("sp2")],
        extrapolate=("constant", "constant"),
    )
    spline.add_to_sbml_model(sbml_model, auto_add=False)

    importer = SbmlImporter(sbml_model)
    model = importer._build_ode_model(
        observation_model=[MeasurementChannel(id_="obs_x", formula="x")],
    )
    model.generate_basic_variables()
    return model


@skip_on_valgrind
def test_spline_static_indices_and_substitution():
    """Splines occur in the model as `AmiciSpline`/`AmiciSplineSensitivity`
    sympy `Function` calls, which must be substituted for `spl`/`sspl`
    symbols so `static_indices()` can classify spline-dependent rows as
    dynamic without falling back to string matching (see the FIXMEs this
    replaces)."""
    model = _build_spline_model(with_spline_dependent_event=False)

    # the spline's own row in `w` is time-varying and must not be static
    w = model.eq("w")
    static_w = set(model.static_indices("w"))
    spline_row = next(
        i for i, sym in enumerate(model.sym("w")) if str(sym) == "u"
    )
    assert spline_row not in static_w

    # no raw spline Function calls should survive substitution anywhere
    for name in ("w", "dwdx", "dwdw", "dwdp"):
        eq = model.eq(name)
        for entry in eq:
            assert "AmiciSpline" not in str(entry)

    # dwdp's entry for the spline row references the sensitivity symbols
    assert set(model.sym("sspl")) & w[spline_row].free_symbols == set()
    dwdp_spline_row = model.eq("dwdp")[spline_row, :]
    assert set(model.sym("sspl")) & dwdp_spline_row.free_symbols


@skip_on_valgrind
def test_spline_derivative_in_event_not_supported():
    """`AmiciSplineDerivative` (a spline's time derivative) has no C++
    codegen support yet. A model where an event trigger depends on a
    time-varying spline must fail loudly at import time instead of
    silently generating C++ that references an undefined `dspl_N`."""
    model = _build_spline_model(with_spline_dependent_event=True)

    with pytest.raises(NotImplementedError, match="time derivative"):
        model.eq("drootdt_total")


@skip_on_valgrind
def test_model_quantity_reserved_name():
    """`t` and the fixed array-parameter names (x, p, k, h, w, y) are
    reserved for ModelQuantity symbols: the JAX exporter has no mangling of
    its own and relies on these being renamed before codegen ever sees
    them (unlike the C++ backend, whose printer mangles every identifier).
    A name that only collides with a C++ keyword/macro (e.g. NULL) is
    handled by that mangling instead and doesn't need to be rejected
    here."""
    FreeParameter(symbol=sp.Symbol("NULL"), name="NULL", value=1.0)

    for name in ("t", "x", "p", "k", "h", "w", "y"):
        with pytest.raises(ValueError, match="Cannot add"):
            FreeParameter(symbol=sp.Symbol(name), name=name, value=1.0)


@skip_on_valgrind
def test_event_trigger_time():
    e = Event(
        symbol=sp.Symbol("event1"),
        name="event name",
        value=amici_time_symbol - 10,
        assignments=sp.Float(1),
        use_values_from_trigger_time=False,
    )
    assert e.triggers_at_fixed_timepoint() is True
    assert e.get_trigger_time() == 10

    # fixed, but multiple timepoints - not (yet) supported
    e = Event(
        symbol=sp.Symbol("event1"),
        name="event name",
        value=sp.sin(amici_time_symbol),
        assignments=sp.Float(1),
        use_values_from_trigger_time=False,
    )
    assert e.triggers_at_fixed_timepoint() is False

    e = Event(
        symbol=sp.Symbol("event1"),
        name="event name",
        value=amici_time_symbol / 2,
        assignments=sp.Float(1),
        use_values_from_trigger_time=False,
    )
    assert e.triggers_at_fixed_timepoint() is True
    assert e.get_trigger_time() == 0

    # parameter-dependent triggers - not (yet) supported
    e = Event(
        symbol=sp.Symbol("event1"),
        name="event name",
        value=amici_time_symbol - sp.Symbol("delay"),
        assignments=sp.Float(1),
        use_values_from_trigger_time=False,
    )
    assert e.triggers_at_fixed_timepoint() is False
