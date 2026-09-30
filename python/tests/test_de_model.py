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


@skip_on_valgrind
def test_event_trigger_time_minmax():
    """``Min``/``Max``-based event triggers (as generated for ``And``/``Or``
    SBML triggers, and for PEtab v2 period-start events) could not
    previously be solved for ``t`` at all (AMICI#3126), since
    ``sympy.solve`` raises ``NotImplementedError`` on ``Min``/``Max``."""
    t = amici_time_symbol
    preeq, e0 = sp.symbols("preeq e0")
    static_syms = {preeq, e0}

    def make_event(value):
        return Event(
            symbol=sp.Symbol("event1"),
            name="event name",
            value=value,
            assignments=sp.Float(1),
            use_values_from_trigger_time=False,
        )

    # the exact trigger from the issue / PEtab v2 period-start events:
    # `And(preeq <= 1/2, e0 >= 1/2, t >= 10)`
    e = make_event(
        sp.Min(sp.Rational(1, 2) - preeq, e0 - sp.Rational(1, 2), t - 10)
    )
    assert e.has_explicit_trigger_times(static_syms)
    (t_root,) = e.get_trigger_times()
    expected = {
        (0, 0): sp.oo,
        (0, 1): 10,
        (1, 0): sp.oo,
        (1, 1): sp.oo,
    }
    for (preeq_v, e0_v), exp in expected.items():
        assert t_root.subs({preeq: preeq_v, e0: e0_v}) == exp

    # a pure gating `Min` (no time-dependence at all, as for the
    # preequilibration-period trigger) has no trigger *time* to compute
    e = make_event(sp.Min(e0 - sp.Rational(1, 2), preeq - sp.Rational(1, 2)))
    assert not e.has_explicit_trigger_times(static_syms)

    # `Or` of purely time-dependent, differently-scaled conditions:
    # triggers at the earliest
    e = make_event(sp.Max(t - 5, 2 * t - 20))
    assert e.has_explicit_trigger_times(static_syms)
    assert e.get_trigger_times() == {5}

    # `Max`/`Or` with a gating (time-independent) leaf is out of scope
    e = make_event(sp.Max(t - 5, e0 - sp.Rational(1, 2)))
    assert not e.has_explicit_trigger_times(static_syms)

    # decreasing leaves (e.g. from `t < 10`) are fine as long as every
    # time-dependent leaf goes the same way: the `Min` starts non-negative
    # and drops below zero as soon as the *first* condition fails
    e = make_event(sp.Min(10 - t, 20 - 2 * t, e0 - sp.Rational(1, 2)))
    assert e.has_explicit_trigger_times(static_syms)
    (t_root,) = e.get_trigger_times()
    assert t_root.subs({e0: 1}) == 10
    assert t_root.subs({e0: 0}) == sp.oo

    # ... but mixed monotonicity is not, since the expression may then
    # cross zero more than once
    e = make_event(sp.Min(t - 5, 10 - t))
    assert not e.has_explicit_trigger_times(static_syms)

    # a trigger and its negated counterpart (AMICI tracks both for
    # persisted Heaviside variables) must resolve to the same time, or
    # `_reorder_events` would separate the pair (AMICI#3126 follow-up)
    trigger = sp.Min(sp.Rational(1, 2) - preeq, e0 - sp.Rational(1, 2), t - 10)
    pos, neg = make_event(trigger), make_event(-trigger)
    assert pos.has_explicit_trigger_times(static_syms)
    assert neg.has_explicit_trigger_times(static_syms)
    assert pos.get_trigger_times() == neg.get_trigger_times()

    # a non-affine leaf is out of scope
    e = make_event(sp.Min(t**2 - 100, e0 - sp.Rational(1, 2)))
    assert not e.has_explicit_trigger_times(static_syms)

    # nested Min-inside-Max (or vice versa) is not (yet) supported
    e = make_event(sp.Max(t - 5, sp.Min(t - 10, e0 - sp.Rational(1, 2))))
    assert not e.has_explicit_trigger_times(static_syms)


def _build_event_model(trigger: str):
    """Build a minimal one-state, one-observable `DEModel` with a single
    event using the given Antimony trigger expression."""
    ant_str = rf"""
    model m
        compartment_ = 1;
        species x = 1;
        y = 0;
        p1 = 10;
        x' = -x;
        at ({trigger}): y = y + 1;
    end
    """
    sbml_str = antimony2sbml(ant_str)
    sbml_model = libsbml.SBMLReader().readSBMLFromString(sbml_str).getModel()
    model = SbmlImporter(sbml_model)._build_ode_model(
        observation_model=[MeasurementChannel(id_="obs_y", formula="y")],
    )
    model.generate_basic_variables()
    return model


@skip_on_valgrind
def test_event_trigger_time_state_dependence_stays_implicit():
    """A trigger comparing `time` directly to a state (`time >= x`) is
    solved by plain `sympy.solve` to a state-dependent expression
    (`solve(t - x, t) == [x]`) -- `x` is a dynamic state, not a static
    parameter, so this must not make the event count as having an explicit
    (precomputable) trigger time anywhere. Before the fix,
    `DEModel._reorder_events` (physical event/Heaviside-array order) and
    the JAX exporter's `iroot`/`eroot`/`ih`/`eh` split (root-detection
    order) disagreed on this, permuting the two orderings relative to each
    other (AMICI#3286)."""
    model = _build_event_model("time >= x")
    (event,) = model.events()

    # the trigger time *is* solved, but still references the state `x`
    (t_root,) = event.get_trigger_times()
    assert model.sym("x")[0] in t_root.free_symbols

    assert model.num_events_solver() == 1
    assert len(model.eq("iroot")) == 1
    assert len(model.eq("eroot")) == 0
    assert model.sym("ih").shape == (1, 1)
    assert model.sym("eh").shape == (0, 0)
    assert model.get_explicit_roots() == []
    assert len(model.get_implicit_roots()) == 1


@skip_on_valgrind
def test_event_trigger_time_purely_static_is_explicit():
    """The static-parameter counterpart of the above: a trigger comparing
    `time` to a genuine (fixed) parameter must be classified explicit
    everywhere, consistently."""
    model = _build_event_model("time >= p1")
    (event,) = model.events()

    assert model.num_events_solver() == 0
    assert len(model.eq("iroot")) == 0
    assert len(model.eq("eroot")) == 1
    assert model.sym("ih").shape == (0, 0)
    assert model.sym("eh").shape == (1, 1)
    assert len(model.get_explicit_roots()) == 1
    assert model.get_implicit_roots() == []


@skip_on_valgrind
def test_event_ordering_consistent_across_classifications():
    """Direct regression test for AMICI#3286: the physical event order
    (`DEModel._reorder_events`, used for `h`/`event_initial_values`/
    `deltax`) must never disagree with the order `iroot`++`eroot`
    reconstructs (used for the JAX backend's root-detection vector) --
    otherwise root crossings get attributed to the wrong Heaviside slot.
    Uses the issue's own reproduction shape: a state-vs-time trigger
    (solvable but not static) alongside a periodic, time-only trigger
    (genuinely unsolvable, agreeing under any criterion)."""
    ant_str = r"""
    model m
        compartment_ = 1;
        species x1 = 20;
        y = 0;
        x1' = -x1;
        ev_state: at (time >= x1): y = y + 1;
        ev_periodic: at (sin(time) > 0.5): y = y + 2;
    end
    """
    sbml_str = antimony2sbml(ant_str)
    sbml_model = libsbml.SBMLReader().readSBMLFromString(sbml_str).getModel()
    model = SbmlImporter(sbml_model)._build_ode_model(
        observation_model=[MeasurementChannel(id_="obs_y", formula="y")],
    )
    model.generate_basic_variables()

    static_syms = model.static_symbols
    iroot_events = [
        e
        for e in model.events()
        if not e.has_explicit_trigger_times(static_syms)
    ]
    eroot_events = [
        e for e in model.events() if e.has_explicit_trigger_times(static_syms)
    ]
    assert [e.get_sym() for e in iroot_events + eroot_events] == [
        e.get_sym() for e in model.events()
    ]


@skip_on_valgrind
def test_has_implicit_event_assignments():
    """`has_implicit_event_assignments` gates JAX export of state-updating
    events. What it must reject is a trigger whose firing time depends on
    the trajectory; a trigger referencing only time and parameters is
    fine, whether or not it can be solved for `t`."""
    for trigger in ("time >= 10", "time >= p1", "sin(time) > 0.5"):
        model = _build_event_model(trigger)
        assert not model.has_implicit_event_assignments(), trigger

    for trigger in ("time >= x", "x > 0.5"):
        model = _build_event_model(trigger)
        assert model.has_implicit_event_assignments(), trigger
