import sympy as sp
from amici._symbolic.de_model_components import Event
from amici.importers.utils import amici_time_symbol
from amici.testing import skip_on_valgrind


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
