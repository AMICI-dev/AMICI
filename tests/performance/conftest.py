"""Pytest configuration for JAX regression tests.

Provides:
- ``--results-dir`` option: path where jax_regression_results.json is written.
- Session-scoped solver/model fixtures shared across tests.
- ``results_collector`` fixture that accumulates metrics and writes JSON at end.
"""

import json
import platform
from pathlib import Path

import diffrax
import jax
import optimistix
import pytest

# ── Pytest option ──────────────────────────────────────────────────────────


def pytest_addoption(parser):
    parser.addoption(
        "--results-dir",
        default=None,
        metavar="DIR",
        help="Directory in which to write jax_regression_results.json",
    )


# ── Session initialisation ────────────────────────────────────────────────


@pytest.fixture(scope="session", autouse=True)
def configure_jax():
    jax.config.update("jax_enable_x64", True)


# ── Default solver configuration ──────────────────────────────────────────

SOLVER_KWARGS = dict(
    solver=diffrax.Kvaerno5(),
    controller=diffrax.PIDController(
        atol=1e-8, rtol=1e-8, pcoeff=0.4, icoeff=0.3, dcoeff=0.0
    ),
    root_finder=optimistix.Newton(atol=1e-12, rtol=1e-12),
    adjoint=diffrax.RecursiveCheckpointAdjoint(),
    steady_state_event=diffrax.steady_state_event(),
    max_steps=2**14,
)


@pytest.fixture(scope="session")
def solver_kwargs():
    return SOLVER_KWARGS


# ── Session-scoped model instances ────────────────────────────────────────


@pytest.fixture(scope="session")
def tier1_models():
    """Instantiate all synthetic Tier 1 models once per session."""
    from tests.performance.synthetic_models.conservation_law import (
        ConservationLaw,
    )
    from tests.performance.synthetic_models.linear_decay import LinearDecay
    from tests.performance.synthetic_models.lotka_volterra import LotkaVolterra
    from tests.performance.synthetic_models.multi_event import MultiEvent
    from tests.performance.synthetic_models.robertson import Robertson
    from tests.performance.synthetic_models.single_event import SingleEvent

    return {
        "LinearDecay": LinearDecay(),
        "Robertson": Robertson(),
        "LotkaVolterra": LotkaVolterra(),
        "ConservationLaw": ConservationLaw(),
        "SingleEvent": SingleEvent(),
        "MultiEvent": MultiEvent(),
    }


# ── Results collector ─────────────────────────────────────────────────────


class ResultsCollector:
    """Accumulates per-test metrics and writes them to a JSON file."""

    def __init__(self):
        self._data: dict = {}

    def add(self, model_id: str, operation: str, metrics: dict) -> None:
        self._data.setdefault(model_id, {})[operation] = metrics

    def dump(self, path: Path) -> None:
        import amici
        import equinox

        payload = {
            "_metadata": {
                "jax_version": jax.__version__,
                "diffrax_version": diffrax.__version__,
                "equinox_version": equinox.__version__,
                "optimistix_version": optimistix.__version__,
                "amici_version": amici.__version__,
                "python_version": platform.python_version(),
            },
            **self._data,
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(payload, f, indent=2)


@pytest.fixture(scope="session")
def results_collector(request):
    collector = ResultsCollector()
    yield collector
    results_dir = request.config.getoption("--results-dir")
    if results_dir:
        collector.dump(Path(results_dir) / "jax_regression_results.json")


# ── Timing helper ─────────────────────────────────────────────────────────
