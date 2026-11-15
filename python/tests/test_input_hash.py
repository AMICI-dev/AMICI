"""Tests for reusing imported models based on hashes of the import inputs."""

import enum
import inspect
import json
import logging
import os
import subprocess
import sys
from types import SimpleNamespace

import amici
import pytest
import sympy as sp
from amici._input_hash import (
    NO_MODEL_REUSE_ENV_VAR,
    NON_HASHED_ARGS,
    UNHASHABLE,
    _normalize,
    _UnhashableError,
    compute_input_hash,
    get_changed_inputs,
)
from amici._model_dir import MARKER_FILENAME
from amici.importers.antimony import antimony2sbml
from amici.importers.sbml import SbmlImporter
from amici.importers.utils import MeasurementChannel, _default_simplify

ANTIMONY_MODEL = "x = 1; x' = -k * x; k = 0.1; c = 2"


class _Stop(Exception):
    pass


def _get_hashed_options(monkeypatch, import_function, *args, **kwargs):
    """Get the options an import function passes to the input hash."""
    import amici._input_hash

    captured = {}

    def compute(model, options, **_):
        captured.update(options)
        raise _Stop

    for module in (
        sys.modules["amici.importers.sbml"],
        sys.modules.get("amici.importers.pysb"),
    ):
        if module is not None:
            monkeypatch.setattr(module, "compute_input_hash", compute)
    monkeypatch.setattr(amici._input_hash, "compute_input_hash", compute)
    with pytest.raises(_Stop):
        import_function(*args, **kwargs)
    return captured


@pytest.mark.parametrize("entry_point", ["sbml2amici", "sbml2jax"])
def test_all_sbml_args_hashed(tmp_path, monkeypatch, entry_point):
    """All arguments that may affect the generated code are hashed."""
    if entry_point == "sbml2jax":
        pytest.importorskip("jax")
    importer = SbmlImporter(antimony2sbml(ANTIMONY_MODEL), from_file=False)
    function = getattr(importer, entry_point)
    options = _get_hashed_options(
        monkeypatch, function, "test_model", tmp_path
    )
    params = set(inspect.signature(function).parameters)
    assert params - NON_HASHED_ARGS <= options.keys()
    assert not NON_HASHED_ARGS & options.keys()
    assert options["entry_point"] == entry_point


@pytest.mark.parametrize("entry_point", ["pysb2amici", "pysb2jax"])
def test_all_pysb_args_hashed(tmp_path, monkeypatch, entry_point):
    """All arguments that may affect the generated code are hashed."""
    pysb = pytest.importorskip("pysb")
    if entry_point == "pysb2jax":
        pytest.importorskip("jax")
    import amici.importers.pysb

    function = getattr(amici.importers.pysb, entry_point)
    model = pysb.Model("test_model", _export=False)
    options = _get_hashed_options(monkeypatch, function, model, tmp_path)
    params = set(inspect.signature(function).parameters)
    assert params - NON_HASHED_ARGS - {"model"} <= options.keys()
    assert not (NON_HASHED_ARGS | {"model"}) & options.keys()
    assert options["model_name"] == "test_model"
    assert options["entry_point"] == entry_point


def test_pysb_export_failure(monkeypatch):
    """PySB models that can't be exported can be imported, not reused."""
    pysb = pytest.importorskip("pysb")
    pysb_export = pytest.importorskip("pysb.export")
    from amici.importers.pysb import _compute_input_hash

    def raise_error(*args, **kwargs):
        raise ValueError("not exportable")

    monkeypatch.delenv(NO_MODEL_REUSE_ENV_VAR, raising=False)
    monkeypatch.setattr(pysb_export, "export", raise_error)
    input_hash, inputs = _compute_input_hash(
        {"model": pysb.Model("test_model", _export=False)},
        "pysb2amici",
        build=False,
    )
    assert input_hash is None
    assert inputs["model"] == UNHASHABLE


def test_normalize(tmp_path):
    x, y = sp.symbols("x y")
    # equal values of different types are distinguished
    assert (
        len({json.dumps(_normalize(value)) for value in (1, 1.0, True)}) == 3
    )
    assert _normalize([1]) == _normalize((1,))
    assert _normalize({"a", "b"}) == _normalize({"b", "a"})
    assert _normalize({x: 1, y: 2}) == _normalize({y: 2, x: 1})
    assert _normalize(x + y) != _normalize(x * y)
    assert _normalize(sp.Symbol("x", real=True)) != _normalize(x)
    assert _normalize(float("nan")) == _normalize(float("nan"))
    json.dumps(_normalize(float("inf")), allow_nan=False)
    assert _normalize(_default_simplify) == {"callable": "default_simplify"}
    tuple_enum = enum.Enum("TupleEnum", {"A": (1, 2)})
    assert json.loads(json.dumps(_normalize(tuple_enum.A))) == _normalize(
        tuple_enum.A
    )
    with pytest.raises(_UnhashableError):
        _normalize(enum.Enum("ObjectEnum", {"A": object()}).A)
    assert _normalize(
        MeasurementChannel("obs", formula="x", sigma=1.0)
    ) != _normalize(MeasurementChannel("obs", formula="x", sigma=2.0))

    # files are represented by their contents
    file = tmp_path / "file.txt"
    file.write_text("a")
    normalized = _normalize(file)
    assert _normalize(file) == normalized
    file.write_text("b")
    assert _normalize(file) != normalized

    with pytest.raises(_UnhashableError):
        _normalize(tmp_path)
    with pytest.raises(_UnhashableError):
        _normalize(tmp_path / "missing")
    with pytest.raises(_UnhashableError):
        _normalize(lambda x: x)
    with pytest.raises(_UnhashableError):
        _normalize(
            MeasurementChannel(
                "obs", formula="x", noise_distribution=lambda s: s
            )
        )
    with pytest.raises(_UnhashableError):
        _normalize(object())


def test_compute_input_hash(monkeypatch):
    monkeypatch.delenv(NO_MODEL_REUSE_ENV_VAR, raising=False)
    options = {"a": 1, "b": {"x", "y"}}
    input_hash, inputs = compute_input_hash(
        "model", options, packages=(), build=False
    )
    assert input_hash is not None
    assert "model" not in json.dumps(inputs["options"])
    assert "build" not in inputs["environment"]
    assert (
        compute_input_hash("model", options, packages=(), build=False)[0]
        == input_hash
    )

    # stored inputs compare equal after reading them back
    _, inputs = compute_input_hash(
        "model",
        {"d": {"a": 1}, "mc": [MeasurementChannel("obs", formula="x")]},
        packages=(),
        build=True,
    )
    assert json.loads(json.dumps(inputs)) == inputs
    assert (
        compute_input_hash("model2", options, packages=(), build=False)[0]
        != input_hash
    )
    assert (
        compute_input_hash("model", options, packages=(), build=True)[0]
        != input_hash
    )
    monkeypatch.setenv("AMICI_CXXFLAGS", "-O0")
    assert (
        compute_input_hash("model", options, packages=(), build=True)[0]
        != compute_input_hash("model", options, packages=(), build=False)[0]
    )

    # the order of fixed parameters doesn't matter, but the order of other
    #  lists does
    def get_hash(**options):
        return compute_input_hash("model", options, packages=(), build=False)[
            0
        ]

    assert get_hash(fixed_parameters=["a", "b"]) == get_hash(
        fixed_parameters=["b", "a"]
    )
    assert get_hash(fixed_parameters=None) != get_hash(fixed_parameters=[])
    assert get_hash(other=["a", "b"]) != get_hash(other=["b", "a"])

    # unhashable inputs
    input_hash, inputs = compute_input_hash(
        "model", {"f": lambda: 1}, packages=(), build=False
    )
    assert input_hash is None
    assert inputs["options"]["f"] == UNHASHABLE
    input_hash, inputs = compute_input_hash(
        None, options, packages=(), build=False
    )
    assert input_hash is None
    assert inputs["model"] == UNHASHABLE

    # disabled reuse
    for value, disabled in (
        ("", False),
        ("0", False),
        ("false", False),
        ("1", True),
        ("True", True),
    ):
        monkeypatch.setenv(NO_MODEL_REUSE_ENV_VAR, value)
        input_hash, _ = compute_input_hash(
            "model", options, packages=(), build=False
        )
        assert (input_hash is None) == disabled


def test_input_hash_is_deterministic(tmp_path):
    """The hash doesn't depend on the process, e.g., on set ordering."""
    script = tmp_path / "hash.py"
    script.write_text(
        "from amici._input_hash import compute_input_hash\n"
        "import sympy as sp\n"
        # dummy symbols, e.g., PySB components, have a process-dependent
        #  index
        "print(compute_input_hash('model', {"
        "'ids': {f'p{i}' for i in range(20)},"
        "'expr': sp.sympify('a + b * c') * sp.Dummy('d', real=True),"
        "}, packages=('python-libsbml',), build=True)[0])\n"
    )
    env = os.environ.copy()
    env.pop(NO_MODEL_REUSE_ENV_VAR, None)
    hashes = {
        subprocess.run(
            [sys.executable, str(script)],
            env=env | {"PYTHONHASHSEED": seed},
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        for seed in ("1", "2")
    }
    assert len(hashes) == 1


def test_get_changed_inputs():
    old = {"model": "a", "options": {"x": 1, "y": 2}, "environment": {}}
    assert get_changed_inputs(old, old) == []
    new = {"model": "b", "options": {"x": 1, "z": 2}, "environment": {}}
    assert get_changed_inputs(old, new) == [
        "model",
        "options.y",
        "options.z",
    ]
    assert get_changed_inputs(None, new) == ["<unknown>"]


def test_sbml_reuse(tmp_path, monkeypatch, caplog):
    monkeypatch.delenv(NO_MODEL_REUSE_ENV_VAR, raising=False)
    model_dir = tmp_path / "model"
    sbml = antimony2sbml(ANTIMONY_MODEL)
    num_imports = 0
    build_ode_model = SbmlImporter._build_ode_model

    def counting_build_ode_model(*args, **kwargs):
        nonlocal num_imports
        num_imports += 1
        return build_ode_model(*args, **kwargs)

    monkeypatch.setattr(
        SbmlImporter, "_build_ode_model", counting_build_ode_model
    )

    def import_model(
        sbml=sbml, skip_if_unchanged=True, compile=False, **kwargs
    ):
        SbmlImporter(sbml, from_file=False).sbml2amici(
            "reuse_test",
            model_dir,
            compile=compile,
            skip_if_unchanged=skip_if_unchanged,
            **kwargs,
        )
        return num_imports

    assert import_model() == 1
    marker = json.loads((model_dir / MARKER_FILENAME).read_text())
    assert marker["input_hash"]
    assert marker["inputs"]["options"]["entry_point"] == "sbml2amici"

    assert import_model() == 1
    # reuse is opt-in
    assert import_model(skip_if_unchanged=False) == 2
    assert import_model() == 2

    # reported according to the verbosity of the import
    with caplog.at_level(logging.INFO):
        assert import_model(fixed_parameters=["k"]) == 3
        assert "Cannot reuse" not in caplog.text
        assert import_model(fixed_parameters=["c"], verbose=logging.INFO) == 4
    assert "Changed import inputs: options.fixed_parameters." in caplog.text
    assert import_model(fixed_parameters=["c"]) == 4

    caplog.clear()
    with caplog.at_level(logging.INFO):
        assert (
            import_model(
                sbml=antimony2sbml(ANTIMONY_MODEL + "; d = 3"),
                verbose=logging.INFO,
            )
            == 5
        )
    assert "Changed import inputs: model, options.fixed_parameters." in (
        caplog.text
    )

    # unhashable inputs are never reused
    assert import_model(simplify=lambda x: x) == 6
    caplog.clear()
    with caplog.at_level(logging.INFO):
        assert import_model(simplify=lambda x: x, verbose=logging.INFO) == 7
    assert "Import options can't be hashed: simplify." in caplog.text

    assert import_model() == 8
    monkeypatch.setenv(NO_MODEL_REUSE_ENV_VAR, "1")
    assert import_model() == 9

    # incomplete models are not reused
    monkeypatch.delenv(NO_MODEL_REUSE_ENV_VAR)
    assert import_model() == 10
    assert import_model() == 10
    marker = json.loads((model_dir / MARKER_FILENAME).read_text())
    (model_dir / MARKER_FILENAME).write_text(
        json.dumps(marker | {"status": "importing"})
    )
    assert import_model() == 11

    # a generated model is not reused if a compiled model is required,
    #  but a compiled model is reused if only code generation is required
    from amici.exporters.sundials import compile

    monkeypatch.setattr(
        compile.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(stdout=b""),
    )
    monkeypatch.setattr(compile, "_check_swig_version", lambda swig: None)
    monkeypatch.setattr(
        amici,
        "import_model_module",
        lambda **_: SimpleNamespace(get_model=lambda: None),
    )
    assert import_model() == 11
    assert import_model(compile=True) == 12
    assert import_model(compile=True) == 12
    assert import_model() == 12

    # a model extension that can't be loaded is rebuilt
    def raise_import_error(**_):
        raise ImportError("broken")

    monkeypatch.setattr(amici, "import_model_module", raise_import_error)
    with pytest.raises(ImportError, match="broken"):
        import_model(compile=True)
    assert num_imports == 13


def test_petab_v1_reuse(tmp_path, monkeypatch):
    """`compile_` controls reuse in PEtab v1 imports."""
    pytest.importorskip("petab.v1")
    from amici.importers.petab.v1 import _petab_import
    from test_model_dir import _petab_v1_problem

    calls = []

    def import_model_sbml(*, skip_if_unchanged, **_):
        calls.append(skip_if_unchanged)

    monkeypatch.setattr(_petab_import, "import_model_sbml", import_model_sbml)
    monkeypatch.setattr(amici, "import_model_module", lambda *args: None)
    problem = _petab_v1_problem()
    for compile_ in (None, True, False):
        # fails on loading the mocked model module
        with pytest.raises(AttributeError):
            _petab_import.import_petab_problem(
                problem, tmp_path, model_name="m", compile_=compile_
            )
    assert calls == [True, False]
