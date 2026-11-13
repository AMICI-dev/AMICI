"""Tests for model directory ownership and completion markers."""

import json
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import amici
import pytest
from amici._model_dir import (
    MARKER_FILENAME,
    MARKER_FORMAT_VERSION,
    ModelDirStatus,
    get_foreign_entries,
    get_model_dir_status,
    is_model_dir_complete,
    prepare_model_dir,
    set_model_dir_status,
)
from amici.importers.antimony import antimony2amici, antimony2sbml

LEGACY_SETUP_PY = '''"""AMICI model package setup"""
setup(
    name="old_model",
    description="AMICI-generated module for model old_model",
)
'''


def _make_files(root: Path, *names: str):
    for name in names:
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)


def _make_model_dir(model_dir: Path, status=ModelDirStatus.COMPILED):
    """Create a fake model directory with a marker."""
    prepare_model_dir(model_dir, "model")
    _make_files(
        model_dir,
        "model.cpp",
        "setup.py",
        "swig/model.i",
        "model/__init__.py",
    )
    set_model_dir_status(model_dir, status)


def test_prepare_new_dir(tmp_path):
    model_dir = tmp_path / "a" / "model"
    prepare_model_dir(model_dir, "model")
    assert [p.name for p in model_dir.iterdir()] == [MARKER_FILENAME]
    assert get_model_dir_status(model_dir) == ModelDirStatus.IMPORTING
    assert not is_model_dir_complete(model_dir, compiled=False)

    # empty directory
    model_dir = tmp_path / "empty"
    model_dir.mkdir()
    prepare_model_dir(model_dir, "model", jax=True)
    assert get_model_dir_status(model_dir) == ModelDirStatus.IMPORTING


def test_prepare_marked_dir(tmp_path):
    """A previous import is removed completely, including subdirectories."""
    model_dir = tmp_path / "model"
    _make_model_dir(model_dir)
    marker = json.loads((model_dir / MARKER_FILENAME).read_text())
    assert marker["entries"] == ["model", "model.cpp", "setup.py", "swig"]
    # unrecorded build directories of a failed compilation
    _make_files(
        model_dir,
        "build_model_ext/CMakeCache.txt",
        "build/x.o",
        "__pycache__/x.pyc",
    )
    # leftover of an interrupted marker update
    _make_files(model_dir, f"{MARKER_FILENAME}.123.tmp")
    assert get_foreign_entries(model_dir) == []
    # entries created later are not recorded by status updates
    _make_files(model_dir, "build_results/x.csv", "notes.txt")
    set_model_dir_status(model_dir, ModelDirStatus.GENERATED)
    marker = json.loads((model_dir / MARKER_FILENAME).read_text())
    assert marker["entries"] == ["model", "model.cpp", "setup.py", "swig"]
    assert get_foreign_entries(model_dir) == ["build_results", "notes.txt"]
    (model_dir / "notes.txt").unlink()
    (model_dir / "build_results" / "x.csv").unlink()
    (model_dir / "build_results").rmdir()

    prepare_model_dir(model_dir, "other_name")
    assert [p.name for p in model_dir.iterdir()] == [MARKER_FILENAME]
    marker = json.loads((model_dir / MARKER_FILENAME).read_text())
    assert marker["model_name"] == "other_name"
    assert marker["status"] == ModelDirStatus.IMPORTING
    assert marker["entries"] == []
    assert marker["amici_version"] == amici.__version__
    assert marker["format_version"] == MARKER_FORMAT_VERSION


@pytest.mark.parametrize(
    "status", [ModelDirStatus.COMPILED, ModelDirStatus.IMPORTING]
)
def test_prepare_marked_dir_with_foreign_files(tmp_path, status):
    """Foreign files prevent removal, also after an interrupted import."""
    model_dir = tmp_path / "model"
    _make_model_dir(model_dir, status)
    _make_files(model_dir, "notes.txt", "data/data.csv")

    with pytest.raises(FileExistsError, match="data, notes.txt"):
        prepare_model_dir(model_dir, "model")
    assert (model_dir / "model.cpp").exists()
    assert (model_dir / "notes.txt").exists()
    assert (model_dir / "data" / "data.csv").exists()
    assert get_model_dir_status(model_dir) == status


def test_prepare_foreign_dir(tmp_path):
    """Directories not created by AMICI are never removed."""
    model_dir = tmp_path / "model"
    _make_files(model_dir, "data.csv")
    with pytest.raises(FileExistsError, match="data.csv"):
        prepare_model_dir(model_dir, "model")

    # a setup.py that just contains the AMICI docstring is not enough
    model_dir = tmp_path / "model2"
    model_dir.mkdir()
    (model_dir / "setup.py").write_text('"""AMICI model package setup"""')
    with pytest.raises(FileExistsError, match="setup.py"):
        prepare_model_dir(model_dir, "model")
    assert (model_dir / "setup.py").exists()

    # invalid marker
    model_dir = tmp_path / "model3"
    _make_files(model_dir, MARKER_FILENAME, "model.cpp")
    with pytest.raises(FileExistsError, match="model.cpp"):
        prepare_model_dir(model_dir, "model")

    # unknown status or backend
    for key, value in [("status", "done"), ("backend", "cvodes")]:
        model_dir = tmp_path / f"model_{key}"
        _make_model_dir(model_dir)
        marker = json.loads((model_dir / MARKER_FILENAME).read_text())
        marker[key] = value
        (model_dir / MARKER_FILENAME).write_text(json.dumps(marker))
        assert get_model_dir_status(model_dir) is None
        with pytest.raises(FileExistsError, match="model, model.cpp"):
            prepare_model_dir(model_dir, "model")

    # marker without format version
    model_dir = tmp_path / "model4"
    _make_model_dir(model_dir)
    marker = json.loads((model_dir / MARKER_FILENAME).read_text())
    del marker["format_version"]
    (model_dir / MARKER_FILENAME).write_text(json.dumps(marker))
    with pytest.raises(FileExistsError, match="model, model.cpp"):
        prepare_model_dir(model_dir, "model")


@pytest.mark.parametrize("platform", ["linux", "win32"])
def test_prepare_permission_error(tmp_path, monkeypatch, platform):
    """Previous imports that cannot be removed, e.g., loaded on Windows."""
    model_dir = tmp_path / "model"
    _make_model_dir(model_dir)

    def rmtree(path, *args, **kwargs):
        raise PermissionError(f"Access denied: '{path}'")

    monkeypatch.setattr(shutil, "rmtree", rmtree)
    monkeypatch.setattr(sys, "platform", platform)
    with pytest.raises(PermissionError, match="Access denied") as exc_info:
        prepare_model_dir(model_dir, "model")
    assert ("still loaded" in str(exc_info.value)) == (platform == "win32")
    # still recognized as created by AMICI
    monkeypatch.undo()
    prepare_model_dir(model_dir, "model")
    assert [p.name for p in model_dir.iterdir()] == [MARKER_FILENAME]


def test_prepare_newer_format(tmp_path):
    """Directories with a marker of a newer format are not removed."""
    model_dir = tmp_path / "model"
    _make_model_dir(model_dir)
    marker = json.loads((model_dir / MARKER_FILENAME).read_text())
    marker["format_version"] = MARKER_FORMAT_VERSION + 1
    (model_dir / MARKER_FILENAME).write_text(json.dumps(marker))

    assert get_foreign_entries(model_dir) == marker["entries"]
    with pytest.raises(FileExistsError, match="newer AMICI version"):
        prepare_model_dir(model_dir, "model")
    assert (model_dir / "model.cpp").exists()
    assert get_model_dir_status(model_dir) is None
    assert not is_model_dir_complete(model_dir)
    # not modified
    set_model_dir_status(model_dir, ModelDirStatus.GENERATED)
    assert json.loads((model_dir / MARKER_FILENAME).read_text()) == marker


@pytest.mark.parametrize("quote", ['"', "'"])
def test_prepare_legacy_dir(tmp_path, quote):
    """Model directories created before the marker file existed."""
    model_dir = tmp_path / "model"
    _make_files(
        model_dir,
        "old_model.cpp",
        "old_model.h",
        "CMakeLists.txt",
        ".gitignore",
        "swig/old_model.i",
        "old_model/_old_model.so",
        "build_model_ext/CMakeCache.txt",
    )
    (model_dir / "setup.py").write_text(
        LEGACY_SETUP_PY.replace(
            'description="', f"description={quote}"
        ).replace('old_model",\n)', f"old_model{quote},\n)")
    )
    _make_files(
        model_dir, "other_model/__init__.py", "build_results/x", "notes.txt"
    )
    assert get_foreign_entries(model_dir) == [
        "build_results",
        "notes.txt",
        "other_model",
    ]
    with pytest.raises(FileExistsError, match="notes.txt, other_model"):
        prepare_model_dir(model_dir, "model")
    assert (model_dir / "old_model" / "_old_model.so").exists()

    shutil.rmtree(model_dir / "build_results")
    shutil.rmtree(model_dir / "other_model")
    (model_dir / "notes.txt").unlink()
    prepare_model_dir(model_dir, "model")
    assert [p.name for p in model_dir.iterdir()] == [MARKER_FILENAME]


@pytest.mark.parametrize("module", ["amici.sim.jax.model", "amici.jax.model"])
def test_prepare_legacy_jax_dir(tmp_path, module):
    model_dir = tmp_path / "model"
    model_dir.mkdir()
    # importing JAXModel alone is not sufficient
    init_py = f"import jax\nfrom {module} import JAXModel, safe_div\n"
    (model_dir / "__init__.py").write_text(init_py)
    with pytest.raises(FileExistsError, match="__init__.py"):
        prepare_model_dir(model_dir, "model", jax=True)

    (model_dir / "__init__.py").write_text(
        f"{init_py}\n\nclass JAXModel_old(JAXModel):\n    pass\n\n\n"
        "Model = JAXModel_old\n"
    )
    _make_files(model_dir, "net1.py", "__pycache__/x.pyc", "data.csv")
    with pytest.raises(FileExistsError, match="data.csv"):
        prepare_model_dir(model_dir, "model", jax=True)

    (model_dir / "data.csv").unlink()
    prepare_model_dir(model_dir, "model", jax=True)
    assert [p.name for p in model_dir.iterdir()] == [MARKER_FILENAME]


@pytest.mark.skipif(sys.platform == "win32", reason="symlinks")
def test_prepare_does_not_follow_symlinks(tmp_path):
    outside = tmp_path / "outside"
    _make_files(outside, "keep.txt")
    model_dir = tmp_path / "model"
    prepare_model_dir(model_dir, "model")
    (model_dir / "model").symlink_to(outside, target_is_directory=True)
    set_model_dir_status(model_dir, ModelDirStatus.GENERATED)

    prepare_model_dir(model_dir, "model")
    assert not (model_dir / "model").exists()
    assert (outside / "keep.txt").exists()

    # not accepted as a generated directory unless recorded
    (model_dir / "swig").symlink_to(outside, target_is_directory=True)
    with pytest.raises(FileExistsError, match="swig"):
        prepare_model_dir(model_dir, "model")


def test_is_model_dir_complete(tmp_path, monkeypatch):
    assert not is_model_dir_complete(tmp_path)
    prepare_model_dir(tmp_path, "model")
    assert not is_model_dir_complete(tmp_path, compiled=False)

    set_model_dir_status(tmp_path, ModelDirStatus.GENERATED)
    assert not is_model_dir_complete(tmp_path)
    assert is_model_dir_complete(tmp_path, compiled=False)
    assert not is_model_dir_complete(tmp_path, jax=True)

    set_model_dir_status(tmp_path, ModelDirStatus.COMPILED)
    assert is_model_dir_complete(tmp_path)
    monkeypatch.setattr(amici, "__version__", "0.0.0")
    assert not is_model_dir_complete(tmp_path)
    monkeypatch.undo()

    prepare_model_dir(tmp_path, "model", jax=True)
    set_model_dir_status(tmp_path, ModelDirStatus.GENERATED)
    assert is_model_dir_complete(tmp_path, jax=True)
    assert not is_model_dir_complete(tmp_path, compiled=False)

    (tmp_path / MARKER_FILENAME).write_text("{")
    assert not is_model_dir_complete(tmp_path, jax=True)


def test_sbml_import(tmp_path, monkeypatch):
    """Marker handling by code generation and compilation."""
    from amici.exporters.sundials import compile

    build_fails = False

    def run(args, cwd, **kwargs):
        (Path(cwd) / "build_model_ext").mkdir(exist_ok=True)
        if build_fails:
            raise subprocess.CalledProcessError(1, args, output=b"failed")
        return SimpleNamespace(stdout=b"")

    monkeypatch.setattr(compile.subprocess, "run", run)
    monkeypatch.setattr(compile, "_check_swig_version", lambda swig: None)
    monkeypatch.setattr(
        amici,
        "import_model_module",
        lambda **_: SimpleNamespace(get_model=lambda: None),
    )
    model_dir = tmp_path / "model"

    def import_model(compile=True):
        antimony2amici(
            "x' = 1",
            model_name="marker_test",
            output_dir=model_dir,
            compile=compile,
        )

    import_model(compile=False)
    assert get_model_dir_status(model_dir) == ModelDirStatus.GENERATED
    marker = json.loads((model_dir / MARKER_FILENAME).read_text())
    assert {"setup.py", "swig", "marker_test"} <= set(marker["entries"])
    # building later
    compile.build_model_extension(model_dir)
    assert get_model_dir_status(model_dir) == ModelDirStatus.COMPILED
    # removed along with the previous import
    _make_files(model_dir, "swig/leftover.i")

    import_model()
    assert not (model_dir / "swig" / "leftover.i").exists()
    assert get_model_dir_status(model_dir) == ModelDirStatus.COMPILED

    # failed build of a previously compiled model
    build_fails = True
    with pytest.raises(subprocess.CalledProcessError):
        compile.build_model_extension(model_dir)
    assert get_model_dir_status(model_dir) == ModelDirStatus.GENERATED
    with pytest.raises(subprocess.CalledProcessError):
        import_model()
    assert get_model_dir_status(model_dir) == ModelDirStatus.GENERATED
    build_fails = False
    # the failed build directory is removed on re-import
    import_model(compile=False)
    assert not (model_dir / "build_model_ext").exists()

    _make_files(model_dir, "notes.txt")
    with pytest.raises(FileExistsError, match="notes.txt"):
        import_model(compile=False)
    assert (model_dir / "setup.py").exists()


def test_jax_import(tmp_path):
    pytest.importorskip("jax")
    from amici.importers.sbml import SbmlImporter

    model_dir = tmp_path / "model"
    importer = SbmlImporter(antimony2sbml("x' = 1"), from_file=False)
    importer.sbml2jax("marker_test", output_dir=model_dir)
    assert is_model_dir_complete(model_dir, jax=True)

    # importing the model creates __pycache__
    amici._module_from_path("marker_test", model_dir / "__init__.py")
    assert (model_dir / "__pycache__").is_dir()
    importer.sbml2jax("marker_test", output_dir=model_dir)
    assert is_model_dir_complete(model_dir, jax=True)

    _make_files(model_dir, "notes.txt")
    with pytest.raises(FileExistsError, match="notes.txt"):
        importer.sbml2jax("marker_test", output_dir=model_dir)
    assert is_model_dir_complete(model_dir, jax=True)


def _petab_v2_problem():
    from petab.v2 import Problem
    from petab.v2.models.sbml_model import SbmlModel

    problem = Problem()
    problem.model = SbmlModel.from_antimony("x = 1; x' = k;")
    problem.add_parameter("k", nominal_value=1.0, estimate=True, lb=0, ub=10)
    problem.add_observable("obs1", "x", noise_formula="1")
    problem.add_measurement("obs1", time=1, measurement=1)
    return problem


def test_petab_v2_import_module(tmp_path, monkeypatch):
    """Re-import into incomplete model directories."""
    pytest.importorskip("petab.v2")
    from amici.importers.petab import PetabImporter

    model_dir = tmp_path / "model"
    importer = PetabImporter(_petab_v2_problem(), output_dir=model_dir)
    imports = []

    def do_import_sbml():
        imports.append(1)
        prepare_model_dir(model_dir, "model")
        set_model_dir_status(model_dir, ModelDirStatus.COMPILED)

    monkeypatch.setattr(importer, "_do_import_sbml", do_import_sbml)
    monkeypatch.setattr(amici, "import_model_module", lambda *args: None)

    importer.import_module()
    assert len(imports) == 1
    # complete model
    importer.import_module()
    assert len(imports) == 1
    importer.import_module(force_import=True)
    assert len(imports) == 2
    # left behind by a failed import
    set_model_dir_status(model_dir, ModelDirStatus.GENERATED)
    importer.import_module()
    assert len(imports) == 3


@pytest.mark.parametrize("force_import", [False, True])
def test_petab_v2_foreign_dir(tmp_path, force_import):
    """Foreign directories are not removed, even with `force_import`."""
    pytest.importorskip("petab.v2")
    from amici.importers.petab import PetabImporter

    model_dir = tmp_path / "model"
    _make_files(model_dir, "data.csv")
    importer = PetabImporter(_petab_v2_problem(), output_dir=model_dir)
    with pytest.raises(FileExistsError, match="data.csv"):
        importer.import_module(force_import=force_import)
    assert (model_dir / "data.csv").exists()


@pytest.mark.parametrize("compile_", [None, True])
def test_petab_v1_foreign_dir(tmp_path, compile_):
    """Foreign directories are not removed, even with `compile_=True`."""
    petab = pytest.importorskip("petab.v1")
    import pandas as pd
    from amici.importers.petab.v1 import import_petab_problem
    from petab.v1.models.sbml_model import SbmlModel

    petab_problem = petab.Problem(
        model=SbmlModel.from_antimony("x = 1; x' = k; k = 1"),
        condition_df=petab.get_condition_df(
            pd.DataFrame({petab.CONDITION_ID: ["c0"]})
        ),
        parameter_df=petab.get_parameter_df(
            pd.DataFrame({petab.PARAMETER_ID: [], petab.ESTIMATE: []})
        ),
        observable_df=petab.get_observable_df(
            pd.DataFrame(
                {
                    petab.OBSERVABLE_ID: ["obs1"],
                    petab.OBSERVABLE_FORMULA: ["x"],
                    petab.NOISE_FORMULA: [1],
                }
            )
        ),
    )
    model_dir = tmp_path / "model"
    _make_files(model_dir, "data.csv")
    with pytest.raises(FileExistsError, match="data.csv"):
        import_petab_problem(petab_problem, model_dir, compile_=compile_)
    assert (model_dir / "data.csv").exists()


def test_marker_written_by_import(sbml_example_presimulation_module):
    model_dir = Path(sbml_example_presimulation_module.__file__).parents[1]
    assert is_model_dir_complete(model_dir)
    assert get_foreign_entries(model_dir) == []
