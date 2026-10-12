"""Miscellaneous AMICI Python interface tests"""

import os
import re
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import amici
import pytest
import sympy as sp
from amici.importers.utils import (
    maximum_bipartite_matching,
    smart_subs_dict,
)
from amici.sim.sundials import (
    ParameterScaling,
    parameter_scaling_from_int_vector,
)
from amici.testing import skip_on_valgrind


def test_parameter_scaling_from_int_vector():
    """Ensure we can generate a ParameterScaling vector from Python"""

    scale_vector = parameter_scaling_from_int_vector(
        [
            ParameterScaling.log10,
            ParameterScaling.ln,
            ParameterScaling.none,
        ]
    )

    assert scale_vector[0] == ParameterScaling.log10
    assert scale_vector[1] == ParameterScaling.ln
    assert scale_vector[2] == ParameterScaling.none


@skip_on_valgrind
def test_hill_function_dwdx():
    """Kinetic laws with Hill functions, may lead to NaNs in the Jacobian
    if involved states are zero if not properly arranged symbolically.
    Test that what we are applying the right sympy simplification."""

    w = sp.Matrix([[sp.sympify("Pow(x1, p1) / (Pow(x1, p1) + a)")]])
    dwdx = w.diff(sp.Symbol("x1"))

    # Verify that without simplification we fail
    with pytest.raises(ZeroDivisionError):
        with sp.evaluate(False):
            res = dwdx.subs({"x1": 0.0})
        _ = str(res)

    # Test that powsimp does the job
    dwdx = dwdx.applyfunc(lambda x: sp.powsimp(x, deep=True))
    with sp.evaluate(False):
        res = dwdx.subs({"x1": 0.0})
    _ = str(res)


@skip_on_valgrind
@pytest.mark.skipif(
    os.environ.get("AMICI_SKIP_CMAKE_TESTS", "") == "TRUE",
    reason="skipping cmake based test",
)
def test_cmake_compilation(sbml_example_presimulation_module):
    """Check that CMake build succeeds for one of the models generated during
    Python tests"""

    source_dir = os.path.dirname(sbml_example_presimulation_module.__path__[0])
    build_dir = f"{source_dir}/build"
    # path hint for amici base installation, in case CMake configuration has
    #  not been exported
    amici_dir = (Path(__file__).parents[2] / "build").absolute()
    cmd = (
        f"set -e; "
        f"cmake -S {source_dir} -B '{build_dir}' "
        f"-DCMAKE_BUILD_TYPE=Debug -DAmici_DIR={amici_dir}; "
        f"cmake --build '{build_dir}' --config Debug"
    )

    try:
        subprocess.run(
            cmd,
            shell=True,
            check=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as e:
        print(e.stdout.decode())
        print(e.stderr.decode())
        raise


@skip_on_valgrind
@pytest.mark.skipif(sys.platform == "win32", reason="requires POSIX shell")
@pytest.mark.filterwarnings("ignore:SWIG version mismatch:RuntimeWarning")
def test_model_rebuild_reruns_swig_on_swig_change(tmp_path, monkeypatch):
    """Check that rebuilding a model reruns SWIG if the SWIG version changed,
    even if the SWIG executable path and the interface file didn't."""
    from amici.exporters.sundials.compile import build_model_extension
    from amici.importers.antimony import antimony2amici

    real_swig = os.environ.get("SWIG") or shutil.which("swig")

    # SWIG wrapper that can report a fake version and logs code generation
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake_version_file = bin_dir / "fake_version"
    swig_runs_file = bin_dir / "swig_runs"
    swig = bin_dir / "swig"
    swig.write_text(
        "#!/bin/sh\n"
        f'if [ "$1" = "-version" ] && [ -f "{fake_version_file}" ]; then\n'
        f'  echo "SWIG Version $(cat "{fake_version_file}")"; exit 0\n'
        "fi\n"
        f'case " $* " in *" -python "*) echo run >> "{swig_runs_file}";; esac\n'
        f'exec "{real_swig}" "$@"\n'
    )
    swig.chmod(0o755)
    monkeypatch.setenv("SWIG", str(swig))

    module_name = "test_model_rebuild_reruns_swig"
    model_dir = tmp_path / module_name
    antimony2amici(
        "x = 1; x' = -x",
        model_name=module_name,
        output_dir=model_dir,
    )
    assert swig_runs_file.read_text().count("run") == 1

    # rebuild with a "different" SWIG version, without regenerating the model
    fake_version_file.write_text("9.9.9")
    build_model_extension(model_dir)
    assert swig_runs_file.read_text().count("run") == 2

    # the CMake cache must not keep the previous SWIG version
    (cmake_cache,) = model_dir.rglob("CMakeCache.txt")
    assert re.search(
        r"^SWIG_VERSION:STRING=9\.9\.9$",
        cmake_cache.read_text(),
        re.MULTILINE,
    )


@skip_on_valgrind
@pytest.mark.parametrize("flatten_first", [True, False])
def test_smart_subs_dict(flatten_first):
    expr_str = "c + d"
    subs_dict = {
        "c": "a + b",
        "d": "c + a",
    }
    expected_default_str = "3*a + 2*b"
    expected_reverse_str = "2*a + b + c"

    expr_sym = sp.sympify(expr_str)
    subs_sym = {sp.sympify(k): sp.sympify(v) for k, v in subs_dict.items()}
    expected_default = sp.sympify(expected_default_str)
    expected_reverse = sp.sympify(expected_reverse_str)

    result_default = smart_subs_dict(
        expr_sym, subs_sym, flatten_first=flatten_first
    )
    result_reverse = smart_subs_dict(
        expr_sym, subs_sym, reverse=False, flatten_first=flatten_first
    )

    assert sp.simplify(result_default - expected_default).is_zero
    assert sp.simplify(result_reverse - expected_reverse).is_zero


@skip_on_valgrind
def test_get_default_argument():
    # no default
    with pytest.raises(ValueError):
        amici._get_default_argument(lambda x: x, "x")

    # non-existant parameter
    with pytest.raises(KeyError):
        amici._get_default_argument(lambda x: x, "y")

    # okay
    assert amici._get_default_argument(lambda x=1: x, "x") == 1


def _write_fake_swig(path: Path) -> Path:
    """Create a fake ``swig`` executable that only reports its version."""
    path.mkdir(parents=True, exist_ok=True)
    swig = path / "swig"
    swig.write_text(
        "#!/bin/sh\n"
        'if [ "$1" = "-version" ]; then echo "SWIG Version 1.2.3"; exit 0; fi\n'
        "exit 1\n"
    )
    swig.chmod(0o755)
    return swig


@skip_on_valgrind
@pytest.mark.skipif(sys.platform == "win32", reason="requires POSIX shell")
def test_model_build_uses_env_swig(tmp_path, monkeypatch):
    """Check that model compilation uses the SWIG of the current Python
    environment, not the first one on PATH."""
    from amici.exporters.sundials import compile as compile_module

    if (env_swig := compile_module._get_env_swig()) is None:
        pytest.skip("no swig in the current Python environment")

    fake_swig = _write_fake_swig(tmp_path / "bin")
    monkeypatch.delenv("SWIG", raising=False)
    monkeypatch.setenv(
        "PATH", f"{fake_swig.parent}{os.pathsep}{os.environ['PATH']}"
    )

    # record the environment of the build command instead of building
    build_envs = []
    orig_run = subprocess.run

    def run(args, **kwargs):
        if args[0] == sys.executable:
            build_envs.append(kwargs["env"])
            return subprocess.CompletedProcess(args, 0, stdout=b"")
        return orig_run(args, **kwargs)

    monkeypatch.setattr(compile_module.subprocess, "run", run)

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        compile_module.build_model_extension(tmp_path)

    assert len(build_envs) == 1
    assert build_envs[0]["SWIG"] == env_swig != str(fake_swig)


@skip_on_valgrind
@pytest.mark.skipif(sys.platform == "win32", reason="requires POSIX shell")
def test_model_build_swig_version_mismatch(tmp_path, monkeypatch):
    """Check that a mismatching SWIG version is reported before
    compilation."""
    from amici.exporters.sundials.compile import build_model_extension

    fake_swig = _write_fake_swig(tmp_path / "bin")
    monkeypatch.setenv("SWIG", str(fake_swig))

    # There is nothing to build, so the build fails after the check
    with (
        pytest.warns(RuntimeWarning, match=r"SWIG version mismatch.*1\.2\.3"),
        pytest.raises(subprocess.CalledProcessError),
    ):
        build_model_extension(tmp_path)


def _assert_valid_matching(edges, matching):
    assert all(right in edges[left] for left, right in matching.items())
    assert len(set(matching.values())) == len(matching)


@pytest.mark.parametrize(
    "edges, expected_size",
    [
        ({}, 0),
        ({"a": []}, 0),
        ({"a": [1], "b": []}, 1),
        # `b` needs `a` to switch from 1 to 2
        ({"a": [1, 2], "b": [1]}, 2),
        # `c` needs `b` to switch to 3, `b` needs `a` to switch to 2
        ({"a": [1, 2], "b": [1, 3], "c": [1]}, 3),
        # two left vertices compete for a single right vertex
        ({"a": [1], "b": [1]}, 1),
        ({"a": [1, 2], "b": [1, 2], "c": [1, 2]}, 2),
    ],
)
def test_maximum_bipartite_matching(edges, expected_size):
    matching = maximum_bipartite_matching(edges)
    _assert_valid_matching(edges, matching)
    assert len(matching) == expected_size


def test_maximum_bipartite_matching_long_augmenting_path():
    """Test that long augmenting paths don't hit the recursion limit."""
    n = 5000
    # left i is first matched to right i; the last left vertex then needs
    # all others to shift by one
    edges = {i: [i, i + 1] for i in range(n - 1)}
    edges[n - 1] = [0]
    matching = maximum_bipartite_matching(edges)
    _assert_valid_matching(edges, matching)
    assert len(matching) == n
