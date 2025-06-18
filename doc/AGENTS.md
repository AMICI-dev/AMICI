# Project Agents.md Guide for OpenAI Codex

This Agents.md file provides comprehensive guidance for OpenAI Codex and other AI agents working with this codebase.

## Project Structure for OpenAI Codex Navigation

AMICI is a python package that uses SWIG to generate python bindings to C++ code. There are also a deprecated matlab interface to the C++ code, which will be removed at some point in the future.

- `/binder`: binder configuration
- `/cmake`: various cmake utility functions
- `/container`: docker configuration
- `/doc`: high level documentation, all API documentation is automatically generated using doxygen/sphinx
- `/include`: C++ header files
- `/matlab`: matlab interface
- `/models`: pre-generated c++ models for testing
- `/python`: python source code
  - `/benchmark`: helper scripts for benchmarking
  - `/sdist`: python package
    - `/amici`: python module
      - `/_codegen` helper functions for C++ code generation
      - `/debugging` helpfer functions for debugging simulation failures
      - `/include` symlink to C++ headers
      - `/jax` code for JAX backend, alternative to C++ backend
      - `/petab` interface for simulating models/problems in the PEtab format
      - `/testing` helpfer functions for tests
    - `/tests`: self contained python package specific tests
- `/scripts`: bash scripts for testing and installing the package
- `/src`: C++ source files
- `/swig`: definition of the SWIG interface
- `/tests`: C++ tests and python tests that require third party resources, base directory contains tests for the SBML testsuite
  - `/benchmark-models` regression tests for the PEtab benchmark collection
  - `/cpp` C++ tests
  - `/generateTestConfg` helper functions to generate python test configurations for pregenerated C++ models
  - `/performace` performance tests for the PEtab benchmark collection
  - `/petab_test_suite` regression tests for the PEtab testsuite

## Contribution Guidelines

Please see instructions in `doc/CONTRIBUTING.md`

# Agent Instructions

To ensure all tools and dependencies are available, activate the virtual environment before running any commands:

```bash
source ./venv/bin/activate
```

This project uses `pre-commit` for linting and `pytest` for tests. Run them on changed files whenever you make modifications.

When running the tests locally, change into the test directory first:

```bash
cd tests/benchmark-models
pytest test_petab_benchmark.py
pytest test_petab_benchmark_jax.py
```

To quickly verify the benchmark tests, you can limit execution to a small model:

```bash
pytest -k Boehm_JProteomeRes2014 test_petab_benchmark.py
```
