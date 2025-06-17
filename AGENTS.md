# Agent Instructions

To ensure all tools and dependencies are available, activate the virtual environment before running any commands:

```bash
source ./venv/bin/activate
```

This project uses `pre-commit` for linting and `pytest` for tests. Run them on changed files whenever you make modifications.

When running the benchmark tests locally, change into the test directory first:

```bash
cd tests/benchmark-models
pytest test_petab_benchmark.py
pytest test_petab_benchmark_jax.py
```

To quickly verify the benchmark tests, you can limit execution to a small model:

```bash
pytest -k Boehm_JProteomeRes2014 test_petab_benchmark.py
```
