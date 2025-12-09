AMICI developer’s guide
=======================

This document contains information for AMICI developers, not too
relevant to regular users.

When starting to work on some issue
-----------------------------------

When starting to work on some GitHub issue, please assign yourself to
let other developers know that you are working on it to avoid duplicate
work. If the respective issue is not completely clear, it is generally a
good idea to ask for clarification before starting to work on it.

If you want to work on something new, please create a
`GitHub issue <https://github.com/AMICI-dev/AMICI/issues>`__ first.

Code contributions
------------------

When making code contributions, please follow our style guide and the
process described below:

-  Check if you agree to release your contribution under
   `AMICI's license conditions <https://github.com/AMICI-dev/AMICI/blob/master/LICENSE.md>`__.
   By opening a pull requests you confirm us that you do agree.

-  Start a new branch from ``main`` (on your fork, or at the main
   repository if you have access)

-  Implement your changes

-  Update ``CHANGELOG.md`` with a short description of your changes.
   This can be omitted for small changes that do not affect the
   functionality of AMICI, e.g. fixing typos or formatting issues,
   private refactoring, ...

-  Submit a pull request to the ``main`` branch

-  Ensure all tests pass

-  When adding new functionality, please also provide test cases (see
   ``tests/cpp/``, ``python/tests/``,
   and `documentation/CI.md <documentation/CI.md>`__)

-  Write meaningful commit messages

-  Run all tests to ensure nothing was broken (`more
   details <documentation/CI.md>`__)

   -  Run ``scripts/buildAll.sh && scripts/run-cpp-tests.sh``.

   -  If you made changes to the Matlab or C++ code and have a Matlab
      license, please also run ``tests/cpp/wrapTestModels.m`` and
      ``tests/testModels.m``

   -  If you made changes to the Python or C++ code, run
      ``make python-tests`` in ``build``

-  When all tests are passing and you think your code is ready to merge,
   request a code review (see also our `code review
   guideline <documentation/code_review_guide.md>`__)

-  Wait for feedback. If you do not receive feedback to your pull
   request within a week, please give us a friendly reminder.

Style/compatibility guide
~~~~~~~~~~~~~~~~~~~~~~~~~

General
^^^^^^^

-  All files and functions should come with file-level and
   function-level documentation.

-  All new functionality should be covered by unit or integration tests.
   Runtime of those tests should be kept as short as possible.

Python
^^^^^^

-  In terms of Python compatibility, we follow numpy's
   `NEP 29 <https://numpy.org/neps/nep-0029-deprecation_policy.html>`__.

-  For the Python code we want to follow
   `PEP8 <https://www.python.org/dev/peps/pep-0008/>`__. Although this
   is not the case for all existing code, any new contributions should
   do so. We use `ruff <https://docs.astral.sh/ruff/>`__ for automated
   code formatting.

   To run ruff as pre-commit hook, install the
   `pre-commit <https://pre-commit.com/>`_ package
   (e.g. ``pip install pre-commit``), and enable AMICI-hooks by running
   ``pre-commit install`` from within the AMICI directory.

-  We use Python `type
   hints <https://docs.python.org/3/library/typing.html>`__ for all
   functions and attributes. We do not include any other redundant
   type annotation in docstrings.

-  We use the `sphinx docstring-style <https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html>`__ for new code.

C++
^^^

-  We use C++20

-  We want to maintain compatibility with g++, clang, and the Intel C++
   compiler

-  For code formatting, we use ``clang-format`` and ``cmake-format``. They can
   be invoked by ``make clang-format cmake-format`` from the CMake build
   directory.

-  For new code, we use `Google's C++ style guide <https://google.github.io/styleguide/cppguide.html>`__ as a reference.


Branches / releases
-------------------

For AMICI, we mostly do `trunk-based development <https://trunkbaseddevelopment.com/>`__.
All new contributions are merged into ``main`` after passing the test suite
and code review. Releases are usually created directly from ``main``.
New releases are created on GitHub and are automatically deployed to
`Zenodo <https://doi.org/10.5281/zenodo.597928>`__ for
archiving and to obtain a digital object identifier (DOI) to make them
citable. Furthermore, our `CI pipeline <documentation/CI.md>`__ will
automatically create and deploy a new release on
`PyPI <https://pypi.org/project/amici/>`__.

We try to keep a clean git history. Therefore, feature pull requests are
squash-merged to ``main``.

Release process
~~~~~~~~~~~~~~~

Releases are created by the maintainer team.

To create a new release, please follow these steps:

1. Ensure that all changes intended for the new release are merged
   into ``main``.

2. Update ``CHANGELOG.md`` with a short description of the changes
   included in the new release. Focus on user-relevant changes.

3. Bump the version number in `version.txt` according to
   `Semantic Versioning <https://semver.org/>`__.

4. Regenerate the test models by running

    ```shell
    python -c "from amici.testing.models import import_test_models; import_test_models()"
    ```

    This ensures that the models can be imported with the new version.

5. Create a new release on GitHub, using the new version number prefixed
   by "v" (e.g., "v0.12.0"). Copy the relevant parts of
   ``CHANGELOG.md`` into the release notes.

6. After creating the release, our GitHub Actions pipeline will automatically
   create and deploy the new release on Zenodo and PyPI.
   Verify that this was successful.

7. Bump the version number in `version.txt` back to a development
   version (e.g., after release "1.1.0", set the version to "1.2.0-dev")
   and commit this change to ``main``.
   This ensures that documentation at https://amici.readthedocs.io/en/latest/
   will show the correct development version and won't be confused with the
   latest release, and that models imported with a development version
   will be marked as such.

In rare cases, it might be necessary to create a hotfix release for a
critical bug in an existing release. In this case, create a new branch
from the respective tag, apply the fix (usually a backport from ``main``),
and follow the same steps as above, starting from step 2.
Merge the updated CHANGELOG and version bump back into ``main`` afterwards.


Further topics
--------------

.. toctree::
   :maxdepth: 2

   Organization of the documentation <README>
   code_review_guide
   CI
   debugging
