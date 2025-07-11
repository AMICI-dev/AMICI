AMICI developer’s guide
=======================

This document contains information for AMICI developers, not too
relevant to regular users.

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


Matlab
^^^^^^

No new Matlab code should be added to AMICI
(see `#2727 <https://github.com/AMICI-dev/AMICI/issues/2727>`__).

Further topics
--------------

.. toctree::
   :maxdepth: 2

   Organization of the documentation <README>
   code_review_guide
   CI
   debugging
