#!/usr/bin/env bash
# Check that the test models checked in under `models/` are up to date with
# what the current code generator produces.
#
# Regenerates all test models in place, then diffs the result against the
# checked-in version. Differences in `get_amici_commit()` are ignored, since
# that value is expected to differ from run to run and does not indicate
# stale generated code. Any other difference -- including files that are no
# longer produced by the generator, or new files it now produces -- fails
# the check.
#
# Requires a working, installed `amici` package (see
# scripts/installAmiciSource.sh).
set -euo pipefail

SCRIPT_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
AMICI_PATH=$(cd "$SCRIPT_PATH/.." && pwd)
cd "$AMICI_PATH"

# This script force-cleans 'models/' on exit (see cleanup() below), which
# would silently destroy any pre-existing uncommitted or untracked content
# there. Refuse to run unless 'models/' already matches HEAD exactly.
if [[ -n "$(git status --porcelain -- models)" ]]; then
    echo "error: 'models/' has uncommitted or untracked changes." >&2
    echo "This script regenerates and force-cleans 'models/'; please commit" >&2
    echo "or stash any changes there first." >&2
    exit 1
fi

# Force-deletes any untracked files under 'models/' left over from the
# regeneration below. Safe only because of the preflight check above.
cleanup() {
    git reset --quiet HEAD -- models
    git checkout --quiet -- models
    git clean --quiet -fd -- models
}
trap cleanup EXIT

# Remove all currently tracked model files first, so that any file the
# generator no longer produces shows up as a deletion, rather than being
# silently left in place.
git rm -r --quiet models

# `compile=False`: we only care about the generated source here, not
# compiled binaries -- nothing downstream of this script needs them.
# Skipping compilation also avoids populating `models/` with build
# byproducts (CMake build tree, compiled `.so`, `__pycache__`, ...), which
# would otherwise get swept up by the force-add below.
python -c "from amici.testing.models import import_test_models; import_test_models(compile=False)"

# `-f` is required: the per-model directories under `models/` are only
# tracked because they were originally force-added -- they (and everything
# under them) are `.gitignore`d (see `models/*` in the top-level
# `.gitignore`, plus each model directory's own `**` `.gitignore`). Since
# `git rm -r` above fully unstaged them, a plain `git add -A` would silently
# skip re-adding the regenerated files as "ignored", leaving every model
# looking like a pure deletion instead of a diff against the regenerated
# content.
git add -A -f -- models

diff_output=$(git diff --cached -I 'return "[0-9a-f]{40}";' -- models)
if [[ -n "$diff_output" ]]; then
    echo ""
    echo "$diff_output"
    echo ""
    echo "The checked-in test models under 'models/' are out of date with" >&2
    echo "the current code generator (diff shown above). Please regenerate" >&2
    echo "them and commit the result:" >&2
    echo "" >&2
    echo '    python -c "from amici.testing.models import import_test_models; import_test_models()"' >&2
    echo "" >&2
    exit 1
fi

echo "Test models under 'models/' are up to date."
