"""
Summarize (un)supported tags from SBML semantic test suite runs.

Expected to be run as part of test result consolidation GHA job.

The list of unsupported tags can be pasted into the SBML test result database
submission form.
"""

import json
from pathlib import Path

# where all result artifacts have been unpacked to
result_dir = Path("combined")

# tags encountered across all tests (from the simulation check only)
encountered_tags: set[str] = set()
# tags for which at least one test passed
supported_tags: set[str] = set()

# test IDs of passed tests, by check
passed_ids: dict[str, set[str]] = {"simulation": set(), "sensitivity": set()}
# failed or skipped tests with error message, by check
failed_or_skipped: dict[str, dict[str, str]] = {
    "simulation": dict(),
    "sensitivity": dict(),
}

for tag_file in result_dir.glob("results_*.json"):
    with open(tag_file) as f:
        cur_tags = json.load(f)
    encountered_tags |= set(cur_tags["encountered_tags"])
    supported_tags |= set(cur_tags["supported_tags"])
    for check in ("simulation", "sensitivity"):
        passed_ids[check] |= set(cur_tags[f"passed_tests_{check}"])
        failed_or_skipped[check] |= cur_tags[f"failed_or_skipped_{check}"]

for check in ("simulation", "sensitivity"):
    num_tests_success = len(passed_ids[check])
    num_tests_total = num_tests_success + len(failed_or_skipped[check])
    frac_tests_passed = num_tests_success / num_tests_total
    print(
        f"[{check}] {num_tests_success}/{num_tests_total} ≈ "
        f"{frac_tests_passed:.2%} tests passed."
    )

# tags for which not a single test passed
unsupported_tags = set(encountered_tags) - set(supported_tags)

print("Supported tags")
print("--------------")
print()
print(",".join(sorted(list(supported_tags))))
print()
print("Unsupported tags")
print("----------------")
print()
print(",".join(sorted(list(unsupported_tags))))
print()
for check in ("simulation", "sensitivity"):
    print(f"Failed or-skipped tests [{check}]")
    print("-----------------------" + "-" * len(check))
    print()
    for test_id in sorted(failed_or_skipped[check]):
        msg = failed_or_skipped[check][test_id]
        print(f"{test_id}: {msg}")
    print()
