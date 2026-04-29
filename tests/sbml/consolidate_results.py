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

# tags encountered across all tests
encountered_tags: set[str] = set()
# tags for which at least one test passed
supported_tags: set[str] = set()

# test IDs of passed tests
passed_ids: set[str] = set()
# failed or skipped tests with error message
failed_or_skipped: dict[str, str] = dict()

for tag_file in result_dir.glob("results_*.json"):
    with open(tag_file) as f:
        cur_tags = json.load(f)
    encountered_tags |= set(cur_tags["encountered_tags"])
    supported_tags |= set(cur_tags["supported_tags"])
    passed_ids |= set(cur_tags["passed_tests"])
    failed_or_skipped |= cur_tags["failed_or_skipped"]

num_tests_success = len(passed_ids)
num_tests_total = num_tests_success + len(failed_or_skipped)
frac_tests_passed = num_tests_success / num_tests_total
print(
    f"{num_tests_success}/{num_tests_total} ≈ "
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
print("Failed or-skipped tests")
print("-----------------------")
print()
for test_id in sorted(failed_or_skipped):
    msg = failed_or_skipped[test_id]
    print(f"{test_id}: {msg}")
print()
