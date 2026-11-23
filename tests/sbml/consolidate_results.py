"""
Summarize (un)supported tags from SBML semantic test suite runs.

Expected to be run as part of test result consolidation GHA job.

The list of unsupported tags can be pasted into the SBML test result database
submission form.
"""

from pathlib import Path

from suite_results import CHECKS, SuiteResults

# where all result artifacts have been unpacked to
result_dir = Path("combined")

results = SuiteResults()
for result_file in result_dir.glob("results_*.json"):
    results.update(SuiteResults.from_json(result_file))

for check in CHECKS:
    num_tests_success = len(results.passed[check])
    num_tests_total = num_tests_success + len(results.failed_or_skipped[check])
    frac_tests_passed = num_tests_success / num_tests_total
    print(
        f"[{check}] {num_tests_success}/{num_tests_total} ≈ "
        f"{frac_tests_passed:.2%} tests passed."
    )

print("Supported tags")
print("--------------")
print()
print(",".join(sorted(results.supported_tags)))
print()
print("Unsupported tags")
print("----------------")
print()
print(",".join(sorted(results.unsupported_tags)))
print()
for check in CHECKS:
    print(f"Failed or-skipped tests [{check}]")
    print("-----------------------" + "-" * len(check))
    print()
    for test_id in sorted(results.failed_or_skipped[check]):
        msg = results.failed_or_skipped[check][test_id]
        print(f"{test_id}: {msg}")
    print()
