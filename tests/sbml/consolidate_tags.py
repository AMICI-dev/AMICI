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
all_tags: set[str] = set()
# tags for which at least one test passed
supported_tags: set[str] = set()

for tag_file in result_dir.glob("tags_*.json"):
    with open(tag_file) as f:
        cur_tags = json.load(f)
    all_tags |= set(cur_tags["all"])
    supported_tags |= set(cur_tags["passed"])

# tags for which not a single test passed
unsupported_tags = set(all_tags) - set(supported_tags)

print("Supported tags")
print("--------------")
print()
print(",".join(sorted(list(supported_tags))))
print()
print("Unsupported tags")
print("----------------")
print()
print(",".join(sorted(list(unsupported_tags))))
