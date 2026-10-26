"""Report the size and the bundled shared libraries of AMICI wheels.

Fails if a wheel contains the same shared library more than once, e.g.
``libscipy_openblas64_`` both in ``amici/_installation/lib`` and, vendored
again by auditwheel or delvewheel, in ``amici.libs``.

Usage::

    python inspect_wheels.py WHEEL [WHEEL ...]

If ``GITHUB_STEP_SUMMARY`` is set, a Markdown summary is appended to it.
"""

import os
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import PurePosixPath

SHARED_LIB_RE = re.compile(r"\.(so(\.\d+)*|dylib|dll|pyd)$")
# hashes that auditwheel / delvewheel / delocate add to vendored libraries
VENDOR_HASH_RE = re.compile(r"-[0-9a-f]{8,}(?=\.)")


def lib_key(name: str) -> str:
    """Library name with version suffix and vendoring hash removed"""
    name = VENDOR_HASH_RE.sub("", name)
    return re.sub(r"\.(so(\.\d+)*|dylib|dll|pyd)$", "", name)


def inspect_wheel(path: str) -> tuple[list[str], bool]:
    lines = []
    with zipfile.ZipFile(path) as zf:
        infos = zf.infolist()
    size = os.path.getsize(path)
    unpacked = sum(i.file_size for i in infos)
    lines.append(
        f"### `{os.path.basename(path)}`\n\n"
        f"{size / 2**20:.1f} MiB, {unpacked / 2**20:.1f} MiB unpacked, "
        f"{len(infos)} files\n"
    )

    libs = defaultdict(list)
    for info in infos:
        name = PurePosixPath(info.filename).name
        if SHARED_LIB_RE.search(name):
            libs[lib_key(name)].append(info)

    lines.append("| shared library | MiB |\n|---|---:|")
    for key in sorted(libs):
        for info in libs[key]:
            lines.append(
                f"| `{info.filename}` | {info.file_size / 2**20:.1f} |"
            )
    lines.append("")

    largest = sorted(infos, key=lambda i: i.file_size, reverse=True)[:10]
    lines.append("<details><summary>Largest files</summary>\n")
    lines.append("| file | MiB |\n|---|---:|")
    for info in largest:
        lines.append(f"| `{info.filename}` | {info.file_size / 2**20:.1f} |")
    lines.append("\n</details>\n")

    duplicates = {k: v for k, v in libs.items() if len(v) > 1}
    for key, infos in duplicates.items():
        lines.append(
            f"**Duplicate library `{key}`:** "
            + ", ".join(f"`{i.filename}`" for i in infos)
            + "\n"
        )
    return lines, bool(duplicates)


def main():
    report = ["## Wheel contents\n"]
    found_duplicates = False
    for path in sys.argv[1:]:
        lines, duplicates = inspect_wheel(path)
        report.extend(lines)
        found_duplicates |= duplicates

    report = "\n".join(report)
    print(report)
    if summary := os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(summary, "a") as f:
            f.write(report + "\n")

    if found_duplicates:
        sys.exit("Found shared libraries that are included more than once.")


if __name__ == "__main__":
    main()
