"""Compare two jax_regression_results.json files and emit a markdown summary.

Usage:
    python compare_jax_results.py <prev.json> <curr.json>

Writes a markdown comparison to stdout.  Exits 0 always (regressions are
warnings, not hard errors).  Designed to be piped into $GITHUB_STEP_SUMMARY.

Regression thresholds:
- Step counts: exact equality required (any change is flagged)
- LLH: relative tolerance 1e-4
- Execution time: 2x relative increase (100 % slower) flagged
"""

import json
import sys
from pathlib import Path


def _load(path: str) -> dict:
    p = Path(path)
    if not p.exists():
        return {}
    with open(p) as f:
        return json.load(f)


def _fmt_pct(ratio: float) -> str:
    pct = (ratio - 1.0) * 100
    sign = "+" if pct >= 0 else ""
    return f"{sign}{pct:.1f}%"


def compare(prev: dict, curr: dict) -> tuple[str, int]:
    """Return (markdown_report, n_warnings)."""
    lines = []
    n_warn = 0

    # Header and metadata
    prev_meta = prev.pop("_metadata", {})
    curr_meta = curr.pop("_metadata", {})

    lines.append("## JAX Regression Comparison\n")
    if prev_meta or curr_meta:
        lines.append("### Dependency versions\n")
        lines.append("| Dep | Previous | Current |")
        lines.append("|-----|---------|---------|")
        all_keys = sorted(set(prev_meta) | set(curr_meta))
        for k in all_keys:
            pv = prev_meta.get(k, "–")
            cv = curr_meta.get(k, "–")
            flag = " ⚠️" if pv != cv else ""
            lines.append(f"| {k} | {pv} | {cv}{flag} |")
        lines.append("")

    dep_changed = any(prev_meta.get(k) != curr_meta.get(k) for k in curr_meta)

    # Per-model comparison
    all_models = sorted(set(prev) | set(curr))
    step_issues = []
    llh_issues = []
    timing_issues = []

    for model_id in all_models:
        prev_m = prev.get(model_id, {})
        curr_m = curr.get(model_id, {})
        all_ops = sorted(set(prev_m) | set(curr_m))

        for op in all_ops:
            prev_op = prev_m.get(op, {})
            curr_op = curr_m.get(op, {})

            tag = f"{model_id}/{op}"

            # Step counts
            for stat_key in ("stats_dyn", "stats_posteq", "stats_preeq"):
                prev_s = prev_op.get(stat_key)
                curr_s = curr_op.get(stat_key)
                if prev_s is None or curr_s is None:
                    continue
                for cnt_key in ("num_accepted_steps", "num_rejected_steps"):
                    pv = prev_s.get(cnt_key)
                    cv = curr_s.get(cnt_key)
                    if pv is None or cv is None:
                        continue
                    if pv != cv:
                        step_issues.append(
                            f"| `{tag}` | `{stat_key}.{cnt_key}` | {pv} | {cv} |"
                        )
                        n_warn += 1

            # LLH
            pv_llh = prev_op.get("llh")
            cv_llh = curr_op.get("llh")
            if pv_llh is not None and cv_llh is not None:
                if abs(cv_llh - pv_llh) > 1e-4 * (abs(pv_llh) + 1e-10):
                    llh_issues.append(
                        f"| `{tag}` | {pv_llh:.6g} | {cv_llh:.6g} | "
                        f"{cv_llh - pv_llh:+.3g} |"
                    )
                    n_warn += 1

            # x_preeq (for preeq operation)
            pv_xp = prev_op.get("x_preeq")
            cv_xp = curr_op.get("x_preeq")
            if pv_xp is not None and cv_xp is not None and pv_xp != cv_xp:
                llh_issues.append(
                    f"| `{tag}` | x_preeq={pv_xp} | x_preeq={cv_xp} | changed |"
                )
                n_warn += 1

            # Timing (2x threshold)
            for t_key in ("t_exec_s",):
                for sub in ("", "fwd.", "adj."):
                    full_key = f"{sub}{t_key}"
                    if "." in full_key:
                        subdict, sk = full_key.split(".", 1)
                        pv_t = (prev_op.get(subdict) or {}).get(sk)
                        cv_t = (curr_op.get(subdict) or {}).get(sk)
                    else:
                        pv_t = prev_op.get(full_key)
                        cv_t = curr_op.get(full_key)
                    if pv_t is None or cv_t is None or pv_t <= 0:
                        continue
                    ratio = cv_t / pv_t
                    if ratio > 2.0:
                        timing_issues.append(
                            f"| `{tag}` | `{full_key}` | {pv_t * 1e3:.1f} ms "
                            f"| {cv_t * 1e3:.1f} ms | {_fmt_pct(ratio)} |"
                        )
                        n_warn += 1

    if step_issues:
        lines.append("### ⚠️ Step count changes (deterministic — investigate)")
        lines.append("")
        lines.append("| Test | Metric | Previous | Current |")
        lines.append("|------|--------|---------|---------|")
        lines.extend(step_issues)
        lines.append("")
    else:
        lines.append("### ✅ Step counts — no changes\n")

    if llh_issues:
        lines.append("### ⚠️ LLH changes (rtol > 1e-4)")
        lines.append("")
        lines.append("| Test | Previous | Current | Delta |")
        lines.append("|------|---------|---------|-------|")
        lines.extend(llh_issues)
        lines.append("")
    else:
        lines.append("### ✅ LLH — within tolerance\n")

    if timing_issues:
        lines.append("### ⚠️ Timing regressions (>2× slower)")
        lines.append("")
        lines.append("| Test | Metric | Previous | Current | Change |")
        lines.append("|------|--------|---------|---------|--------|")
        lines.extend(timing_issues)
        lines.append("")
    else:
        lines.append("### ✅ Timing — no significant regressions\n")

    if n_warn == 0:
        lines.insert(2, "> ✅ No regressions detected.\n")
    else:
        lines.insert(
            2,
            f"> ⚠️ **{n_warn} potential regression(s) detected** — see details below.\n",
        )

    if dep_changed:
        lines.insert(3, "> ℹ️ Dependency versions changed between runs.\n")

    return "\n".join(lines), n_warn


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <prev.json> <curr.json>", file=sys.stderr)
        sys.exit(1)

    prev = _load(sys.argv[1])
    curr = _load(sys.argv[2])

    if not prev:
        print("## JAX Regression Comparison\n")
        print("> ℹ️ No previous results available — this is the first run.\n")
        return

    report, n_warn = compare(prev, curr)
    print(report)


if __name__ == "__main__":
    main()
