#!/usr/bin/env python3
"""Validate key connectivity totals against values reported in the manuscript."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = ROOT / "config" / "expected_summary_counts.csv"
OBSERVED = ROOT / "outputs" / "connectivity_summary_all.csv"


def main() -> None:
    if not OBSERVED.exists():
        raise SystemExit(
            f"Missing {OBSERVED}. Run scripts/01_synapse_counts.py first."
        )

    exp = pd.read_csv(EXPECTED)
    obs = pd.read_csv(OBSERVED)
    obs["weightHP"] = pd.to_numeric(obs["weightHP"], errors="coerce")

    rows = []
    failures = 0
    for _, e in exp.iterrows():
        hit = obs[
            (obs["dataset"] == e["dataset"])
            & (obs["side"] == e["side"])
            & (obs["direction"] == e["direction"])
            & (obs["group"] == e["group"])
        ]
        observed = int(hit["weightHP"].sum()) if len(hit) else 0
        expected = int(e["expected_weightHP"])
        ok = observed == expected
        failures += int(not ok)
        rows.append({
            "dataset": e["dataset"],
            "side": e["side"],
            "direction": e["direction"],
            "group": e["group"],
            "expected_weightHP": expected,
            "observed_weightHP": observed,
            "match": ok,
        })

    report = pd.DataFrame(rows)
    report.to_csv(ROOT / "outputs" / "validation_report.csv", index=False)
    print(report.to_string(index=False))

    if failures:
        print(
            f"\nVALIDATION FAILED: {failures} key value(s) did not match. "
            "Audit annotation/grouping rules before public release.",
            file=sys.stderr,
        )
        raise SystemExit(1)

    print("\nVALIDATION PASSED: all configured key totals matched.")


if __name__ == "__main__":
    main()
