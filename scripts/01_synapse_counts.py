#!/usr/bin/env python3
"""Reproduce SIFa high-precision connectivity summaries from neuPrint.

This script consolidates the final logic used in an originally interactive
custom-Cypher workflow. It is intentionally transparent: raw queried
connections are saved before any grouping is applied.

Authentication
--------------
Set NEUPRINT_APPLICATION_CREDENTIALS in the environment. Never commit a token.
"""

from __future__ import annotations

import re
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd
from neuprint import Client

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"
OUTPUTS = ROOT / "outputs"
SERVER = "neuprint.janelia.org"

DATASETS = ("hemibrain:v1.2.1", "male-cns:v1.0")

# Fields are queried explicitly so that the raw CSV preserves the labels used
# for grouping and permits later audit/reclassification.
PARTNER_FIELDS = [
    "bodyId",
    "type",
    "instance",
    "hemibrainType",
    "class",
    "subclass",
    "superclass",
    "somaSide",
    "rootSide",
]


def cypher_value(var: str, prop: str, alias: str) -> str:
    return f"{var}.`{prop}` AS `{alias}`"


def query_connections(client: Client, sifa_ids: list[int], direction: str) -> pd.DataFrame:
    """Return all high-precision connections involving the supplied SIFa IDs."""
    if direction == "SIFa_to_target":
        match = "MATCH (s:Neuron)-[c:ConnectsTo]->(t:Neuron)"
        where = "WHERE s.bodyId IN $sifa_ids AND coalesce(c.weightHP, 0) > 0"
        s_var, t_var = "s", "t"
    elif direction == "target_to_SIFa":
        match = "MATCH (t:Neuron)-[c:ConnectsTo]->(s:Neuron)"
        where = "WHERE s.bodyId IN $sifa_ids AND coalesce(c.weightHP, 0) > 0"
        s_var, t_var = "s", "t"
    else:
        raise ValueError(direction)

    returns = [
        f"{s_var}.bodyId AS sifa_body_id",
        f"{t_var}.bodyId AS partner_body_id",
        "c.weight AS weight",
        "c.weightHP AS weightHP",
    ]
    for field in PARTNER_FIELDS[1:]:
        returns.append(cypher_value(t_var, field, f"partner_{field}"))

    q = "\n".join([
        match,
        where,
        "RETURN " + ",\n       ".join(returns),
        "ORDER BY sifa_body_id, partner_body_id",
    ])

    # neuprint-python's Client.fetch_custom does not provide parameter binding,
    # so insert a validated numeric list only.
    numeric = ", ".join(str(int(x)) for x in sifa_ids)
    q = q.replace("$sifa_ids", f"[{numeric}]")
    df = client.fetch_custom(q)
    df["direction"] = direction
    return df


def asciiish(text: object) -> str:
    if text is None or (isinstance(text, float) and np.isnan(text)):
        return ""
    s = str(text)
    s = s.replace("α", "alpha").replace("β", "beta").replace("γ", "gamma")
    s = s.replace("′", "'").replace("’", "'").replace("`", "'")
    s = unicodedata.normalize("NFKD", s)
    s = "".join(ch for ch in s if not unicodedata.combining(ch))
    return s.lower()


def partner_text(row: pd.Series) -> str:
    fields = [
        row.get("partner_type"),
        row.get("partner_instance"),
        row.get("partner_hemibrainType"),
        row.get("partner_class"),
        row.get("partner_subclass"),
        row.get("partner_superclass"),
    ]
    return " | ".join(asciiish(x) for x in fields if asciiish(x))


def infer_side(row: pd.Series, dataset: str) -> str:
    """Infer L/R side from metadata; Hemibrain defaults to reference/right side."""
    if dataset.startswith("hemibrain"):
        return "R"

    vals = [row.get("partner_somaSide"), row.get("partner_rootSide"), row.get("partner_instance")]
    txt = " | ".join(asciiish(x) for x in vals if asciiish(x))
    if re.search(r"\b(lhs|left)\b|_l\b", txt):
        return "L"
    if re.search(r"\b(rhs|right)\b|_r\b", txt):
        return "R"
    return "unknown"


def classify_partner(row: pd.Series) -> str:
    """Map partner annotations to the manuscript-level groups.

    These rules are intentionally explicit and should be audited against the
    final Supplementary Tables before repository release.
    """
    text = partner_text(row)
    compact = re.sub(r"[^a-z0-9']+", "", text)

    # Clock-neuron groups.
    if "5ths-lnv" in text or "5th s-lnv" in text or "lnd6" in text:
        return "fifth_s_LNv_LNd6_associated"
    if "l-lnv" in text or "llnv" in compact:
        return "PDF_positive_LNv"
    if ("s-lnv" in text or "slnv" in compact) and "5th" not in text:
        return "PDF_positive_LNv"
    if re.search(r"\blnd[1-5]\b", compact) or "lnd1-5" in text:
        return "LNd1_5"

    # Kenyon-cell classes. Hemibrain-style type labels include KCab, KCg,
    # and KCa'b'; MaleCNS provides cross-dataset hemibrainType annotations.
    if "kca'b'" in compact or "kcalpha'beta'" in compact or "kcalphaprimebetaprime" in compact:
        return "KC_alpha_prime_beta_prime"
    if compact.startswith("kcab") or "kcalphabeta" in compact:
        return "KC_alpha_beta"
    if compact.startswith("kcg") or "kcgamma" in compact:
        return "KC_gamma"

    # Other MB-associated classes used in Supplementary Table 4.
    if "mbon" in compact:
        return "MBON"
    if "dpm" in compact:
        return "DPM"
    if "apl" in compact:
        return "APL"
    # Typical fly DAN types are PAM/PPL; annotation class may also explicitly
    # contain dopaminergic/DAN.
    if re.search(r"\b(pam|ppl)[a-z0-9]", compact) or "dopaminergic" in text or "dan" in text:
        return "DAN"

    return "other"


def make_summary(raw: pd.DataFrame, dataset: str) -> pd.DataFrame:
    df = raw.copy()
    df["group"] = df.apply(classify_partner, axis=1)
    df["side"] = df.apply(lambda r: infer_side(r, dataset), axis=1)
    df["weightHP"] = pd.to_numeric(df["weightHP"], errors="coerce").fillna(0).astype(int)

    summary = (
        df.groupby(["direction", "side", "group"], dropna=False, as_index=False)["weightHP"]
        .sum()
        .sort_values(["direction", "side", "group"])
    )
    summary.insert(0, "dataset", dataset)

    # Add bilateral totals for MaleCNS to match the KC totals reported in text.
    if dataset.startswith("male-cns"):
        lr = summary[summary["side"].isin(["L", "R"])].copy()
        bilateral = (
            lr.groupby(["dataset", "direction", "group"], as_index=False)["weightHP"]
            .sum()
        )
        bilateral.insert(2, "side", "L+R")
        summary = pd.concat([summary, bilateral], ignore_index=True)
        summary = summary.sort_values(["direction", "side", "group"]).reset_index(drop=True)

    return df, summary


def main() -> None:
    OUTPUTS.mkdir(exist_ok=True)
    ids = pd.read_csv(CONFIG / "neuron_ids.csv")

    all_summaries = []
    meta_rows = []

    for dataset in DATASETS:
        sifa_ids = ids.loc[ids["dataset"] == dataset, "body_id"].astype(int).tolist()
        if len(sifa_ids) != 4:
            raise RuntimeError(f"Expected four SIFa IDs for {dataset}; found {len(sifa_ids)}")

        client = Client(SERVER, dataset=dataset, progress=False)
        meta = client.fetch_custom(
            "MATCH (m:Meta) RETURN m.postHPThreshold AS postHPThreshold, "
            "m.preHPThreshold AS preHPThreshold"
        )
        if not meta.empty:
            meta_rows.append({
                "dataset": dataset,
                "postHPThreshold": meta.iloc[0].get("postHPThreshold"),
                "preHPThreshold": meta.iloc[0].get("preHPThreshold"),
            })

        parts = [
            query_connections(client, sifa_ids, "SIFa_to_target"),
            query_connections(client, sifa_ids, "target_to_SIFa"),
        ]
        raw = pd.concat(parts, ignore_index=True)
        raw["dataset"] = dataset

        classified, summary = make_summary(raw, dataset)
        safe = dataset.replace(":", "_")
        classified.to_csv(OUTPUTS / f"raw_connections_{safe}.csv", index=False)
        summary.to_csv(OUTPUTS / f"grouped_summary_{safe}.csv", index=False)
        all_summaries.append(summary)

    combined = pd.concat(all_summaries, ignore_index=True)
    combined.to_csv(OUTPUTS / "connectivity_summary_all.csv", index=False)
    pd.DataFrame(meta_rows).to_csv(OUTPUTS / "neuprint_hp_thresholds.csv", index=False)

    # Convenience exports corresponding to the two supplementary tables.
    t3_groups = {"PDF_positive_LNv", "LNd1_5", "fifth_s_LNv_LNd6_associated"}
    t4_groups = {"KC_alpha_beta", "KC_gamma", "KC_alpha_prime_beta_prime", "DAN", "MBON", "DPM", "APL"}
    combined[combined["group"].isin(t3_groups)].to_csv(
        OUTPUTS / "supplementary_table3_reproduced.csv", index=False
    )
    combined[combined["group"].isin(t4_groups)].to_csv(
        OUTPUTS / "supplementary_table4_reproduced.csv", index=False
    )

    print(f"Wrote outputs to: {OUTPUTS}")
    print("Next: python scripts/03_validate_expected_counts.py")


if __name__ == "__main__":
    main()
