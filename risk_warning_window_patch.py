#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
risk_warning_window_patch.py
============================

Patch an individualized AD risk-timing analysis script to add GI-aware warning
window outputs.

This utility modifies an existing risk-timing script by:
1. fixing an off-by-one error in the 3-year risk-acceleration window loop;
2. adding GI burden and GI disease indicators to the individual-level risk-curve
   output table;
3. creating GI burden groups and any-GI indicators;
4. replacing the subgroup high-risk summary with richer warning-window summaries;
5. saving a backward-compatible alias for older downstream plotting scripts.

The script is intentionally conservative: it checks that each target code block
is found exactly once before writing the patched output file.

Example
-------
python risk_warning_window_patch.py \
    --input scripts/risk_timing_stratification.py \
    --output scripts/risk_timing_stratification_warning_windows.py

To inspect whether the patch can be applied without writing a file:

python risk_warning_window_patch.py \
    --input scripts/risk_timing_stratification.py \
    --dry-run
"""

from __future__ import annotations

import argparse
import difflib
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple


@dataclass
class PatchResult:
    """Summary of one patch operation."""

    name: str
    status: str
    replacements: int = 0
    message: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Patch an AD risk-timing script to add GI-aware warning-window "
            "outputs and subgroup summaries."
        )
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Path to the original risk-timing Python script.",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="",
        help=(
            "Path for the patched script. If omitted, a new file with suffix "
            "'_warning_windows.py' is created next to the input script."
        ),
    )
    parser.add_argument(
        "--backup",
        action="store_true",
        help="Create a .bak copy of the output path before overwriting it.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwriting an existing output file.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Check and preview the patch without writing the output file.",
    )
    parser.add_argument(
        "--show-diff",
        action="store_true",
        help="Print a unified diff between the original and patched script.",
    )
    return parser.parse_args()


def fail(message: str) -> None:
    raise SystemExit(f"[ERROR] {message}")


def default_output_path(input_path: Path) -> Path:
    return input_path.with_name(input_path.stem + "_warning_windows.py")


def replace_exact_once(text: str, old: str, new: str, patch_name: str) -> Tuple[str, PatchResult]:
    count = text.count(old)
    if count == 1:
        return text.replace(old, new, 1), PatchResult(patch_name, "patched", 1)
    if count == 0 and new in text:
        return text, PatchResult(patch_name, "already_patched", 0)
    fail(f"{patch_name}: expected exactly one target occurrence, found {count}.")


def replace_regex_once(text: str, pattern: re.Pattern[str], replacement: str, patch_name: str) -> Tuple[str, PatchResult]:
    new_text, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        fail(f"{patch_name}: expected exactly one regex replacement, made {count}.")
    return new_text, PatchResult(patch_name, "patched", count)


def patch_acceleration_window(text: str) -> Tuple[str, PatchResult]:
    old_loop = "for s in range(len(years) - RISK_ACCEL_WINDOW_WIDTH):"
    new_loop = "for s in range(len(years) - RISK_ACCEL_WINDOW_WIDTH + 1):"
    return replace_exact_once(text, old_loop, new_loop, "risk_acceleration_window_off_by_one")


def patch_curves_output_block(text: str) -> Tuple[str, PatchResult]:
    pattern = re.compile(
        r"    curves_all = derive_risk_timing\(risk_all, thresholds\).*?\n"
        r"    save_table\(curves_all, \"phaseC_primary_risk_curves_all\.csv\"\)",
        flags=re.DOTALL,
    )

    replacement = r'''    # ------------------------------------------------------------------
    # Individual-level risk curves and risk-timing metrics.
    # ------------------------------------------------------------------
    # Keep GI disease status and GI burden variables in the individual-level
    # output so that AD warning windows can be summarized by gastrointestinal
    # disease burden and disease category.
    gi_exposure_candidates = [
        "exp_diverticular",
        "exp_other_functional_intestinal",
        "exp_ibd",
        "exp_other_chronic_intestinal",
        "exp_ibs",
        "exp_malabsorption",
        "exp_appendiceal",
        "exp_anorectal",
    ]

    extra_cols = [
        "sex",
        "age_at_baseline",
        "apoe4_carrier",
        "ad_prs_std",
        "n_grouped_gi",
    ] + gi_exposure_candidates
    extra_cols = [c for c in extra_cols if c in analysis.columns]

    curves_all = derive_risk_timing(risk_all, thresholds).reset_index(drop=True).join(
        analysis[["eid", time_col, event_col] + extra_cols].reset_index(drop=True)
    )

    gi_exposure_cols = [c for c in gi_exposure_candidates if c in curves_all.columns]
    for c in gi_exposure_cols:
        curves_all[c] = ensure_binary(curves_all[c])

    # If n_grouped_gi is not available from preprocessing, reconstruct it from
    # retained grouped GI exposure indicators.
    if "n_grouped_gi" not in curves_all.columns and len(gi_exposure_cols) > 0:
        curves_all["n_grouped_gi"] = curves_all[gi_exposure_cols].apply(
            lambda row: pd.to_numeric(row, errors="coerce").sum(min_count=1),
            axis=1,
        )

    if len(gi_exposure_cols) > 0:
        curves_all["any_gi"] = (
            curves_all[gi_exposure_cols].max(axis=1, skipna=True) == 1
        ).astype(float)

    if "n_grouped_gi" in curves_all.columns:
        curves_all["gi_burden_group"] = pd.cut(
            safe_numeric(curves_all["n_grouped_gi"]),
            bins=[-0.1, 0.5, 1.5, 100],
            labels=["No GI disease", "Single GI disease", "Multiple GI diseases"],
        )

    if "age_at_baseline" in curves_all.columns:
        age = safe_numeric(curves_all["age_at_baseline"])
        curves_all["age_group"] = pd.cut(
            age,
            bins=[0, 55, 65, 75, 200],
            labels=["<55", "55-64", "65-74", "75+"],
        )

    risk_cols = [f"risk_{y}y" for y in RISK_YEARS]
    annual_cols = [f"annual_inc_{y}y" for y in RISK_YEARS]
    timing_cols = [
        "high_risk_5y",
        "high_risk_10y",
        "current_high_risk",
        "warning_onset_year",
        "peak_window_start",
        "peak_window_end",
        "peak_window_gain",
        "risk_accel_window_start",
        "risk_accel_window_end",
        "risk_accel_gain",
    ]
    background_cols = [
        "sex",
        "age_at_baseline",
        "age_group",
        "apoe4_carrier",
        "ad_prs_std",
        "n_grouped_gi",
        "gi_burden_group",
        "any_gi",
    ] + gi_exposure_cols
    background_cols = [c for c in background_cols if c in curves_all.columns]

    keep_cols = ["eid", time_col, event_col] + background_cols + risk_cols + annual_cols + timing_cols
    curves_all = curves_all[[c for c in keep_cols if c in curves_all.columns]]
    save_table(curves_all, "phaseC_primary_risk_curves_all.csv")'''

    return replace_regex_once(text, pattern, replacement, "gi_aware_primary_risk_curve_output")


def patch_subgroup_summary_block(text: str) -> Tuple[str, PatchResult]:
    pattern = re.compile(
        r"    # Subgroup high-risk summaries by selected risk-background variables\..*?\n"
        r"    save_table\(pd\.DataFrame\(subgroup_rows\), \"phaseC_subgroup_high_risk_summary\.csv\"\)",
        flags=re.DOTALL,
    )

    replacement = r'''    # ------------------------------------------------------------------
    # Subgroup summaries for individualized warning-window metrics.
    # ------------------------------------------------------------------
    # These summaries support Results text and supplementary figures/tables on
    # individualized AD warning windows by age, APOE and GI disease burden.
    def _nanmedian_or_na(x):
        x = safe_numeric(pd.Series(x)).dropna()
        return float(np.nanmedian(x)) if len(x) else np.nan

    def _nanq_or_na(x, q):
        x = safe_numeric(pd.Series(x)).dropna()
        return float(np.nanquantile(x, q)) if len(x) else np.nan

    def _add_timing_summary(rows, subgroup_var, subgroup, ss):
        if ss is None or len(ss) == 0:
            return
        rows.append({
            "subgroup_var": subgroup_var,
            "subgroup": str(subgroup),
            "n": int(len(ss)),
            "n_event": int(ss[event_col].sum()) if event_col in ss.columns else np.nan,
            "event_rate": float(ss[event_col].mean()) if event_col in ss.columns and len(ss) else np.nan,
            "prop_current_high_risk": float(ss["current_high_risk"].mean()) if "current_high_risk" in ss.columns else np.nan,
            "median_risk_5y": _nanmedian_or_na(ss["risk_5y"]) if "risk_5y" in ss.columns else np.nan,
            "median_risk_10y": _nanmedian_or_na(ss["risk_10y"]) if "risk_10y" in ss.columns else np.nan,
            "median_risk_15y": _nanmedian_or_na(ss["risk_15y"]) if "risk_15y" in ss.columns else np.nan,
            "median_risk_20y": _nanmedian_or_na(ss["risk_20y"]) if "risk_20y" in ss.columns else np.nan,
            "median_warning_onset_year": _nanmedian_or_na(ss["warning_onset_year"]) if "warning_onset_year" in ss.columns else np.nan,
            "q25_warning_onset_year": _nanq_or_na(ss["warning_onset_year"], 0.25) if "warning_onset_year" in ss.columns else np.nan,
            "q75_warning_onset_year": _nanq_or_na(ss["warning_onset_year"], 0.75) if "warning_onset_year" in ss.columns else np.nan,
            "median_peak_window_start": _nanmedian_or_na(ss["peak_window_start"]) if "peak_window_start" in ss.columns else np.nan,
            "median_peak_window_end": _nanmedian_or_na(ss["peak_window_end"]) if "peak_window_end" in ss.columns else np.nan,
            "median_peak_window_gain": _nanmedian_or_na(ss["peak_window_gain"]) if "peak_window_gain" in ss.columns else np.nan,
            "median_risk_accel_window_start": _nanmedian_or_na(ss["risk_accel_window_start"]) if "risk_accel_window_start" in ss.columns else np.nan,
            "median_risk_accel_window_end": _nanmedian_or_na(ss["risk_accel_window_end"]) if "risk_accel_window_end" in ss.columns else np.nan,
            "median_risk_accel_gain": _nanmedian_or_na(ss["risk_accel_gain"]) if "risk_accel_gain" in ss.columns else np.nan,
        })

    subgroup_rows = []

    for gcol in ["sex", "age_group", "apoe4_carrier", "gi_burden_group", "any_gi"]:
        if gcol in curves_all.columns:
            for g, ss in curves_all.groupby(gcol, observed=False):
                if pd.isna(g) or len(ss) == 0:
                    continue
                _add_timing_summary(subgroup_rows, gcol, g, ss)

    gi_abbr = {
        "exp_diverticular": "DIV",
        "exp_other_functional_intestinal": "OFID",
        "exp_ibd": "IBD",
        "exp_other_chronic_intestinal": "OCID",
        "exp_ibs": "IBS",
        "exp_malabsorption": "MAL",
        "exp_appendiceal": "APD",
        "exp_anorectal": "ARD",
    }
    for exp_col in [c for c in gi_abbr if c in curves_all.columns]:
        for level, ss in curves_all.groupby(exp_col):
            if pd.isna(level) or len(ss) == 0:
                continue
            label = f"{gi_abbr[exp_col]}+" if float(level) == 1.0 else f"{gi_abbr[exp_col]}-"
            _add_timing_summary(subgroup_rows, f"GI disease: {gi_abbr[exp_col]}", label, ss)

    warning_summary = pd.DataFrame(subgroup_rows)
    save_table(warning_summary, "phaseC_warning_window_subgroup_summary.csv")

    # Backward-compatible alias used by previous figure scripts.
    save_table(warning_summary, "phaseC_subgroup_high_risk_summary.csv")'''

    return replace_regex_once(text, pattern, replacement, "gi_aware_warning_window_subgroup_summary")


def patch_report_preview(text: str) -> Tuple[str, PatchResult]:
    if "Warning-window subgroup summary preview" in text:
        return text, PatchResult("warning_window_report_preview", "already_patched", 0)

    old = (
        '        "High-risk summary:",\n'
        '        hr_summary.to_string(index=False) if len(hr_summary) > 0 else "No high-risk summary.",\n'
    )
    new = (
        '        "High-risk summary:",\n'
        '        hr_summary.to_string(index=False) if len(hr_summary) > 0 else "No high-risk summary.",\n'
        '        "",\n'
        '        "Warning-window subgroup summary preview:",\n'
        '        warning_summary.head(60).to_string(index=False) if len(warning_summary) > 0 else "No warning-window subgroup summary.",\n'
    )

    count = text.count(old)
    if count == 1:
        return text.replace(old, new, 1), PatchResult("warning_window_report_preview", "patched", 1)
    if count == 0:
        return text, PatchResult(
            "warning_window_report_preview",
            "skipped",
            0,
            "High-risk report preview block was not found; core patch was still applied.",
        )
    fail(f"warning_window_report_preview: expected zero or one target occurrence, found {count}.")


def apply_patch(text: str) -> Tuple[str, List[PatchResult]]:
    results: List[PatchResult] = []
    patched = text

    for patch_fn in [
        patch_acceleration_window,
        patch_curves_output_block,
        patch_subgroup_summary_block,
        patch_report_preview,
    ]:
        patched, result = patch_fn(patched)
        results.append(result)

    if patched == text:
        fail("No changes were made; the input may already be fully patched or the patterns changed.")
    return patched, results


def print_results(results: List[PatchResult]) -> None:
    print("\nPatch summary")
    print("-" * 72)
    for r in results:
        msg = f" | {r.message}" if r.message else ""
        print(f"{r.name}: {r.status}; replacements={r.replacements}{msg}")
    print("-" * 72)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).expanduser().resolve()
    if not input_path.exists():
        fail(f"Input script not found: {input_path}")
    if input_path.suffix != ".py":
        fail(f"Input file should be a Python script: {input_path}")

    output_path = Path(args.output).expanduser().resolve() if args.output else default_output_path(input_path)

    if output_path.exists() and not args.overwrite and not args.dry_run:
        fail(
            f"Output file already exists: {output_path}. "
            "Use --overwrite to replace it or choose a different --output path."
        )

    original = input_path.read_text(encoding="utf-8")
    patched, results = apply_patch(original)
    print_results(results)

    if args.show_diff or args.dry_run:
        diff = difflib.unified_diff(
            original.splitlines(),
            patched.splitlines(),
            fromfile=str(input_path),
            tofile=str(output_path),
            lineterm="",
        )
        print("\n".join(diff))

    if args.dry_run:
        print("\n[DRY-RUN] Patch can be applied. No file was written.")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and args.backup:
        backup_path = output_path.with_suffix(output_path.suffix + ".bak")
        shutil.copy2(output_path, backup_path)
        print(f"[BACKUP] {backup_path}")

    output_path.write_text(patched, encoding="utf-8", newline="\n")
    print(f"[SAVE] Patched script written to: {output_path}")


if __name__ == "__main__":
    main()
