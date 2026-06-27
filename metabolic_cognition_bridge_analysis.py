#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Metabolic-cognition bridge analysis.

This script performs marker-level pathway screening for:

    gastrointestinal disease -> metabolic marker -> cognition phenotype

It produces exact path-level tables for Sankey/alluvial plots and supporting
analyses linking cognition phenotypes to incident Alzheimer's disease (AD) and
all-cause dementia.

Inputs
------
1. Master participant-level table with GI exposures, AD/dementia outcomes,
   covariates, and genetics.
2. Prepared cognition phenotype table.
3. Prepared metabolic marker table.
4. Optional metabolic marker label map with columns: feature,label.

Outputs
-------
results/metabolic_cognition_bridge_analysis/
  tables/
    stage1_gi_to_marker_all.csv
    stage2_marker_to_cognition_all.csv
    marker_cognition_bridge_paths_all.csv
    marker_cognition_bridge_paths_significant.csv
    sankey_source_paths.csv
    gi_marker_cognition_mediation_all.csv
    support_cognition_to_ad_cox.csv
    support_cognition_to_dementia_cox.csv
    marker_cognition_summary_table.csv
  figures/
  qc/
  run_metadata.json
  analysis_report.txt
  run.log

Example
-------
python metabolic_cognition_bridge_analysis.py \
    --master data/master_preprocessed_with_genetics.csv \
    --cognition data/cognition_prepared.csv \
    --metabolic data/metabolic_prepared.csv \
    --outdir results/metabolic_cognition_bridge_analysis
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import sys
from dataclasses import dataclass
from math import erfc
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from lifelines import CoxPHFitter
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression


# =============================================================================
# Constants
# =============================================================================

MAIN_EXPOSURES = [
    "exp_diverticular",
    "exp_other_functional_intestinal",
    "exp_ibd",
    "exp_other_chronic_intestinal",
    "exp_ibs",
]
SECONDARY_EXPOSURES = ["exp_appendiceal", "exp_malabsorption"]
ALL_EXPOSURES_ORDERED = MAIN_EXPOSURES + SECONDARY_EXPOSURES

GROUPED_GI_MAP = {
    "exp_ibd": ["exp_k50", "exp_k51"],
    "exp_ibs": ["exp_k58"],
    "exp_other_functional_intestinal": ["exp_k59"],
    "exp_diverticular": ["exp_k57"],
    "exp_anorectal": ["exp_k60", "exp_k61", "exp_k62", "exp_k64"],
    "exp_malabsorption": ["exp_k90"],
    "exp_other_chronic_intestinal": ["exp_k52", "exp_k55", "exp_k56", "exp_k63"],
    "exp_appendiceal": ["exp_k35", "exp_k36", "exp_k37", "exp_k38"],
}

EXPOSURE_LABELS = {
    "exp_diverticular": "Diverticular disease",
    "exp_other_functional_intestinal": "Other functional intestinal disorders",
    "exp_ibd": "Inflammatory bowel disease",
    "exp_other_chronic_intestinal": "Other chronic intestinal disease",
    "exp_ibs": "Irritable bowel syndrome",
    "exp_appendiceal": "Appendiceal disease",
    "exp_malabsorption": "Malabsorption",
}

EXPOSURE_ABBR = {
    "exp_diverticular": "DIV",
    "exp_other_functional_intestinal": "OFID",
    "exp_other_chronic_intestinal": "OCID",
    "exp_ibd": "IBD",
    "exp_ibs": "IBS",
    "exp_appendiceal": "APD",
    "exp_malabsorption": "MAL",
}

COGNITION_DOMAIN_CANDIDATES = {
    "cognition_speed_pc1": [
        "reaction_time_mean",
        "symbol_digit_entry_time_mean",
        "symbol_digit_entry_time_median",
        "trail1_duration",
        "trail2_duration",
        "pm_complete_time_sum",
        "pm_complete_time_mean",
    ],
    "cognition_accuracy_pc1": [
        "fluid_intelligence",
        "prospective_memory",
        "numeric_memory_maxdigits",
        "symbol_digit_accuracy",
        "pm_cols_max",
        "pm_rows_max",
        "pm_correct_sum",
        "symbol_digit_correct",
        "puzzles_solved",
    ],
    "cognition_error_pc1": [
        "trail1_errors",
        "trail2_errors",
        "pm_incorrect_sum",
    ],
}


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class Config:
    master: Path
    cognition: Path
    metabolic: Path
    outdir: Path
    label_map: Optional[Path]
    min_total_n: int
    min_exposed: int
    min_unexposed: int
    min_event_support: int
    fdr_alpha_stage1: float
    fdr_alpha_stage2: float
    p_alpha_path: float
    max_markers_per_exposure_layer: int
    max_mediation_paths: int
    med_boot: int
    support_horizon: int
    cox_penalizer: float
    seed: int
    resume_from_stage_tables: bool
    force_rerun_screens: bool
    stop_after_paths: bool
    no_figures: bool

    @property
    def table_dir(self) -> Path:
        return self.outdir / "tables"

    @property
    def figure_dir(self) -> Path:
        return self.outdir / "figures"

    @property
    def qc_dir(self) -> Path:
        return self.outdir / "qc"


def parse_args(argv: Optional[Sequence[str]] = None) -> Config:
    parser = argparse.ArgumentParser(
        description=(
            "Marker-level GI -> metabolic marker -> cognition bridge analysis. "
            "Generates path-level source tables for Sankey/alluvial plots and "
            "supporting cognition-to-AD/dementia analyses."
        )
    )
    parser.add_argument("--master", required=True, type=Path, help="Master preprocessed participant table.")
    parser.add_argument("--cognition", required=True, type=Path, help="Prepared cognition phenotype table.")
    parser.add_argument("--metabolic", required=True, type=Path, help="Prepared metabolic marker table.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    parser.add_argument("--label-map", default=None, type=Path, help="Optional CSV with columns feature,label.")
    parser.add_argument("--min-total-n", type=int, default=800, help="Minimum complete-case N for each regression.")
    parser.add_argument("--min-exposed", type=int, default=15, help="Minimum exposed N for GI models.")
    parser.add_argument("--min-unexposed", type=int, default=15, help="Minimum unexposed N for GI models.")
    parser.add_argument("--min-event-support", type=int, default=20, help="Minimum AD/dementia events for support models.")
    parser.add_argument("--fdr-alpha-stage1", type=float, default=0.05, help="FDR threshold for GI -> marker screening.")
    parser.add_argument("--fdr-alpha-stage2", type=float, default=0.05, help="FDR threshold for marker -> cognition screening.")
    parser.add_argument("--p-alpha-path", type=float, default=0.05, help="Nominal P threshold for exploratory path support.")
    parser.add_argument(
        "--max-markers-per-exposure-layer",
        type=int,
        default=300,
        help="Safety cap for top stage-1 markers passed to stage 2 per exposure/layer/adjustment. Set <=0 to disable.",
    )
    parser.add_argument(
        "--max-mediation-paths",
        type=int,
        default=300,
        help="Maximum selected paths for bootstrap g-computation mediation. Set 0 to skip.",
    )
    parser.add_argument("--med-boot", type=int, default=200, help="Bootstrap iterations for selected mediation paths.")
    parser.add_argument("--support-horizon", type=int, default=10, help="AD horizon in years for support logistic analysis.")
    parser.add_argument("--cox-penalizer", type=float, default=0.01, help="L2 penalizer used in Cox models.")
    parser.add_argument("--seed", type=int, default=20260510, help="Random seed.")
    parser.add_argument(
        "--resume-from-stage-tables",
        action="store_true",
        help="Reuse existing stage1_gi_to_marker_all.csv and stage2_marker_to_cognition_all.csv if present.",
    )
    parser.add_argument(
        "--force-rerun-screens",
        action="store_true",
        help="Force rerunning stage 1 and stage 2 even if cached stage tables exist.",
    )
    parser.add_argument(
        "--stop-after-paths",
        action="store_true",
        help="Stop after creating path and Sankey source tables.",
    )
    parser.add_argument("--no-figures", action="store_true", help="Skip generating figures.")

    a = parser.parse_args(argv)
    return Config(
        master=a.master,
        cognition=a.cognition,
        metabolic=a.metabolic,
        outdir=a.outdir,
        label_map=a.label_map,
        min_total_n=a.min_total_n,
        min_exposed=a.min_exposed,
        min_unexposed=a.min_unexposed,
        min_event_support=a.min_event_support,
        fdr_alpha_stage1=a.fdr_alpha_stage1,
        fdr_alpha_stage2=a.fdr_alpha_stage2,
        p_alpha_path=a.p_alpha_path,
        max_markers_per_exposure_layer=a.max_markers_per_exposure_layer,
        max_mediation_paths=a.max_mediation_paths,
        med_boot=a.med_boot,
        support_horizon=a.support_horizon,
        cox_penalizer=a.cox_penalizer,
        seed=a.seed,
        resume_from_stage_tables=a.resume_from_stage_tables,
        force_rerun_screens=a.force_rerun_screens,
        stop_after_paths=a.stop_after_paths,
        no_figures=a.no_figures,
    )


# =============================================================================
# Generic utilities
# =============================================================================

def setup_outputs(cfg: Config) -> None:
    for d in [cfg.outdir, cfg.table_dir, cfg.figure_dir, cfg.qc_dir]:
        d.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(cfg.outdir / "run.log", mode="w", encoding="utf-8"),
        ],
    )


def stable_int_seed(*parts: object, modulo: int = 100000) -> int:
    text = "|".join(str(p) for p in parts)
    digest = hashlib.md5(text.encode("utf-8")).hexdigest()
    return int(digest[:12], 16) % modulo


def safe_numeric(x: object) -> pd.Series:
    return pd.to_numeric(x, errors="coerce")


def ensure_binary(x: object) -> pd.Series:
    z = safe_numeric(x)
    return z.where(z.isin([0, 1]), np.nan)


def winsorize(x: object, q: Tuple[float, float] = (0.005, 0.995)) -> pd.Series:
    z = safe_numeric(x).copy()
    lo, hi = z.quantile(q[0]), z.quantile(q[1])
    if pd.notna(lo) and pd.notna(hi):
        z[z < lo] = lo
        z[z > hi] = hi
    return z


def zscore(x: object) -> pd.Series:
    z = safe_numeric(x)
    mu = z.mean(skipna=True)
    sd = z.std(skipna=True)
    if pd.isna(sd) or sd == 0:
        return pd.Series(np.nan, index=z.index)
    return (z - mu) / sd


def fdr_bh(pvals: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvals), dtype=float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    idx = np.where(ok)[0]
    pv = p[idx]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(ranked)
    q_ranked = ranked * m / (np.arange(m) + 1)
    q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
    q_ranked = np.clip(q_ranked, 0, 1)
    q[idx[order]] = q_ranked
    return q


def add_fdr_by_group(
    df: pd.DataFrame,
    p_col: str = "p_value",
    group_cols: Optional[List[str]] = None,
    fdr_col: str = "fdr",
) -> pd.DataFrame:
    if df is None or len(df) == 0:
        return pd.DataFrame() if df is None else df
    out = df.copy()
    out[fdr_col] = np.nan
    if group_cols is None:
        out[fdr_col] = fdr_bh(out[p_col].values)
    else:
        for _, idx in out.groupby(group_cols, dropna=False).groups.items():
            out.loc[idx, fdr_col] = fdr_bh(out.loc[idx, p_col].values)
    return out


def clean_model_frame(df: pd.DataFrame, protected: Optional[Sequence[str]] = None) -> pd.DataFrame:
    protected_set = set(protected or [])
    d = df.copy().replace([np.inf, -np.inf], np.nan)
    for c in list(d.columns):
        if c in protected_set:
            continue
        if d[c].isna().all():
            d = d.drop(columns=[c])
            continue
        if d[c].dropna().nunique() <= 1:
            d = d.drop(columns=[c])
    for c in d.columns:
        if c in protected_set:
            continue
        if d[c].isna().any():
            med = d[c].median(skipna=True)
            d[c] = d[c].fillna(0.0 if pd.isna(med) else med)
    return d


def save_table(df: Optional[pd.DataFrame], path: Path) -> None:
    out = pd.DataFrame() if df is None else df
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    logging.info("Saved %s rows=%s", path, len(out))


def infer_effect_direction(x: float) -> str:
    if not np.isfinite(x):
        return "NA"
    if x > 0:
        return "Positive"
    if x < 0:
        return "Negative"
    return "Zero"


def short_label(x: object, width: int = 45) -> str:
    s = str(x).replace("_", " ")
    return s if len(s) <= width else s[: width - 3] + "..."


# =============================================================================
# Exposure and feature helpers
# =============================================================================

def duration_col_from_exposure(exp: str) -> str:
    return exp.replace("exp_", "dur_") + "_to_baseline"


def date_col_candidates_from_exposure(exp: str) -> List[str]:
    suffix = exp.replace("exp_", "")
    return [f"date_{suffix}", f"dt_{suffix}", exp.replace("exp_", "date_"), exp.replace("exp_", "dt_")]


def get_date_col(df: pd.DataFrame, exp: str) -> Optional[str]:
    for c in date_col_candidates_from_exposure(exp):
        if c in df.columns:
            return c
    return None


def max_binary_across(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    cols = [c for c in cols if c in df.columns]
    if not cols:
        return pd.Series(np.nan, index=df.index)
    tmp = pd.concat([safe_numeric(df[c]) for c in cols], axis=1)
    out = tmp.max(axis=1, skipna=True)
    out[tmp.isna().all(axis=1)] = np.nan
    return out


def min_numeric_across(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    cols = [c for c in cols if c in df.columns]
    if not cols:
        return pd.Series(np.nan, index=df.index)
    tmp = pd.concat([safe_numeric(df[c]) for c in cols], axis=1)
    out = tmp.min(axis=1, skipna=True)
    out[tmp.isna().all(axis=1)] = np.nan
    return out


def min_date_across(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    cols = [c for c in cols if c in df.columns]
    if not cols:
        return pd.Series(pd.NaT, index=df.index)
    tmp = pd.concat([pd.to_datetime(df[c], errors="coerce") for c in cols], axis=1)
    return tmp.min(axis=1)


def add_grouped_exposures(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    for grouped, comps in GROUPED_GI_MAP.items():
        if grouped not in d.columns:
            d[grouped] = max_binary_across(d, comps)
        dur = duration_col_from_exposure(grouped)
        if dur not in d.columns:
            d[dur] = min_numeric_across(d, [duration_col_from_exposure(c) for c in comps])
        dcol = f"date_{grouped.replace('exp_', '')}"
        if dcol not in d.columns:
            date_cols = [get_date_col(d, comp) for comp in comps]
            d[dcol] = min_date_across(d, [c for c in date_cols if c])
    return d


def apply_temporality_filter(df: pd.DataFrame, exposure: str, date_col: str = "cognition_date") -> pd.DataFrame:
    """Keep unexposed and exposed participants with GI date before/equal outcome date when GI date exists."""
    d = df.copy()
    if exposure not in d.columns or date_col not in d.columns:
        return d
    exp = ensure_binary(d[exposure])
    dc = get_date_col(d, exposure)
    if dc is None:
        return d
    exp_dt = pd.to_datetime(d[dc], errors="coerce")
    out_dt = pd.to_datetime(d[date_col], errors="coerce")
    keep = (exp == 0) | ((exp == 1) & exp_dt.notna() & out_dt.notna() & (exp_dt <= out_dt))
    return d.loc[keep.fillna(False)].copy()


def classify_metabolic_feature(feature: str) -> str:
    s = str(feature)
    if s.startswith("p30"):
        return "blood_biochemistry"
    if s.startswith("p23"):
        return "nmr_metabolomics"
    return "other"


def feature_label_map(path: Optional[Path]) -> Dict[str, str]:
    if path and path.exists():
        df = pd.read_csv(path)
        df.columns = [c.lower() for c in df.columns]
        if "feature" in df.columns and "label" in df.columns:
            return dict(zip(df["feature"].astype(str), df["label"].astype(str)))
    return {}


def label_feature(feature: str, fmap: Dict[str, str]) -> str:
    return fmap.get(str(feature), str(feature))


# =============================================================================
# Data preparation
# =============================================================================

def load_table(path: Path, name: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"{name} file does not exist: {path}")
    logging.info("Loading %s: %s", name, path)
    return pd.read_csv(path, low_memory=False)


def parse_dates_inplace(df: pd.DataFrame) -> None:
    date_like = {"time0", "ad_date", "dementia_date", "cognition_date", "censor_date"}
    for c in df.columns:
        if c in date_like or c.startswith("date_") or c.startswith("dt_"):
            df[c] = pd.to_datetime(df[c], errors="coerce")


def prepare_covariates(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    alias_pairs = {
        "age_z": "age_at_baseline",
        "bmi_z": "bmi_i0",
        "whr_z": "whr_i0",
        "townsend_z": "townsend_i0",
        "smoke_pack_z": "smoke_pack_i0",
        "alcohol_freq_z": "alcohol_freq",
        "sleep_hours_z": "sleep_hours",
        "sedentary_z": "total_sedentary_hours",
        "activity_z": "activity_score_days",
        "diet_quality_z": "diet_quality_proxy",
        "parental_dementia_count_z": "parental_dementia_count",
        "apoe_e4_count_z": "apoe_e4_count",
    }
    for newc, oldc in alias_pairs.items():
        if newc not in d.columns and oldc in d.columns:
            d[newc] = zscore(winsorize(d[oldc]))
    if "ad_prs_std_z" not in d.columns and "ad_prs_std" in d.columns:
        d["ad_prs_std_z"] = zscore(winsorize(d["ad_prs_std"]))
    if "ad_prs_enh_z" not in d.columns and "ad_prs_enhanced_std" in d.columns:
        d["ad_prs_enh_z"] = zscore(winsorize(d["ad_prs_enhanced_std"]))

    categorical_aliases = [
        ("sex_model", "sex"),
        ("ethnic_model", "ethnicity_5cat"),
        ("edu_model", "edu_level"),
        ("health_model", "health_4cat"),
        ("smoke_status_model", "smoke_status_3cat"),
        ("apoe_e4_carrier_model", "apoe4_carrier"),
    ]
    for newc, oldc in categorical_aliases:
        if newc not in d.columns and oldc in d.columns:
            d[newc] = safe_numeric(d[oldc])

    binary_cols = [
        "hx_diabetes", "hx_hypertension", "hx_hyperlipidemia", "hx_chd_cvd", "hx_stroke_tia",
        "hx_depression_anxiety", "hx_sleep_disorder", "hx_parkinson_other_nd",
        "med_statin", "med_antihypertensive", "med_antidiabetic", "med_ppi_gi_drug",
        "med_laxative", "med_antidepressant", "med_steroid_immunosuppressive",
        "parental_dementia_any", "father_dementia", "mother_dementia",
        "bowel_cancer_screening", "any_current_medication", "apoe_e4_carrier_model",
        "sleep_short", "sleep_long", "insomnia_any",
    ]
    for c in binary_cols:
        if c in d.columns:
            d[c] = ensure_binary(d[c]).fillna(0)
    return d


def build_covariate_blocks(analysis: pd.DataFrame) -> Dict[str, List[str]]:
    base = [c for c in ["age_z", "sex_model", "ethnic_model", "edu_model", "townsend_z"] if c in analysis.columns]
    lifestyle = [c for c in [
        "bmi_z", "whr_z", "smoke_status_model", "smoke_pack_z",
        "alcohol_freq_z", "sleep_hours_z", "sleep_short", "sleep_long", "insomnia_any",
        "sedentary_z", "activity_z", "diet_quality_z",
    ] if c in analysis.columns]
    comorbidity = [c for c in [
        "hx_diabetes", "hx_hypertension", "hx_hyperlipidemia", "hx_chd_cvd",
        "hx_stroke_tia", "hx_depression_anxiety", "hx_sleep_disorder", "hx_parkinson_other_nd",
    ] if c in analysis.columns]
    medication = [c for c in [
        "med_statin", "med_antihypertensive", "med_antidiabetic", "med_ppi_gi_drug",
        "med_laxative", "med_antidepressant", "med_steroid_immunosuppressive", "any_current_medication",
    ] if c in analysis.columns]
    family_genetic = [c for c in [
        "parental_dementia_any", "parental_dementia_count_z",
        "apoe_e4_carrier_model", "apoe_e4_count_z", "ad_prs_std_z", "ad_prs_enh_z",
    ] if c in analysis.columns]
    detection = [c for c in ["bowel_cancer_screening"] if c in analysis.columns]

    return {
        "base": list(dict.fromkeys(base)),
        "lifestyle": list(dict.fromkeys(base + lifestyle)),
        "full": list(dict.fromkeys(base + lifestyle + comorbidity + medication + family_genetic)),
        "detection": list(dict.fromkeys(base + lifestyle + comorbidity + medication + family_genetic + detection)),
    }


def standardize_covariates(analysis: pd.DataFrame, covar_blocks: Dict[str, List[str]]) -> pd.DataFrame:
    d = analysis.copy()
    all_covars = list(dict.fromkeys(sum(covar_blocks.values(), [])))
    for c in all_covars:
        x = safe_numeric(d[c])
        vals = set(pd.unique(x.dropna()))
        if len(vals) <= 2 and vals.issubset({0, 1}):
            d[c] = x.fillna(0)
        else:
            d[c] = zscore(winsorize(x)).fillna(0)
    return d


def load_and_prepare_data(cfg: Config) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, List[str]]]:
    master = load_table(cfg.master, "master")
    cognition = load_table(cfg.cognition, "cognition")
    metabolic = load_table(cfg.metabolic, "metabolic")

    for df in [master, cognition]:
        parse_dates_inplace(df)

    master = add_grouped_exposures(master)
    for df in [master, cognition, metabolic]:
        if "eid" in df.columns:
            df["eid"] = df["eid"].astype(str).str.strip()

    if "eid" not in master.columns or "eid" not in cognition.columns or "eid" not in metabolic.columns:
        raise ValueError("master, cognition, and metabolic tables must all contain 'eid'.")

    cog2 = cognition.drop(columns=["time0"], errors="ignore").copy()
    metab2 = metabolic.drop(columns=["time0"], errors="ignore").copy()
    analysis = master.merge(cog2, on="eid", how="inner", suffixes=("", "_cog"))
    analysis = analysis.merge(metab2, on="eid", how="inner", suffixes=("", "_metab"))

    if "time0" not in analysis.columns:
        raise ValueError("The merged analysis table requires 'time0' from master.")
    analysis["time0"] = pd.to_datetime(analysis["time0"], errors="coerce")
    if "cognition_date" not in analysis.columns:
        analysis["cognition_date"] = analysis["time0"]
    analysis["cognition_date"] = pd.to_datetime(analysis["cognition_date"], errors="coerce")
    analysis = analysis.loc[
        analysis["time0"].notna()
        & analysis["cognition_date"].notna()
        & (analysis["cognition_date"] >= analysis["time0"])
    ].copy()

    for c in ["ad_event", "dementia_event"]:
        if c in analysis.columns:
            analysis[c] = ensure_binary(analysis[c])
    for c in ["ad_time_years", "dementia_time_years"]:
        if c in analysis.columns:
            analysis[c] = safe_numeric(analysis[c])

    if "ad_date" in analysis.columns:
        analysis["ad_date"] = pd.to_datetime(analysis["ad_date"], errors="coerce")
        ad_event = analysis.get("ad_event", pd.Series(0, index=analysis.index))
        mask_ok = (ad_event != 1) | analysis["ad_date"].isna() | (analysis["cognition_date"] < analysis["ad_date"])
        analysis = analysis.loc[mask_ok].copy()

    analysis = prepare_covariates(analysis)
    covar_blocks = build_covariate_blocks(analysis)
    analysis = standardize_covariates(analysis, covar_blocks)
    return analysis, cognition, metabolic, covar_blocks


# =============================================================================
# Statistical engines
# =============================================================================

def stage1_gi_to_marker(
    df: pd.DataFrame,
    exposure: str,
    features: Sequence[str],
    covars: Sequence[str],
    layer: str,
    adjustment: str,
    fmap: Dict[str, str],
    cfg: Config,
) -> pd.DataFrame:
    rows = []
    covars = [c for c in covars if c in df.columns]
    for feat in features:
        if feat not in df.columns:
            continue
        sub = df[[exposure, feat] + covars].copy()
        sub[exposure] = safe_numeric(sub[exposure])
        sub[feat] = safe_numeric(sub[feat])
        sub = sub.loc[sub[exposure].notna() & sub[feat].notna()].copy()
        if len(sub) < cfg.min_total_n:
            continue
        n_exp = int((sub[exposure] == 1).sum())
        n_unexp = int((sub[exposure] == 0).sum())
        if n_exp < cfg.min_exposed or n_unexp < cfg.min_unexposed:
            continue
        for c in covars:
            sub[c] = safe_numeric(sub[c]).fillna(0)
        x0 = clean_model_frame(sub[[exposure] + covars], protected=[])
        if exposure not in x0.columns or x0[exposure].nunique() <= 1:
            continue
        x = sm.add_constant(x0, has_constant="add")
        y = sub.loc[x.index, feat].values
        try:
            fit = sm.OLS(y, x).fit(cov_type="HC3")
            rows.append({
                "analysis": "GI_to_marker",
                "adjustment": adjustment,
                "exposure": exposure,
                "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                "exposure_abbr": EXPOSURE_ABBR.get(exposure, exposure),
                "marker": feat,
                "marker_label": label_feature(feat, fmap),
                "marker_layer": layer,
                "n": int(len(x)),
                "n_exposed": n_exp,
                "n_unexposed": n_unexp,
                "alpha": float(fit.params[exposure]),
                "alpha_se": float(fit.bse[exposure]),
                "alpha_ci_low": float(fit.params[exposure] - 1.96 * fit.bse[exposure]),
                "alpha_ci_high": float(fit.params[exposure] + 1.96 * fit.bse[exposure]),
                "alpha_p": float(fit.pvalues[exposure]),
                "alpha_direction": infer_effect_direction(float(fit.params[exposure])),
            })
        except Exception as exc:
            logging.debug("Stage1 failed exposure=%s marker=%s: %s", exposure, feat, exc)
            continue
    out = pd.DataFrame(rows)
    if len(out) > 0:
        out = add_fdr_by_group(
            out,
            p_col="alpha_p",
            group_cols=["adjustment", "exposure", "marker_layer"],
            fdr_col="alpha_fdr",
        )
    return out


def stage2_marker_to_cognition(
    df: pd.DataFrame,
    exposure: str,
    markers: Sequence[str],
    outcomes: Sequence[str],
    covars: Sequence[str],
    layer: str,
    cognition_layer: str,
    adjustment: str,
    fmap: Dict[str, str],
    cfg: Config,
) -> pd.DataFrame:
    rows = []
    covars = [c for c in covars if c in df.columns]
    for marker in markers:
        if marker not in df.columns:
            continue
        for outcome in outcomes:
            if outcome not in df.columns:
                continue
            need = list(dict.fromkeys([outcome, marker, exposure] + covars))
            sub = df[need].copy()
            sub[outcome] = safe_numeric(sub[outcome])
            sub[marker] = safe_numeric(sub[marker])
            sub[exposure] = safe_numeric(sub[exposure])
            sub = sub.loc[sub[outcome].notna() & sub[marker].notna() & sub[exposure].notna()].copy()
            if len(sub) < cfg.min_total_n:
                continue
            for c in covars:
                sub[c] = safe_numeric(sub[c]).fillna(0)
            x0 = clean_model_frame(sub[[marker, exposure] + covars], protected=[])
            if marker not in x0.columns:
                continue
            x = sm.add_constant(x0, has_constant="add")
            y = sub.loc[x.index, outcome].values
            try:
                fit = sm.OLS(y, x).fit(cov_type="HC3")
                rows.append({
                    "analysis": "marker_to_cognition",
                    "adjustment": adjustment,
                    "exposure": exposure,
                    "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                    "exposure_abbr": EXPOSURE_ABBR.get(exposure, exposure),
                    "marker": marker,
                    "marker_label": label_feature(marker, fmap),
                    "marker_layer": layer,
                    "cognition": outcome,
                    "cognition_label": outcome,
                    "cognition_layer": cognition_layer,
                    "n": int(len(x)),
                    "beta": float(fit.params[marker]),
                    "beta_se": float(fit.bse[marker]),
                    "beta_ci_low": float(fit.params[marker] - 1.96 * fit.bse[marker]),
                    "beta_ci_high": float(fit.params[marker] + 1.96 * fit.bse[marker]),
                    "beta_p": float(fit.pvalues[marker]),
                    "beta_direction": infer_effect_direction(float(fit.params[marker])),
                })
            except Exception as exc:
                logging.debug("Stage2 failed marker=%s outcome=%s: %s", marker, outcome, exc)
                continue
    out = pd.DataFrame(rows)
    if len(out) > 0:
        out = add_fdr_by_group(
            out,
            p_col="beta_p",
            group_cols=["adjustment", "exposure", "marker_layer", "cognition_layer"],
            fdr_col="beta_fdr",
        )
    return out


def build_path_table(stage1_all: pd.DataFrame, stage2_all: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    if len(stage1_all) == 0 or len(stage2_all) == 0:
        return pd.DataFrame()
    s1_cols = [
        "adjustment", "exposure", "exposure_label", "exposure_abbr", "marker", "marker_label", "marker_layer",
        "n", "n_exposed", "n_unexposed", "alpha", "alpha_se", "alpha_p", "alpha_fdr", "alpha_direction",
    ]
    s2_cols = [
        "adjustment", "exposure", "marker", "marker_layer", "cognition", "cognition_label", "cognition_layer",
        "n", "beta", "beta_se", "beta_p", "beta_fdr", "beta_direction",
    ]
    s1 = stage1_all[[c for c in s1_cols if c in stage1_all.columns]].copy()
    s2 = stage2_all[[c for c in s2_cols if c in stage2_all.columns]].copy()
    paths = s1.merge(s2, on=["adjustment", "exposure", "marker", "marker_layer"], how="inner", suffixes=("_stage1", "_stage2"))
    if len(paths) == 0:
        return paths

    paths["indirect_effect"] = paths["alpha"] * paths["beta"]
    paths["abs_indirect_effect"] = paths["indirect_effect"].abs()
    paths["sobel_se"] = np.sqrt((paths["beta"] ** 2) * (paths["alpha_se"] ** 2) + (paths["alpha"] ** 2) * (paths["beta_se"] ** 2))
    paths["sobel_z"] = paths["indirect_effect"] / paths["sobel_se"].replace(0, np.nan)
    z_abs = pd.to_numeric(paths["sobel_z"], errors="coerce").abs().to_numpy(dtype=float)
    paths["sobel_p"] = np.array([erfc(z / np.sqrt(2.0)) if np.isfinite(z) else np.nan for z in z_abs], dtype=float)
    paths["sobel_fdr"] = fdr_bh(paths["sobel_p"].values)
    paths["path_direction"] = paths["indirect_effect"].apply(infer_effect_direction)
    paths["stage1_fdr_sig"] = paths["alpha_fdr"] < cfg.fdr_alpha_stage1
    paths["stage2_fdr_sig"] = paths["beta_fdr"] < cfg.fdr_alpha_stage2
    paths["stage1_nominal"] = paths["alpha_p"] < cfg.p_alpha_path
    paths["stage2_nominal"] = paths["beta_p"] < cfg.p_alpha_path
    paths["path_fdr_sig"] = paths["stage1_fdr_sig"] & paths["stage2_fdr_sig"]
    paths["path_nominal"] = paths["stage1_nominal"] & paths["stage2_nominal"]
    paths["path_sobel_sig"] = paths["sobel_p"] < cfg.p_alpha_path
    paths["path_support_level"] = np.where(
        paths["path_fdr_sig"] & (paths["sobel_fdr"] < 0.05),
        "stage1+stage2 FDR and Sobel FDR",
        np.where(
            paths["path_fdr_sig"],
            "stage1+stage2 FDR",
            np.where(paths["path_nominal"], "stage1+stage2 nominal", "exploratory"),
        ),
    )
    paths["sankey_weight"] = paths["abs_indirect_effect"].replace(0, np.nan)
    paths["sankey_weight"] = paths["sankey_weight"].fillna(paths["sankey_weight"].median(skipna=True)).fillna(1e-6)
    return paths.sort_values(["path_fdr_sig", "path_nominal", "abs_indirect_effect"], ascending=[False, False, False])


def mediation_continuous_outcome(
    df: pd.DataFrame,
    exposure: str,
    mediator: str,
    outcome: str,
    covars: Sequence[str],
    adjustment: str,
    cfg: Config,
    n_boot: int,
    seed: int,
) -> Optional[Dict[str, float]]:
    covars = [c for c in covars if c in df.columns]
    need = list(dict.fromkeys([exposure, mediator, outcome] + covars))
    dd = df[need].copy()
    if exposure not in dd.columns or mediator not in dd.columns or outcome not in dd.columns:
        return None
    dd[exposure] = safe_numeric(dd[exposure])
    dd[mediator] = safe_numeric(dd[mediator])
    dd[outcome] = safe_numeric(dd[outcome])
    dd = dd.loc[dd[exposure].notna() & dd[mediator].notna() & dd[outcome].notna()].copy()
    if len(dd) < cfg.min_total_n:
        return None
    if int((dd[exposure] == 1).sum()) < cfg.min_exposed or int((dd[exposure] == 0).sum()) < cfg.min_unexposed:
        return None
    for c in covars:
        dd[c] = safe_numeric(dd[c]).fillna(0)

    def one_pass(dat: pd.DataFrame) -> Tuple[float, float, float]:
        xm0 = clean_model_frame(dat[[exposure] + covars].copy(), protected=[])
        if exposure not in xm0.columns:
            raise ValueError("exposure dropped from mediator model")
        xm = sm.add_constant(xm0, has_constant="add")
        fit_m = sm.OLS(dat.loc[xm.index, mediator].values, xm).fit()

        dat2 = dat.loc[xm.index].copy()
        xy0 = clean_model_frame(dat2[[exposure, mediator] + covars].copy(), protected=[])
        if exposure not in xy0.columns or mediator not in xy0.columns:
            raise ValueError("exposure or mediator dropped from outcome model")
        xy = sm.add_constant(xy0, has_constant="add")
        fit_y = sm.OLS(dat2.loc[xy.index, outcome].values, xy).fit()

        x1m = dat2.loc[xy.index, [exposure] + [c for c in covars if c in xm0.columns]].copy()
        x0m = x1m.copy()
        x1m[exposure] = 1
        x0m[exposure] = 0
        x1m = sm.add_constant(x1m, has_constant="add").reindex(columns=xm.columns, fill_value=0)
        x0m = sm.add_constant(x0m, has_constant="add").reindex(columns=xm.columns, fill_value=0)
        m1 = fit_m.predict(x1m)
        m0 = fit_m.predict(x0m)

        def pred(xval: int, mvec: np.ndarray) -> np.ndarray:
            tmp = dat2.loc[xy.index, [exposure, mediator] + [c for c in covars if c in xy0.columns]].copy()
            tmp[exposure] = xval
            tmp[mediator] = mvec
            tmp = sm.add_constant(tmp, has_constant="add").reindex(columns=xy.columns, fill_value=0)
            return fit_y.predict(tmp)

        y11 = pred(1, m1)
        y10 = pred(1, m0)
        y00 = pred(0, m0)
        return float(np.mean(y11 - y00)), float(np.mean(y11 - y10)), float(np.mean(y10 - y00))

    try:
        te, nie, nde = one_pass(dd)
    except Exception:
        return None

    rng = np.random.default_rng(seed)
    boot = []
    n = len(dd)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            boot.append(one_pass(dd.iloc[idx].reset_index(drop=True)))
        except Exception:
            pass
    if len(boot) < max(30, int(0.25 * n_boot)):
        return None
    arr = np.asarray(boot)
    return {
        "adjustment": adjustment,
        "n": int(len(dd)),
        "te": te,
        "te_ci_low": float(np.quantile(arr[:, 0], 0.025)),
        "te_ci_high": float(np.quantile(arr[:, 0], 0.975)),
        "nie": nie,
        "nie_ci_low": float(np.quantile(arr[:, 1], 0.025)),
        "nie_ci_high": float(np.quantile(arr[:, 1], 0.975)),
        "nde": nde,
        "nde_ci_low": float(np.quantile(arr[:, 2], 0.025)),
        "nde_ci_high": float(np.quantile(arr[:, 2], 0.975)),
        "prop_mediated": float(nie / te) if np.isfinite(te) and abs(te) > 1e-12 else np.nan,
    }


def build_horizon_binary_outcome(df: pd.DataFrame, time_col: str, event_col: str, horizon: int) -> Tuple[pd.Series, pd.Series]:
    timev = safe_numeric(df[time_col])
    event = ensure_binary(df[event_col])
    event_by_h = (event == 1) & (timev <= horizon)
    observed_past_h = timev >= horizon
    eligible = event_by_h | observed_past_h
    y = pd.Series(np.nan, index=df.index, dtype=float)
    y.loc[eligible] = 0
    y.loc[event_by_h] = 1
    return eligible, y


def support_cognition_to_event_cox(
    df: pd.DataFrame,
    features: Sequence[str],
    covars: Sequence[str],
    event_col: str,
    time_col: str,
    outcome_name: str,
    cfg: Config,
) -> pd.DataFrame:
    rows = []
    if event_col not in df.columns or time_col not in df.columns:
        return pd.DataFrame()
    covars = [c for c in covars if c in df.columns]
    for feat in features:
        if feat not in df.columns:
            continue
        need = list(dict.fromkeys([time_col, event_col, feat] + covars))
        sub = df[need].copy()
        sub[time_col] = safe_numeric(sub[time_col])
        sub[event_col] = ensure_binary(sub[event_col])
        sub[feat] = safe_numeric(sub[feat])
        sub = sub.loc[sub[time_col].notna() & (sub[time_col] > 0) & sub[event_col].notna() & sub[feat].notna()].copy()
        if len(sub) < cfg.min_total_n or sub[event_col].sum() < cfg.min_event_support:
            continue
        for c in covars:
            sub[c] = safe_numeric(sub[c]).fillna(0)
        model_df = sub.rename(columns={time_col: "time_years", event_col: "event"}).copy()
        model_df = clean_model_frame(model_df, protected=["time_years", "event"])
        if feat not in model_df.columns:
            continue
        try:
            cph = CoxPHFitter(penalizer=cfg.cox_penalizer)
            cph.fit(model_df, duration_col="time_years", event_col="event")
            s = cph.summary.loc[feat]
            rows.append({
                "analysis": "support_cognition_to_event_cox",
                "outcome": outcome_name,
                "feature": feat,
                "n": int(len(model_df)),
                "n_event": int(model_df["event"].sum()),
                "coef": float(s["coef"]),
                "hr": float(np.exp(s["coef"])),
                "ci_low": float(np.exp(s["coef lower 95%"])),
                "ci_high": float(np.exp(s["coef upper 95%"])),
                "p_value": float(s["p"]),
            })
        except Exception as exc:
            logging.debug("Support Cox failed feature=%s outcome=%s: %s", feat, outcome_name, exc)
            continue
    out = pd.DataFrame(rows)
    if len(out) > 0:
        out["fdr"] = fdr_bh(out["p_value"].values)
    return out


def support_cognition_to_ad_horizon_logit(
    df: pd.DataFrame,
    features: Sequence[str],
    covars: Sequence[str],
    cfg: Config,
) -> pd.DataFrame:
    rows = []
    if "ad_event" not in df.columns or "ad_time_years" not in df.columns:
        return pd.DataFrame()
    eligible, y = build_horizon_binary_outcome(df, "ad_time_years", "ad_event", cfg.support_horizon)
    base = df.loc[eligible].copy()
    y_col = f"ad_{cfg.support_horizon}y"
    base[y_col] = y.loc[base.index]
    base = base.loc[base[y_col].notna()].copy()
    base[y_col] = base[y_col].astype(int)
    if base[y_col].sum() < cfg.min_event_support:
        return pd.DataFrame()
    covars = [c for c in covars if c in base.columns]
    for feat in features:
        if feat not in base.columns:
            continue
        need = list(dict.fromkeys([y_col, feat] + covars))
        sub = base[need].copy()
        sub[feat] = safe_numeric(sub[feat])
        sub = sub.loc[sub[feat].notna()].copy()
        if len(sub) < cfg.min_total_n or sub[y_col].sum() < cfg.min_event_support:
            continue
        for c in covars:
            sub[c] = safe_numeric(sub[c]).fillna(0)
        x = clean_model_frame(sub[[feat] + covars], protected=[])
        if feat not in x.columns:
            continue
        yv = sub.loc[x.index, y_col].astype(int).values
        try:
            model = LogisticRegression(max_iter=800, solver="lbfgs")
            model.fit(x.values, yv)
            coef = float(model.coef_[0][list(x.columns).index(feat)])
            rows.append({
                "analysis": "support_cognition_to_ad_horizon_logit",
                "horizon_years": cfg.support_horizon,
                "feature": feat,
                "n": int(len(x)),
                "n_event": int(np.sum(yv)),
                "or": float(np.exp(coef)),
                "coef": coef,
            })
        except Exception as exc:
            logging.debug("Support logit failed feature=%s: %s", feat, exc)
            continue
    return pd.DataFrame(rows)


# =============================================================================
# Feature construction
# =============================================================================

def define_metabolic_features(analysis: pd.DataFrame, metabolic: pd.DataFrame) -> Tuple[List[str], List[str], List[str], List[Tuple[str, List[str]]]]:
    meta_features = [c for c in metabolic.columns if c not in {"eid", "time0"} and c in analysis.columns]
    for c in meta_features:
        analysis[c] = zscore(winsorize(analysis[c]))
    blood_features = [c for c in meta_features if classify_metabolic_feature(c) == "blood_biochemistry"]
    nmr_features = [c for c in meta_features if classify_metabolic_feature(c) == "nmr_metabolomics"]
    other_features = [c for c in meta_features if classify_metabolic_feature(c) == "other"]
    layer_sets = [("blood_biochemistry", blood_features), ("nmr_metabolomics", nmr_features)]
    if other_features:
        layer_sets.append(("other_metabolic", other_features))
    return meta_features, blood_features, nmr_features, layer_sets


def define_cognition_features(analysis: pd.DataFrame, cognition: pd.DataFrame, seed: int) -> Tuple[List[str], List[str], List[str], List[str], List[Tuple[str, List[str]]]]:
    cog_meta_cols = {"eid", "time0", "cognition_date", "cognition_instance_used"}
    cognition_cols = [c for c in cognition.columns if c not in cog_meta_cols and c in analysis.columns]
    cognition_cols = [c for c in cognition_cols if safe_numeric(analysis[c]).notna().sum() > 0]
    for c in cognition_cols:
        analysis[c] = zscore(winsorize(analysis[c]))

    added_composites: List[str] = []
    for comp_name, feats in COGNITION_DOMAIN_CANDIDATES.items():
        present = [f for f in feats if f in analysis.columns]
        if len(present) >= 2:
            x = analysis[present].copy()
            for c in present:
                x[c] = safe_numeric(x[c])
                med = x[c].median(skipna=True)
                x[c] = x[c].fillna(0.0 if pd.isna(med) else med)
            try:
                pca = PCA(n_components=1, random_state=seed)
                analysis[comp_name] = zscore(pd.Series(pca.fit_transform(x.values).ravel(), index=analysis.index))
                if comp_name not in cognition_cols:
                    cognition_cols.append(comp_name)
                added_composites.append(comp_name)
            except Exception as exc:
                logging.warning("Failed to construct cognition composite %s: %s", comp_name, exc)

    primary_scores = [c for c in ["cognition_speed_pc1", "cognition_accuracy_pc1", "cognition_error_pc1"] if c in analysis.columns]
    raw_features = [c for c in cognition_cols if c not in primary_scores]
    cognition_layers = [("primary_cognition_scores", primary_scores), ("raw_cognition_features", raw_features)]
    return cognition_cols, primary_scores, raw_features, added_composites, cognition_layers


# =============================================================================
# Main analysis workflow
# =============================================================================

def run_stage_screens(
    analysis: pd.DataFrame,
    layer_sets: Sequence[Tuple[str, List[str]]],
    cognition_layers: Sequence[Tuple[str, List[str]]],
    covar_blocks: Dict[str, List[str]],
    fmap: Dict[str, str],
    cfg: Config,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    stage1_cache = cfg.table_dir / "stage1_gi_to_marker_all.csv"
    stage2_cache = cfg.table_dir / "stage2_marker_to_cognition_all.csv"

    if cfg.resume_from_stage_tables and not cfg.force_rerun_screens and stage1_cache.exists() and stage2_cache.exists():
        logging.info("Loading cached stage tables.")
        stage1_all = pd.read_csv(stage1_cache, low_memory=False)
        stage2_all = pd.read_csv(stage2_cache, low_memory=False)
        summary_df = reconstruct_stage1_summary(stage1_all, cfg)
        return stage1_all, stage2_all, summary_df

    if cfg.resume_from_stage_tables:
        logging.warning("Cached stage tables requested but missing; rerunning stage screens.")

    main_exposures = [x for x in MAIN_EXPOSURES if x in analysis.columns]
    secondary_exposures = [x for x in SECONDARY_EXPOSURES if x in analysis.columns]
    adjustments = ["base", "full", "detection"]
    all_stage1: List[pd.DataFrame] = []
    all_stage2: List[pd.DataFrame] = []
    summary_rows: List[Dict[str, object]] = []

    for exposure_set, exposures in [("main", main_exposures), ("secondary", secondary_exposures)]:
        for exposure in exposures:
            d_exp = apply_temporality_filter(analysis, exposure, date_col="cognition_date")
            exp_vals = safe_numeric(d_exp[exposure])
            n_exp = int((exp_vals == 1).sum())
            n_unexp = int((exp_vals == 0).sum())
            if n_exp < cfg.min_exposed or n_unexp < cfg.min_unexposed:
                logging.info("Skipping %s: exposed=%s, unexposed=%s", exposure, n_exp, n_unexp)
                continue
            logging.info("Exposure %s/%s: exposed=%s, unexposed=%s, n=%s", exposure_set, exposure, n_exp, n_unexp, len(d_exp))

            for marker_layer, marker_features in layer_sets:
                if not marker_features:
                    continue
                logging.info("Marker layer %s: %s markers", marker_layer, len(marker_features))
                stage1_by_adj: Dict[str, pd.DataFrame] = {}

                for adj in adjustments:
                    s1 = stage1_gi_to_marker(d_exp, exposure, marker_features, covar_blocks[adj], marker_layer, adj, fmap, cfg)
                    if len(s1) > 0:
                        s1["exposure_set"] = exposure_set
                        all_stage1.append(s1)
                        save_table(s1, cfg.table_dir / f"stage1_gi_to_marker_{adj}_{exposure}_{marker_layer}.csv")
                    stage1_by_adj[adj] = s1

                for adj in adjustments:
                    s1 = stage1_by_adj.get(adj, pd.DataFrame())
                    if len(s1) == 0:
                        continue
                    marker_pool = s1.loc[s1["alpha_p"].notna()].sort_values("alpha_p")
                    markers = marker_pool.loc[marker_pool["alpha_fdr"] < cfg.fdr_alpha_stage1, "marker"].dropna().astype(str).tolist()
                    if not markers:
                        markers = marker_pool.loc[marker_pool["alpha_p"] < cfg.p_alpha_path, "marker"].dropna().astype(str).tolist()
                    if not markers:
                        continue
                    if cfg.max_markers_per_exposure_layer and cfg.max_markers_per_exposure_layer > 0:
                        markers = markers[: cfg.max_markers_per_exposure_layer]

                    for cog_layer, outcomes in cognition_layers:
                        if not outcomes:
                            continue
                        s2 = stage2_marker_to_cognition(
                            d_exp,
                            exposure,
                            markers,
                            outcomes,
                            covar_blocks[adj],
                            marker_layer,
                            cog_layer,
                            adj,
                            fmap,
                            cfg,
                        )
                        if len(s2) > 0:
                            s2["exposure_set"] = exposure_set
                            all_stage2.append(s2)
                            save_table(s2, cfg.table_dir / f"stage2_marker_to_cognition_{adj}_{exposure}_{marker_layer}_{cog_layer}.csv")

                s1_base = stage1_by_adj.get("base", pd.DataFrame())
                summary_rows.append({
                    "exposure_set": exposure_set,
                    "exposure": exposure,
                    "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                    "marker_layer": marker_layer,
                    "n_stage1_base_tested": int(len(s1_base)) if len(s1_base) else 0,
                    "n_stage1_base_fdr_sig": int((s1_base["alpha_fdr"] < cfg.fdr_alpha_stage1).sum()) if len(s1_base) else 0,
                    "n_stage1_base_nominal": int((s1_base["alpha_p"] < cfg.p_alpha_path).sum()) if len(s1_base) else 0,
                })

    stage1_all = pd.concat(all_stage1, ignore_index=True) if all_stage1 else pd.DataFrame()
    stage2_all = pd.concat(all_stage2, ignore_index=True) if all_stage2 else pd.DataFrame()
    save_table(stage1_all, stage1_cache)
    save_table(stage2_all, stage2_cache)
    return stage1_all, stage2_all, pd.DataFrame(summary_rows)


def reconstruct_stage1_summary(stage1_all: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    if len(stage1_all) == 0:
        return pd.DataFrame()
    try:
        base_s1 = stage1_all.loc[stage1_all["adjustment"].astype(str) == "base"].copy()
        if len(base_s1) == 0:
            return pd.DataFrame()
        return base_s1.groupby(["exposure_set", "exposure", "exposure_label", "marker_layer"], dropna=False).agg(
            n_stage1_base_tested=("marker", "size"),
            n_stage1_base_fdr_sig=("alpha_fdr", lambda x: int((pd.to_numeric(x, errors="coerce") < cfg.fdr_alpha_stage1).sum())),
            n_stage1_base_nominal=("alpha_p", lambda x: int((pd.to_numeric(x, errors="coerce") < cfg.p_alpha_path).sum())),
        ).reset_index()
    except Exception as exc:
        logging.warning("Failed to reconstruct stage 1 summary from cached table: %s", exc)
        return pd.DataFrame()


def select_display_paths(path_all: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    if len(path_all) == 0:
        return pd.DataFrame()
    path_sig = path_all.loc[
        (path_all["adjustment"] == "base")
        & ((path_all["path_fdr_sig"] == True) | ((path_all["path_nominal"] == True) & (path_all["path_sobel_sig"] == True)))
    ].copy()
    if len(path_sig) == 0:
        path_sig = path_all.loc[path_all["adjustment"] == "base"].sort_values("abs_indirect_effect", ascending=False).head(300).copy()
        path_sig["path_support_level"] = "top exploratory by abs indirect effect"
    return path_sig


def build_sankey_source(path_sig: pd.DataFrame) -> pd.DataFrame:
    if len(path_sig) == 0:
        return pd.DataFrame()
    sankey = path_sig.sort_values("abs_indirect_effect", ascending=False).copy()
    balanced = []
    for (_, _), g in sankey.groupby(["exposure", "cognition_layer"], dropna=False):
        balanced.append(g.head(25))
    sankey = pd.concat(balanced, ignore_index=True) if balanced else pd.DataFrame()
    sankey = sankey.sort_values("abs_indirect_effect", ascending=False).head(250)
    sankey = sankey.rename(columns={
        "exposure_abbr": "gi",
        "cognition": "cognition_outcome",
        "indirect_effect": "bridge_effect",
        "abs_indirect_effect": "abs_bridge_effect",
        "sobel_p": "p_indirect",
        "sobel_fdr": "fdr_indirect",
    })
    keep_cols = [
        "gi", "exposure", "exposure_label", "marker", "marker_label", "marker_layer",
        "cognition_outcome", "cognition_label", "cognition_layer",
        "alpha", "alpha_p", "alpha_fdr", "beta", "beta_p", "beta_fdr",
        "bridge_effect", "abs_bridge_effect", "p_indirect", "fdr_indirect",
        "path_support_level", "path_direction", "sankey_weight",
    ]
    return sankey[[c for c in keep_cols if c in sankey.columns]]


def summarize_paths(path_all: pd.DataFrame, summary_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if len(path_all) > 0:
        path_summary = path_all.groupby(["adjustment", "exposure", "exposure_label", "marker_layer", "cognition_layer"], dropna=False).agg(
            n_paths_tested=("marker", "size"),
            n_unique_markers=("marker", "nunique"),
            n_unique_cognition=("cognition", "nunique"),
            n_stage1_stage2_fdr_paths=("path_fdr_sig", "sum"),
            n_nominal_paths=("path_nominal", "sum"),
            sum_abs_indirect_effect=("abs_indirect_effect", "sum"),
        ).reset_index()
    else:
        path_summary = pd.DataFrame()

    summary_out = summary_df.copy()
    if len(path_summary) > 0 and len(summary_out) > 0:
        base_path_sum = path_summary.loc[path_summary["adjustment"] == "base"].copy()
        summary_out = summary_out.merge(base_path_sum, on=["exposure", "exposure_label", "marker_layer"], how="left")
    return path_summary, summary_out


def run_bootstrap_mediation(
    analysis: pd.DataFrame,
    path_sig: pd.DataFrame,
    covar_blocks: Dict[str, List[str]],
    fmap: Dict[str, str],
    cfg: Config,
) -> pd.DataFrame:
    if not cfg.max_mediation_paths or cfg.max_mediation_paths <= 0 or len(path_sig) == 0:
        return pd.DataFrame()
    med_candidates = path_sig.loc[path_sig["adjustment"] == "base"].sort_values("abs_indirect_effect", ascending=False).head(cfg.max_mediation_paths)
    logging.info("Running bootstrap g-computation for %s selected marker-level paths.", len(med_candidates))
    rows = []
    for _, r in med_candidates.reset_index(drop=True).iterrows():
        exposure = str(r["exposure"])
        marker = str(r["marker"])
        outcome = str(r["cognition"])
        d_exp = apply_temporality_filter(analysis, exposure, date_col="cognition_date")
        for adj in ["base", "full"]:
            seed = cfg.seed + stable_int_seed(exposure, marker, outcome, adj)
            res = mediation_continuous_outcome(
                d_exp,
                exposure,
                marker,
                outcome,
                covar_blocks[adj],
                adjustment=adj,
                cfg=cfg,
                n_boot=cfg.med_boot,
                seed=seed,
            )
            if res is None:
                continue
            res.update({
                "analysis": "GI_to_marker_to_cognition_gcomp",
                "exposure": exposure,
                "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                "exposure_abbr": EXPOSURE_ABBR.get(exposure, exposure),
                "mediator": marker,
                "mediator_label": label_feature(marker, fmap),
                "marker_layer": r.get("marker_layer"),
                "outcome": outcome,
                "cognition_layer": r.get("cognition_layer"),
                "screen_alpha": r.get("alpha"),
                "screen_beta": r.get("beta"),
                "screen_product": r.get("indirect_effect"),
                "screen_sobel_p": r.get("sobel_p"),
            })
            rows.append(res)
    return pd.DataFrame(rows)


def draw_hr_forest(df: pd.DataFrame, title: str, outfile: Path, label_col: str = "feature", top_n: int = 30) -> None:
    if df is None or len(df) == 0:
        return
    dd = df.loc[np.isfinite(df["hr"]) & np.isfinite(df["ci_low"]) & np.isfinite(df["ci_high"])].copy()
    if len(dd) == 0:
        return
    dd = dd.sort_values("p_value").head(top_n).sort_values("hr")
    y = np.arange(len(dd))
    plt.figure(figsize=(9, max(4.5, 0.36 * len(dd) + 1.5)))
    plt.errorbar(dd["hr"], y, xerr=[dd["hr"] - dd["ci_low"], dd["ci_high"] - dd["hr"]], fmt="o", capsize=3)
    plt.axvline(1.0, linestyle="--", linewidth=1)
    plt.xscale("log")
    plt.yticks(y, dd[label_col].tolist())
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()


def draw_summary_figures(path_summary: pd.DataFrame, path_sig: pd.DataFrame, cfg: Config) -> None:
    if cfg.no_figures:
        return
    try:
        if len(path_summary) > 0:
            for val in ["n_stage1_stage2_fdr_paths", "n_nominal_paths", "n_unique_markers", "n_unique_cognition"]:
                dd = path_summary.loc[path_summary["adjustment"] == "base"].copy()
                if val not in dd.columns or len(dd) == 0:
                    continue
                dd["xlab"] = dd["exposure_label"].astype(str) + "\n" + dd["marker_layer"].astype(str) + "\n" + dd["cognition_layer"].astype(str)
                dd = dd.sort_values(val, ascending=False).head(40)
                plt.figure(figsize=(max(10, 0.32 * len(dd)), 5.2))
                plt.bar(np.arange(len(dd)), dd[val].fillna(0).values)
                plt.xticks(np.arange(len(dd)), dd["xlab"], rotation=70, ha="right", fontsize=7)
                plt.ylabel(val)
                plt.title(f"Marker-level path summary: {val}")
                plt.tight_layout()
                plt.savefig(cfg.figure_dir / f"summary_{val}.png", dpi=300, bbox_inches="tight")
                plt.close()
        if len(path_sig) > 0:
            dd = path_sig.sort_values("abs_indirect_effect", ascending=False).head(30).copy()
            dd["label"] = (
                dd["exposure_abbr"].astype(str)
                + " -> "
                + dd["marker_label"].astype(str).map(lambda x: short_label(x, 22))
                + " -> "
                + dd["cognition"].astype(str).map(lambda x: short_label(x, 22))
            )
            dd = dd.sort_values("indirect_effect")
            y = np.arange(len(dd))
            plt.figure(figsize=(12, max(5, 0.36 * len(dd) + 1.5)))
            plt.scatter(dd["indirect_effect"], y, s=45, alpha=0.85)
            plt.axvline(0, linestyle="--", linewidth=1)
            plt.yticks(y, dd["label"], fontsize=7)
            plt.xlabel("alpha × beta product")
            plt.title("Top marker-level GI -> metabolic marker -> cognition paths")
            plt.tight_layout()
            plt.savefig(cfg.figure_dir / "top_marker_cognition_paths_product_forest.png", dpi=300, bbox_inches="tight")
            plt.close()
    except Exception as exc:
        logging.warning("Summary plotting failed: %s", exc)


def write_qc_tables(
    cfg: Config,
    analysis: pd.DataFrame,
    meta_features: List[str],
    blood_features: List[str],
    nmr_features: List[str],
    other_meta_features: List[str],
    cognition_cols: List[str],
    primary_scores: List[str],
    added_composites: List[str],
    main_exposures: List[str],
    secondary_exposures: List[str],
    covar_blocks: Dict[str, List[str]],
    raw_cognition_features: List[str],
) -> None:
    save_table(pd.DataFrame([{
        "n_analysis": int(len(analysis)),
        "n_ad_events_support": int(analysis["ad_event"].sum()) if "ad_event" in analysis.columns else np.nan,
        "n_dementia_events_support": int(analysis["dementia_event"].sum()) if "dementia_event" in analysis.columns else np.nan,
        "n_metabolic_features": len(meta_features),
        "blood_features": len(blood_features),
        "nmr_features": len(nmr_features),
        "other_metabolic_features": len(other_meta_features),
        "n_cognition_features": len(cognition_cols),
        "n_primary_cognition_scores": len(primary_scores),
        "added_composites": ";".join(added_composites),
        "main_exposures": ";".join(main_exposures),
        "secondary_exposures": ";".join(secondary_exposures),
    }]), cfg.table_dir / "marker_cognition_cohort_summary.csv")

    save_table(pd.DataFrame([
        {"adjustment": k, "n_covariates": len(v), "covariates": ";".join(v)} for k, v in covar_blocks.items()
    ]), cfg.qc_dir / "marker_cognition_covariate_blocks.csv")

    save_table(pd.DataFrame([
        {"cognition_layer": "primary_scores", "n_features": len(primary_scores), "features": ";".join(primary_scores)},
        {"cognition_layer": "raw_features", "n_features": len(raw_cognition_features), "features": ";".join(raw_cognition_features)},
    ]), cfg.qc_dir / "marker_cognition_feature_qc.csv")


def run_analysis(cfg: Config) -> None:
    setup_outputs(cfg)
    np.random.seed(cfg.seed)

    logging.info("Metabolic-cognition bridge analysis")
    logging.info("Master: %s", cfg.master)
    logging.info("Cognition: %s", cfg.cognition)
    logging.info("Metabolic: %s", cfg.metabolic)
    logging.info("Outdir: %s", cfg.outdir)

    analysis, cognition, metabolic, covar_blocks = load_and_prepare_data(cfg)
    fmap = feature_label_map(cfg.label_map)

    meta_features, blood_features, nmr_features, layer_sets = define_metabolic_features(analysis, metabolic)
    other_meta_features = [c for c in meta_features if classify_metabolic_feature(c) == "other"]
    cognition_cols, primary_scores, raw_cognition_features, added_composites, cognition_layers = define_cognition_features(analysis, cognition, cfg.seed)

    main_exposures = [x for x in MAIN_EXPOSURES if x in analysis.columns]
    secondary_exposures = [x for x in SECONDARY_EXPOSURES if x in analysis.columns]

    write_qc_tables(
        cfg,
        analysis,
        meta_features,
        blood_features,
        nmr_features,
        other_meta_features,
        cognition_cols,
        primary_scores,
        added_composites,
        main_exposures,
        secondary_exposures,
        covar_blocks,
        raw_cognition_features,
    )

    stage1_all, stage2_all, stage1_summary = run_stage_screens(analysis, layer_sets, cognition_layers, covar_blocks, fmap, cfg)

    path_all = build_path_table(stage1_all, stage2_all, cfg)
    save_table(path_all, cfg.table_dir / "marker_cognition_bridge_paths_all.csv")

    path_sig = select_display_paths(path_all, cfg)
    save_table(path_sig, cfg.table_dir / "marker_cognition_bridge_paths_significant.csv")

    sankey_src = build_sankey_source(path_sig)
    save_table(sankey_src, cfg.table_dir / "sankey_source_paths.csv")

    path_summary, summary_df = summarize_paths(path_all, stage1_summary)
    save_table(path_summary, cfg.table_dir / "marker_cognition_path_summary_by_exposure.csv")
    save_table(summary_df, cfg.table_dir / "marker_cognition_summary_table.csv")

    if cfg.stop_after_paths:
        write_metadata_and_report(
            cfg,
            analysis,
            meta_features,
            blood_features,
            nmr_features,
            other_meta_features,
            cognition_cols,
            primary_scores,
            added_composites,
            main_exposures,
            secondary_exposures,
            covar_blocks,
            stage1_all,
            stage2_all,
            path_all,
            path_sig,
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            summary_df,
        )
        logging.info("Stop-after-paths requested. Core path tables are complete.")
        return

    mediation_all = run_bootstrap_mediation(analysis, path_sig, covar_blocks, fmap, cfg)
    save_table(mediation_all, cfg.table_dir / "gi_marker_cognition_mediation_all.csv")

    support_features = list(dict.fromkeys(primary_scores + raw_cognition_features))
    support_ad = support_cognition_to_event_cox(analysis, support_features, covar_blocks["full"], "ad_event", "ad_time_years", "strict_ad", cfg)
    support_dem = support_cognition_to_event_cox(analysis, support_features, covar_blocks["full"], "dementia_event", "dementia_time_years", "all_cause_dementia", cfg)
    support_logit = support_cognition_to_ad_horizon_logit(analysis, support_features, covar_blocks["full"], cfg)

    save_table(support_ad, cfg.table_dir / "support_cognition_to_ad_cox.csv")
    save_table(support_dem, cfg.table_dir / "support_cognition_to_dementia_cox.csv")
    save_table(support_logit, cfg.table_dir / f"support_cognition_to_ad_logit_{cfg.support_horizon}y.csv")

    if not cfg.no_figures:
        if len(support_ad) > 0:
            draw_hr_forest(support_ad, "Support only: cognition -> AD", cfg.figure_dir / "support_cognition_to_ad_hr_forest.png")
        if len(support_dem) > 0:
            draw_hr_forest(support_dem, "Support only: cognition -> all-cause dementia", cfg.figure_dir / "support_cognition_to_dementia_hr_forest.png")
        draw_summary_figures(path_summary, path_sig, cfg)

    write_metadata_and_report(
        cfg,
        analysis,
        meta_features,
        blood_features,
        nmr_features,
        other_meta_features,
        cognition_cols,
        primary_scores,
        added_composites,
        main_exposures,
        secondary_exposures,
        covar_blocks,
        stage1_all,
        stage2_all,
        path_all,
        path_sig,
        mediation_all,
        support_ad,
        support_dem,
        support_logit,
        summary_df,
    )

    logging.info("Analysis complete. Read first: %s", cfg.outdir / "analysis_report.txt")


def write_metadata_and_report(
    cfg: Config,
    analysis: pd.DataFrame,
    meta_features: List[str],
    blood_features: List[str],
    nmr_features: List[str],
    other_meta_features: List[str],
    cognition_cols: List[str],
    primary_scores: List[str],
    added_composites: List[str],
    main_exposures: List[str],
    secondary_exposures: List[str],
    covar_blocks: Dict[str, List[str]],
    stage1_all: pd.DataFrame,
    stage2_all: pd.DataFrame,
    path_all: pd.DataFrame,
    path_sig: pd.DataFrame,
    mediation_all: pd.DataFrame,
    support_ad: pd.DataFrame,
    support_dem: pd.DataFrame,
    support_logit: pd.DataFrame,
    summary_df: pd.DataFrame,
) -> None:
    meta = {
        "script": "metabolic_cognition_bridge_analysis.py",
        "master_path": str(cfg.master),
        "cognition_path": str(cfg.cognition),
        "metabolic_path": str(cfg.metabolic),
        "outdir": str(cfg.outdir),
        "design": {
            "main_question": "Which GI disease connects to which metabolic marker and which cognition phenotype?",
            "stage1": "GI -> metabolic marker using robust OLS",
            "stage2": "metabolic marker -> cognition phenotype using robust OLS adjusted for GI exposure",
            "path_effect": "alpha * beta product with Sobel screening statistic",
            "optional_mediation": "continuous-outcome g-computation for selected top paths",
            "support_only": "cognition -> AD/dementia Cox/logistic support",
        },
        "thresholds": {
            "stage1_fdr": cfg.fdr_alpha_stage1,
            "stage2_fdr": cfg.fdr_alpha_stage2,
            "path_nominal_p": cfg.p_alpha_path,
            "min_total_n": cfg.min_total_n,
            "min_exposed": cfg.min_exposed,
            "min_unexposed": cfg.min_unexposed,
            "min_event_support": cfg.min_event_support,
        },
        "n_analysis": int(len(analysis)),
        "n_metabolic_features": len(meta_features),
        "blood_features": len(blood_features),
        "nmr_features": len(nmr_features),
        "other_metabolic_features": len(other_meta_features),
        "n_cognition_features": len(cognition_cols),
        "primary_cognition_scores": primary_scores,
        "added_composites": added_composites,
        "main_exposures": main_exposures,
        "secondary_exposures": secondary_exposures,
        "covariate_blocks": covar_blocks,
        "key_standardization_changes": [
            "Removed personal hard-coded default paths; all inputs are command-line arguments.",
            "Renamed outputs to generic GitHub-friendly names without phase-specific prefixes.",
            "Wrapped workflow into parse_args/run_analysis/main for import safety.",
            "Replaced Python hash-based bootstrap seeds with stable md5-based seeds.",
            "Kept temporality filtering: GI exposure date must be before/equal cognition date when available.",
        ],
    }
    with open(cfg.outdir / "run_metadata.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)

    report = [
        "Metabolic-cognition bridge analysis finished.",
        "",
        "Main target:",
        "  GI disease -> metabolic marker -> cognition phenotype",
        "",
        f"Analysis cohort n = {len(analysis)}",
        f"Metabolic markers = {len(meta_features)}; blood = {len(blood_features)}; NMR = {len(nmr_features)}; other = {len(other_meta_features)}",
        f"Cognition features = {len(cognition_cols)}",
        f"Main exposures = {main_exposures}",
        f"Secondary exposures = {secondary_exposures}",
        "",
        "Core output tables:",
        "  tables/stage1_gi_to_marker_all.csv",
        "  tables/stage2_marker_to_cognition_all.csv",
        "  tables/marker_cognition_bridge_paths_all.csv",
        "  tables/marker_cognition_bridge_paths_significant.csv",
        "  tables/sankey_source_paths.csv",
        "  tables/gi_marker_cognition_mediation_all.csv",
        "",
        "Interpretation:",
        "  The path-level table gives exact GI | marker | cognition triplets for Sankey plotting.",
        "  indirect_effect = alpha(GI->marker) * beta(marker->cognition | GI,covariates).",
        "  This is a marker-level mechanistic screening layer, not proof of causal mediation.",
        "",
        "Summary:",
        summary_df.to_string(index=False) if len(summary_df) > 0 else "No summary rows produced.",
        "",
        f"Stage1 rows = {len(stage1_all)}",
        f"Stage2 rows = {len(stage2_all)}",
        f"All path rows = {len(path_all)}",
        f"Significant/display path rows = {len(path_sig)}",
        f"Bootstrap mediation rows = {len(mediation_all)}",
        f"Support AD Cox rows = {len(support_ad)}",
        f"Support dementia Cox rows = {len(support_dem)}",
        f"Support AD logistic rows = {len(support_logit)}",
    ]
    with open(cfg.outdir / "analysis_report.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(report))


# =============================================================================
# Entrypoint
# =============================================================================

def main(argv: Optional[Sequence[str]] = None) -> None:
    cfg = parse_args(argv)
    run_analysis(cfg)


if __name__ == "__main__":
    main()
