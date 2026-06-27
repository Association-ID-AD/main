#!/usr/bin/env python3
"""
Metabolic-neural bridge analysis for gastrointestinal disease and Alzheimer-related outcomes.

This script identifies marker-level bridge paths of the form:

    GI disease -> metabolic marker -> brain imaging-derived phenotype (IDP)

It also performs supportive screens for:

    GI disease -> brain IDP
    brain IDP -> incident AD / all-cause dementia

The script is designed for repository use. It avoids hard-coded personal paths and
uses command-line arguments for all inputs and outputs.

Example
-------
python metabolic_neural_bridge_analysis.py \
    --master data/master_preprocessed_with_genetics.csv \
    --brain data/brain_idp_prepared.csv \
    --metabolic data/metabolic_prepared.csv \
    --metabolic-bridge-dir results/metabolic_bridge_analysis \
    --outdir results/metabolic_neural_bridge_analysis

Optional IDP panel file:
python metabolic_neural_bridge_analysis.py \
    --master data/master_preprocessed_with_genetics.csv \
    --brain data/brain_idp_prepared.csv \
    --metabolic data/metabolic_prepared.csv \
    --metabolic-bridge-dir results/metabolic_bridge_analysis \
    --outdir results/metabolic_neural_bridge_analysis \
    --panel-file data/selected_idps.txt
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import time
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm
from lifelines import CoxPHFitter
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression

warnings.filterwarnings("default")


# =============================================================================
# Constants
# =============================================================================

RANDOM_SEED = 20260519

MAIN_EXPOSURES = [
    "exp_diverticular",
    "exp_other_functional_intestinal",
    "exp_ibd",
    "exp_other_chronic_intestinal",
    "exp_ibs",
]
SECONDARY_EXPOSURES = ["exp_appendiceal", "exp_malabsorption"]
ALL_EXPOSURES = MAIN_EXPOSURES + SECONDARY_EXPOSURES

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
    "exp_anorectal": "Anorectal disease",
}

EXPOSURE_ABBR = {
    "exp_diverticular": "DIV",
    "exp_other_functional_intestinal": "OFID",
    "exp_other_chronic_intestinal": "OCID",
    "exp_ibd": "IBD",
    "exp_ibs": "IBS",
    "exp_appendiceal": "APD",
    "exp_malabsorption": "MAL",
    "exp_anorectal": "ARD",
}
ABBR_TO_EXPOSURE = {v: k for k, v in EXPOSURE_ABBR.items()}

BRAIN_ID_COLUMNS = {"eid", "time0", "imaging_date", "repeat_imaging_date"}
METABOLIC_ID_COLUMNS = {"eid", "time0", "metabolic_date", "assessment_date"}


# =============================================================================
# Configuration
# =============================================================================


@dataclass
class AnalysisConfig:
    master: Path
    brain: Path
    metabolic: Path
    outdir: Path
    metabolic_bridge_dir: Optional[Path] = None
    panel_file: Optional[Path] = None
    min_total_n: int = 1200
    min_exposed: int = 20
    min_unexposed: int = 20
    min_event_support: int = 20
    fdr_alpha: float = 0.05
    marker_to_idp_p_screen: float = 0.01
    marker_to_idp_fdr_screen: float = 0.10
    max_markers_per_exposure_total: int = 50
    max_markers_per_layer: int = 30
    max_idps_per_exposure: int = 60
    max_idps_per_exposure_fallback: int = 40
    max_total_bridge_tests_per_exposure: int = 3000
    support_horizon: int = 10
    cox_penalizer: float = 0.01
    min_feature_nonmissing: int = 100
    seed: int = RANDOM_SEED

    @property
    def tables_dir(self) -> Path:
        return self.outdir / "tables"

    @property
    def qc_dir(self) -> Path:
        return self.outdir / "qc"

    @property
    def figures_dir(self) -> Path:
        return self.outdir / "figures"


# =============================================================================
# Basic utilities
# =============================================================================


def setup_output_dirs(config: AnalysisConfig) -> None:
    for directory in [config.outdir, config.tables_dir, config.qc_dir, config.figures_dir]:
        directory.mkdir(parents=True, exist_ok=True)


def setup_logging(outdir: Path) -> None:
    log_path = outdir / "run.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_path, mode="w", encoding="utf-8")],
    )


def safe_numeric(x: pd.Series | np.ndarray | object) -> pd.Series:
    return pd.to_numeric(x, errors="coerce")


def ensure_binary(x: pd.Series | np.ndarray | object) -> pd.Series:
    z = safe_numeric(x)
    return z.where(z.isin([0, 1]), np.nan)


def winsorize(x: pd.Series, q: Tuple[float, float] = (0.005, 0.995)) -> pd.Series:
    z = safe_numeric(x).copy()
    lo, hi = z.quantile(q[0]), z.quantile(q[1])
    if pd.notna(lo) and pd.notna(hi):
        z[z < lo] = lo
        z[z > hi] = hi
    return z


def zscore(x: pd.Series) -> pd.Series:
    z = safe_numeric(x)
    mu = z.mean(skipna=True)
    sd = z.std(skipna=True)
    if pd.isna(sd) or sd == 0:
        return pd.Series(np.nan, index=z.index)
    return (z - mu) / sd


def fdr_bh(pvals: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
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


def save_table(df: pd.DataFrame, name: str, directory: Path) -> None:
    if isinstance(df, pd.DataFrame):
        df.to_csv(directory / name, index=False)


def read_csv_if_exists(path: Optional[Path]) -> pd.DataFrame:
    if path is not None and path.exists():
        logging.info("Reading %s", path)
        return pd.read_csv(path, low_memory=False)
    logging.warning("Missing optional file: %s", path)
    return pd.DataFrame()


def stable_int_seed(*parts: object, modulo: int = 1_000_000) -> int:
    text = "::".join(str(p) for p in parts)
    digest = hashlib.md5(text.encode("utf-8")).hexdigest()
    return int(digest[:12], 16) % modulo


def clean_model_frame(df: pd.DataFrame, protected: Optional[Iterable[str]] = None) -> pd.DataFrame:
    protected_set = set(protected or [])
    d = df.copy().replace([np.inf, -np.inf], np.nan)

    for col in list(d.columns):
        if col in protected_set:
            continue
        if d[col].isna().all():
            d = d.drop(columns=[col])
            continue
        if d[col].dropna().nunique() <= 1:
            d = d.drop(columns=[col])

    for col in d.columns:
        if col in protected_set:
            continue
        if d[col].isna().any():
            med = d[col].median(skipna=True)
            d[col] = d[col].fillna(0.0 if pd.isna(med) else med)
    return d


# =============================================================================
# Exposure and date utilities
# =============================================================================


def duration_col_from_exposure(exposure: str) -> str:
    return exposure.replace("exp_", "dur_") + "_to_baseline"


def date_col_candidates_from_exposure(exposure: str) -> List[str]:
    suffix = exposure.replace("exp_", "")
    return [f"date_{suffix}", f"dt_{suffix}", exposure.replace("exp_", "date_"), exposure.replace("exp_", "dt_")]


def get_date_col(df: pd.DataFrame, exposure: str) -> Optional[str]:
    for col in date_col_candidates_from_exposure(exposure):
        if col in df.columns:
            return col
    return None


def max_binary_across(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    available = [c for c in cols if c in df.columns]
    if not available:
        return pd.Series(np.nan, index=df.index)
    tmp = pd.concat([safe_numeric(df[c]) for c in available], axis=1)
    out = tmp.max(axis=1, skipna=True)
    out[tmp.isna().all(axis=1)] = np.nan
    return out


def min_numeric_across(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    available = [c for c in cols if c in df.columns]
    if not available:
        return pd.Series(np.nan, index=df.index)
    tmp = pd.concat([safe_numeric(df[c]) for c in available], axis=1)
    out = tmp.min(axis=1, skipna=True)
    out[tmp.isna().all(axis=1)] = np.nan
    return out


def min_date_across(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    available = [c for c in cols if c in df.columns]
    if not available:
        return pd.Series(pd.NaT, index=df.index)
    tmp = pd.concat([pd.to_datetime(df[c], errors="coerce") for c in available], axis=1)
    return tmp.min(axis=1)


def add_grouped_exposures(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    for grouped, components in GROUPED_GI_MAP.items():
        if grouped not in d.columns:
            d[grouped] = max_binary_across(d, components)

        duration_col = duration_col_from_exposure(grouped)
        if duration_col not in d.columns:
            d[duration_col] = min_numeric_across(d, [duration_col_from_exposure(c) for c in components])

        date_col = f"date_{grouped.replace('exp_', '')}"
        if date_col not in d.columns:
            component_date_cols = [get_date_col(d, component) for component in components]
            d[date_col] = min_date_across(d, [c for c in component_date_cols if c])
    return d


def apply_temporality_filter(df: pd.DataFrame, exposure: str, imaging_col: str = "imaging_date") -> pd.DataFrame:
    d = df.copy()
    if exposure not in d.columns or imaging_col not in d.columns:
        return d
    date_col = get_date_col(d, exposure)
    if date_col is None:
        return d

    exposure_status = ensure_binary(d[exposure])
    exposure_date = pd.to_datetime(d[date_col], errors="coerce")
    imaging_date = pd.to_datetime(d[imaging_col], errors="coerce")

    keep = (exposure_status == 0) | (
        (exposure_status == 1) & exposure_date.notna() & imaging_date.notna() & (exposure_date < imaging_date)
    )
    keep = keep.fillna(False)
    return d.loc[keep].copy()


def standard_exposure_from_any(value: object) -> object:
    s = str(value)
    if s in ALL_EXPOSURES:
        return s
    if s in ABBR_TO_EXPOSURE:
        return ABBR_TO_EXPOSURE[s]
    low = s.lower()
    if "divert" in low:
        return "exp_diverticular"
    if "functional" in low:
        return "exp_other_functional_intestinal"
    if "chronic" in low:
        return "exp_other_chronic_intestinal"
    if "ibd" in low or "inflammatory bowel" in low:
        return "exp_ibd"
    if "ibs" in low or "irritable" in low:
        return "exp_ibs"
    if "append" in low:
        return "exp_appendiceal"
    if "malabs" in low:
        return "exp_malabsorption"
    return np.nan


# =============================================================================
# Feature utilities
# =============================================================================


def parse_panel_file(path: Optional[Path]) -> List[str]:
    if path is None or not path.exists():
        return []
    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path)
        for col in ["feature", "idp"]:
            if col in df.columns:
                return df[col].dropna().astype(str).tolist()
        return df.iloc[:, 0].dropna().astype(str).tolist()
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def classify_metabolic_marker(feature: object, label: Optional[object] = None, layer: Optional[object] = None) -> str:
    feature_str = str(feature).lower()
    label_str = "" if label is None else str(label).lower()
    layer_str = "" if layer is None else str(layer).lower()
    text = f"{feature_str} {label_str} {layer_str}"

    if any(k in text for k in [
        "crp", "glyca", "glycoprotein", "white blood", "red blood", "platelet",
        "neutrophil", "lymphocyte", "monocyte", "haemoglobin", "hemoglobin",
    ]):
        return "Inflammation / hematology"
    if any(k in text for k in ["glucose", "hba1c", "glycaemic", "glycemic", "insulin"]):
        return "Glucose / glycaemic"
    if any(k in text for k in [
        "albumin", "creatinine", "urea", "cystatin", "alt", "ast", "ggt",
        "bilirubin", "liver", "renal", "kidney",
    ]):
        return "Liver / renal"
    if any(k in text for k in [
        "alanine", "glutamine", "glycine", "histidine", "isoleucine", "leucine",
        "valine", "phenylalanine", "tyrosine", "amino",
    ]):
        return "Amino acids"
    if any(k in text for k in ["fatty acid", "omega", "linoleic", "dha", "pufa", "mufa", "saturated"]):
        return "Fatty acids"
    if any(k in text for k in [
        "cholesterol", "triglycer", "hdl", "ldl", "vldl", "apolipoprotein",
        "lipoprotein", "phospholipid", "sphingomyelin",
    ]):
        return "Lipid / lipoprotein"
    if str(feature).startswith("p23"):
        return "NMR metabolomics"
    if str(feature).startswith("p30"):
        return "Blood biochemistry"
    return "Other metabolites"


def ci_crosses_zero(low: float, high: float) -> bool:
    if not np.isfinite(low) or not np.isfinite(high):
        return False
    return low <= 0 <= high


# =============================================================================
# Outcome and model helpers
# =============================================================================


def build_horizon_binary_outcome(
    df: pd.DataFrame,
    time_col: str,
    event_col: str,
    horizon: int | float,
) -> Tuple[pd.Series, pd.Series]:
    time_value = safe_numeric(df[time_col])
    event = ensure_binary(df[event_col])
    event_by_horizon = (event == 1) & (time_value <= horizon)
    observed_past_horizon = time_value >= horizon
    eligible = event_by_horizon | observed_past_horizon

    y = pd.Series(np.nan, index=df.index, dtype=float)
    y.loc[eligible] = 0
    y.loc[event_by_horizon] = 1
    return eligible, y


def model_ols(
    df: pd.DataFrame,
    outcome: str,
    predictors: Sequence[str],
    covariates: Sequence[str],
    min_total_n: int,
    robust: bool = True,
) -> Optional[Tuple[sm.regression.linear_model.RegressionResultsWrapper, pd.DataFrame]]:
    cols = [outcome] + list(predictors) + list(covariates)
    cols = [c for c in cols if c in df.columns]
    if outcome not in cols or any(p not in cols for p in predictors):
        return None

    d = df[cols].copy().replace([np.inf, -np.inf], np.nan)
    for col in cols:
        d[col] = safe_numeric(d[col])
    d = d.dropna(subset=[outcome] + list(predictors))
    if len(d) < min_total_n:
        return None

    for col in covariates:
        if col in d.columns and d[col].isna().any():
            med = d[col].median(skipna=True)
            d[col] = d[col].fillna(0.0 if pd.isna(med) else med)

    x_cols = list(predictors) + [c for c in covariates if c in d.columns]
    x = clean_model_frame(d[x_cols], protected=predictors)
    if any(p not in x.columns for p in predictors):
        return None

    x = sm.add_constant(x, has_constant="add")
    y = d.loc[x.index, outcome]
    try:
        fit = sm.OLS(y.values, x).fit(cov_type="HC3" if robust else "nonrobust")
        return fit, d.loc[x.index].copy()
    except Exception:
        return None


def extract_ols_coef(fit: sm.regression.linear_model.RegressionResultsWrapper, term: str) -> Optional[Dict[str, float]]:
    if fit is None or term not in fit.params.index:
        return None
    beta = float(fit.params[term])
    se = float(fit.bse[term])
    p = float(fit.pvalues[term])
    return {
        "beta": beta,
        "se": se,
        "ci_low": beta - 1.96 * se,
        "ci_high": beta + 1.96 * se,
        "p_value": p,
    }


# =============================================================================
# Data loading and covariates
# =============================================================================


def convert_date_columns(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    date_like = {"time0", "ad_date", "dementia_date", "imaging_date", "repeat_imaging_date", "censor_date"}
    for col in d.columns:
        if col in date_like or col.startswith("date_") or col.startswith("dt_"):
            d[col] = pd.to_datetime(d[col], errors="coerce")
    return d


def prepare_covariates(df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, List[str]]]:
    d = df.copy()

    continuous_aliases = {
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
    for new_col, old_col in continuous_aliases.items():
        if new_col not in d.columns and old_col in d.columns:
            d[new_col] = zscore(winsorize(d[old_col]))

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
    for new_col, old_col in categorical_aliases:
        if new_col not in d.columns and old_col in d.columns:
            d[new_col] = safe_numeric(d[old_col])

    binary_cols = [
        "hx_diabetes", "hx_hypertension", "hx_hyperlipidemia", "hx_chd_cvd", "hx_stroke_tia",
        "hx_depression_anxiety", "hx_sleep_disorder", "hx_parkinson_other_nd", "med_statin",
        "med_antihypertensive", "med_antidiabetic", "med_ppi_gi_drug", "med_laxative",
        "med_antidepressant", "med_steroid_immunosuppressive", "parental_dementia_any",
        "father_dementia", "mother_dementia", "bowel_cancer_screening", "any_current_medication",
        "apoe_e4_carrier_model", "sleep_short", "sleep_long", "insomnia_any",
    ]
    for col in binary_cols:
        if col in d.columns:
            d[col] = ensure_binary(d[col]).fillna(0)

    base = [c for c in ["age_z", "sex_model", "ethnic_model", "edu_model", "townsend_z"] if c in d.columns]
    lifestyle = [c for c in [
        "bmi_z", "whr_z", "smoke_status_model", "smoke_pack_z", "alcohol_freq_z",
        "sleep_hours_z", "sleep_short", "sleep_long", "insomnia_any", "sedentary_z",
        "activity_z", "diet_quality_z",
    ] if c in d.columns]
    comorbidity = [c for c in [
        "hx_diabetes", "hx_hypertension", "hx_hyperlipidemia", "hx_chd_cvd", "hx_stroke_tia",
        "hx_depression_anxiety", "hx_sleep_disorder", "hx_parkinson_other_nd",
    ] if c in d.columns]
    medication = [c for c in [
        "med_statin", "med_antihypertensive", "med_antidiabetic", "med_ppi_gi_drug",
        "med_laxative", "med_antidepressant", "med_steroid_immunosuppressive", "any_current_medication",
    ] if c in d.columns]
    family_genetic = [c for c in [
        "parental_dementia_any", "parental_dementia_count_z", "apoe_e4_carrier_model",
        "apoe_e4_count_z", "ad_prs_std_z", "ad_prs_enh_z",
    ] if c in d.columns]
    detection = [c for c in ["bowel_cancer_screening"] if c in d.columns]

    covariate_blocks = {
        "base": list(dict.fromkeys(base)),
        "full": list(dict.fromkeys(base + lifestyle + comorbidity + medication + family_genetic)),
        "detection": list(dict.fromkeys(base + lifestyle + comorbidity + medication + family_genetic + detection)),
    }

    for cols in covariate_blocks.values():
        for col in cols:
            x = safe_numeric(d[col])
            vals = set(pd.unique(x.dropna()))
            if len(vals) <= 2 and vals.issubset({0, 1}):
                d[col] = x.fillna(0)
            else:
                d[col] = zscore(winsorize(x)).fillna(0)

    return d, covariate_blocks


def load_and_prepare_data(config: AnalysisConfig) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, List[str]], List[str], List[str], List[str]]:
    logging.info("Loading master, brain, and metabolic tables")
    master = convert_date_columns(pd.read_csv(config.master, low_memory=False))
    brain = convert_date_columns(pd.read_csv(config.brain, low_memory=False))
    metabolic = convert_date_columns(pd.read_csv(config.metabolic, low_memory=False))

    if "eid" not in master.columns or "eid" not in brain.columns or "eid" not in metabolic.columns:
        raise ValueError("master, brain, and metabolic tables must all contain an 'eid' column.")

    master = add_grouped_exposures(master)
    brain2 = brain.drop(columns=["time0"], errors="ignore").copy()
    metabolic2 = metabolic.drop(columns=["time0"], errors="ignore").copy()

    analysis = master.merge(brain2, on="eid", how="inner", suffixes=("", "_brain"))
    analysis = analysis.merge(metabolic2, on="eid", how="left", suffixes=("", "_metabolic"))

    if "time0" not in analysis.columns or "imaging_date" not in analysis.columns:
        raise ValueError("The merged analysis table requires 'time0' and 'imaging_date'.")

    analysis["time0"] = pd.to_datetime(analysis["time0"], errors="coerce")
    analysis["imaging_date"] = pd.to_datetime(analysis["imaging_date"], errors="coerce")
    analysis = analysis.loc[
        analysis["time0"].notna()
        & analysis["imaging_date"].notna()
        & (analysis["imaging_date"] > analysis["time0"])
    ].copy()

    for col in ["ad_event", "dementia_event"]:
        if col in analysis.columns:
            analysis[col] = ensure_binary(analysis[col])
    for col in ["ad_time_years", "dementia_time_years"]:
        if col in analysis.columns:
            analysis[col] = safe_numeric(analysis[col])

    analysis, covariate_blocks = prepare_covariates(analysis)

    brain_features = [c for c in brain.columns if c not in BRAIN_ID_COLUMNS and c in analysis.columns]
    brain_features = [c for c in brain_features if safe_numeric(analysis[c]).notna().sum() >= config.min_feature_nonmissing]

    metabolic_features = [c for c in metabolic.columns if c not in METABOLIC_ID_COLUMNS and c in analysis.columns]
    metabolic_features = [c for c in metabolic_features if safe_numeric(analysis[c]).notna().sum() >= config.min_feature_nonmissing]

    logging.info("Z-scoring %d brain IDPs and %d metabolic markers", len(brain_features), len(metabolic_features))
    for col in brain_features:
        analysis[col] = zscore(winsorize(analysis[col]))
    for col in metabolic_features:
        analysis[col] = zscore(winsorize(analysis[col]))

    available_exposures = [exp for exp in ALL_EXPOSURES if exp in analysis.columns]
    panel_idps = [x for x in parse_panel_file(config.panel_file) if x in brain_features]
    panel_idps = list(dict.fromkeys(panel_idps))

    return analysis, brain, metabolic, covariate_blocks, brain_features, metabolic_features, panel_idps


# =============================================================================
# Phase B1 candidate standardization
# =============================================================================


def find_metabolic_bridge_table(config: AnalysisConfig, name: str) -> Optional[Path]:
    if config.metabolic_bridge_dir is None:
        return None
    candidates = [
        config.metabolic_bridge_dir / "tables" / name,
        config.metabolic_bridge_dir / name,
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


def standardise_metabolic_bridge_candidates(
    candidate_raw: pd.DataFrame,
    bridge_raw: pd.DataFrame,
    metabolic_features: Sequence[str],
    fdr_alpha: float,
) -> pd.DataFrame:
    rows: List[pd.DataFrame] = []

    if len(candidate_raw) > 0:
        c = candidate_raw.copy()
        exp_col = "exposure_abbr" if "exposure_abbr" in c.columns else ("exposure" if "exposure" in c.columns else None)
        feat_col = "feature" if "feature" in c.columns else ("mediator" if "mediator" in c.columns else None)
        if exp_col and feat_col:
            cand = pd.DataFrame()
            cand["exposure"] = c[exp_col].apply(standard_exposure_from_any)
            cand["feature"] = c[feat_col].astype(str)
            cand["feature_label"] = c["feature_label"].astype(str) if "feature_label" in c.columns else cand["feature"]
            cand["layer"] = c["layer"].astype(str) if "layer" in c.columns else cand["feature"].apply(
                lambda z: "nmr_metabolomics" if str(z).startswith("p23") else ("blood_biochemistry" if str(z).startswith("p30") else "metabolic")
            )
            cand["stage1_beta"] = safe_numeric(c["beta"]) if "beta" in c.columns else np.nan
            cand["stage1_p"] = safe_numeric(c["p_value"]) if "p_value" in c.columns else np.nan
            cand["stage1_fdr"] = safe_numeric(c["fdr"]) if "fdr" in c.columns else np.nan
            cand["marker_to_ad_hr"] = safe_numeric(c["hr"]) if "hr" in c.columns else np.nan
            cand["marker_to_ad_fdr"] = safe_numeric(c["fdr_stage2cox"]) if "fdr_stage2cox" in c.columns else np.nan
            if "candidate_sig" in c.columns:
                cand["candidate_sig"] = c["candidate_sig"].astype(str).str.lower().isin(["true", "1", "yes"])
            else:
                cand["candidate_sig"] = (cand["stage1_fdr"] < fdr_alpha) | (cand["stage1_p"] < 0.001)
            rows.append(cand)

    if len(bridge_raw) > 0:
        b = bridge_raw.copy()
        exp_col = "exposure_abbr" if "exposure_abbr" in b.columns else ("exposure" if "exposure" in b.columns else None)
        feat_col = "feature" if "feature" in b.columns else ("mediator" if "mediator" in b.columns else None)
        if exp_col and feat_col:
            bm = pd.DataFrame()
            bm["exposure"] = b[exp_col].apply(standard_exposure_from_any)
            bm["feature"] = b[feat_col].astype(str)
            if "mediator_label" in b.columns:
                bm["feature_label"] = b["mediator_label"].astype(str)
            elif "feature_label" in b.columns:
                bm["feature_label"] = b["feature_label"].astype(str)
            else:
                bm["feature_label"] = bm["feature"]
            bm["layer"] = b["layer"].astype(str) if "layer" in b.columns else bm["feature"].apply(
                lambda z: "nmr_metabolomics" if str(z).startswith("p23") else ("blood_biochemistry" if str(z).startswith("p30") else "metabolic")
            )
            bm["bridge_nie_ad"] = safe_numeric(b["nie"]) if "nie" in b.columns else np.nan
            bm["bridge_nie_ci_low_ad"] = safe_numeric(b["nie_ci_low"]) if "nie_ci_low" in b.columns else np.nan
            bm["bridge_nie_ci_high_ad"] = safe_numeric(b["nie_ci_high"]) if "nie_ci_high" in b.columns else np.nan
            bm["bridge_sig_ad"] = (
                ~bm.apply(lambda r: ci_crosses_zero(r["bridge_nie_ci_low_ad"], r["bridge_nie_ci_high_ad"]), axis=1)
                & bm["bridge_nie_ci_low_ad"].notna()
                & bm["bridge_nie_ci_high_ad"].notna()
            )
            rows.append(bm)

    if not rows:
        return pd.DataFrame()

    all_candidates = pd.concat(rows, ignore_index=True, sort=False).dropna(subset=["exposure", "feature"])
    all_candidates = all_candidates[all_candidates["feature"].isin(metabolic_features)].copy()
    if len(all_candidates) == 0:
        return all_candidates

    def first_nonmissing(series: pd.Series) -> object:
        series = series.dropna()
        return series.iloc[0] if len(series) else np.nan

    agg = all_candidates.groupby(["exposure", "feature"], as_index=False).agg({
        "feature_label": first_nonmissing,
        "layer": first_nonmissing,
        "stage1_beta": "first",
        "stage1_p": "first",
        "stage1_fdr": "first",
        "marker_to_ad_hr": "first",
        "marker_to_ad_fdr": "first",
        "candidate_sig": "max",
        "bridge_nie_ad": "first",
        "bridge_nie_ci_low_ad": "first",
        "bridge_nie_ci_high_ad": "first",
        "bridge_sig_ad": "max",
    })

    agg["marker_category"] = [
        classify_metabolic_marker(f, label, layer)
        for f, label, layer in zip(agg["feature"], agg["feature_label"], agg["layer"])
    ]
    agg["priority_score"] = 0.0
    agg["priority_score"] += np.where(agg["bridge_sig_ad"].fillna(False), 100.0, 0.0)
    agg["priority_score"] += np.where(agg["candidate_sig"].fillna(False), 25.0, 0.0)
    agg["priority_score"] += np.nan_to_num(np.abs(agg["bridge_nie_ad"].astype(float)), nan=0.0) * 200.0
    agg["priority_score"] += np.nan_to_num(-np.log10(np.maximum(agg["stage1_fdr"].astype(float), 1e-300)), nan=0.0)
    return agg


def fallback_gi_to_marker_screen(
    analysis: pd.DataFrame,
    available_exposures: Sequence[str],
    metabolic_features: Sequence[str],
    covariate_blocks: Dict[str, List[str]],
    config: AnalysisConfig,
) -> pd.DataFrame:
    logging.warning("No metabolic bridge candidate table found. Running fallback GI-to-marker screen.")
    rows = []
    for exposure in available_exposures:
        exposure_df = apply_temporality_filter(analysis, exposure, "imaging_date")
        if int((safe_numeric(exposure_df[exposure]) == 1).sum()) < config.min_exposed:
            continue
        for feature in metabolic_features:
            result = model_ols(exposure_df, feature, [exposure], covariate_blocks["base"], config.min_total_n, robust=True)
            if result is None:
                continue
            coef = extract_ols_coef(result[0], exposure)
            if coef is None:
                continue
            rows.append({
                "exposure": exposure,
                "feature": feature,
                "feature_label": feature,
                "layer": "nmr_metabolomics" if feature.startswith("p23") else ("blood_biochemistry" if feature.startswith("p30") else "metabolic"),
                "stage1_beta": coef["beta"],
                "stage1_p": coef["p_value"],
                "stage1_fdr": np.nan,
                "candidate_sig": coef["p_value"] < 0.001,
                "bridge_sig_ad": False,
                "bridge_nie_ad": np.nan,
                "priority_score": -np.log10(max(coef["p_value"], 1e-300)),
            })
    out = pd.DataFrame(rows)
    if len(out) > 0:
        out["marker_category"] = [
            classify_metabolic_marker(f, label, layer)
            for f, label, layer in zip(out["feature"], out["feature_label"], out["layer"])
        ]
    return out


def select_marker_candidates(marker_candidates: pd.DataFrame, available_exposures: Sequence[str], config: AnalysisConfig) -> pd.DataFrame:
    selected = []
    for exposure in available_exposures:
        sub = marker_candidates[marker_candidates["exposure"] == exposure].copy()
        if len(sub) == 0:
            continue
        primary = sub[(sub.get("bridge_sig_ad", False) == True) | (sub.get("candidate_sig", False) == True)].copy()
        if len(primary) == 0:
            primary = sub.copy()
        chosen_parts = []
        for _, layer_df in primary.groupby("layer"):
            chosen_parts.append(layer_df.sort_values("priority_score", ascending=False).head(config.max_markers_per_layer))
        chosen = pd.concat(chosen_parts, ignore_index=True) if chosen_parts else primary
        chosen = chosen.sort_values("priority_score", ascending=False).head(config.max_markers_per_exposure_total)
        selected.append(chosen)
    return pd.concat(selected, ignore_index=True) if selected else pd.DataFrame()


# =============================================================================
# Analysis steps
# =============================================================================


def run_direct_gi_to_idp_screen(
    analysis: pd.DataFrame,
    available_exposures: Sequence[str],
    brain_features: Sequence[str],
    covariate_blocks: Dict[str, List[str]],
    config: AnalysisConfig,
) -> pd.DataFrame:
    logging.info("Running direct GI-to-IDP screen")
    all_results = []
    for exposure in available_exposures:
        exposure_df = apply_temporality_filter(analysis, exposure, "imaging_date")
        exposure_vector = ensure_binary(exposure_df[exposure])
        n_exposed = int((exposure_vector == 1).sum())
        n_unexposed = int((exposure_vector == 0).sum())
        if n_exposed < config.min_exposed or n_unexposed < config.min_unexposed:
            logging.info("Skipping %s: exposed=%d, unexposed=%d", exposure, n_exposed, n_unexposed)
            continue

        rows = []
        for idp in brain_features:
            result = model_ols(exposure_df, idp, [exposure], covariate_blocks["full"], config.min_total_n, robust=True)
            if result is None:
                continue
            coef = extract_ols_coef(result[0], exposure)
            if coef is None:
                continue
            rows.append({
                "exposure": exposure,
                "exposure_abbr": EXPOSURE_ABBR.get(exposure, exposure),
                "idp": idp,
                "idp_label": idp,
                "n": int(len(result[1])),
                "n_exposed": n_exposed,
                "n_unexposed": n_unexposed,
                "beta_gi_idp": coef["beta"],
                "se_gi_idp": coef["se"],
                "ci_low_gi_idp": coef["ci_low"],
                "ci_high_gi_idp": coef["ci_high"],
                "p_gi_idp": coef["p_value"],
            })
        result_df = pd.DataFrame(rows)
        if len(result_df) > 0:
            result_df["fdr_gi_idp"] = fdr_bh(result_df["p_gi_idp"].values)
            save_table(result_df, f"direct_gi_to_idp_screen_full_{exposure}.csv", config.tables_dir)
            all_results.append(result_df)
    out = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
    save_table(out, "direct_gi_to_idp_screen_full.csv", config.tables_dir)
    return out


def select_idp_candidates(gi_idp_all: pd.DataFrame, available_exposures: Sequence[str], config: AnalysisConfig) -> pd.DataFrame:
    selected = []
    for exposure in available_exposures:
        sub = gi_idp_all[gi_idp_all["exposure"] == exposure].copy()
        if len(sub) == 0:
            continue
        sig = sub[sub["fdr_gi_idp"] < config.fdr_alpha].copy()
        if len(sig) == 0:
            sig = sub.sort_values(["p_gi_idp", "beta_gi_idp"], ascending=[True, False]).head(config.max_idps_per_exposure_fallback).copy()
            sig["idp_candidate_source"] = "top_nominal_fallback"
        else:
            sig = sig.sort_values(["fdr_gi_idp", "p_gi_idp", "beta_gi_idp"], ascending=[True, True, False]).head(config.max_idps_per_exposure).copy()
            sig["idp_candidate_source"] = "FDR_significant"
        selected.append(sig)
    return pd.concat(selected, ignore_index=True) if selected else pd.DataFrame()


def write_idp_mapping_template(idp_selected: pd.DataFrame, config: AnalysisConfig) -> None:
    if len(idp_selected) > 0:
        template = idp_selected[["idp", "idp_label"]].drop_duplicates().copy()
    else:
        template = pd.DataFrame(columns=["idp", "idp_label"])
    for col in ["manual_region_label", "manual_brain_category", "manual_hemisphere", "manual_slice_x", "manual_slice_y", "manual_slice_z"]:
        template[col] = ""
    save_table(template, "idp_mapping_template.csv", config.tables_dir)


def run_marker_level_bridge(
    analysis: pd.DataFrame,
    available_exposures: Sequence[str],
    marker_selected: pd.DataFrame,
    idp_selected: pd.DataFrame,
    covariate_blocks: Dict[str, List[str]],
    config: AnalysisConfig,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    logging.info("Running marker-level GI -> marker -> IDP bridge decomposition")
    bridge_rows = []
    screen_rows = []

    for exposure in available_exposures:
        exp_markers = marker_selected[marker_selected["exposure"] == exposure].copy()
        exp_idps = idp_selected[idp_selected["exposure"] == exposure].copy()
        if len(exp_markers) == 0 or len(exp_idps) == 0:
            logging.info("Skipping bridge for %s: markers=%d, idps=%d", exposure, len(exp_markers), len(exp_idps))
            continue

        marker_list = exp_markers["feature"].dropna().astype(str).tolist()
        idp_list = exp_idps["idp"].dropna().astype(str).tolist()
        if len(marker_list) * len(idp_list) > config.max_total_bridge_tests_per_exposure:
            n_markers = min(len(marker_list), max(10, int(config.max_total_bridge_tests_per_exposure / max(1, len(idp_list)))))
            marker_list = exp_markers.sort_values("priority_score", ascending=False)["feature"].head(n_markers).astype(str).tolist()
            if len(marker_list) * len(idp_list) > config.max_total_bridge_tests_per_exposure:
                n_idps = max(10, int(config.max_total_bridge_tests_per_exposure / max(1, len(marker_list))))
                idp_list = exp_idps.sort_values(["fdr_gi_idp", "p_gi_idp"], ascending=[True, True])["idp"].head(n_idps).astype(str).tolist()

        exposure_df = apply_temporality_filter(analysis, exposure, "imaging_date")
        exposure_vector = ensure_binary(exposure_df[exposure])
        if int((exposure_vector == 1).sum()) < config.min_exposed or int((exposure_vector == 0).sum()) < config.min_unexposed:
            continue

        alpha_map: Dict[str, Dict[str, float]] = {}
        for marker in marker_list:
            result_alpha = model_ols(exposure_df, marker, [exposure], covariate_blocks["full"], config.min_total_n, robust=True)
            if result_alpha is None:
                continue
            coef_alpha = extract_ols_coef(result_alpha[0], exposure)
            if coef_alpha is not None:
                alpha_map[marker] = coef_alpha

        for marker in marker_list:
            if marker not in alpha_map:
                continue
            marker_row = exp_markers[exp_markers["feature"] == marker].iloc[0].to_dict()
            alpha = alpha_map[marker]["beta"]
            se_alpha = alpha_map[marker]["se"]
            p_alpha = alpha_map[marker]["p_value"]

            for idp in idp_list:
                result_beta = model_ols(exposure_df, idp, [exposure, marker], covariate_blocks["full"], config.min_total_n, robust=True)
                if result_beta is None:
                    continue
                coef_beta = extract_ols_coef(result_beta[0], marker)
                coef_direct = extract_ols_coef(result_beta[0], exposure)
                if coef_beta is None:
                    continue

                beta = coef_beta["beta"]
                se_beta = coef_beta["se"]
                indirect = alpha * beta
                se_indirect = np.sqrt((beta ** 2) * (se_alpha ** 2) + (alpha ** 2) * (se_beta ** 2))
                z_indirect = indirect / se_indirect if np.isfinite(se_indirect) and se_indirect > 0 else np.nan
                p_indirect = 2 * (1 - norm.cdf(abs(z_indirect))) if np.isfinite(z_indirect) else np.nan

                gi_idp = exp_idps[exp_idps["idp"] == idp]
                row = {
                    "analysis": "marker_level_gi_marker_idp_bridge",
                    "exposure": exposure,
                    "exposure_abbr": EXPOSURE_ABBR.get(exposure, exposure),
                    "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                    "marker": marker,
                    "marker_label": marker_row.get("feature_label", marker),
                    "marker_layer": marker_row.get("layer", ""),
                    "marker_category": marker_row.get("marker_category", classify_metabolic_marker(marker)),
                    "idp": idp,
                    "idp_label": idp,
                    "n": int(len(result_beta[1])),
                    "alpha_gi_marker": alpha,
                    "se_alpha_gi_marker": se_alpha,
                    "p_alpha_gi_marker": p_alpha,
                    "beta_marker_idp_adj_gi": beta,
                    "se_beta_marker_idp_adj_gi": se_beta,
                    "p_beta_marker_idp_adj_gi": coef_beta["p_value"],
                    "direct_gi_idp_adj_marker_beta": coef_direct["beta"] if coef_direct else np.nan,
                    "direct_gi_idp_adj_marker_p": coef_direct["p_value"] if coef_direct else np.nan,
                    "direct_gi_idp_total_beta": gi_idp["beta_gi_idp"].iloc[0] if len(gi_idp) else np.nan,
                    "direct_gi_idp_total_p": gi_idp["p_gi_idp"].iloc[0] if len(gi_idp) else np.nan,
                    "direct_gi_idp_total_fdr": gi_idp["fdr_gi_idp"].iloc[0] if len(gi_idp) else np.nan,
                    "indirect_effect_alpha_beta": indirect,
                    "se_indirect_delta": se_indirect,
                    "z_indirect_delta": z_indirect,
                    "p_indirect_delta": p_indirect,
                    "bridge_direction": "positive" if indirect > 0 else "negative",
                    "ad_bridge_nie_from_metabolic_bridge": marker_row.get("bridge_nie_ad", np.nan),
                    "ad_bridge_sig_from_metabolic_bridge": bool(marker_row.get("bridge_sig_ad", False)),
                    "marker_candidate_sig_from_metabolic_bridge": bool(marker_row.get("candidate_sig", False)),
                    "priority_score_marker": marker_row.get("priority_score", np.nan),
                }
                bridge_rows.append(row)
                screen_rows.append({
                    "exposure": exposure,
                    "exposure_abbr": EXPOSURE_ABBR.get(exposure, exposure),
                    "marker": marker,
                    "marker_label": marker_row.get("feature_label", marker),
                    "idp": idp,
                    "idp_label": idp,
                    "beta_marker_idp_adj_gi": beta,
                    "p_beta_marker_idp_adj_gi": coef_beta["p_value"],
                    "n": int(len(result_beta[1])),
                })

    screen_all = pd.DataFrame(screen_rows)
    if len(screen_all) > 0:
        screen_all["fdr_marker_idp_by_exposure"] = np.nan
        for exposure, idx in screen_all.groupby("exposure").groups.items():
            screen_all.loc[idx, "fdr_marker_idp_by_exposure"] = fdr_bh(screen_all.loc[idx, "p_beta_marker_idp_adj_gi"].values)
    save_table(screen_all, "marker_to_idp_screen.csv", config.tables_dir)

    bridge_all = pd.DataFrame(bridge_rows)
    if len(bridge_all) > 0:
        bridge_all["fdr_indirect_global"] = fdr_bh(bridge_all["p_indirect_delta"].values)
        bridge_all["fdr_indirect_by_exposure"] = np.nan
        bridge_all["fdr_marker_idp_by_exposure"] = np.nan
        for exposure, idx in bridge_all.groupby("exposure").groups.items():
            bridge_all.loc[idx, "fdr_indirect_by_exposure"] = fdr_bh(bridge_all.loc[idx, "p_indirect_delta"].values)
            bridge_all.loc[idx, "fdr_marker_idp_by_exposure"] = fdr_bh(bridge_all.loc[idx, "p_beta_marker_idp_adj_gi"].values)
        bridge_all["bridge_supported"] = (
            (bridge_all["p_alpha_gi_marker"] < 0.05)
            & (
                (bridge_all["p_beta_marker_idp_adj_gi"] < config.marker_to_idp_p_screen)
                | (bridge_all["fdr_marker_idp_by_exposure"] < config.marker_to_idp_fdr_screen)
            )
            & (
                (bridge_all["p_indirect_delta"] < 0.05)
                | (bridge_all["fdr_indirect_by_exposure"] < config.fdr_alpha)
            )
        )
    save_table(bridge_all, "marker_idp_bridge_effects.csv", config.tables_dir)

    bridge_sig = bridge_all[bridge_all["bridge_supported"]].copy() if len(bridge_all) > 0 else pd.DataFrame()
    save_table(bridge_sig, "marker_idp_bridge_effects_significant.csv", config.tables_dir)

    if len(bridge_all) > 0:
        figure_paths = bridge_all.copy()
        figure_paths["rank_key"] = np.where(figure_paths["bridge_supported"], 0, 1)
        figure_paths = figure_paths.sort_values(
            ["rank_key", "fdr_indirect_by_exposure", "p_indirect_delta", "direct_gi_idp_total_fdr"],
            ascending=True,
        ).head(150)
    else:
        figure_paths = pd.DataFrame()
    save_table(figure_paths, "figure_source_paths.csv", config.tables_dir)

    brain_map = build_brain_map_source(bridge_all)
    save_table(brain_map, "brain_map_source_data.csv", config.tables_dir)

    return bridge_all, bridge_sig, screen_all, brain_map


def build_brain_map_source(bridge_all: pd.DataFrame) -> pd.DataFrame:
    if len(bridge_all) == 0:
        return pd.DataFrame()
    brain_map = bridge_all[bridge_all["bridge_supported"]].copy()
    if len(brain_map) == 0:
        brain_map = bridge_all.sort_values("p_indirect_delta").head(200).copy()
        brain_map["map_source"] = "top_nominal_bridge"
    else:
        brain_map["map_source"] = "supported_bridge"
    out = brain_map.groupby(["exposure", "exposure_abbr", "idp", "idp_label"], as_index=False).agg(
        n_markers=("marker", "nunique"),
        sum_abs_indirect=("indirect_effect_alpha_beta", lambda x: np.sum(np.abs(x))),
        signed_sum_indirect=("indirect_effect_alpha_beta", "sum"),
        min_p_indirect=("p_indirect_delta", "min"),
        min_fdr_indirect_by_exposure=("fdr_indirect_by_exposure", "min"),
        top_marker=("marker_label", lambda x: "|".join(list(pd.Series(x).astype(str).head(5)))),
        map_source=("map_source", "first"),
    )
    out["effect_direction"] = np.where(out["signed_sum_indirect"] >= 0, "higher_IDP_bridge", "lower_IDP_bridge")
    return out


# =============================================================================
# IDP support analyses
# =============================================================================


def support_idp_to_event_cox(
    df: pd.DataFrame,
    phenotypes: Sequence[str],
    covariates: Sequence[str],
    event_col: str,
    time_col: str,
    outcome_name: str,
    config: AnalysisConfig,
) -> pd.DataFrame:
    rows = []
    if event_col not in df.columns or time_col not in df.columns:
        return pd.DataFrame()

    for phenotype in phenotypes:
        needed = [time_col, event_col, phenotype] + list(covariates)
        sub = df[[c for c in needed if c in df.columns]].copy()
        if phenotype not in sub.columns:
            continue
        sub[time_col] = safe_numeric(sub[time_col])
        sub[event_col] = ensure_binary(sub[event_col])
        sub[phenotype] = safe_numeric(sub[phenotype])
        sub = sub.loc[sub[time_col].notna() & (sub[time_col] > 0) & sub[event_col].notna() & sub[phenotype].notna()].copy()
        if len(sub) < config.min_total_n or sub[event_col].sum() < config.min_event_support:
            continue
        for col in covariates:
            if col in sub.columns:
                sub[col] = safe_numeric(sub[col]).fillna(0)

        model_df = sub.rename(columns={time_col: "time_years", event_col: "event"}).copy()
        model_df = clean_model_frame(model_df, protected=["time_years", "event", phenotype])
        if phenotype not in model_df.columns:
            continue
        try:
            cph = CoxPHFitter(penalizer=config.cox_penalizer)
            cph.fit(model_df, duration_col="time_years", event_col="event", show_progress=False)
            if phenotype not in cph.summary.index:
                continue
            s = cph.summary.loc[phenotype]
            rows.append({
                "analysis": "brain_idp_to_event_cox",
                "outcome": outcome_name,
                "idp": phenotype,
                "idp_label": phenotype,
                "n": int(len(model_df)),
                "n_event": int(model_df["event"].sum()),
                "coef": float(s["coef"]),
                "hr": float(np.exp(s["coef"])),
                "ci_low": float(np.exp(s["coef lower 95%"])),
                "ci_high": float(np.exp(s["coef upper 95%"])),
                "p_value": float(s["p"]),
            })
        except Exception:
            continue

    out = pd.DataFrame(rows)
    if len(out) > 0:
        out["fdr"] = fdr_bh(out["p_value"].values)
        out["ad_direction"] = np.where(out["coef"] > 0, "higher_event_risk", "lower_event_risk")
    return out


def support_idp_to_ad_horizon_logit(
    df: pd.DataFrame,
    phenotypes: Sequence[str],
    covariates: Sequence[str],
    config: AnalysisConfig,
) -> pd.DataFrame:
    rows = []
    if "ad_event" not in df.columns or "ad_time_years" not in df.columns:
        return pd.DataFrame()

    eligible, y = build_horizon_binary_outcome(df, "ad_time_years", "ad_event", config.support_horizon)
    base = df.loc[eligible].copy()
    outcome_col = f"ad_{config.support_horizon}y"
    base[outcome_col] = y.loc[base.index]
    base = base.loc[base[outcome_col].notna()].copy()
    base[outcome_col] = base[outcome_col].astype(int)
    if base[outcome_col].sum() < config.min_event_support:
        return pd.DataFrame()

    for phenotype in phenotypes:
        needed = [outcome_col, phenotype] + list(covariates)
        sub = base[[c for c in needed if c in base.columns]].copy()
        if phenotype not in sub.columns:
            continue
        sub[phenotype] = safe_numeric(sub[phenotype])
        sub = sub.loc[sub[phenotype].notna()].copy()
        if len(sub) < config.min_total_n or sub[outcome_col].sum() < config.min_event_support:
            continue
        for col in covariates:
            if col in sub.columns:
                sub[col] = safe_numeric(sub[col]).fillna(0)

        x = clean_model_frame(sub[[phenotype] + [c for c in covariates if c in sub.columns]], protected=[phenotype])
        if phenotype not in x.columns:
            continue
        yv = sub.loc[x.index, outcome_col].astype(int).values
        try:
            model = LogisticRegression(max_iter=800, solver="lbfgs")
            model.fit(x.values, yv)
            coef = float(model.coef_[0][list(x.columns).index(phenotype)])
            rows.append({
                "analysis": "brain_idp_to_ad_horizon_logit",
                "horizon_years": config.support_horizon,
                "idp": phenotype,
                "idp_label": phenotype,
                "n": int(len(x)),
                "n_event": int(np.sum(yv)),
                "coef": coef,
                "or": float(np.exp(coef)),
            })
        except Exception:
            continue
    return pd.DataFrame(rows)


def run_idp_support_analysis(
    analysis: pd.DataFrame,
    brain_features: Sequence[str],
    panel_idps: Sequence[str],
    idp_selected: pd.DataFrame,
    bridge_all: pd.DataFrame,
    covariate_blocks: Dict[str, List[str]],
    config: AnalysisConfig,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    logging.info("Running IDP-to-AD/dementia support analyses")
    support_idps = list(dict.fromkeys(brain_features))
    if len(support_idps) > 1000:
        candidate_idps: List[str] = []
        if len(idp_selected) > 0:
            candidate_idps.extend(idp_selected["idp"].dropna().astype(str).tolist())
        if len(bridge_all) > 0:
            candidate_idps.extend(bridge_all.sort_values("p_indirect_delta").head(300)["idp"].dropna().astype(str).tolist())
        candidate_idps.extend(panel_idps)
        support_idps = list(dict.fromkeys([x for x in candidate_idps if x in brain_features]))

    idp_ad = support_idp_to_event_cox(
        analysis, support_idps, covariate_blocks["full"], "ad_event", "ad_time_years", "strict_ad", config
    )
    idp_dementia = support_idp_to_event_cox(
        analysis, support_idps, covariate_blocks["full"], "dementia_event", "dementia_time_years", "all_cause_dementia", config
    )
    idp_logit = support_idp_to_ad_horizon_logit(analysis, support_idps, covariate_blocks["full"], config)

    save_table(idp_ad, "idp_to_ad_cox_all.csv", config.tables_dir)
    save_table(idp_dementia, "idp_to_dementia_cox_all.csv", config.tables_dir)
    save_table(idp_logit, f"idp_to_ad_logit_{config.support_horizon}y.csv", config.tables_dir)

    if len(idp_ad) > 0:
        ad_map = idp_ad.sort_values(["fdr", "p_value"]).copy()
        ad_map["is_fdr_sig"] = ad_map["fdr"] < config.fdr_alpha
        save_table(ad_map.head(300), "ad_related_idp_top_for_brain_map.csv", config.tables_dir)

    if len(bridge_all) > 0 and len(idp_ad) > 0:
        overlay = bridge_all.merge(
            idp_ad[["idp", "hr", "ci_low", "ci_high", "p_value", "fdr", "ad_direction"]].rename(columns={
                "hr": "idp_to_ad_hr",
                "ci_low": "idp_to_ad_ci_low",
                "ci_high": "idp_to_ad_ci_high",
                "p_value": "idp_to_ad_p",
                "fdr": "idp_to_ad_fdr",
                "ad_direction": "idp_to_ad_direction",
            }),
            on="idp",
            how="left",
        )
        overlay["path_has_ad_idp_support"] = (overlay["idp_to_ad_p"] < 0.05) | (overlay["idp_to_ad_fdr"] < config.fdr_alpha)
        overlay = overlay.sort_values(
            ["bridge_supported", "path_has_ad_idp_support", "p_indirect_delta"],
            ascending=[False, False, True],
        )
    else:
        overlay = pd.DataFrame()
    save_table(overlay, "marker_idp_bridge_with_ad_overlay.csv", config.tables_dir)

    return idp_ad, idp_dementia, idp_logit, overlay


# =============================================================================
# Reporting
# =============================================================================


def write_qc_tables(
    analysis: pd.DataFrame,
    covariate_blocks: Dict[str, List[str]],
    brain_features: Sequence[str],
    metabolic_features: Sequence[str],
    available_exposures: Sequence[str],
    config: AnalysisConfig,
) -> None:
    save_table(
        pd.DataFrame([
            {"adjustment": name, "n_covariates": len(cols), "covariates": ";".join(cols)}
            for name, cols in covariate_blocks.items()
        ]),
        "covariate_blocks.csv",
        config.qc_dir,
    )
    rows = []
    for exposure in available_exposures:
        x = ensure_binary(analysis[exposure])
        rows.append({
            "exposure": exposure,
            "exposure_abbr": EXPOSURE_ABBR.get(exposure, exposure),
            "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
            "n_exposed": int((x == 1).sum()),
            "n_unexposed": int((x == 0).sum()),
            "missing_rate": float(x.isna().mean()),
            "date_col": get_date_col(analysis, exposure) or "",
        })
    save_table(pd.DataFrame(rows), "exposure_qc.csv", config.qc_dir)
    save_table(
        pd.DataFrame([{
            "n_analysis": int(len(analysis)),
            "n_brain_features": int(len(brain_features)),
            "n_metabolic_features": int(len(metabolic_features)),
            "n_available_exposures": int(len(available_exposures)),
            "available_exposures": ";".join(available_exposures),
            "n_ad_events": int(ensure_binary(analysis["ad_event"]).sum(skipna=True)) if "ad_event" in analysis else np.nan,
            "n_dementia_events": int(ensure_binary(analysis["dementia_event"]).sum(skipna=True)) if "dementia_event" in analysis else np.nan,
        }]),
        "cohort_qc.csv",
        config.qc_dir,
    )


def build_summary_by_exposure(
    available_exposures: Sequence[str],
    marker_selected: pd.DataFrame,
    idp_selected: pd.DataFrame,
    bridge_all: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for exposure in available_exposures:
        if len(bridge_all) > 0:
            supported_mask = (bridge_all["exposure"] == exposure) & bridge_all["bridge_supported"]
            n_bridge_tests = int((bridge_all["exposure"] == exposure).sum())
            n_supported = int(supported_mask.sum())
            n_supported_markers = int(bridge_all.loc[supported_mask, "marker"].nunique())
            n_supported_idps = int(bridge_all.loc[supported_mask, "idp"].nunique())
        else:
            n_bridge_tests = n_supported = n_supported_markers = n_supported_idps = 0
        rows.append({
            "exposure": exposure,
            "exposure_abbr": EXPOSURE_ABBR.get(exposure, exposure),
            "n_marker_candidates": int((marker_selected["exposure"] == exposure).sum()) if len(marker_selected) else 0,
            "n_idp_candidates": int((idp_selected["exposure"] == exposure).sum()) if len(idp_selected) else 0,
            "n_bridge_tests": n_bridge_tests,
            "n_supported_bridge_paths": n_supported,
            "n_unique_supported_markers": n_supported_markers,
            "n_unique_supported_idps": n_supported_idps,
        })
    return pd.DataFrame(rows)


def write_metadata_and_report(
    config: AnalysisConfig,
    elapsed_minutes: float,
    analysis: pd.DataFrame,
    brain_features: Sequence[str],
    metabolic_features: Sequence[str],
    available_exposures: Sequence[str],
    marker_selected: pd.DataFrame,
    idp_selected: pd.DataFrame,
    bridge_all: pd.DataFrame,
    idp_ad: pd.DataFrame,
    idp_dementia: pd.DataFrame,
    summary: pd.DataFrame,
) -> None:
    metadata = {
        "script": "metabolic_neural_bridge_analysis.py",
        "inputs": {
            "master": str(config.master),
            "brain": str(config.brain),
            "metabolic": str(config.metabolic),
            "metabolic_bridge_dir": str(config.metabolic_bridge_dir) if config.metabolic_bridge_dir else None,
            "panel_file": str(config.panel_file) if config.panel_file else None,
        },
        "outdir": str(config.outdir),
        "design": {
            "main_question": "Which GI disease is linked to which brain IDP through which individual blood/NMR marker?",
            "step_1": "GI -> brain IDP direct screen",
            "step_2": "metabolic-bridge-derived individual marker candidates",
            "step_3": "marker-level product-of-coefficients bridge: GI -> marker -> IDP",
            "step_4": "IDP -> AD/dementia support Cox screen",
            "step_5": "AD overlay onto marker-level GI-marker-IDP bridge paths",
        },
        "interpretation": "Marker-level metabolic-neural bridge evidence; causal interpretation should remain cautious because metabolic marker and imaging timing may be cross-sectional.",
        "parameters": asdict(config) | {
            "master": str(config.master),
            "brain": str(config.brain),
            "metabolic": str(config.metabolic),
            "outdir": str(config.outdir),
            "metabolic_bridge_dir": str(config.metabolic_bridge_dir) if config.metabolic_bridge_dir else None,
            "panel_file": str(config.panel_file) if config.panel_file else None,
        },
        "summary": {
            "n_analysis": int(len(analysis)),
            "n_brain_features": int(len(brain_features)),
            "n_metabolic_features": int(len(metabolic_features)),
            "available_exposures": list(available_exposures),
            "n_marker_candidates_selected": int(len(marker_selected)) if len(marker_selected) else 0,
            "n_idp_candidates_selected": int(len(idp_selected)) if len(idp_selected) else 0,
            "n_bridge_tests": int(len(bridge_all)) if len(bridge_all) else 0,
            "n_supported_bridge_paths": int(bridge_all["bridge_supported"].sum()) if len(bridge_all) else 0,
            "n_idp_to_ad_tests": int(len(idp_ad)) if len(idp_ad) else 0,
            "n_idp_to_ad_fdr_sig": int((idp_ad["fdr"] < config.fdr_alpha).sum()) if len(idp_ad) else 0,
            "n_idp_to_dementia_tests": int(len(idp_dementia)) if len(idp_dementia) else 0,
            "elapsed_minutes": elapsed_minutes,
        },
    }
    with open(config.outdir / "run_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2, default=str)

    report = []
    report.append("Metabolic-neural marker-level bridge analysis finished.")
    report.append("=" * 70)
    report.append(f"Analysis cohort n = {len(analysis)}")
    report.append(f"Brain IDPs tested = {len(brain_features)}")
    report.append(f"Metabolic markers available = {len(metabolic_features)}")
    report.append(f"Marker candidates selected = {len(marker_selected) if len(marker_selected) else 0}")
    report.append(f"IDP candidates selected = {len(idp_selected) if len(idp_selected) else 0}")
    report.append(f"Bridge tests = {len(bridge_all) if len(bridge_all) else 0}")
    report.append(f"Supported bridge paths = {int(bridge_all['bridge_supported'].sum()) if len(bridge_all) else 0}")
    report.append(f"IDP -> AD Cox tests = {len(idp_ad) if len(idp_ad) else 0}")
    report.append(f"IDP -> AD FDR significant = {int((idp_ad['fdr'] < config.fdr_alpha).sum()) if len(idp_ad) else 0}")
    report.append(f"IDP -> AD nominal P < 0.05 = {int((idp_ad['p_value'] < 0.05).sum()) if len(idp_ad) else 0}")
    report.append(f"Elapsed minutes = {elapsed_minutes}")
    report.append("")
    report.append("Primary outputs:")
    report.append("  tables/marker_idp_bridge_effects.csv")
    report.append("  tables/marker_idp_bridge_effects_significant.csv")
    report.append("  tables/figure_source_paths.csv")
    report.append("  tables/brain_map_source_data.csv")
    report.append("  tables/idp_to_ad_cox_all.csv")
    report.append("  tables/ad_related_idp_top_for_brain_map.csv")
    report.append("")
    report.append("Summary by exposure:")
    report.append(summary.to_string(index=False) if len(summary) else "No exposure summary rows.")
    with open(config.outdir / "analysis_report.txt", "w", encoding="utf-8") as handle:
        handle.write("\n".join(report))


# =============================================================================
# Main workflow
# =============================================================================


def run_analysis(config: AnalysisConfig) -> None:
    np.random.seed(config.seed)
    setup_output_dirs(config)
    setup_logging(config.outdir)
    start_time = time.time()

    analysis, brain, metabolic, covariate_blocks, brain_features, metabolic_features, panel_idps = load_and_prepare_data(config)
    available_exposures = [exp for exp in ALL_EXPOSURES if exp in analysis.columns]
    logging.info("Analysis N=%d", len(analysis))
    logging.info("Available exposures=%s", available_exposures)
    logging.info("Brain features=%d; metabolic features=%d", len(brain_features), len(metabolic_features))

    write_qc_tables(analysis, covariate_blocks, brain_features, metabolic_features, available_exposures, config)

    candidate_file = find_metabolic_bridge_table(config, "candidate_markers_all.csv")
    bridge_marker_file = find_metabolic_bridge_table(config, "bridge_marker_all.csv")
    candidate_raw = read_csv_if_exists(candidate_file)
    bridge_raw = read_csv_if_exists(bridge_marker_file)

    marker_candidates = standardise_metabolic_bridge_candidates(candidate_raw, bridge_raw, metabolic_features, config.fdr_alpha)
    if len(marker_candidates) == 0:
        marker_candidates = fallback_gi_to_marker_screen(analysis, available_exposures, metabolic_features, covariate_blocks, config)
        if len(marker_candidates) == 0:
            logging.warning("No marker candidates were available or generated. Downstream bridge analysis may be empty.")

    marker_selected = select_marker_candidates(marker_candidates, available_exposures, config)
    save_table(marker_candidates, "all_marker_candidates_from_metabolic_bridge_standardised.csv", config.tables_dir)
    save_table(marker_selected, "marker_candidates_from_metabolic_bridge.csv", config.tables_dir)

    gi_idp_all = run_direct_gi_to_idp_screen(analysis, available_exposures, brain_features, covariate_blocks, config)
    idp_selected = select_idp_candidates(gi_idp_all, available_exposures, config)
    save_table(idp_selected, "idp_candidates_by_exposure.csv", config.tables_dir)
    write_idp_mapping_template(idp_selected, config)

    bridge_all, bridge_sig, marker_idp_screen, brain_map = run_marker_level_bridge(
        analysis, available_exposures, marker_selected, idp_selected, covariate_blocks, config
    )

    idp_ad, idp_dementia, idp_logit, overlay = run_idp_support_analysis(
        analysis, brain_features, panel_idps, idp_selected, bridge_all, covariate_blocks, config
    )

    summary = build_summary_by_exposure(available_exposures, marker_selected, idp_selected, bridge_all)
    save_table(summary, "marker_level_bridge_summary_by_exposure.csv", config.tables_dir)

    cohort_summary = pd.DataFrame([{
        "n_analysis": int(len(analysis)),
        "n_brain_features": int(len(brain_features)),
        "n_metabolic_features": int(len(metabolic_features)),
        "n_available_exposures": int(len(available_exposures)),
        "available_exposures": ";".join(available_exposures),
        "n_marker_candidates_total": int(len(marker_selected)) if len(marker_selected) else 0,
        "n_idp_candidates_total": int(len(idp_selected)) if len(idp_selected) else 0,
        "n_bridge_tests_total": int(len(bridge_all)) if len(bridge_all) else 0,
        "n_supported_bridge_paths_total": int(bridge_all["bridge_supported"].sum()) if len(bridge_all) else 0,
        "n_idp_to_ad_tests": int(len(idp_ad)) if len(idp_ad) else 0,
        "n_idp_to_ad_fdr_sig": int((idp_ad["fdr"] < config.fdr_alpha).sum()) if len(idp_ad) else 0,
        "n_idp_to_ad_nominal": int((idp_ad["p_value"] < 0.05).sum()) if len(idp_ad) else 0,
        "n_idp_to_dementia_tests": int(len(idp_dementia)) if len(idp_dementia) else 0,
        "n_idp_to_dementia_fdr_sig": int((idp_dementia["fdr"] < config.fdr_alpha).sum()) if len(idp_dementia) else 0,
        "elapsed_minutes": round((time.time() - start_time) / 60, 2),
    }])
    save_table(cohort_summary, "marker_level_bridge_cohort_summary.csv", config.tables_dir)

    elapsed_minutes = round((time.time() - start_time) / 60, 2)
    write_metadata_and_report(
        config=config,
        elapsed_minutes=elapsed_minutes,
        analysis=analysis,
        brain_features=brain_features,
        metabolic_features=metabolic_features,
        available_exposures=available_exposures,
        marker_selected=marker_selected,
        idp_selected=idp_selected,
        bridge_all=bridge_all,
        idp_ad=idp_ad,
        idp_dementia=idp_dementia,
        summary=summary,
    )

    logging.info("Done. Outputs saved to %s", config.outdir)
    logging.info("Read first: %s", config.outdir / "analysis_report.txt")


def parse_args() -> AnalysisConfig:
    parser = argparse.ArgumentParser(
        description="Marker-level metabolic-neural bridge analysis for GI disease and AD-related outcomes."
    )
    parser.add_argument("--master", required=True, type=Path, help="Participant-level master CSV file.")
    parser.add_argument("--brain", required=True, type=Path, help="Prepared brain IDP CSV file.")
    parser.add_argument("--metabolic", required=True, type=Path, help="Prepared blood/NMR metabolic marker CSV file.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    parser.add_argument(
        "--metabolic-bridge-dir",
        type=Path,
        default=None,
        help="Directory containing metabolic bridge outputs, usually with tables/candidate_markers_all.csv and tables/bridge_marker_all.csv.",
    )
    parser.add_argument("--panel-file", type=Path, default=None, help="Optional text/CSV file listing IDPs to force into support analyses.")
    parser.add_argument("--min-total-n", type=int, default=1200, help="Minimum sample size for model fitting.")
    parser.add_argument("--min-exposed", type=int, default=20, help="Minimum exposed sample size per exposure.")
    parser.add_argument("--min-unexposed", type=int, default=20, help="Minimum unexposed sample size per exposure.")
    parser.add_argument("--min-event-support", type=int, default=20, help="Minimum events for IDP-to-event support models.")
    parser.add_argument("--fdr-alpha", type=float, default=0.05, help="FDR threshold for supported evidence.")
    parser.add_argument("--marker-to-idp-p-screen", type=float, default=0.01, help="Nominal marker-to-IDP P-value threshold.")
    parser.add_argument("--marker-to-idp-fdr-screen", type=float, default=0.10, help="Marker-to-IDP FDR threshold within exposure.")
    parser.add_argument("--max-markers-per-exposure-total", type=int, default=50, help="Maximum selected metabolic markers per exposure.")
    parser.add_argument("--max-markers-per-layer", type=int, default=30, help="Maximum selected metabolic markers per layer and exposure.")
    parser.add_argument("--max-idps-per-exposure", type=int, default=60, help="Maximum selected IDPs per exposure when FDR-significant IDPs exist.")
    parser.add_argument("--max-idps-per-exposure-fallback", type=int, default=40, help="Maximum fallback IDPs per exposure when no FDR-significant IDPs exist.")
    parser.add_argument("--max-total-bridge-tests-per-exposure", type=int, default=3000, help="Maximum marker-IDP bridge tests per exposure.")
    parser.add_argument("--support-horizon", type=int, default=10, help="AD horizon in years for logistic support analysis.")
    parser.add_argument("--cox-penalizer", type=float, default=0.01, help="L2 penalizer used in Cox models.")
    parser.add_argument("--min-feature-nonmissing", type=int, default=100, help="Minimum non-missing values for retaining brain/metabolic features.")
    parser.add_argument("--seed", type=int, default=RANDOM_SEED, help="Random seed.")

    args = parser.parse_args()
    return AnalysisConfig(
        master=args.master,
        brain=args.brain,
        metabolic=args.metabolic,
        outdir=args.outdir,
        metabolic_bridge_dir=args.metabolic_bridge_dir,
        panel_file=args.panel_file,
        min_total_n=args.min_total_n,
        min_exposed=args.min_exposed,
        min_unexposed=args.min_unexposed,
        min_event_support=args.min_event_support,
        fdr_alpha=args.fdr_alpha,
        marker_to_idp_p_screen=args.marker_to_idp_p_screen,
        marker_to_idp_fdr_screen=args.marker_to_idp_fdr_screen,
        max_markers_per_exposure_total=args.max_markers_per_exposure_total,
        max_markers_per_layer=args.max_markers_per_layer,
        max_idps_per_exposure=args.max_idps_per_exposure,
        max_idps_per_exposure_fallback=args.max_idps_per_exposure_fallback,
        max_total_bridge_tests_per_exposure=args.max_total_bridge_tests_per_exposure,
        support_horizon=args.support_horizon,
        cox_penalizer=args.cox_penalizer,
        min_feature_nonmissing=args.min_feature_nonmissing,
        seed=args.seed,
    )


def main() -> None:
    config = parse_args()
    run_analysis(config)


if __name__ == "__main__":
    main()
