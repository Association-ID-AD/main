#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gi_ad_analysis.py

Gastrointestinal disease exposure and Alzheimer-related outcome analysis.

This script performs participant-level association and risk analyses for a
longitudinal GI-AD study. It supports:

1. construction of grouped gastrointestinal disease exposures;
2. construction or reuse of MCI-or-AD composite outcomes;
3. staged Cox proportional hazards models;
4. baseline logistic support analysis;
5. optional DeepHit absolute-risk modelling;
6. reverse-causation, exposure-recency, detection-bias, subgroup and
   interaction analyses;
7. QC tables and run metadata for reproducibility.

Example
-------
python gi_ad_analysis.py \
    --input data/master_preprocessed_with_genetics.csv \
    --outdir results/gi_ad_analysis \
    --run-deephit

Minimal run without DeepHit
---------------------------
python gi_ad_analysis.py \
    --input data/master_preprocessed_with_genetics.csv \
    --outdir results/gi_ad_analysis \
    --skip-deephit
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import random
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

try:
    import torch
    import torchtuples as tt
    from pycox.evaluation import EvalSurv
    from pycox.models import DeepHitSingle
except Exception:  # pragma: no cover - optional dependency checked at runtime
    torch = None
    tt = None
    EvalSurv = None
    DeepHitSingle = None

warnings.filterwarnings("default")


# -----------------------------------------------------------------------------
# Study definitions
# -----------------------------------------------------------------------------

GROUPED_GI_MAP: Dict[str, List[str]] = {
    "exp_ibd": ["exp_k50", "exp_k51"],
    "exp_ibs": ["exp_k58"],
    "exp_other_functional_intestinal": ["exp_k59"],
    "exp_diverticular": ["exp_k57"],
    "exp_anorectal": ["exp_k60", "exp_k61", "exp_k62", "exp_k64"],
    "exp_malabsorption": ["exp_k90"],
    "exp_other_chronic_intestinal": ["exp_k52", "exp_k55", "exp_k56", "exp_k63"],
    "exp_appendiceal": ["exp_k35", "exp_k36", "exp_k37", "exp_k38"],
}

GROUPED_LABEL: Dict[str, str] = {
    "exp_ibd": "IBD (K50+K51)",
    "exp_ibs": "IBS (K58)",
    "exp_other_functional_intestinal": "Other functional intestinal disorders (K59)",
    "exp_diverticular": "Diverticular disease (K57)",
    "exp_anorectal": "Anorectal disease (K60+K61+K62+K64)",
    "exp_malabsorption": "Malabsorption (K90)",
    "exp_other_chronic_intestinal": "Other chronic intestinal disease (K52+K55+K56+K63)",
    "exp_appendiceal": "Appendiceal disease (K35-K38)",
}

DISEASE_LABEL: Dict[str, str] = {
    "exp_k26": "K26 Duodenal ulcer",
    "exp_k27": "K27 Peptic ulcer, site unspecified",
    "exp_k29": "K29 Gastritis and duodenitis",
    "exp_k35": "K35 Acute appendicitis",
    "exp_k36": "K36 Other appendicitis",
    "exp_k37": "K37 Appendicitis, unspecified",
    "exp_k38": "K38 Other diseases of appendix",
    "exp_k50": "K50 Crohn's disease",
    "exp_k51": "K51 Ulcerative colitis",
    "exp_k52": "K52 Other noninfective gastroenteritis/colitis",
    "exp_k55": "K55 Vascular disorders of intestine",
    "exp_k56": "K56 Paralytic ileus / intestinal obstruction",
    "exp_k57": "K57 Diverticular disease",
    "exp_k58": "K58 Irritable bowel syndrome",
    "exp_k59": "K59 Other functional intestinal disorders",
    "exp_k60": "K60 Fissure / fistula of anal / rectal region",
    "exp_k61": "K61 Abscess of anal / rectal region",
    "exp_k62": "K62 Other diseases of anus / rectum",
    "exp_k63": "K63 Other diseases of intestine",
    "exp_k64": "K64 Hemorrhoids / perianal venous thrombosis",
    "exp_k90": "K90 Intestinal malabsorption",
}

MAIN_GROUPED_EXPOSURES = [
    "exp_diverticular",
    "exp_other_functional_intestinal",
    "exp_ibd",
    "exp_other_chronic_intestinal",
    "exp_ibs",
]
SECONDARY_GROUPED_EXPOSURES = ["exp_appendiceal", "exp_malabsorption"]
SUPPLEMENTARY_GROUPED_EXPOSURES = ["exp_anorectal"]
DEEPHIT_EXPOSURES = MAIN_GROUPED_EXPOSURES + SECONDARY_GROUPED_EXPOSURES

SURVIVAL_OUTCOMES: Dict[str, Dict[str, str]] = {
    "strict_ad": {
        "event_col": "ad_event",
        "time_col": "ad_time_years",
        "label": "Incident AD",
        "role": "primary",
    },
    "all_cause_dementia": {
        "event_col": "dementia_event",
        "time_col": "dementia_time_years",
        "label": "All-cause dementia",
        "role": "secondary",
    },
    "cognition_defined_mci": {
        "event_col": "mci_event",
        "time_col": "mci_time_years",
        "label": "Cognition-defined MCI-like status",
        "role": "secondary/support",
    },
    "incident_mci_or_ad": {
        "event_col": "mci_ad_event",
        "time_col": "mci_ad_time_years",
        "label": "Incident MCI-or-AD composite",
        "role": "secondary/support composite",
    },
}

LOGISTIC_OUTCOMES: Dict[str, Dict[str, str]] = {
    "baseline_broad_mci_ad": {
        "event_col": "baseline_broad_mci_ad",
        "label": "Baseline broad MCI/AD",
    }
}

DEEPHIT_OUTCOMES = [
    "strict_ad",
    "all_cause_dementia",
    "cognition_defined_mci",
    "incident_mci_or_ad",
]


@dataclass(frozen=True)
class AnalysisConfig:
    input_file: Path
    output_dir: Path
    random_seed: int = 20260509
    min_total_n: int = 200
    min_events: int = 40
    min_exposed: int = 30
    min_unexposed: int = 30
    cox_penalizer: float = 0.01
    run_deephit: bool = False
    deephit_max_n: Optional[int] = 120_000
    num_durations: int = 50
    batch_size: int = 512
    epochs: int = 120
    patience: int = 15
    learning_rate: float = 1e-3
    horizons: Tuple[int, ...] = (5, 10, 15, 20)
    trajectory_grid_start: int = 1
    trajectory_grid_end: int = 20
    n_boot: int = 200
    reverse_causation_lags: Tuple[int, ...] = (1, 2, 5)
    exposure_recency_years: Tuple[int, ...] = (1, 2, 5)

    @property
    def trajectory_grid(self) -> np.ndarray:
        return np.arange(self.trajectory_grid_start, self.trajectory_grid_end + 1)


# -----------------------------------------------------------------------------
# Generic utilities
# -----------------------------------------------------------------------------


def setup_logging(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / "run.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, encoding="utf-8")],
    )


def set_random_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    if torch is not None:
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)


def safe_numeric(x: object) -> pd.Series:
    return pd.to_numeric(x, errors="coerce")


def ensure_binary(x: object) -> pd.Series:
    s = safe_numeric(x)
    return s.where(s.isin([0, 1]), np.nan)


def zscore(x: object) -> pd.Series:
    s = safe_numeric(x)
    mu = s.mean(skipna=True)
    sd = s.std(skipna=True)
    if pd.isna(sd) or sd == 0:
        return pd.Series(np.nan, index=s.index)
    return (s - mu) / sd


def fill_numeric(x: object) -> pd.Series:
    s = safe_numeric(x)
    med = s.median(skipna=True)
    if pd.isna(med):
        med = 0.0
    return s.fillna(med)


def missing_indicator(x: object) -> pd.Series:
    return pd.isna(x).astype(int)


def exposure_label(col: str) -> str:
    return GROUPED_LABEL.get(col, DISEASE_LABEL.get(col, col))


def duration_col_from_exposure(exp_col: str) -> str:
    return exp_col.replace("exp_", "dur_") + "_to_baseline"


def exposure_date_candidates(exp_col: str) -> List[str]:
    suffix = exp_col.replace("exp_", "")
    return [
        f"date_{suffix}",
        f"dt_{suffix}",
        exp_col.replace("exp_", "date_"),
        exp_col.replace("exp_", "dt_"),
    ]


def get_exposure_date_col(df: pd.DataFrame, exp_col: str) -> Optional[str]:
    for col in exposure_date_candidates(exp_col):
        if col in df.columns:
            return col
    return None


def first_existing_col(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    for col in candidates:
        if col in df.columns:
            return col
    return None


def safe_concat(dfs: Iterable[pd.DataFrame]) -> pd.DataFrame:
    valid = [x for x in dfs if isinstance(x, pd.DataFrame) and not x.empty]
    if not valid:
        return pd.DataFrame()
    return pd.concat(valid, axis=0, ignore_index=True)


def bh_fdr(p_values: Sequence[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
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
    q[idx[order]] = np.clip(q_ranked, 0, 1)
    return q


def add_fdr(df: pd.DataFrame, group_cols: Sequence[str]) -> pd.DataFrame:
    if df is None or df.empty or "p_value" not in df.columns:
        return df

    out = df.copy()
    out["q_value"] = np.nan
    if all(col in out.columns for col in group_cols):
        for _, idx in out.groupby(list(group_cols), dropna=False).groups.items():
            out.loc[idx, "q_value"] = bh_fdr(out.loc[idx, "p_value"].values)
    else:
        out["q_value"] = bh_fdr(out["p_value"].values)
    return out


# -----------------------------------------------------------------------------
# Data loading and preprocessing
# -----------------------------------------------------------------------------


def load_data(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    df = pd.read_csv(path, low_memory=False)
    for col in df.columns:
        if (
            col in [
                "time0",
                "ad_date",
                "dementia_date",
                "censor_date",
                "death_date",
                "lost_followup_date",
                "mci_date",
                "mci_first_date",
                "dt_mci",
                "dt_mci_pre",
                "date_mci",
            ]
            or col.startswith("date_")
            or col.startswith("dt_")
        ):
            df[col] = pd.to_datetime(df[col], errors="coerce")
    return df


def max_binary_across(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    valid_cols = [col for col in cols if col in df.columns]
    if not valid_cols:
        return pd.Series(np.nan, index=df.index)

    tmp = pd.concat([safe_numeric(df[col]) for col in valid_cols], axis=1)
    out = tmp.max(axis=1, skipna=True)
    out[tmp.isna().all(axis=1)] = np.nan
    return out


def min_duration_across(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    valid_cols = [col for col in cols if col in df.columns]
    if not valid_cols:
        return pd.Series(np.nan, index=df.index)

    tmp = pd.concat([safe_numeric(df[col]) for col in valid_cols], axis=1)
    return tmp.min(axis=1, skipna=True)


def earliest_datetime_across(df: pd.DataFrame, cols: Sequence[Optional[str]]) -> pd.Series:
    valid_cols = [col for col in cols if col is not None and col in df.columns]
    if not valid_cols:
        return pd.Series(pd.NaT, index=df.index)

    tmp = pd.concat([pd.to_datetime(df[col], errors="coerce") for col in valid_cols], axis=1)
    return tmp.min(axis=1, skipna=True)


def add_grouped_exposures(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()

    for grouped_exp, components in GROUPED_GI_MAP.items():
        d[grouped_exp] = max_binary_across(d, components)
        d[duration_col_from_exposure(grouped_exp)] = min_duration_across(
            d, [duration_col_from_exposure(c) for c in components]
        )
        d[f"date_{grouped_exp.replace('exp_', '')}"] = earliest_datetime_across(
            d, [get_exposure_date_col(d, c) for c in components]
        )

    grouped_cols = [col for col in GROUPED_GI_MAP if col in d.columns]
    if grouped_cols:
        tmp = d[grouped_cols].apply(pd.to_numeric, errors="coerce")
        d["n_grouped_gi"] = tmp.sum(axis=1, skipna=True)
        d.loc[tmp.isna().all(axis=1), "n_grouped_gi"] = np.nan
    else:
        d["n_grouped_gi"] = np.nan

    return d


def build_mci_ad_composite(df: pd.DataFrame) -> pd.DataFrame:
    """Build incident MCI-or-AD composite if it is not already present.

    If integrated preprocessing has already created mci_ad_event and
    mci_ad_time_years, those columns are reused without overwriting.
    """
    d = df.copy()

    if "mci_ad_event" in d.columns and "mci_ad_time_years" in d.columns:
        return d

    if "time0" not in d.columns:
        d["mci_ad_event"] = np.nan
        d["mci_ad_time_years"] = np.nan
        d["mci_ad_date"] = pd.NaT
        return d

    time0 = pd.to_datetime(d["time0"], errors="coerce")
    ad_event = ensure_binary(d.get("ad_event", np.nan)).fillna(0)
    ad_time = safe_numeric(d.get("ad_time_years", np.nan))

    ad_date_col = first_existing_col(d, ["ad_date", "dt_ad", "date_ad"])
    ad_date = (
        pd.to_datetime(d[ad_date_col], errors="coerce")
        if ad_date_col is not None
        else pd.Series(pd.NaT, index=d.index)
    )

    mci_date_col = first_existing_col(
        d, ["mci_date", "dt_mci", "dt_mci_pre", "date_mci", "mci_first_date"]
    )
    if mci_date_col is None:
        d["mci_ad_event"] = np.nan
        d["mci_ad_time_years"] = np.nan
        d["mci_ad_date"] = pd.NaT
        return d

    mci_date = pd.to_datetime(d[mci_date_col], errors="coerce")
    mci_time = (mci_date - time0).dt.days / 365.25
    mci_event = mci_time.notna() & (mci_time >= 0)

    # Avoid zero duration in survival models for baseline MCI-like observations.
    mci_time_adj = mci_time.copy()
    mci_time_adj[mci_event & (mci_time_adj <= 0)] = 1e-4

    d["mci_ad_event"] = ad_event.astype(float)
    d["mci_ad_time_years"] = ad_time.copy()

    valid_ad_follow = ad_time.notna() & (ad_time > 0)
    use_mci = mci_event & (valid_ad_follow | (ad_event == 0))
    use_mci &= (~valid_ad_follow) | (mci_time_adj <= ad_time + 1e-10) | (ad_event == 0)

    d.loc[use_mci, "mci_ad_event"] = 1.0
    d.loc[use_mci, "mci_ad_time_years"] = mci_time_adj.loc[use_mci]

    comp_date = pd.Series(pd.NaT, index=d.index, dtype="datetime64[ns]")
    comp_date.loc[use_mci] = mci_date.loc[use_mci]
    ad_only = (ad_event == 1) & ~use_mci
    comp_date.loc[ad_only] = ad_date.loc[ad_only]
    d["mci_ad_date"] = comp_date
    return d


def add_baseline_broad_mci_ad_if_possible(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    if "baseline_broad_mci_ad" in d.columns:
        return d

    candidates = []
    for col in ["MCI_pre", "AD_pre", "mci_pre", "ad_pre", "baseline_mci", "baseline_ad"]:
        if col in d.columns:
            candidates.append(ensure_binary(d[col]).fillna(0))

    if candidates:
        tmp = pd.concat(candidates, axis=1)
        d["baseline_broad_mci_ad"] = tmp.max(axis=1).astype(float)
    return d


def prepare_dataset(df: pd.DataFrame) -> pd.DataFrame:
    d = add_grouped_exposures(df)
    d = build_mci_ad_composite(d)
    d = add_baseline_broad_mci_ad_if_possible(d)

    if "sex" in d.columns:
        d["sex_model"] = fill_numeric(d["sex"])
        d["sex_missing"] = missing_indicator(d["sex"])

    if "age_at_baseline" in d.columns:
        d["age_at_baseline_raw"] = safe_numeric(d["age_at_baseline"])
        d["age_z"] = fill_numeric(zscore(d["age_at_baseline_raw"]))
        d["age_missing"] = missing_indicator(d["age_at_baseline_raw"])

    categorical_map = {
        "ethnicity_5cat": "ethnic_model",
        "edu_level": "edu_model",
        "health_4cat": "health_model",
    }
    for raw_col, out_col in categorical_map.items():
        if raw_col in d.columns:
            d[out_col] = fill_numeric(d[raw_col])
            d[out_col.replace("_model", "_missing")] = missing_indicator(d[raw_col])

    continuous_map = {
        "townsend_i0": "townsend_z",
        "bmi_i0": "bmi_z",
        "whr_i0": "whr_z",
        "sbp_i0": "sbp_z",
        "dbp_i0": "dbp_z",
        "smoke_pack_i0": "smoke_pack_z",
        "smoke_pack_proxy": "smoke_pack_proxy_z",
        "alcohol_freq": "alcohol_freq_z",
        "sleep_hours": "sleep_hours_z",
        "total_sedentary_hours": "sedentary_z",
        "activity_score_days": "activity_z",
        "diet_quality_proxy": "diet_quality_z",
        "parental_dementia_count": "parental_dementia_count_z",
        "ad_prs_std": "ad_prs_std_z",
        "ad_prs_enhanced_std": "ad_prs_enh_z",
        "p22009": "ad_prs_std_z",
        "p22020": "ad_prs_enh_z",
        "apoe_e4_count": "apoe_e4_count_z",
    }
    for raw_col, out_col in continuous_map.items():
        if raw_col in d.columns and out_col not in d.columns:
            d[out_col] = fill_numeric(zscore(d[raw_col]))
            d[out_col.replace("_z", "_missing")] = missing_indicator(d[raw_col])

    if "smoke_status_3cat" in d.columns:
        d["smoke_status_model"] = fill_numeric(d["smoke_status_3cat"])
        d["smoke_status_missing"] = missing_indicator(d["smoke_status_3cat"])
    elif "smoking_status" in d.columns:
        d["smoke_status_model"] = fill_numeric(d["smoking_status"])
        d["smoke_status_missing"] = missing_indicator(d["smoking_status"])

    binary_cols = [
        "sleep_short",
        "sleep_long",
        "insomnia_any",
        "hx_diabetes",
        "hx_hypertension",
        "hx_hyperlipidemia",
        "hx_chd_cvd",
        "hx_stroke_tia",
        "hx_depression_anxiety",
        "hx_sleep_disorder",
        "hx_parkinson_other_nd",
        "med_statin",
        "med_antihypertensive",
        "med_antidiabetic",
        "med_ppi_gi_drug",
        "med_laxative",
        "med_antidepressant",
        "med_antipsychotic",
        "med_steroid_immunosuppressive",
        "parental_dementia_any",
        "father_dementia",
        "mother_dementia",
        "apoe4_carrier",
        "bowel_cancer_screening",
        "any_current_medication",
    ]
    for col in binary_cols:
        if col in d.columns:
            d[col] = ensure_binary(d[col]).fillna(0)

    if "apoe4_carrier" in d.columns:
        d["apoe_e4_carrier_model"] = ensure_binary(d["apoe4_carrier"]).fillna(0)

    if "bmi_i0" in d.columns:
        d["obesity"] = (safe_numeric(d["bmi_i0"]) >= 30).astype(float)

    return d


# -----------------------------------------------------------------------------
# Covariate blocks
# -----------------------------------------------------------------------------


def unique_existing(df: pd.DataFrame, cols: Sequence[str]) -> List[str]:
    return list(dict.fromkeys([col for col in cols if col in df.columns]))


def model_blocks(df: pd.DataFrame) -> Dict[str, List[str]]:
    m0 = unique_existing(df, ["sex_model", "sex_missing", "age_z", "age_missing"])

    m1 = unique_existing(
        df,
        m0
        + [
            "ethnic_model",
            "ethnic_missing",
            "edu_model",
            "edu_missing",
            "townsend_z",
            "townsend_missing",
            "health_model",
            "health_missing",
        ],
    )

    m2 = unique_existing(
        df,
        m1
        + [
            "bmi_z",
            "bmi_missing",
            "whr_z",
            "whr_missing",
            "sbp_z",
            "sbp_missing",
            "dbp_z",
            "dbp_missing",
            "smoke_status_model",
            "smoke_status_missing",
            "smoke_pack_z",
            "smoke_pack_missing",
            "smoke_pack_proxy_z",
            "smoke_pack_proxy_missing",
            "alcohol_freq_z",
            "alcohol_freq_missing",
            "sleep_hours_z",
            "sleep_hours_missing",
            "sleep_short",
            "sleep_long",
            "insomnia_any",
            "sedentary_z",
            "sedentary_missing",
            "activity_z",
            "activity_missing",
            "diet_quality_z",
            "diet_quality_missing",
        ],
    )

    m3 = unique_existing(
        df,
        m2
        + [
            "hx_diabetes",
            "hx_hypertension",
            "hx_hyperlipidemia",
            "hx_chd_cvd",
            "hx_stroke_tia",
            "hx_depression_anxiety",
            "hx_sleep_disorder",
            "hx_parkinson_other_nd",
        ],
    )

    m4 = unique_existing(
        df,
        m3
        + [
            "med_statin",
            "med_antihypertensive",
            "med_antidiabetic",
            "med_ppi_gi_drug",
            "med_laxative",
            "med_antidepressant",
            "med_antipsychotic",
            "med_steroid_immunosuppressive",
            "any_current_medication",
        ],
    )

    m5 = unique_existing(
        df,
        m4
        + [
            "parental_dementia_any",
            "parental_dementia_count_z",
            "parental_dementia_count_missing",
            "apoe_e4_carrier_model",
            "apoe_e4_count_z",
            "apoe_e4_count_missing",
            "ad_prs_std_z",
            "ad_prs_std_missing",
            "ad_prs_enh_z",
            "ad_prs_enh_missing",
        ],
    )

    m6 = unique_existing(df, m5 + ["bowel_cancer_screening"])

    return {
        "Model0_age_sex": m0,
        "Model1_socioeconomic": m1,
        "Model2_lifestyle": m2,
        "Model3_comorbidity": m3,
        "Model4_medication": m4,
        "Model5_family_genetic": m5,
        "Model6_detection_bias": m6,
    }


# -----------------------------------------------------------------------------
# Model-ready subsets and cleaning
# -----------------------------------------------------------------------------


def passes_minimum_counts(
    dd: pd.DataFrame,
    exp_col: str,
    event_col: str,
    cfg: AnalysisConfig,
) -> bool:
    if len(dd) < cfg.min_total_n:
        return False
    if safe_numeric(dd[event_col]).sum(skipna=True) < cfg.min_events:
        return False
    if dd[exp_col].nunique(dropna=True) < 2:
        return False
    if (dd[exp_col] == 1).sum() < cfg.min_exposed:
        return False
    if (dd[exp_col] == 0).sum() < cfg.min_unexposed:
        return False
    return True


def survival_subset(
    df: pd.DataFrame,
    exp_col: str,
    event_col: str,
    time_col: str,
    covariates: Optional[Sequence[str]],
    cfg: AnalysisConfig,
    include_duration: bool = False,
) -> Optional[pd.DataFrame]:
    """Create an exposure-specific survival subset while retaining covariates."""
    if exp_col not in df.columns or event_col not in df.columns or time_col not in df.columns:
        return None

    keep = ["eid", exp_col, event_col, time_col] + [c for c in (covariates or []) if c in df.columns]
    if include_duration and duration_col_from_exposure(exp_col) in df.columns:
        keep.append(duration_col_from_exposure(exp_col))

    dd = df[list(dict.fromkeys(keep))].copy()
    dd[exp_col] = safe_numeric(dd[exp_col])
    dd[event_col] = ensure_binary(dd[event_col])
    dd[time_col] = safe_numeric(dd[time_col])
    dd = dd.loc[
        dd[exp_col].notna()
        & dd[event_col].notna()
        & dd[time_col].notna()
        & (dd[time_col] > 0)
    ].copy()

    if not passes_minimum_counts(dd, exp_col, event_col, cfg):
        return None
    return dd


def logistic_subset(
    df: pd.DataFrame,
    exp_col: str,
    event_col: str,
    covariates: Optional[Sequence[str]],
    cfg: AnalysisConfig,
) -> Optional[pd.DataFrame]:
    """Create exposure-specific logistic subset while retaining covariates."""
    if exp_col not in df.columns or event_col not in df.columns:
        return None

    keep = ["eid", exp_col, event_col] + [c for c in (covariates or []) if c in df.columns]
    dd = df[list(dict.fromkeys(keep))].copy()
    dd[exp_col] = safe_numeric(dd[exp_col])
    dd[event_col] = ensure_binary(dd[event_col])
    dd = dd.loc[dd[exp_col].notna() & dd[event_col].notna()].copy()

    if not passes_minimum_counts(dd, exp_col, event_col, cfg):
        return None
    return dd


def clean_model_frame(
    model_df: pd.DataFrame,
    event_col: str,
    duration_col: Optional[str] = None,
    protected_cols: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    z = model_df.copy().replace([np.inf, -np.inf], np.nan)
    protected = set(protected_cols or []) | {event_col}
    if duration_col is not None:
        protected.add(duration_col)

    drop_cols = []
    for col in z.columns:
        if col in protected:
            continue
        if z[col].isna().all() or z[col].dropna().nunique() <= 1:
            drop_cols.append(col)
    if drop_cols:
        z = z.drop(columns=drop_cols)

    for col in z.columns:
        if col in protected:
            continue
        if z[col].isna().any():
            median_value = z[col].median(skipna=True)
            if pd.isna(median_value):
                median_value = 0.0
            z[col] = z[col].fillna(median_value)

    if duration_col is not None:
        z = z.loc[z[duration_col].notna() & (z[duration_col] > 0) & z[event_col].notna()].copy()
    else:
        z = z.loc[z[event_col].notna()].copy()
    return z


# -----------------------------------------------------------------------------
# Cox and logistic models
# -----------------------------------------------------------------------------


def run_cox(
    dd: pd.DataFrame,
    exp_col: str,
    event_col: str,
    time_col: str,
    covariates: Sequence[str],
    model_name: str,
    outcome_name: str,
    analysis_name: str,
    exposure_group: str,
    cfg: AnalysisConfig,
) -> Dict[str, object]:
    covariates = [c for c in covariates if c in dd.columns and c != exp_col]
    model_df = dd[[time_col, event_col, exp_col] + covariates].copy()
    for col in model_df.columns:
        model_df[col] = safe_numeric(model_df[col])

    model_df = model_df.rename(columns={time_col: "time_years", event_col: "event"})
    model_df = clean_model_frame(
        model_df,
        event_col="event",
        duration_col="time_years",
        protected_cols=["time_years", "event", exp_col],
    )

    if exp_col not in model_df.columns or model_df[exp_col].nunique() < 2:
        raise ValueError("Exposure invalid after model-frame cleaning")

    cph = CoxPHFitter(penalizer=cfg.cox_penalizer)
    cph.fit(model_df, duration_col="time_years", event_col="event")
    summary = cph.summary.loc[exp_col]

    covariates_used = [c for c in model_df.columns if c not in ["time_years", "event", exp_col]]
    return {
        "analysis_name": analysis_name,
        "outcome": outcome_name,
        "model_name": model_name,
        "exposure_group": exposure_group,
        "exposure": exp_col,
        "exposure_label": exposure_label(exp_col),
        "n_total_after_clean": int(len(model_df)),
        "n_event_after_clean": int(model_df["event"].sum()),
        "n_exposed_after_clean": int((model_df[exp_col] == 1).sum()),
        "n_unexposed_after_clean": int((model_df[exp_col] == 0).sum()),
        "n_covariates_used": int(len(covariates_used)),
        "covariates_used": ";".join(covariates_used),
        "hr": float(np.exp(summary["coef"])),
        "ci_low": float(np.exp(summary["coef lower 95%"])),
        "ci_high": float(np.exp(summary["coef upper 95%"])),
        "p_value": float(summary["p"]),
        "method": "penalized_CoxPH",
        "cox_penalizer": cfg.cox_penalizer,
    }


def run_logistic(
    dd: pd.DataFrame,
    exp_col: str,
    event_col: str,
    covariates: Sequence[str],
    model_name: str,
    outcome_name: str,
    analysis_name: str,
    exposure_group: str,
) -> Dict[str, object]:
    covariates = [c for c in covariates if c in dd.columns and c != exp_col]
    model_df = dd[[event_col, exp_col] + covariates].copy()
    for col in model_df.columns:
        model_df[col] = safe_numeric(model_df[col])

    model_df = clean_model_frame(
        model_df,
        event_col=event_col,
        duration_col=None,
        protected_cols=[event_col, exp_col],
    )

    if exp_col not in model_df.columns or model_df[exp_col].nunique() < 2:
        raise ValueError("Exposure invalid after model-frame cleaning")

    x = model_df.drop(columns=[event_col]).values
    y = model_df[event_col].astype(int).values
    if len(np.unique(y)) < 2:
        raise ValueError("Logistic outcome has only one class")

    model = LogisticRegression(max_iter=1000, solver="lbfgs")
    model.fit(x, y)
    exposure_position = list(model_df.drop(columns=[event_col]).columns).index(exp_col)
    coef = float(model.coef_[0][exposure_position])
    covariates_used = [c for c in model_df.columns if c not in [event_col, exp_col]]

    return {
        "analysis_name": analysis_name,
        "outcome": outcome_name,
        "model_name": model_name,
        "exposure_group": exposure_group,
        "exposure": exp_col,
        "exposure_label": exposure_label(exp_col),
        "n_total": int(len(model_df)),
        "n_event": int(model_df[event_col].sum()),
        "n_covariates_used": int(len(covariates_used)),
        "covariates_used": ";".join(covariates_used),
        "or": float(np.exp(coef)),
        "coef": coef,
        "method": "sklearn_LogisticRegression_support",
    }


# -----------------------------------------------------------------------------
# DeepHit supplementary absolute-risk modelling
# -----------------------------------------------------------------------------


def ensure_deephit_available() -> None:
    missing = []
    if torch is None:
        missing.append("torch")
    if tt is None:
        missing.append("torchtuples")
    if DeepHitSingle is None or EvalSurv is None:
        missing.append("pycox")
    if missing:
        raise ImportError(
            "DeepHit requires optional dependencies that are not available: "
            + ", ".join(sorted(set(missing)))
            + ". Install them or run with --skip-deephit."
        )


def sample_for_deephit(
    df_model: pd.DataFrame,
    event_col: str,
    exp_col: str,
    max_n: Optional[int],
    seed: int,
) -> pd.DataFrame:
    if max_n is None or len(df_model) <= max_n:
        return df_model.copy()

    rng = np.random.default_rng(seed)
    mandatory = np.unique(
        np.concatenate(
            [
                df_model.index[df_model[event_col] == 1].to_numpy(),
                df_model.index[df_model[exp_col] == 1].to_numpy(),
            ]
        )
    )
    remaining = np.setdiff1d(df_model.index.to_numpy(), mandatory)
    n_remaining = max(0, max_n - len(mandatory))

    if n_remaining > 0 and len(remaining) > 0:
        selected = np.unique(
            np.concatenate(
                [mandatory, rng.choice(remaining, size=min(n_remaining, len(remaining)), replace=False)]
            )
        )
    else:
        selected = mandatory[:max_n]
    return df_model.loc[selected].copy()


def split_and_scale_survival(
    df_model: pd.DataFrame,
    event_col: str,
    time_col: str,
    seed: int,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    feature_cols = [c for c in df_model.columns if c not in ["eid", event_col, time_col]]

    stratify = df_model[event_col].astype(int) if df_model[event_col].nunique() > 1 else None
    train_val, test = train_test_split(
        df_model,
        test_size=0.20,
        random_state=seed,
        stratify=stratify,
    )

    stratify_train = train_val[event_col].astype(int) if train_val[event_col].nunique() > 1 else None
    train, valid = train_test_split(
        train_val,
        test_size=0.20,
        random_state=seed,
        stratify=stratify_train,
    )

    x_train = train[feature_cols].copy()
    x_valid = valid[feature_cols].copy()
    x_test = test[feature_cols].copy()

    numeric_cols = []
    for col in feature_cols:
        values = set(pd.unique(x_train[col].dropna()))
        if not (len(values) <= 2 and values.issubset({0, 1})):
            numeric_cols.append(col)

    if numeric_cols:
        scaler = StandardScaler()
        x_train[numeric_cols] = scaler.fit_transform(x_train[numeric_cols])
        x_valid[numeric_cols] = scaler.transform(x_valid[numeric_cols])
        x_test[numeric_cols] = scaler.transform(x_test[numeric_cols])

    return train, valid, test, x_train, x_valid, x_test


def fit_deephit_model(
    x_train: pd.DataFrame,
    train: pd.DataFrame,
    x_valid: pd.DataFrame,
    valid: pd.DataFrame,
    event_col: str,
    time_col: str,
    cfg: AnalysisConfig,
):
    ensure_deephit_available()

    labtrans = DeepHitSingle.label_transform(cfg.num_durations)
    y_train = labtrans.fit_transform(
        np.maximum(safe_numeric(train[time_col]).values.astype("float32"), 1e-4),
        ensure_binary(train[event_col]).fillna(0).values.astype("int64"),
    )
    y_valid = labtrans.transform(
        np.maximum(safe_numeric(valid[time_col]).values.astype("float32"), 1e-4),
        ensure_binary(valid[event_col]).fillna(0).values.astype("int64"),
    )

    net = tt.practical.MLPVanilla(
        in_features=x_train.shape[1],
        num_nodes=[128, 64],
        out_features=labtrans.out_features,
        batch_norm=True,
        dropout=0.10,
    )
    model = DeepHitSingle(
        net,
        tt.optim.Adam(lr=cfg.learning_rate),
        alpha=0.2,
        sigma=0.1,
        duration_index=labtrans.cuts,
    )
    model.fit(
        x_train.astype("float32").values,
        y_train,
        batch_size=cfg.batch_size,
        epochs=cfg.epochs,
        callbacks=[tt.callbacks.EarlyStopping(patience=cfg.patience)],
        verbose=False,
        val_data=(x_valid.astype("float32").values, y_valid),
    )
    return model


def risk_from_survival(surv_df: pd.DataFrame, horizon: int) -> np.ndarray:
    idx = np.abs(surv_df.index.values.astype(float) - horizon).argmin()
    return 1.0 - surv_df.iloc[idx].values


def horizon_metrics(
    test: pd.DataFrame,
    predicted_risk: np.ndarray,
    event_col: str,
    time_col: str,
    horizon: int,
) -> Dict[str, object]:
    time = safe_numeric(test[time_col])
    event = ensure_binary(test[event_col]).fillna(0)

    is_event_by_horizon = (event == 1) & (time <= horizon)
    eligible = is_event_by_horizon | (time >= horizon)

    if eligible.sum() < 100 or is_event_by_horizon.sum() < 5:
        return {
            "auc": np.nan,
            "brier": np.nan,
            "n_eval": int(eligible.sum()),
            "n_event": int(is_event_by_horizon.sum()),
        }

    y = is_event_by_horizon.loc[eligible].astype(int).values
    p = np.clip(predicted_risk[eligible.values], 1e-6, 1 - 1e-6)
    return {
        "auc": float(roc_auc_score(y, p)) if len(np.unique(y)) > 1 else np.nan,
        "brier": float(brier_score_loss(y, p)),
        "n_eval": int(eligible.sum()),
        "n_event": int(is_event_by_horizon.sum()),
    }


def bootstrap_risk_contrast(
    risk_exposed: np.ndarray,
    risk_unexposed: np.ndarray,
    n_boot: int,
    seed: int,
) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    n = len(risk_exposed)
    risk_difference = np.empty(n_boot)
    risk_ratio = np.empty(n_boot)

    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        mean_1 = np.mean(risk_exposed[idx])
        mean_0 = np.mean(risk_unexposed[idx])
        risk_difference[b] = mean_1 - mean_0
        risk_ratio[b] = np.nan if mean_0 <= 0 else mean_1 / mean_0

    risk_ratio_valid = risk_ratio[np.isfinite(risk_ratio)]
    return {
        "risk_difference_ci_low": float(np.quantile(risk_difference, 0.025)),
        "risk_difference_ci_high": float(np.quantile(risk_difference, 0.975)),
        "p_value_rd": float(np.clip(2 * min(np.mean(risk_difference <= 0), np.mean(risk_difference >= 0)), 0, 1)),
        "risk_ratio_ci_low": float(np.quantile(risk_ratio_valid, 0.025)) if len(risk_ratio_valid) > 0 else np.nan,
        "risk_ratio_ci_high": float(np.quantile(risk_ratio_valid, 0.975)) if len(risk_ratio_valid) > 0 else np.nan,
        "p_value_rr": float(np.clip(2 * min(np.mean(risk_ratio_valid <= 1), np.mean(risk_ratio_valid >= 1)), 0, 1))
        if len(risk_ratio_valid) > 0
        else np.nan,
    }


def run_deephit(
    df: pd.DataFrame,
    exp_col: str,
    event_col: str,
    time_col: str,
    covariates: Sequence[str],
    outcome_name: str,
    cfg: AnalysisConfig,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ensure_deephit_available()

    keep = list(dict.fromkeys(["eid", event_col, time_col, exp_col] + [c for c in covariates if c in df.columns]))
    z = df[keep].copy()
    for col in z.columns:
        if col == "eid":
            continue
        z[col] = safe_numeric(z[col])
        median_value = z[col].median(skipna=True)
        if pd.isna(median_value):
            median_value = 0.0
        z[col] = z[col].fillna(median_value)

    z = z.loc[(z[time_col] > 0) & z[event_col].isin([0, 1]) & z[exp_col].isin([0, 1])].copy()
    if not passes_minimum_counts(z, exp_col, event_col, cfg):
        raise ValueError("Invalid DeepHit data after filtering")

    z = sample_for_deephit(z, event_col, exp_col, cfg.deephit_max_n, cfg.random_seed)
    train, valid, test, x_train, x_valid, x_test = split_and_scale_survival(
        z, event_col, time_col, cfg.random_seed
    )

    model = fit_deephit_model(x_train, train, x_valid, valid, event_col, time_col, cfg)
    surv_obs = model.predict_surv_df(x_test.astype("float32").values)
    concordance_td = float(
        EvalSurv(
            surv_obs,
            test[time_col].values,
            test[event_col].values,
            censor_surv="km",
        ).concordance_td("antolini")
    )

    x_exposed = x_test.copy()
    x_unexposed = x_test.copy()
    x_exposed[exp_col] = 1
    x_unexposed[exp_col] = 0

    surv_exposed = model.predict_surv_df(x_exposed.astype("float32").values)
    surv_unexposed = model.predict_surv_df(x_unexposed.astype("float32").values)

    assoc_rows = []
    metrics_rows = []
    trajectory_rows = []

    for horizon in cfg.horizons:
        risk_obs = risk_from_survival(surv_obs, horizon)
        risk_1 = risk_from_survival(surv_exposed, horizon)
        risk_0 = risk_from_survival(surv_unexposed, horizon)

        metrics = horizon_metrics(test, risk_obs, event_col, time_col, horizon)
        boot = bootstrap_risk_contrast(
            risk_1,
            risk_0,
            n_boot=cfg.n_boot,
            seed=cfg.random_seed + int(horizon * 100),
        )

        mean_1 = float(np.mean(risk_1))
        mean_0 = float(np.mean(risk_0))
        assoc_rows.append(
            {
                "outcome": outcome_name,
                "exposure": exp_col,
                "exposure_label": exposure_label(exp_col),
                "horizon_years": horizon,
                "risk_exposed": mean_1,
                "risk_unexposed": mean_0,
                "risk_difference": mean_1 - mean_0,
                "risk_ratio": mean_1 / mean_0 if mean_0 > 0 else np.nan,
                **boot,
                "concordance_td": concordance_td,
                "n_training_source": int(len(z)),
                "method": "DeepHit_model_based_risk_contrast",
            }
        )
        metrics_rows.append(
            {
                "outcome": outcome_name,
                "exposure": exp_col,
                "horizon_years": horizon,
                "n_eval": metrics["n_eval"],
                "n_event": metrics["n_event"],
                "auc": metrics["auc"],
                "brier": metrics["brier"],
                "concordance_td": concordance_td,
                "method": "DeepHit",
            }
        )

    for horizon in cfg.trajectory_grid:
        risk_1 = risk_from_survival(surv_exposed, int(horizon))
        risk_0 = risk_from_survival(surv_unexposed, int(horizon))
        mean_1 = float(np.mean(risk_1))
        mean_0 = float(np.mean(risk_0))
        trajectory_rows.append(
            {
                "outcome": outcome_name,
                "exposure": exp_col,
                "exposure_label": exposure_label(exp_col),
                "horizon_years": int(horizon),
                "risk_exposed": mean_1,
                "risk_unexposed": mean_0,
                "risk_difference": mean_1 - mean_0,
                "risk_ratio": mean_1 / mean_0 if mean_0 > 0 else np.nan,
            }
        )

    subject_risk = test[["eid", time_col, event_col, exp_col]].copy()
    for horizon in cfg.horizons:
        subject_risk[f"risk_observed_{horizon}y"] = risk_from_survival(surv_obs, horizon)
        subject_risk[f"risk_cf_exposed_{horizon}y"] = risk_from_survival(surv_exposed, horizon)
        subject_risk[f"risk_cf_unexposed_{horizon}y"] = risk_from_survival(surv_unexposed, horizon)

    return (
        pd.DataFrame(assoc_rows),
        pd.DataFrame(metrics_rows),
        subject_risk,
        pd.DataFrame(trajectory_rows),
    )


# -----------------------------------------------------------------------------
# Analysis runners
# -----------------------------------------------------------------------------


def available_survival_outcomes(df: pd.DataFrame, cfg: AnalysisConfig) -> Dict[str, Dict[str, str]]:
    out = {}
    for key, spec in SURVIVAL_OUTCOMES.items():
        event_col = spec["event_col"]
        time_col = spec["time_col"]
        if event_col not in df.columns or time_col not in df.columns:
            continue
        n_event = ensure_binary(df[event_col]).sum(skipna=True)
        n_time = safe_numeric(df[time_col]).notna().sum()
        if n_event >= cfg.min_events and n_time >= cfg.min_total_n:
            out[key] = spec
    return out


def disease_level_exposures(df: pd.DataFrame) -> List[str]:
    exposures = [col for col in df.columns if col.startswith("exp_k")]

    def sort_key(x: str) -> int:
        suffix = x.replace("exp_k", "")
        return int(suffix) if suffix.isdigit() else 999

    return sorted(exposures, key=sort_key)


def exposure_lists(df: pd.DataFrame) -> Dict[str, List[str]]:
    return {
        "main_grouped": [e for e in MAIN_GROUPED_EXPOSURES if e in df.columns],
        "secondary_grouped": [e for e in SECONDARY_GROUPED_EXPOSURES if e in df.columns],
        "supplementary_grouped": [e for e in SUPPLEMENTARY_GROUPED_EXPOSURES if e in df.columns],
        "disease": disease_level_exposures(df),
    }


def run_cox_grid(
    df: pd.DataFrame,
    outcome_name: str,
    exposure_groups: Mapping[str, Sequence[str]],
    blocks: Mapping[str, Sequence[str]],
    analysis_name: str,
    cfg: AnalysisConfig,
    model_names: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    spec = SURVIVAL_OUTCOMES[outcome_name]
    event_col = spec["event_col"]
    time_col = spec["time_col"]

    selected_models = list(model_names) if model_names is not None else list(blocks.keys())
    max_covariates = list(dict.fromkeys([c for m in selected_models for c in blocks[m]]))
    rows = []

    for exposure_group, exposures in exposure_groups.items():
        for exp_col in exposures:
            dd = survival_subset(df, exp_col, event_col, time_col, max_covariates, cfg)
            if dd is None:
                continue
            logging.info(
                "[COX] %s/%s/%s/%s: n=%d, events=%d",
                analysis_name,
                outcome_name,
                exposure_group,
                exp_col,
                len(dd),
                int(dd[event_col].sum()),
            )
            for model_name in selected_models:
                try:
                    rows.append(
                        run_cox(
                            dd,
                            exp_col,
                            event_col,
                            time_col,
                            blocks[model_name],
                            model_name,
                            outcome_name,
                            analysis_name,
                            exposure_group,
                            cfg,
                        )
                    )
                except Exception as exc:
                    logging.warning(
                        "Cox failed: %s/%s/%s/%s | %s",
                        analysis_name,
                        outcome_name,
                        exp_col,
                        model_name,
                        exc,
                    )

    return add_fdr(
        pd.DataFrame(rows),
        ["analysis_name", "outcome", "model_name", "exposure_group"],
    )


def run_logistic_grid(
    df: pd.DataFrame,
    outcome_name: str,
    exposure_groups: Mapping[str, Sequence[str]],
    blocks: Mapping[str, Sequence[str]],
    cfg: AnalysisConfig,
) -> pd.DataFrame:
    event_col = LOGISTIC_OUTCOMES[outcome_name]["event_col"]
    max_covariates = list(dict.fromkeys([c for covs in blocks.values() for c in covs]))
    rows = []

    for exposure_group, exposures in exposure_groups.items():
        for exp_col in exposures:
            dd = logistic_subset(df, exp_col, event_col, max_covariates, cfg)
            if dd is None:
                continue
            for model_name, covariates in blocks.items():
                try:
                    rows.append(
                        run_logistic(
                            dd,
                            exp_col,
                            event_col,
                            covariates,
                            model_name,
                            outcome_name,
                            "main",
                            exposure_group,
                        )
                    )
                except Exception as exc:
                    logging.warning(
                        "Logistic failed: %s/%s/%s | %s",
                        outcome_name,
                        exp_col,
                        model_name,
                        exc,
                    )

    return pd.DataFrame(rows)


def reverse_lag_filter(df: pd.DataFrame, event_col: str, time_col: str, lag_years: int) -> pd.DataFrame:
    event = ensure_binary(df[event_col]).fillna(0)
    time = safe_numeric(df[time_col])
    return df.loc[~((event == 1) & (time <= lag_years))].copy()


def exposure_recency_filter(df: pd.DataFrame, exp_col: str, recency_years: int) -> pd.DataFrame:
    dur_col = duration_col_from_exposure(exp_col)
    if dur_col not in df.columns:
        return df.copy()
    exposure = safe_numeric(df[exp_col])
    duration = safe_numeric(df[dur_col])
    return df.loc[~((exposure == 1) & duration.notna() & (duration < recency_years))].copy()


def age_group(df: pd.DataFrame) -> pd.Series:
    if "age_at_baseline_raw" in df.columns:
        age = safe_numeric(df["age_at_baseline_raw"])
    elif "age_at_baseline" in df.columns:
        age = safe_numeric(df["age_at_baseline"])
    else:
        return pd.Series(np.nan, index=df.index, dtype="object")

    out = pd.Series(np.nan, index=df.index, dtype="object")
    out.loc[age < 60] = "<60"
    out.loc[(age >= 60) & (age < 70)] = "60-69"
    out.loc[age >= 70] = ">=70"
    return out


def subgroup_masks(df: pd.DataFrame) -> Dict[str, Dict[str, pd.Series]]:
    masks: Dict[str, Dict[str, pd.Series]] = {}

    if "sex_model" in df.columns:
        sex = safe_numeric(df["sex_model"])
        masks["sex"] = {"female": sex == 0, "male": sex == 1}

    if "apoe_e4_carrier_model" in df.columns:
        apoe = ensure_binary(df["apoe_e4_carrier_model"]).fillna(0)
        masks["apoe4"] = {"noncarrier": apoe == 0, "carrier": apoe == 1}

    if "parental_dementia_any" in df.columns:
        fam = ensure_binary(df["parental_dementia_any"]).fillna(0)
        masks["parental_dementia"] = {
            "no_parental_dementia": fam == 0,
            "parental_dementia": fam == 1,
        }

    ag = age_group(df)
    if ag.notna().any():
        masks["age_group"] = {"<60": ag == "<60", "60-69": ag == "60-69", ">=70": ag == ">=70"}

    binary_map = {
        "diabetes": "hx_diabetes",
        "depression_anxiety": "hx_depression_anxiety",
        "obesity": "obesity",
        "hyperlipidemia": "hx_hyperlipidemia",
        "hypertension": "hx_hypertension",
        "chd_cvd": "hx_chd_cvd",
        "stroke_tia": "hx_stroke_tia",
        "sleep_disorder": "hx_sleep_disorder",
        "bowel_screening": "bowel_cancer_screening",
    }
    for subgroup_name, col in binary_map.items():
        if col in df.columns:
            x = ensure_binary(df[col]).fillna(0)
            masks[subgroup_name] = {f"{subgroup_name}_0": x == 0, f"{subgroup_name}_1": x == 1}

    if "n_grouped_gi" in df.columns:
        n_gi = safe_numeric(df["n_grouped_gi"])
        masks["gi_multimorbidity"] = {"single_gi": n_gi == 1, "multi_gi": n_gi >= 2}

    return masks


def interaction_cox(
    df: pd.DataFrame,
    exp_col: str,
    modifier_col: str,
    event_col: str,
    time_col: str,
    covariates: Sequence[str],
    outcome_name: str,
    cfg: AnalysisConfig,
) -> Optional[Dict[str, object]]:
    if modifier_col not in df.columns:
        return None

    dd = survival_subset(df, exp_col, event_col, time_col, list(covariates) + [modifier_col], cfg)
    if dd is None:
        return None

    dd[modifier_col] = ensure_binary(dd[modifier_col]).fillna(0)
    if dd[modifier_col].nunique() < 2:
        return None

    interaction_col = f"{exp_col}_X_{modifier_col}"
    dd[interaction_col] = safe_numeric(dd[exp_col]) * safe_numeric(dd[modifier_col])

    model_covariates = [
        c for c in covariates if c in dd.columns and c not in [exp_col, modifier_col, interaction_col]
    ]
    model_df = dd[[time_col, event_col, exp_col, modifier_col, interaction_col] + model_covariates].copy()
    for col in model_df.columns:
        model_df[col] = safe_numeric(model_df[col])

    model_df = model_df.rename(columns={time_col: "time_years", event_col: "event"})
    model_df = clean_model_frame(
        model_df,
        event_col="event",
        duration_col="time_years",
        protected_cols=["time_years", "event", exp_col, modifier_col, interaction_col],
    )

    if interaction_col not in model_df.columns or model_df[interaction_col].nunique() <= 1:
        return None

    cph = CoxPHFitter(penalizer=cfg.cox_penalizer)
    cph.fit(model_df, duration_col="time_years", event_col="event")
    summary = cph.summary.loc[interaction_col]

    return {
        "outcome": outcome_name,
        "exposure": exp_col,
        "exposure_label": exposure_label(exp_col),
        "modifier": modifier_col,
        "interaction_term": interaction_col,
        "n_total": int(len(model_df)),
        "n_event": int(model_df["event"].sum()),
        "interaction_hr": float(np.exp(summary["coef"])),
        "ci_low": float(np.exp(summary["coef lower 95%"])),
        "ci_high": float(np.exp(summary["coef upper 95%"])),
        "p_value": float(summary["p"]),
        "method": "Cox_interaction",
    }


# -----------------------------------------------------------------------------
# Output helpers
# -----------------------------------------------------------------------------


def make_output_dirs(output_dir: Path) -> Dict[str, Path]:
    dirs = {
        "main": output_dir / "main",
        "sensitivity": output_dir / "sensitivity",
        "subgroup": output_dir / "subgroup",
        "deephit": output_dir / "deephit",
        "qc": output_dir / "qc",
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def write_qc_tables(
    df: pd.DataFrame,
    outcomes: Mapping[str, Mapping[str, str]],
    exposures: Mapping[str, Sequence[str]],
    blocks: Mapping[str, Sequence[str]],
    output_dir: Path,
) -> None:
    qc_dir = output_dir / "qc"

    pd.DataFrame(
        [
            {"model_name": name, "n_covariates": len(covs), "covariates": ";".join(covs)}
            for name, covs in blocks.items()
        ]
    ).to_csv(qc_dir / "model_blocks.csv", index=False)

    outcome_rows = []
    for key, spec in SURVIVAL_OUTCOMES.items():
        event_col = spec["event_col"]
        time_col = spec["time_col"]
        outcome_rows.append(
            {
                "outcome": key,
                "available": key in outcomes,
                "event_col": event_col,
                "time_col": time_col,
                "n_event": int(ensure_binary(df[event_col]).sum(skipna=True)) if event_col in df else np.nan,
                "n_time_nonmissing": int(safe_numeric(df[time_col]).notna().sum()) if time_col in df else np.nan,
                "role": spec.get("role", ""),
            }
        )
    pd.DataFrame(outcome_rows).to_csv(qc_dir / "outcome_qc.csv", index=False)

    exposure_rows = []
    for exposure_group, exp_cols in exposures.items():
        for exp_col in exp_cols:
            if exp_col in df.columns:
                x = safe_numeric(df[exp_col])
                exposure_rows.append(
                    {
                        "exposure_group": exposure_group,
                        "exposure": exp_col,
                        "label": exposure_label(exp_col),
                        "available": True,
                        "n_exposed": int((x == 1).sum()),
                        "n_unexposed": int((x == 0).sum()),
                        "missing_rate": float(x.isna().mean()),
                    }
                )
            else:
                exposure_rows.append(
                    {
                        "exposure_group": exposure_group,
                        "exposure": exp_col,
                        "label": exposure_label(exp_col),
                        "available": False,
                        "n_exposed": np.nan,
                        "n_unexposed": np.nan,
                        "missing_rate": np.nan,
                    }
                )
    pd.DataFrame(exposure_rows).to_csv(qc_dir / "exposure_qc.csv", index=False)


def write_metadata(
    cfg: AnalysisConfig,
    df: pd.DataFrame,
    outcomes: Mapping[str, Mapping[str, str]],
    exposures: Mapping[str, Sequence[str]],
    blocks: Mapping[str, Sequence[str]],
) -> None:
    metadata = asdict(cfg)
    metadata["input_file"] = str(cfg.input_file)
    metadata["output_dir"] = str(cfg.output_dir)
    metadata["n_total"] = int(len(df))
    metadata["survival_outcomes_available"] = list(outcomes.keys())
    metadata["exposure_lists"] = {key: list(value) for key, value in exposures.items()}
    metadata["model_blocks"] = {key: list(value) for key, value in blocks.items()}
    metadata["key_implementation_notes"] = [
        "Exposure-specific model subsets retain all requested covariates before Cox, logistic and DeepHit fitting.",
        "Both date_* and dt_* exposure-date naming conventions are supported.",
        "Precomputed mci_ad_event and mci_ad_time_years are reused if available.",
        "MCI-or-AD composite construction preserves censoring time for non-events when constructed internally.",
        "DeepHit is optional and is treated as supplementary absolute-risk modelling, not primary etiological inference.",
        "Cox models use a small L2 penalizer for numerical stability; the penalizer value is recorded in outputs.",
    ]

    with open(cfg.output_dir / "run_metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2, default=str)


def write_report(
    cfg: AnalysisConfig,
    df: pd.DataFrame,
    outcomes: Mapping[str, Mapping[str, str]],
    exposures: Mapping[str, Sequence[str]],
    blocks: Mapping[str, Sequence[str]],
) -> None:
    lines = [
        "GI-AD association analysis completed.",
        "",
        f"Input file: {cfg.input_file}",
        f"Output directory: {cfg.output_dir}",
        f"Analysis N: {len(df)}",
        f"Available survival outcomes: {list(outcomes.keys())}",
        f"Run DeepHit: {cfg.run_deephit}",
        "",
        "Main analyses:",
        " - Staged Cox models from age/sex adjustment to family-genetic full adjustment.",
        " - Main, secondary, supplementary grouped GI exposures and disease-level exposures.",
        " - Baseline broad MCI/AD logistic support analysis if the outcome is available.",
        "",
        "Sensitivity analyses:",
        " - Reverse-causation lag exclusion.",
        " - Recent-exposure exclusion.",
        " - Bowel cancer screening adjustment/stratification if available.",
        "",
        "Subgroup and interaction analyses:",
        " - Sex, age, APOE4, parental dementia, clinical/lifestyle subgroups and selected interactions.",
        "",
        "Covariate blocks:",
    ]
    for name, covariates in blocks.items():
        lines.append(f" - {name}: {len(covariates)} covariates")

    with open(cfg.output_dir / "analysis_report.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------


def run_analysis(cfg: AnalysisConfig) -> None:
    setup_logging(cfg.output_dir)
    set_random_seed(cfg.random_seed)
    output_dirs = make_output_dirs(cfg.output_dir)

    logging.info("Loading input file: %s", cfg.input_file)
    raw_df = load_data(cfg.input_file)
    df = prepare_dataset(raw_df)

    blocks = model_blocks(df)
    outcomes = available_survival_outcomes(df, cfg)
    exposures = exposure_lists(df)
    grouped_available = [
        e
        for e in MAIN_GROUPED_EXPOSURES + SECONDARY_GROUPED_EXPOSURES + SUPPLEMENTARY_GROUPED_EXPOSURES
        if e in df.columns
    ]

    logging.info("N=%d", len(df))
    logging.info("Available outcomes=%s", list(outcomes.keys()))
    logging.info("Available grouped exposures=%s", grouped_available)
    logging.info("Disease-level exposures=%d", len(exposures.get("disease", [])))
    logging.info("Model blocks=%s", {k: len(v) for k, v in blocks.items()})

    write_qc_tables(df, outcomes, exposures, blocks, cfg.output_dir)

    # Main Cox analyses
    all_cox = []
    for outcome_name in outcomes:
        result = run_cox_grid(df, outcome_name, exposures, blocks, "main", cfg)
        if not result.empty:
            result.to_csv(output_dirs["main"] / f"main_{outcome_name}_cox_results.csv", index=False)
            all_cox.append(result)

    all_cox_df = safe_concat(all_cox)
    if not all_cox_df.empty:
        all_cox_df.to_csv(output_dirs["main"] / "main_cox_all.csv", index=False)
        all_cox_df.loc[all_cox_df["model_name"] == "Model5_family_genetic"].to_csv(
            output_dirs["main"] / "main_cox_model5_summary.csv",
            index=False,
        )

    # Baseline logistic support analyses
    logistic_results = []
    for outcome_name, spec in LOGISTIC_OUTCOMES.items():
        event_col = spec["event_col"]
        if event_col in df.columns and ensure_binary(df[event_col]).sum(skipna=True) >= cfg.min_events:
            result = run_logistic_grid(df, outcome_name, exposures, blocks, cfg)
            if not result.empty:
                result.to_csv(output_dirs["main"] / f"main_{outcome_name}_logistic_results.csv", index=False)
                logistic_results.append(result)

    all_logistic_df = safe_concat(logistic_results)
    if not all_logistic_df.empty:
        all_logistic_df.to_csv(output_dirs["main"] / "main_logistic_all.csv", index=False)

    # DeepHit supplementary analyses
    if cfg.run_deephit:
        ensure_deephit_available()
        deephit_assoc = []
        deephit_metrics = []
        deephit_trajectories = []
        deephit_covariates = blocks.get("Model5_family_genetic", list(blocks.values())[-1])

        for outcome_name in [x for x in DEEPHIT_OUTCOMES if x in outcomes]:
            event_col = outcomes[outcome_name]["event_col"]
            time_col = outcomes[outcome_name]["time_col"]
            for exp_col in [e for e in DEEPHIT_EXPOSURES if e in df.columns]:
                dd = survival_subset(df, exp_col, event_col, time_col, deephit_covariates, cfg)
                if dd is None:
                    continue
                try:
                    logging.info("[DeepHit] %s/%s", outcome_name, exp_col)
                    assoc, metrics, subject_risk, traj = run_deephit(
                        dd, exp_col, event_col, time_col, deephit_covariates, outcome_name, cfg
                    )
                    assoc["analysis_name"] = "main"
                    metrics["analysis_name"] = "main"
                    traj["analysis_name"] = "main"
                    deephit_assoc.append(assoc)
                    deephit_metrics.append(metrics)
                    deephit_trajectories.append(traj)
                    subject_risk.to_csv(
                        output_dirs["deephit"] / f"deephit_subject_risk_{outcome_name}_{exp_col}.csv",
                        index=False,
                    )
                except Exception as exc:
                    logging.warning("DeepHit failed: %s/%s | %s", outcome_name, exp_col, exc)

        if deephit_assoc:
            safe_concat(deephit_assoc).to_csv(output_dirs["deephit"] / "deephit_risk_contrast_all.csv", index=False)
        if deephit_metrics:
            safe_concat(deephit_metrics).to_csv(output_dirs["deephit"] / "deephit_performance_all.csv", index=False)
        if deephit_trajectories:
            safe_concat(deephit_trajectories).to_csv(output_dirs["deephit"] / "deephit_trajectories_all.csv", index=False)

    # Reverse-causation sensitivity: Model5 only
    reverse_results = []
    for lag_years in cfg.reverse_causation_lags:
        for outcome_name, spec in outcomes.items():
            filtered = reverse_lag_filter(df, spec["event_col"], spec["time_col"], lag_years)
            result = run_cox_grid(
                filtered,
                outcome_name,
                exposures,
                blocks,
                f"reverse_lag_{lag_years}y",
                cfg,
                model_names=["Model5_family_genetic"],
            )
            reverse_results.append(result)

    reverse_df = safe_concat(reverse_results)
    if not reverse_df.empty:
        reverse_df.to_csv(output_dirs["sensitivity"] / "reverse_causation_cox_all.csv", index=False)

    # Exposure-recency sensitivity: Model5 only
    recency_results = []
    for recency_years in cfg.exposure_recency_years:
        for outcome_name in outcomes:
            for exposure_group, exp_cols in exposures.items():
                for exp_col in exp_cols:
                    filtered = exposure_recency_filter(df, exp_col, recency_years)
                    result = run_cox_grid(
                        filtered,
                        outcome_name,
                        {exposure_group: [exp_col]},
                        blocks,
                        f"exposure_recency_excl_{recency_years}y",
                        cfg,
                        model_names=["Model5_family_genetic"],
                    )
                    recency_results.append(result)

    recency_df = safe_concat(recency_results)
    if not recency_df.empty:
        recency_df.to_csv(output_dirs["sensitivity"] / "exposure_recency_cox_all.csv", index=False)

    # Detection-bias sensitivity
    if "bowel_cancer_screening" in df.columns and grouped_available:
        detection_results = []
        grouped_map = {"grouped": grouped_available}
        for outcome_name in outcomes:
            detection_results.append(
                run_cox_grid(
                    df,
                    outcome_name,
                    grouped_map,
                    blocks,
                    "detection_bias_adjusted",
                    cfg,
                    model_names=["Model6_detection_bias"],
                )
            )
            for level, mask in {
                "screening_0": df["bowel_cancer_screening"] == 0,
                "screening_1": df["bowel_cancer_screening"] == 1,
            }.items():
                if mask.sum() >= cfg.min_total_n:
                    detection_results.append(
                        run_cox_grid(
                            df.loc[mask].copy(),
                            outcome_name,
                            grouped_map,
                            blocks,
                            f"detection_bias_stratified_{level}",
                            cfg,
                            model_names=["Model5_family_genetic"],
                        )
                    )
        detection_df = safe_concat(detection_results)
        if not detection_df.empty:
            detection_df.to_csv(
                output_dirs["sensitivity"] / "detection_bias_sensitivity_cox_all.csv",
                index=False,
            )

    # Subgroup and interaction analyses for strict AD
    subgroup_rows = []
    interaction_rows = []
    if "strict_ad" in outcomes and grouped_available:
        event_col = outcomes["strict_ad"]["event_col"]
        time_col = outcomes["strict_ad"]["time_col"]
        full_covariates = blocks["Model5_family_genetic"]
        masks = subgroup_masks(df)

        for subgroup_name, levels in masks.items():
            for level_name, mask in levels.items():
                if mask.sum() < cfg.min_total_n:
                    continue
                for exp_col in grouped_available:
                    dd = survival_subset(
                        df.loc[mask].copy(),
                        exp_col,
                        event_col,
                        time_col,
                        full_covariates,
                        cfg,
                    )
                    if dd is None:
                        continue
                    try:
                        result = run_cox(
                            dd,
                            exp_col,
                            event_col,
                            time_col,
                            full_covariates,
                            "Model5_family_genetic",
                            "strict_ad",
                            "subgroup",
                            "grouped",
                            cfg,
                        )
                        result.update({"subgroup": subgroup_name, "subgroup_level": level_name})
                        subgroup_rows.append(result)
                    except Exception as exc:
                        logging.warning(
                            "Subgroup Cox failed: %s/%s/%s | %s",
                            subgroup_name,
                            level_name,
                            exp_col,
                            exc,
                        )

        modifier_map = {
            "apoe4": "apoe_e4_carrier_model",
            "parental_dementia": "parental_dementia_any",
            "sex": "sex_model",
            "diabetes": "hx_diabetes",
            "depression_anxiety": "hx_depression_anxiety",
            "sleep_disorder": "hx_sleep_disorder",
            "bowel_screening": "bowel_cancer_screening",
        }
        for modifier_name, modifier_col in modifier_map.items():
            if modifier_col not in df.columns:
                continue
            for exp_col in grouped_available:
                try:
                    result = interaction_cox(
                        df,
                        exp_col,
                        modifier_col,
                        event_col,
                        time_col,
                        full_covariates,
                        "strict_ad",
                        cfg,
                    )
                    if result is not None:
                        result["modifier_name"] = modifier_name
                        interaction_rows.append(result)
                except Exception as exc:
                    logging.warning("Interaction failed: %s x %s | %s", exp_col, modifier_col, exc)

    subgroup_df = add_fdr(pd.DataFrame(subgroup_rows), ["subgroup", "subgroup_level"])
    if not subgroup_df.empty:
        subgroup_df.to_csv(output_dirs["subgroup"] / "subgroup_strict_ad_grouped_cox.csv", index=False)

    interaction_df = add_fdr(pd.DataFrame(interaction_rows), ["modifier_name"])
    if not interaction_df.empty:
        interaction_df.to_csv(output_dirs["subgroup"] / "interaction_strict_ad_grouped_cox.csv", index=False)

    write_metadata(cfg, df, outcomes, exposures, blocks)
    write_report(cfg, df, outcomes, exposures, blocks)
    logging.info("Analysis completed. Read first: %s", cfg.output_dir / "analysis_report.txt")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def parse_args() -> AnalysisConfig:
    parser = argparse.ArgumentParser(
        description="GI disease exposure and Alzheimer-related outcome association analysis."
    )
    parser.add_argument("--input", required=True, type=Path, help="Path to preprocessed participant-level CSV.")
    parser.add_argument("--outdir", required=True, type=Path, help="Directory for analysis outputs.")
    parser.add_argument("--seed", type=int, default=20260509, help="Random seed.")
    parser.add_argument("--min-total-n", type=int, default=200, help="Minimum sample size for each model.")
    parser.add_argument("--min-events", type=int, default=40, help="Minimum event count for each model.")
    parser.add_argument("--min-exposed", type=int, default=30, help="Minimum exposed count for each model.")
    parser.add_argument("--min-unexposed", type=int, default=30, help="Minimum unexposed count for each model.")
    parser.add_argument("--cox-penalizer", type=float, default=0.01, help="L2 penalizer for CoxPHFitter.")

    deephit = parser.add_mutually_exclusive_group()
    deephit.add_argument("--run-deephit", action="store_true", help="Run supplementary DeepHit analyses.")
    deephit.add_argument("--skip-deephit", action="store_true", help="Skip supplementary DeepHit analyses.")

    parser.add_argument("--deephit-max-n", type=int, default=120000, help="Maximum sample size for each DeepHit fit. Use -1 for full sample.")
    parser.add_argument("--num-durations", type=int, default=50, help="Number of duration bins for DeepHit.")
    parser.add_argument("--batch-size", type=int, default=512, help="DeepHit batch size.")
    parser.add_argument("--epochs", type=int, default=120, help="DeepHit maximum epochs.")
    parser.add_argument("--patience", type=int, default=15, help="DeepHit early-stopping patience.")
    parser.add_argument("--learning-rate", type=float, default=1e-3, help="DeepHit learning rate.")
    parser.add_argument("--horizons", type=int, nargs="+", default=[5, 10, 15, 20], help="Risk horizons in years.")
    parser.add_argument("--trajectory-start", type=int, default=1, help="Start year for risk trajectory.")
    parser.add_argument("--trajectory-end", type=int, default=20, help="End year for risk trajectory.")
    parser.add_argument("--n-boot", type=int, default=200, help="Bootstrap iterations for DeepHit risk contrasts.")
    parser.add_argument("--reverse-lags", type=int, nargs="+", default=[1, 2, 5], help="Reverse-causation lag exclusions in years.")
    parser.add_argument("--recency-years", type=int, nargs="+", default=[1, 2, 5], help="Recent-exposure exclusions in years.")

    args = parser.parse_args()
    run_deephit = bool(args.run_deephit) and not bool(args.skip_deephit)
    deephit_max_n = None if args.deephit_max_n == -1 else args.deephit_max_n

    return AnalysisConfig(
        input_file=args.input,
        output_dir=args.outdir,
        random_seed=args.seed,
        min_total_n=args.min_total_n,
        min_events=args.min_events,
        min_exposed=args.min_exposed,
        min_unexposed=args.min_unexposed,
        cox_penalizer=args.cox_penalizer,
        run_deephit=run_deephit,
        deephit_max_n=deephit_max_n,
        num_durations=args.num_durations,
        batch_size=args.batch_size,
        epochs=args.epochs,
        patience=args.patience,
        learning_rate=args.learning_rate,
        horizons=tuple(args.horizons),
        trajectory_grid_start=args.trajectory_start,
        trajectory_grid_end=args.trajectory_end,
        n_boot=args.n_boot,
        reverse_causation_lags=tuple(args.reverse_lags),
        exposure_recency_years=tuple(args.recency_years),
    )


def main() -> None:
    cfg = parse_args()
    run_analysis(cfg)


if __name__ == "__main__":
    main()