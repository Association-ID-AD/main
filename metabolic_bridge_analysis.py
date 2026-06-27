#!/usr/bin/env python3
"""
Metabolic bridge analysis for gastrointestinal disease and Alzheimer disease risk.

This script evaluates whether systemic metabolic markers provide pathway-consistent
bridge evidence between gastrointestinal disease exposure and subsequent AD risk.
It performs:

1. grouped gastrointestinal exposure construction;
2. metabolic feature classification into blood-biochemistry and NMR layers;
3. Stage 1 screening: GI exposure -> metabolic marker association by robust OLS;
4. Stage 2 screening: metabolic marker -> AD/dementia risk by CoxPH;
5. horizon-specific logistic support analyses;
6. restricted g-computation bridge decomposition for selected markers;
7. PCA-based metabolic module bridge analysis;
8. QC tables, summary tables, figures, metadata, and a text report.

The bridge analysis is intended as pathway-consistent evidence, not definitive
causal mediation.

Example
-------
python metabolic_bridge_analysis.py \
    --master data/master_preprocessed_with_genetics.csv \
    --metabolic data/metabolic_prepared.csv \
    --outdir results/metabolic_bridge_analysis \
    --label-map data/metabolic_feature_labels.csv

For a faster test run:
python metabolic_bridge_analysis.py \
    --master data/master_preprocessed_with_genetics.csv \
    --metabolic data/metabolic_prepared.csv \
    --outdir results/metabolic_bridge_analysis_test \
    --max-features-per-layer 50 \
    --bootstrap 50
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from lifelines import CoxPHFitter
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression


# -----------------------------------------------------------------------------
# Defaults and constants
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

EXPOSURE_LABELS: Dict[str, str] = {
    "exp_diverticular": "Diverticular disease",
    "exp_other_functional_intestinal": "Other functional intestinal disorders",
    "exp_ibd": "Inflammatory bowel disease",
    "exp_other_chronic_intestinal": "Other chronic intestinal disease",
    "exp_ibs": "Irritable bowel syndrome",
    "exp_appendiceal": "Appendiceal disease",
    "exp_malabsorption": "Malabsorption",
    "exp_anorectal": "Anorectal disease",
}

MAIN_EXPOSURES: List[str] = [
    "exp_diverticular",
    "exp_other_functional_intestinal",
    "exp_ibd",
    "exp_other_chronic_intestinal",
    "exp_ibs",
]
SECONDARY_EXPOSURES: List[str] = ["exp_appendiceal", "exp_malabsorption"]

OUTCOMES: Dict[str, Dict[str, str]] = {
    "strict_ad": {"time": "ad_time_years", "event": "ad_event"},
    "all_cause_dementia": {"time": "dementia_time_years", "event": "dementia_event"},
}


@dataclass
class AnalysisConfig:
    master_path: str
    metabolic_path: str
    outdir: str
    label_map_path: Optional[str] = None
    seed: int = 20260510
    min_total_n: int = 2000
    min_events: int = 80
    min_exposed: int = 50
    min_unexposed: int = 50
    fdr_alpha_stage1: float = 0.05
    fdr_alpha_stage2: float = 0.05
    primary_horizon: int = 10
    support_horizons: Tuple[int, ...] = (5, 10, 15)
    bootstrap: int = 200
    max_markers_for_pca: int = 30
    min_markers_for_pca: int = 2
    cox_penalizer: float = 0.01
    blood_prefix: str = "p30"
    nmr_prefix: str = "p23"
    max_features_per_layer: Optional[int] = None
    make_figures: bool = True


# -----------------------------------------------------------------------------
# Generic utilities
# -----------------------------------------------------------------------------


def configure_logging(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / "run.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_path, mode="w", encoding="utf-8")],
    )


def stable_seed(*parts: object, base_seed: int = 0) -> int:
    text = "::".join(str(x) for x in parts)
    digest = hashlib.md5(text.encode("utf-8")).hexdigest()
    return (int(digest[:8], 16) + int(base_seed)) % (2**32 - 1)


def safe_numeric(x: object) -> pd.Series:
    return pd.to_numeric(x, errors="coerce")


def ensure_binary(x: object) -> pd.Series:
    y = safe_numeric(x)
    return y.where(y.isin([0, 1]), np.nan)


def zscore(x: object) -> pd.Series:
    y = safe_numeric(x)
    mu = y.mean(skipna=True)
    sd = y.std(skipna=True)
    if pd.isna(sd) or sd == 0:
        return pd.Series(np.nan, index=y.index)
    return (y - mu) / sd


def winsorize(x: object, q: Tuple[float, float] = (0.005, 0.995)) -> pd.Series:
    y = safe_numeric(x).copy()
    lo, hi = y.quantile(q[0]), y.quantile(q[1])
    if pd.notna(lo) and pd.notna(hi) and lo < hi:
        y = y.clip(lower=lo, upper=hi)
    return y


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
    qr = ranked * m / (np.arange(m) + 1)
    qr = np.minimum.accumulate(qr[::-1])[::-1]
    q[idx[order]] = np.clip(qr, 0, 1)
    return q


def safe_concat(dfs: Iterable[pd.DataFrame]) -> pd.DataFrame:
    valid = [x for x in dfs if isinstance(x, pd.DataFrame) and len(x) > 0]
    return pd.concat(valid, ignore_index=True) if valid else pd.DataFrame()


def save_table(df: pd.DataFrame, path: Path) -> None:
    if isinstance(df, pd.DataFrame) and len(df) > 0:
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, index=False)


def load_feature_label_map(path: Optional[str]) -> Dict[str, str]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        logging.warning("Label map not found: %s", p)
        return {}
    df = pd.read_csv(p)
    lower = {c.lower(): c for c in df.columns}
    if "feature" not in lower or "label" not in lower:
        logging.warning("Label map must contain columns named feature and label: %s", p)
        return {}
    return dict(zip(df[lower["feature"]].astype(str), df[lower["label"]].astype(str)))


def label_feature(feature: str, label_map: Mapping[str, str]) -> str:
    return label_map.get(feature, feature)


def duration_col(exp: str) -> str:
    return exp.replace("exp_", "dur_") + "_to_baseline"


def date_candidates(exp: str) -> List[str]:
    suffix = exp.replace("exp_", "")
    return [f"date_{suffix}", f"dt_{suffix}", exp.replace("exp_", "date_"), exp.replace("exp_", "dt_")]


def first_existing_col(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    for col in candidates:
        if col in df.columns:
            return col
    return None


def max_binary(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    cols = [c for c in cols if c in df.columns]
    if not cols:
        return pd.Series(np.nan, index=df.index)
    tmp = pd.concat([safe_numeric(df[c]) for c in cols], axis=1)
    out = tmp.max(axis=1, skipna=True)
    out[tmp.isna().all(axis=1)] = np.nan
    return out


def min_numeric(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    cols = [c for c in cols if c in df.columns]
    if not cols:
        return pd.Series(np.nan, index=df.index)
    tmp = pd.concat([safe_numeric(df[c]) for c in cols], axis=1)
    out = tmp.min(axis=1, skipna=True)
    out[tmp.isna().all(axis=1)] = np.nan
    return out


def earliest_date(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    cols = [c for c in cols if c in df.columns]
    if not cols:
        return pd.Series(pd.NaT, index=df.index)
    tmp = pd.concat([pd.to_datetime(df[c], errors="coerce") for c in cols], axis=1)
    return tmp.min(axis=1)


def add_grouped_exposures(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    for grouped_exp, components in GROUPED_GI_MAP.items():
        if grouped_exp not in d.columns:
            d[grouped_exp] = max_binary(d, components)
        dur = duration_col(grouped_exp)
        if dur not in d.columns:
            d[dur] = min_numeric(d, [duration_col(c) for c in components])
        date_col = f"date_{grouped_exp.replace('exp_', '')}"
        if date_col not in d.columns:
            component_dates = []
            for comp in components:
                col = first_existing_col(d, date_candidates(comp))
                if col:
                    component_dates.append(col)
            if component_dates:
                d[date_col] = earliest_date(d, component_dates)
    return d


def parse_horizons(text: str) -> Tuple[int, ...]:
    vals = []
    for part in text.split(","):
        part = part.strip()
        if part:
            vals.append(int(part))
    return tuple(vals)


# -----------------------------------------------------------------------------
# Loading and preparation
# -----------------------------------------------------------------------------


def load_tables(config: AnalysisConfig) -> Tuple[pd.DataFrame, pd.DataFrame]:
    master_path = Path(config.master_path)
    metabolic_path = Path(config.metabolic_path)
    if not master_path.exists():
        raise FileNotFoundError(f"Master file not found: {master_path}")
    if not metabolic_path.exists():
        raise FileNotFoundError(f"Metabolic file not found: {metabolic_path}")

    logging.info("Loading master table: %s", master_path)
    master = pd.read_csv(master_path, low_memory=False)
    logging.info("Loading metabolic table: %s", metabolic_path)
    metabolic = pd.read_csv(metabolic_path, low_memory=False)

    if "eid" not in master.columns or "eid" not in metabolic.columns:
        raise ValueError("Both master and metabolic tables must contain an 'eid' column.")

    for col in master.columns:
        if col in {"time0", "ad_date", "dementia_date", "censor_date", "death_date", "lost_followup_date"} or col.startswith("date_") or col.startswith("dt_"):
            master[col] = pd.to_datetime(master[col], errors="coerce")

    if "time0" in metabolic.columns:
        metabolic = metabolic.drop(columns=["time0"])

    return master, metabolic


def classify_metabolic_feature(feature: str, config: AnalysisConfig) -> str:
    if feature.startswith(config.blood_prefix):
        return "blood_biochemistry"
    if feature.startswith(config.nmr_prefix):
        return "nmr_metabolomics"
    return "other"


def prepare_analysis_dataframe(master: pd.DataFrame, metabolic: pd.DataFrame) -> pd.DataFrame:
    master = add_grouped_exposures(master)
    analysis = master.merge(metabolic, on="eid", how="inner")

    required = ["ad_event", "ad_time_years"]
    missing = [c for c in required if c not in analysis.columns]
    if missing:
        raise ValueError(f"Missing required AD outcome columns: {missing}")

    analysis["ad_event"] = ensure_binary(analysis["ad_event"])
    analysis["ad_time_years"] = safe_numeric(analysis["ad_time_years"])
    analysis = analysis.loc[analysis["ad_time_years"].notna() & (analysis["ad_time_years"] > 0)].copy()

    if "ad_date" in analysis.columns and "time0" in analysis.columns:
        analysis = analysis.loc[analysis["ad_date"].isna() | (analysis["ad_date"] > analysis["time0"])].copy()

    if "dementia_event" in analysis.columns:
        analysis["dementia_event"] = ensure_binary(analysis["dementia_event"])
    if "dementia_time_years" in analysis.columns:
        analysis["dementia_time_years"] = safe_numeric(analysis["dementia_time_years"])

    return analysis


def prepare_covariates(df: pd.DataFrame) -> pd.DataFrame:
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
    for new, old in continuous_aliases.items():
        if new not in d.columns and old in d.columns:
            d[new] = zscore(winsorize(d[old]))

    if "ad_prs_std_z" not in d.columns and "ad_prs_std" in d.columns:
        d["ad_prs_std_z"] = zscore(winsorize(d["ad_prs_std"]))
    if "ad_prs_enh_z" not in d.columns and "ad_prs_enhanced_std" in d.columns:
        d["ad_prs_enh_z"] = zscore(winsorize(d["ad_prs_enhanced_std"]))

    categorical_aliases = {
        "sex_model": "sex",
        "ethnic_model": "ethnicity_5cat",
        "edu_model": "edu_level",
        "health_model": "health_4cat",
        "smoke_status_model": "smoke_status_3cat",
        "apoe_e4_carrier_model": "apoe4_carrier",
    }
    for new, old in categorical_aliases.items():
        if new not in d.columns and old in d.columns:
            d[new] = safe_numeric(d[old])

    binary_cols = [
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
        "med_steroid_immunosuppressive",
        "parental_dementia_any",
        "bowel_cancer_screening",
        "any_current_medication",
        "apoe_e4_carrier_model",
        "sleep_short",
        "sleep_long",
        "insomnia_any",
    ]
    for col in binary_cols:
        if col in d.columns:
            d[col] = ensure_binary(d[col]).fillna(0)

    return d


def build_covariate_blocks(df: pd.DataFrame) -> Dict[str, List[str]]:
    base = [c for c in ["age_z", "sex_model", "ethnic_model", "edu_model", "townsend_z"] if c in df.columns]
    lifestyle = [
        c
        for c in [
            "bmi_z",
            "whr_z",
            "smoke_status_model",
            "smoke_pack_z",
            "alcohol_freq_z",
            "sleep_hours_z",
            "sleep_short",
            "sleep_long",
            "insomnia_any",
            "sedentary_z",
            "activity_z",
            "diet_quality_z",
        ]
        if c in df.columns
    ]
    comorbidity = [
        c
        for c in [
            "hx_diabetes",
            "hx_hypertension",
            "hx_hyperlipidemia",
            "hx_chd_cvd",
            "hx_stroke_tia",
            "hx_depression_anxiety",
            "hx_sleep_disorder",
            "hx_parkinson_other_nd",
        ]
        if c in df.columns
    ]
    medications = [
        c
        for c in [
            "med_statin",
            "med_antihypertensive",
            "med_antidiabetic",
            "med_ppi_gi_drug",
            "med_laxative",
            "med_antidepressant",
            "med_steroid_immunosuppressive",
            "any_current_medication",
        ]
        if c in df.columns
    ]
    family_genetic = [
        c
        for c in [
            "parental_dementia_any",
            "parental_dementia_count_z",
            "apoe_e4_carrier_model",
            "apoe_e4_count_z",
            "ad_prs_std_z",
            "ad_prs_enh_z",
        ]
        if c in df.columns
    ]
    detection = [c for c in ["bowel_cancer_screening"] if c in df.columns]

    blocks = {
        "base": list(dict.fromkeys(base)),
        "lifestyle": list(dict.fromkeys(base + lifestyle)),
        "full": list(dict.fromkeys(base + lifestyle + comorbidity + medications + family_genetic)),
        "detection": list(dict.fromkeys(base + lifestyle + comorbidity + medications + family_genetic + detection)),
    }
    return blocks


def standardize_covariates(df: pd.DataFrame, blocks: Mapping[str, Sequence[str]]) -> pd.DataFrame:
    d = df.copy()
    covariates = list(dict.fromkeys([c for cols in blocks.values() for c in cols]))
    for col in covariates:
        if col not in d.columns:
            continue
        x = safe_numeric(d[col])
        values = set(pd.unique(x.dropna()))
        if len(values) <= 2 and values.issubset({0, 1}):
            d[col] = x.fillna(0)
        else:
            d[col] = zscore(winsorize(x)).fillna(0)
    return d


def get_metabolic_layers(metabolic: pd.DataFrame, config: AnalysisConfig) -> Dict[str, List[str]]:
    features = [c for c in metabolic.columns if c != "eid"]
    layers: Dict[str, List[str]] = {"blood_biochemistry": [], "nmr_metabolomics": [], "other": []}
    for feature in features:
        layers[classify_metabolic_feature(feature, config)].append(feature)

    if config.max_features_per_layer is not None:
        for layer in layers:
            layers[layer] = layers[layer][: config.max_features_per_layer]

    return layers


def clean_frame(df: pd.DataFrame, protected: Sequence[str]) -> pd.DataFrame:
    d = df.copy().replace([np.inf, -np.inf], np.nan)
    protected_set = set(protected)
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


def build_horizon_binary_outcome(
    df: pd.DataFrame, time_col: str, event_col: str, horizon: int
) -> Tuple[pd.Series, pd.Series]:
    time = safe_numeric(df[time_col])
    event = ensure_binary(df[event_col])
    event_by_horizon = (event == 1) & (time <= horizon)
    observed_past_horizon = time >= horizon
    eligible = event_by_horizon | observed_past_horizon
    y = pd.Series(np.nan, index=df.index, dtype=float)
    y.loc[eligible] = 0
    y.loc[event_by_horizon] = 1
    return eligible, y


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------


def draw_volcano(df: pd.DataFrame, title: str, outfile: Path) -> None:
    if df is None or len(df) == 0:
        return
    required = {"beta", "p_value", "feature", "feature_label"}
    if not required.issubset(df.columns):
        return
    dd = df.loc[np.isfinite(df["beta"]) & np.isfinite(df["p_value"])].copy()
    if len(dd) == 0:
        return
    dd["neglog10p"] = -np.log10(dd["p_value"].clip(lower=1e-300))
    plt.figure(figsize=(7.5, 5.5))
    plt.scatter(dd["beta"], dd["neglog10p"], s=12, alpha=0.55)
    plt.axhline(-np.log10(0.05), linestyle="--", linewidth=1)
    plt.xlabel("Beta")
    plt.ylabel("-log10(P)")
    plt.title(title)
    for _, row in dd.sort_values("p_value").head(25).iterrows():
        plt.text(row["beta"], row["neglog10p"], str(row["feature_label"]), fontsize=6)
    plt.tight_layout()
    outfile.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()


def draw_bridge(df: pd.DataFrame, title: str, outfile: Path) -> None:
    if df is None or len(df) == 0:
        return
    required = {"nie", "nie_ci_low", "nie_ci_high", "mediator_label"}
    if not required.issubset(df.columns):
        return
    dd = df.loc[np.isfinite(df["nie"])].copy()
    if len(dd) == 0:
        return
    dd["abs_nie"] = dd["nie"].abs()
    dd = dd.sort_values("abs_nie", ascending=False).head(20).sort_values("nie")
    y = np.arange(len(dd))
    plt.figure(figsize=(9, max(4.5, 0.35 * len(dd) + 1.5)))
    plt.errorbar(
        dd["nie"],
        y,
        xerr=[dd["nie"] - dd["nie_ci_low"], dd["nie_ci_high"] - dd["nie"]],
        fmt="o",
        capsize=3,
    )
    plt.axvline(0, linestyle="--", linewidth=1)
    plt.yticks(y, dd["mediator_label"].tolist())
    plt.xlabel("Natural indirect effect on 10-year risk scale")
    plt.title(title)
    plt.tight_layout()
    outfile.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()


def draw_summary_heatmaps(summary_df: pd.DataFrame, figdir: Path) -> None:
    if summary_df is None or len(summary_df) == 0:
        return
    for value_col in [
        "n_stage1_base_sig",
        "n_stage2_base_sig",
        "n_bridge_candidates",
        "n_stable_full_candidates",
    ]:
        if value_col not in summary_df.columns:
            continue
        try:
            mat = summary_df.pivot(index="exposure_label", columns="layer", values=value_col)
            plt.figure(figsize=(7, max(4, 0.35 * len(mat) + 1.5)))
            plt.imshow(mat.values, aspect="auto")
            plt.xticks(range(mat.shape[1]), mat.columns, rotation=30, ha="right")
            plt.yticks(range(mat.shape[0]), mat.index)
            for i in range(mat.shape[0]):
                for j in range(mat.shape[1]):
                    val = mat.values[i, j]
                    plt.text(j, i, str(int(val)) if np.isfinite(val) else "NA", ha="center", va="center")
            plt.title(value_col)
            plt.colorbar()
            plt.tight_layout()
            plt.savefig(figdir / f"summary_{value_col}.png", dpi=300, bbox_inches="tight")
            plt.close()
        except Exception as exc:
            logging.warning("Failed to draw summary heatmap for %s: %s", value_col, exc)


# -----------------------------------------------------------------------------
# Model engines
# -----------------------------------------------------------------------------


def stage1_gi_to_marker(
    df: pd.DataFrame,
    exposure: str,
    features: Sequence[str],
    covariates: Sequence[str],
    layer: str,
    adjustment: str,
    label_map: Mapping[str, str],
    config: AnalysisConfig,
) -> pd.DataFrame:
    rows = []
    covariates = [c for c in covariates if c in df.columns]

    for feature in features:
        if feature not in df.columns or exposure not in df.columns:
            continue
        sub = df[[exposure, feature] + covariates].copy()
        sub[exposure] = safe_numeric(sub[exposure])
        sub[feature] = safe_numeric(sub[feature])
        sub = sub.loc[sub[exposure].notna() & sub[feature].notna()].copy()

        if len(sub) < config.min_total_n:
            continue
        n_exposed = int((sub[exposure] == 1).sum())
        n_unexposed = int((sub[exposure] == 0).sum())
        if n_exposed < config.min_exposed or n_unexposed < config.min_unexposed:
            continue

        x = clean_frame(sub[[exposure] + covariates], protected=[])
        if exposure not in x.columns or x[exposure].nunique() < 2:
            continue
        x = sm.add_constant(x, has_constant="add")
        y = sub.loc[x.index, feature].values

        try:
            fit = sm.OLS(y, x).fit(cov_type="HC3")
            beta = float(fit.params[exposure])
            rows.append(
                {
                    "analysis": "GI_to_marker",
                    "adjustment": adjustment,
                    "exposure": exposure,
                    "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                    "feature": feature,
                    "feature_label": label_feature(feature, label_map),
                    "layer": layer,
                    "n": int(len(x)),
                    "n_exposed": n_exposed,
                    "n_unexposed": n_unexposed,
                    "beta": beta,
                    "se": float(fit.bse[exposure]),
                    "p_value": float(fit.pvalues[exposure]),
                    "direction": "positive" if beta > 0 else "negative",
                }
            )
        except Exception as exc:
            logging.debug("Stage 1 failed: %s | %s | %s | %s", exposure, feature, adjustment, exc)

    out = pd.DataFrame(rows)
    if len(out) > 0:
        out["fdr"] = fdr_bh(out["p_value"].values)
    return out


def stage2_marker_to_outcome_cox(
    df: pd.DataFrame,
    outcome: str,
    exposure: str,
    features: Sequence[str],
    covariates: Sequence[str],
    layer: str,
    adjustment: str,
    label_map: Mapping[str, str],
    config: AnalysisConfig,
) -> pd.DataFrame:
    if outcome not in OUTCOMES:
        return pd.DataFrame()
    time_col = OUTCOMES[outcome]["time"]
    event_col = OUTCOMES[outcome]["event"]
    if time_col not in df.columns or event_col not in df.columns:
        return pd.DataFrame()

    rows = []
    covariates = [c for c in covariates if c in df.columns]

    for feature in features:
        if feature not in df.columns:
            continue
        keep = [time_col, event_col, exposure, feature] + covariates
        keep = list(dict.fromkeys([c for c in keep if c in df.columns]))
        sub = df[keep].copy()
        sub[time_col] = safe_numeric(sub[time_col])
        sub[event_col] = ensure_binary(sub[event_col])
        sub[feature] = safe_numeric(sub[feature])
        if exposure in sub.columns:
            sub[exposure] = safe_numeric(sub[exposure])
        sub = sub.loc[
            sub[time_col].notna()
            & (sub[time_col] > 0)
            & sub[event_col].notna()
            & sub[feature].notna()
        ].copy()
        if len(sub) < config.min_total_n or sub[event_col].sum() < config.min_events:
            continue

        model_df = sub.rename(columns={time_col: "time_years", event_col: "event"})
        model_df = clean_frame(model_df, protected=["time_years", "event"])
        if feature not in model_df.columns:
            continue

        try:
            cph = CoxPHFitter(penalizer=config.cox_penalizer)
            cph.fit(model_df, duration_col="time_years", event_col="event")
            summ = cph.summary.loc[feature]
            rows.append(
                {
                    "analysis": "marker_to_outcome_cox",
                    "outcome": outcome,
                    "adjustment": adjustment,
                    "exposure": exposure,
                    "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                    "feature": feature,
                    "feature_label": label_feature(feature, label_map),
                    "layer": layer,
                    "n": int(len(model_df)),
                    "n_event": int(model_df["event"].sum()),
                    "hr": float(np.exp(summ["coef"])),
                    "ci_low": float(np.exp(summ["coef lower 95%"])),
                    "ci_high": float(np.exp(summ["coef upper 95%"])),
                    "p_value": float(summ["p"]),
                }
            )
        except Exception as exc:
            logging.debug("Stage 2 Cox failed: %s | %s | %s | %s", outcome, feature, adjustment, exc)

    out = pd.DataFrame(rows)
    if len(out) > 0:
        out["fdr"] = fdr_bh(out["p_value"].values)
    return out


def stage2_marker_to_ad_horizon_logit(
    df: pd.DataFrame,
    exposure: str,
    features: Sequence[str],
    covariates: Sequence[str],
    layer: str,
    adjustment: str,
    horizon: int,
    label_map: Mapping[str, str],
    config: AnalysisConfig,
) -> pd.DataFrame:
    if "ad_time_years" not in df.columns or "ad_event" not in df.columns:
        return pd.DataFrame()
    eligible, y = build_horizon_binary_outcome(df, "ad_time_years", "ad_event", horizon)
    base = df.loc[eligible].copy()
    outcome_col = f"ad_{horizon}y"
    base[outcome_col] = y.loc[base.index]
    base = base.loc[base[outcome_col].notna()].copy()
    if len(base) == 0:
        return pd.DataFrame()
    base[outcome_col] = base[outcome_col].astype(int)
    if base[outcome_col].sum() < config.min_events:
        return pd.DataFrame()

    rows = []
    covariates = [c for c in covariates if c in base.columns]

    for feature in features:
        if feature not in base.columns:
            continue
        keep = [outcome_col, exposure, feature] + covariates
        keep = list(dict.fromkeys([c for c in keep if c in base.columns]))
        sub = base[keep].copy()
        sub[feature] = safe_numeric(sub[feature])
        sub = sub.loc[sub[feature].notna()].copy()
        if len(sub) < config.min_total_n or sub[outcome_col].sum() < config.min_events:
            continue

        x = clean_frame(sub[[feature, exposure] + covariates], protected=[])
        if feature not in x.columns:
            continue
        yv = sub.loc[x.index, outcome_col].astype(int).values
        if len(np.unique(yv)) < 2:
            continue
        try:
            lr = LogisticRegression(max_iter=1000, solver="lbfgs")
            lr.fit(x.values, yv)
            coef = float(lr.coef_[0][list(x.columns).index(feature)])
            rows.append(
                {
                    "analysis": "marker_to_ad_horizon_logit",
                    "horizon_years": horizon,
                    "adjustment": adjustment,
                    "exposure": exposure,
                    "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                    "feature": feature,
                    "feature_label": label_feature(feature, label_map),
                    "layer": layer,
                    "n": int(len(x)),
                    "n_event": int(np.sum(yv)),
                    "or": float(np.exp(coef)),
                    "coef": coef,
                }
            )
        except Exception as exc:
            logging.debug("Stage 2 horizon logit failed: %s | %s | %s", horizon, feature, exc)

    return pd.DataFrame(rows)


def gcomp_bridge(
    df: pd.DataFrame,
    exposure: str,
    mediator: str,
    covariates: Sequence[str],
    horizon: int,
    n_boot: int,
    seed: int,
    config: AnalysisConfig,
) -> Optional[Dict[str, float]]:
    if "ad_time_years" not in df.columns or "ad_event" not in df.columns:
        return None
    eligible, y = build_horizon_binary_outcome(df, "ad_time_years", "ad_event", horizon)
    covariates = [c for c in covariates if c in df.columns]
    keep = [exposure, mediator] + covariates
    keep = list(dict.fromkeys([c for c in keep if c in df.columns]))
    if exposure not in keep or mediator not in keep:
        return None
    base = df.loc[eligible, keep].copy()
    base[f"y_{horizon}"] = y.loc[eligible].astype(int).values
    base[exposure] = safe_numeric(base[exposure])
    base[mediator] = safe_numeric(base[mediator])
    base = base.loc[base[exposure].notna() & base[mediator].notna()].copy()
    if len(base) < config.min_total_n or base[f"y_{horizon}"].sum() < config.min_events:
        return None

    for col in covariates:
        base[col] = safe_numeric(base[col]).fillna(0)

    def one_sample(dd: pd.DataFrame) -> Tuple[float, float, float]:
        xm = clean_frame(dd[[exposure] + covariates], protected=[])
        if exposure not in xm.columns:
            raise ValueError("Exposure removed from mediator model.")
        xm = sm.add_constant(xm, has_constant="add")
        fit_m = sm.OLS(dd.loc[xm.index, mediator].values, xm).fit()

        ddo = dd.loc[xm.index].copy()
        xy = clean_frame(ddo[[exposure, mediator] + covariates], protected=[])
        if exposure not in xy.columns or mediator not in xy.columns:
            raise ValueError("Exposure or mediator removed from outcome model.")
        yv = ddo.loc[xy.index, f"y_{horizon}"].astype(int).values
        if len(np.unique(yv)) < 2:
            raise ValueError("Single outcome class.")

        fit_y = LogisticRegression(max_iter=1000, solver="lbfgs")
        fit_y.fit(xy.values, yv)
        y_cols = xy.columns.tolist()

        x_base = ddo.loc[xy.index, [exposure] + covariates].copy()
        x1 = x_base.copy()
        x0 = x_base.copy()
        x1[exposure] = 1
        x0[exposure] = 0
        x1 = sm.add_constant(x1, has_constant="add").reindex(columns=xm.columns, fill_value=0)
        x0 = sm.add_constant(x0, has_constant="add").reindex(columns=xm.columns, fill_value=0)
        m1 = fit_m.predict(x1)
        m0 = fit_m.predict(x0)

        def pred(xval: int, mval: np.ndarray) -> np.ndarray:
            xx = ddo.loc[xy.index, [exposure, mediator] + covariates].copy()
            xx[exposure] = xval
            xx[mediator] = mval
            xx = xx.reindex(columns=y_cols, fill_value=0)
            return fit_y.predict_proba(xx.values)[:, 1]

        p11 = pred(1, m1)
        p10 = pred(1, m0)
        p00 = pred(0, m0)
        total_effect = float(np.mean(p11 - p00))
        natural_indirect = float(np.mean(p11 - p10))
        natural_direct = float(np.mean(p10 - p00))
        return total_effect, natural_indirect, natural_direct

    try:
        te, nie, nde = one_sample(base)
    except Exception as exc:
        logging.debug("G-computation failed for %s | %s: %s", exposure, mediator, exc)
        return None

    rng = np.random.default_rng(seed)
    boot = []
    n = len(base)
    for _ in range(n_boot):
        try:
            idx = rng.integers(0, n, size=n)
            boot.append(one_sample(base.iloc[idx].reset_index(drop=True)))
        except Exception:
            continue

    min_success = max(30, int(0.25 * n_boot))
    if len(boot) < min_success:
        return None

    arr = np.asarray(boot)
    te_b, nie_b, nde_b = arr[:, 0], arr[:, 1], arr[:, 2]
    return {
        "horizon_years": horizon,
        "te": te,
        "te_ci_low": float(np.quantile(te_b, 0.025)),
        "te_ci_high": float(np.quantile(te_b, 0.975)),
        "nie": nie,
        "nie_ci_low": float(np.quantile(nie_b, 0.025)),
        "nie_ci_high": float(np.quantile(nie_b, 0.975)),
        "nde": nde,
        "nde_ci_low": float(np.quantile(nde_b, 0.025)),
        "nde_ci_high": float(np.quantile(nde_b, 0.975)),
        "prop_mediated": float(nie / te) if np.isfinite(te) and abs(te) > 1e-12 else np.nan,
        "n_boot_success": int(len(boot)),
    }


# -----------------------------------------------------------------------------
# Main analysis workflow
# -----------------------------------------------------------------------------


def write_qc_tables(
    analysis: pd.DataFrame,
    layer_sets: Mapping[str, Sequence[str]],
    main_exposures: Sequence[str],
    secondary_exposures: Sequence[str],
    covariate_blocks: Mapping[str, Sequence[str]],
    outdir: Path,
) -> None:
    tables = outdir / "tables"
    qc = outdir / "qc"

    cohort_row = {
        "n_analysis": int(len(analysis)),
        "n_ad_events": int(ensure_binary(analysis["ad_event"]).sum(skipna=True)) if "ad_event" in analysis else np.nan,
        "n_dementia_events": int(ensure_binary(analysis["dementia_event"]).sum(skipna=True))
        if "dementia_event" in analysis
        else np.nan,
        "blood_features": len(layer_sets.get("blood_biochemistry", [])),
        "nmr_features": len(layer_sets.get("nmr_metabolomics", [])),
        "other_features": len(layer_sets.get("other", [])),
        "main_exposures": ";".join(main_exposures),
        "secondary_exposures": ";".join(secondary_exposures),
    }
    save_table(pd.DataFrame([cohort_row]), tables / "cohort_summary.csv")

    block_rows = [
        {"adjustment": name, "n_covariates": len(cols), "covariates": ";".join(cols)}
        for name, cols in covariate_blocks.items()
    ]
    save_table(pd.DataFrame(block_rows), qc / "covariate_blocks.csv")

    exposure_rows = []
    for group, exposures in {"main": main_exposures, "secondary": secondary_exposures}.items():
        for exposure in exposures:
            x = safe_numeric(analysis[exposure])
            exposure_rows.append(
                {
                    "exposure_set": group,
                    "exposure": exposure,
                    "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                    "n_exposed": int((x == 1).sum()),
                    "n_unexposed": int((x == 0).sum()),
                    "missing_rate": float(x.isna().mean()),
                }
            )
    save_table(pd.DataFrame(exposure_rows), qc / "exposure_qc.csv")


def run_analysis(config: AnalysisConfig) -> None:
    np.random.seed(config.seed)
    warnings.filterwarnings("default")

    outdir = Path(config.outdir)
    tables_dir = outdir / "tables"
    figures_dir = outdir / "figures"
    qc_dir = outdir / "qc"
    for path in [outdir, tables_dir, figures_dir, qc_dir]:
        path.mkdir(parents=True, exist_ok=True)
    configure_logging(outdir)

    logging.info("Starting metabolic bridge analysis")
    logging.info("Configuration: %s", asdict(config))

    label_map = load_feature_label_map(config.label_map_path)
    master, metabolic = load_tables(config)
    analysis = prepare_analysis_dataframe(master, metabolic)
    analysis = prepare_covariates(analysis)
    covariate_blocks = build_covariate_blocks(analysis)
    analysis = standardize_covariates(analysis, covariate_blocks)

    layer_sets = get_metabolic_layers(metabolic, config)
    layer_sets = {k: v for k, v in layer_sets.items() if k in ["blood_biochemistry", "nmr_metabolomics"]}
    main_exposures = [e for e in MAIN_EXPOSURES if e in analysis.columns]
    secondary_exposures = [e for e in SECONDARY_EXPOSURES if e in analysis.columns]

    logging.info("Analysis N: %d", len(analysis))
    logging.info("Main exposures: %s", main_exposures)
    logging.info("Secondary exposures: %s", secondary_exposures)
    logging.info("Metabolic layers: %s", {k: len(v) for k, v in layer_sets.items()})
    logging.info("Covariate blocks: %s", {k: len(v) for k, v in covariate_blocks.items()})

    write_qc_tables(analysis, layer_sets, main_exposures, secondary_exposures, covariate_blocks, outdir)

    all_stage1: List[pd.DataFrame] = []
    all_stage2_cox: List[pd.DataFrame] = []
    all_stage2_logit: List[pd.DataFrame] = []
    all_candidates: List[pd.DataFrame] = []
    all_bridge_markers: List[pd.DataFrame] = []
    all_bridge_modules: List[pd.DataFrame] = []
    all_module_scores: List[pd.DataFrame] = []
    manifest_rows: List[Dict[str, object]] = []

    exposure_groups = {"main": main_exposures, "secondary": secondary_exposures}

    for exposure_set, exposures in exposure_groups.items():
        for exposure in exposures:
            x = safe_numeric(analysis[exposure])
            n_exposed = int((x == 1).sum())
            n_unexposed = int((x == 0).sum())
            if n_exposed < config.min_exposed or n_unexposed < config.min_unexposed:
                logging.info("Skipping %s: exposed=%d, unexposed=%d", exposure, n_exposed, n_unexposed)
                continue
            logging.info("Exposure: %s/%s | exposed=%d", exposure_set, exposure, n_exposed)

            for layer, features in layer_sets.items():
                if not features:
                    continue
                logging.info("Layer: %s | features=%d", layer, len(features))
                stage1_by_adjustment: Dict[str, pd.DataFrame] = {}
                stage2_by_adjustment: Dict[str, pd.DataFrame] = {}

                for adjustment in ["base", "full", "detection"]:
                    if adjustment == "detection" and len(covariate_blocks.get("detection", [])) == len(covariate_blocks.get("full", [])):
                        continue
                    covariates = covariate_blocks[adjustment]

                    s1 = stage1_gi_to_marker(
                        analysis, exposure, features, covariates, layer, adjustment, label_map, config
                    )
                    stage1_by_adjustment[adjustment] = s1
                    if len(s1) > 0:
                        s1["exposure_set"] = exposure_set
                        all_stage1.append(s1)
                        save_table(s1, tables_dir / f"screen1_{adjustment}_{exposure}_{layer}.csv")
                        if config.make_figures and adjustment == "base":
                            draw_volcano(
                                s1,
                                f"Stage 1: {EXPOSURE_LABELS.get(exposure, exposure)} | {layer}",
                                figures_dir / f"volcano_stage1_{exposure}_{layer}.png",
                            )

                    significant_features = (
                        s1.loc[s1["fdr"] < config.fdr_alpha_stage1, "feature"].tolist()
                        if len(s1) > 0 and "fdr" in s1.columns
                        else []
                    )
                    if not significant_features:
                        continue

                    s2 = stage2_marker_to_outcome_cox(
                        analysis,
                        "strict_ad",
                        exposure,
                        significant_features,
                        covariates,
                        layer,
                        adjustment,
                        label_map,
                        config,
                    )
                    stage2_by_adjustment[adjustment] = s2
                    if len(s2) > 0:
                        s2["exposure_set"] = exposure_set
                        all_stage2_cox.append(s2)
                        save_table(s2, tables_dir / f"screen2_cox_{adjustment}_{exposure}_{layer}.csv")

                    for horizon in config.support_horizons:
                        sh = stage2_marker_to_ad_horizon_logit(
                            analysis,
                            exposure,
                            significant_features,
                            covariates,
                            layer,
                            adjustment,
                            horizon,
                            label_map,
                            config,
                        )
                        if len(sh) > 0:
                            sh["exposure_set"] = exposure_set
                            all_stage2_logit.append(sh)
                            save_table(
                                sh,
                                tables_dir / f"screen2_logit_{horizon}y_{adjustment}_{exposure}_{layer}.csv",
                            )

                # Support outcome: all-cause dementia, base adjustment only.
                s1_base = stage1_by_adjustment.get("base", pd.DataFrame())
                if (
                    "dementia_event" in analysis.columns
                    and "dementia_time_years" in analysis.columns
                    and len(s1_base) > 0
                ):
                    dementia_features = s1_base.loc[
                        s1_base["fdr"] < config.fdr_alpha_stage1, "feature"
                    ].tolist()
                    dementia_s2 = stage2_marker_to_outcome_cox(
                        analysis,
                        "all_cause_dementia",
                        exposure,
                        dementia_features,
                        covariate_blocks["base"],
                        layer,
                        "base_support_dementia",
                        label_map,
                        config,
                    )
                    if len(dementia_s2) > 0:
                        dementia_s2["exposure_set"] = exposure_set
                        all_stage2_cox.append(dementia_s2)
                        save_table(
                            dementia_s2,
                            tables_dir / f"screen2_cox_base_support_dementia_{exposure}_{layer}.csv",
                        )

                s2_base = stage2_by_adjustment.get("base", pd.DataFrame())
                if len(s1_base) == 0 or len(s2_base) == 0:
                    continue

                s1_sig = s1_base.loc[s1_base["fdr"] < config.fdr_alpha_stage1].copy()
                s2_sig = s2_base.loc[s2_base["fdr"] < config.fdr_alpha_stage2].copy()
                candidates = sorted(set(s1_sig["feature"]).intersection(set(s2_sig["feature"])))

                s1_full = stage1_by_adjustment.get("full", pd.DataFrame())
                s2_full = stage2_by_adjustment.get("full", pd.DataFrame())
                stable1 = set(s1_full.loc[s1_full["fdr"] < config.fdr_alpha_stage1, "feature"]) if len(s1_full) else set()
                stable2 = set(s2_full.loc[s2_full["fdr"] < config.fdr_alpha_stage2, "feature"]) if len(s2_full) else set()
                stable_candidates = set(candidates).intersection(stable1).intersection(stable2)

                manifest_rows.append(
                    {
                        "exposure_set": exposure_set,
                        "exposure": exposure,
                        "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                        "layer": layer,
                        "n_stage1_base_sig": int(len(s1_sig)),
                        "n_stage2_base_sig": int(len(s2_sig)),
                        "n_bridge_candidates": int(len(candidates)),
                        "n_stable_full_candidates": int(len(stable_candidates)),
                    }
                )

                if not candidates:
                    continue

                candidate_df = s1_sig.merge(
                    s2_sig[["feature", "hr", "ci_low", "ci_high", "p_value", "fdr"]],
                    on="feature",
                    how="inner",
                    suffixes=("_stage1", "_stage2cox"),
                )
                candidate_df["exposure_set"] = exposure_set
                candidate_df["layer"] = layer
                candidate_df["stable_stage1_full"] = candidate_df["feature"].isin(stable1)
                candidate_df["stable_stage2_full"] = candidate_df["feature"].isin(stable2)
                candidate_df["stable_full_candidate"] = (
                    candidate_df["stable_stage1_full"] & candidate_df["stable_stage2_full"]
                )
                all_candidates.append(candidate_df)
                save_table(candidate_df, tables_dir / f"candidate_markers_primary_{exposure}_{layer}.csv")

                bridge_specs = [
                    ("base", candidates, covariate_blocks["base"]),
                    (
                        "full_stable_only",
                        candidate_df.loc[candidate_df["stable_full_candidate"], "feature"].tolist(),
                        covariate_blocks["full"],
                    ),
                ]
                for adjustment, bridge_features, covariates in bridge_specs:
                    if not bridge_features:
                        continue
                    rows = []
                    for feature in bridge_features:
                        result = gcomp_bridge(
                            analysis,
                            exposure,
                            feature,
                            covariates,
                            config.primary_horizon,
                            config.bootstrap,
                            stable_seed(exposure, feature, adjustment, base_seed=config.seed),
                            config,
                        )
                        if result is None:
                            continue
                        result.update(
                            {
                                "analysis": "marker_bridge_gcomp",
                                "adjustment": adjustment,
                                "exposure_set": exposure_set,
                                "exposure": exposure,
                                "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                                "layer": layer,
                                "mediator": feature,
                                "mediator_label": label_feature(feature, label_map),
                                "bridge_type": "marker",
                            }
                        )
                        rows.append(result)
                    if rows:
                        bridge_df = pd.DataFrame(rows)
                        all_bridge_markers.append(bridge_df)
                        save_table(bridge_df, tables_dir / f"bridge_marker_{adjustment}_{exposure}_{layer}.csv")
                        if config.make_figures:
                            draw_bridge(
                                bridge_df,
                                f"Marker bridge effects ({adjustment}): {EXPOSURE_LABELS.get(exposure, exposure)} | {layer}",
                                figures_dir / f"bridge_marker_{adjustment}_{exposure}_{layer}.png",
                            )

                pca_features = candidates[: config.max_markers_for_pca]
                if len(pca_features) >= config.min_markers_for_pca:
                    x_pca = analysis[pca_features].apply(safe_numeric)
                    x_pca = x_pca.apply(lambda col: zscore(winsorize(col))).fillna(0)
                    pca = PCA(n_components=min(3, len(pca_features)), random_state=config.seed)
                    components = pca.fit_transform(x_pca.values)
                    module_name = f"module_pc1__{exposure}__{layer}"
                    score_df = pd.DataFrame({"eid": analysis["eid"].values, module_name: components[:, 0]})
                    all_module_scores.append(score_df)

                    temp = analysis.copy()
                    temp[module_name] = components[:, 0]
                    for adjustment in ["base", "full"]:
                        result = gcomp_bridge(
                            temp,
                            exposure,
                            module_name,
                            covariate_blocks[adjustment],
                            config.primary_horizon,
                            config.bootstrap,
                            stable_seed(module_name, adjustment, base_seed=config.seed),
                            config,
                        )
                        if result is None:
                            continue
                        result.update(
                            {
                                "analysis": "module_bridge_gcomp",
                                "adjustment": adjustment,
                                "exposure_set": exposure_set,
                                "exposure": exposure,
                                "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                                "layer": layer,
                                "mediator": module_name,
                                "mediator_label": f"{layer} module PC1",
                                "bridge_type": "module_pc1",
                                "n_markers_in_module": int(len(pca_features)),
                                "explained_variance_ratio_pc1": float(pca.explained_variance_ratio_[0]),
                                "marker_members": "|".join(pca_features),
                            }
                        )
                        module_df = pd.DataFrame([result])
                        all_bridge_modules.append(module_df)
                        save_table(module_df, tables_dir / f"bridge_module_{adjustment}_{exposure}_{layer}.csv")

    screen1_all = safe_concat(all_stage1)
    screen2_cox_all = safe_concat(all_stage2_cox)
    screen2_logit_all = safe_concat(all_stage2_logit)
    candidate_all = safe_concat(all_candidates)
    bridge_marker_all = safe_concat(all_bridge_markers)
    bridge_module_all = safe_concat(all_bridge_modules)
    manifest_df = pd.DataFrame(manifest_rows)

    save_table(screen1_all, tables_dir / "screen1_all.csv")
    save_table(screen2_cox_all, tables_dir / "screen2_cox_all.csv")
    save_table(screen2_logit_all, tables_dir / "screen2_logit_all.csv")
    save_table(candidate_all, tables_dir / "candidate_markers_all.csv")
    save_table(bridge_marker_all, tables_dir / "bridge_marker_all.csv")
    save_table(bridge_module_all, tables_dir / "bridge_module_all.csv")
    save_table(manifest_df, tables_dir / "candidate_manifest.csv")

    if all_module_scores:
        module_scores = all_module_scores[0]
        for score in all_module_scores[1:]:
            module_scores = module_scores.merge(score, on="eid", how="outer")
        save_table(module_scores, tables_dir / "metabolic_module_scores_for_downstream.csv")

    summary_rows = []
    for exposure in main_exposures + secondary_exposures:
        for layer in ["blood_biochemistry", "nmr_metabolomics"]:
            row = {
                "exposure": exposure,
                "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                "layer": layer,
            }
            if len(screen1_all):
                row["n_stage1_base_sig"] = int(
                    (
                        (screen1_all["exposure"] == exposure)
                        & (screen1_all["layer"] == layer)
                        & (screen1_all["adjustment"] == "base")
                        & (screen1_all["fdr"] < config.fdr_alpha_stage1)
                    ).sum()
                )
            else:
                row["n_stage1_base_sig"] = 0
            if len(screen2_cox_all):
                row["n_stage2_base_sig"] = int(
                    (
                        (screen2_cox_all["exposure"] == exposure)
                        & (screen2_cox_all["layer"] == layer)
                        & (screen2_cox_all["adjustment"] == "base")
                        & (screen2_cox_all["outcome"] == "strict_ad")
                        & (screen2_cox_all["fdr"] < config.fdr_alpha_stage2)
                    ).sum()
                )
            else:
                row["n_stage2_base_sig"] = 0
            if len(candidate_all):
                row["n_bridge_candidates"] = int(
                    ((candidate_all["exposure"] == exposure) & (candidate_all["layer"] == layer)).sum()
                )
                row["n_stable_full_candidates"] = int(
                    (
                        (candidate_all["exposure"] == exposure)
                        & (candidate_all["layer"] == layer)
                        & (candidate_all["stable_full_candidate"] == True)
                    ).sum()
                )
            else:
                row["n_bridge_candidates"] = 0
                row["n_stable_full_candidates"] = 0
            if len(bridge_marker_all):
                row["n_marker_bridge_base"] = int(
                    (
                        (bridge_marker_all["exposure"] == exposure)
                        & (bridge_marker_all["layer"] == layer)
                        & (bridge_marker_all["adjustment"] == "base")
                    ).sum()
                )
            else:
                row["n_marker_bridge_base"] = 0
            if len(bridge_module_all):
                row["n_module_bridge_base"] = int(
                    (
                        (bridge_module_all["exposure"] == exposure)
                        & (bridge_module_all["layer"] == layer)
                        & (bridge_module_all["adjustment"] == "base")
                    ).sum()
                )
            else:
                row["n_module_bridge_base"] = 0
            summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    save_table(summary_df, tables_dir / "summary_table.csv")
    if config.make_figures:
        draw_summary_heatmaps(summary_df, figures_dir)

    metadata = {
        "script": "metabolic_bridge_analysis.py",
        "config": asdict(config),
        "n_analysis": int(len(analysis)),
        "n_ad_events": int(ensure_binary(analysis["ad_event"]).sum(skipna=True)) if "ad_event" in analysis else None,
        "n_dementia_events": int(ensure_binary(analysis["dementia_event"]).sum(skipna=True))
        if "dementia_event" in analysis
        else None,
        "metabolic_layers": {k: len(v) for k, v in layer_sets.items()},
        "main_exposures": main_exposures,
        "secondary_exposures": secondary_exposures,
        "covariate_blocks": covariate_blocks,
        "primary_candidate_rule": "Base-adjusted Stage 1 GI->marker FDR<threshold and base-adjusted Stage 2 marker->AD Cox FDR<threshold.",
        "stability_rule": "Candidate also significant in full-adjusted Stage 1 and Stage 2 models.",
        "bridge_method": f"Restricted g-computation with {config.primary_horizon}-year logistic AD outcome.",
        "interpretation": "Pathway-consistent metabolic bridge evidence, not definitive causal mediation.",
        "key_standardization_changes": [
            "Removed local hard-coded default paths and replaced them with command-line arguments.",
            "Moved top-level execution into a main workflow for import-safe execution.",
            "Added deterministic stable seeds for bootstrap tasks instead of Python's randomized hash().",
            "Added structured output folders, logging, metadata, and QC tables.",
            "Kept CoxPH penalization explicit in configuration and metadata.",
        ],
    }
    with open(outdir / "run_metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    report_lines = [
        "Metabolic bridge analysis completed.",
        "",
        f"Input master table: {config.master_path}",
        f"Input metabolic table: {config.metabolic_path}",
        f"Output directory: {config.outdir}",
        f"Analysis cohort n = {len(analysis)}",
        f"AD events = {int(ensure_binary(analysis['ad_event']).sum(skipna=True)) if 'ad_event' in analysis else 'NA'}",
        f"Main exposures = {main_exposures}",
        f"Secondary exposures = {secondary_exposures}",
        f"Metabolic layers = { {k: len(v) for k, v in layer_sets.items()} }",
        "",
        "Design:",
        " - Primary candidate selection: base-adjusted GI->marker and marker->AD Cox, both FDR below threshold.",
        " - Full-adjusted models are used for stability and sensitivity.",
        " - Bowel cancer screening is used as detection-bias sensitivity when available.",
        f" - Bridge decomposition uses restricted g-computation at the {config.primary_horizon}-year AD horizon.",
        " - All-cause dementia is used as a support outcome when available.",
        "",
        "Summary table:",
        summary_df.to_string(index=False) if len(summary_df) else "No summary rows produced.",
    ]
    with open(outdir / "analysis_report.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    logging.info("Analysis completed. Read first: %s", outdir / "analysis_report.txt")
    logging.info("Summary table: %s", tables_dir / "summary_table.csv")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run metabolic bridge analysis for GI disease exposures and AD risk.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--master", required=True, help="Participant-level preprocessed master CSV file.")
    parser.add_argument("--metabolic", required=True, help="Metabolic feature CSV file with eid and marker columns.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--label-map", default=None, help="Optional CSV with columns feature,label.")
    parser.add_argument("--seed", type=int, default=20260510, help="Random seed.")
    parser.add_argument("--min-total-n", type=int, default=2000, help="Minimum sample size per model.")
    parser.add_argument("--min-events", type=int, default=80, help="Minimum number of events per outcome model.")
    parser.add_argument("--min-exposed", type=int, default=50, help="Minimum number of exposed participants.")
    parser.add_argument("--min-unexposed", type=int, default=50, help="Minimum number of unexposed participants.")
    parser.add_argument("--fdr-alpha-stage1", type=float, default=0.05, help="FDR threshold for GI-to-marker screening.")
    parser.add_argument("--fdr-alpha-stage2", type=float, default=0.05, help="FDR threshold for marker-to-outcome screening.")
    parser.add_argument("--primary-horizon", type=int, default=10, help="Primary AD horizon for bridge decomposition.")
    parser.add_argument("--support-horizons", default="5,10,15", help="Comma-separated horizons for logistic support analyses.")
    parser.add_argument("--bootstrap", type=int, default=200, help="Number of bootstrap replicates for bridge decomposition.")
    parser.add_argument("--max-markers-for-pca", type=int, default=30, help="Maximum candidate markers used for PCA module construction.")
    parser.add_argument("--min-markers-for-pca", type=int, default=2, help="Minimum candidate markers required for PCA module construction.")
    parser.add_argument("--cox-penalizer", type=float, default=0.01, help="L2 penalizer used in lifelines CoxPHFitter.")
    parser.add_argument("--blood-prefix", default="p30", help="Column prefix used to identify blood-biochemistry features.")
    parser.add_argument("--nmr-prefix", default="p23", help="Column prefix used to identify NMR-metabolomics features.")
    parser.add_argument("--max-features-per-layer", type=int, default=None, help="Optional cap for quick testing.")
    parser.add_argument("--no-figures", action="store_true", help="Disable figure generation.")
    return parser


def config_from_args(args: argparse.Namespace) -> AnalysisConfig:
    return AnalysisConfig(
        master_path=args.master,
        metabolic_path=args.metabolic,
        outdir=args.outdir,
        label_map_path=args.label_map,
        seed=args.seed,
        min_total_n=args.min_total_n,
        min_events=args.min_events,
        min_exposed=args.min_exposed,
        min_unexposed=args.min_unexposed,
        fdr_alpha_stage1=args.fdr_alpha_stage1,
        fdr_alpha_stage2=args.fdr_alpha_stage2,
        primary_horizon=args.primary_horizon,
        support_horizons=parse_horizons(args.support_horizons),
        bootstrap=args.bootstrap,
        max_markers_for_pca=args.max_markers_for_pca,
        min_markers_for_pca=args.min_markers_for_pca,
        cox_penalizer=args.cox_penalizer,
        blood_prefix=args.blood_prefix,
        nmr_prefix=args.nmr_prefix,
        max_features_per_layer=args.max_features_per_layer,
        make_figures=not args.no_figures,
    )


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    config = config_from_args(args)
    run_analysis(config)


if __name__ == "__main__":
    main()
