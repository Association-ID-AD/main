#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
genetic_interaction_analysis.py
================================

Genetic susceptibility and effect-modification analysis for GI disease and
AD/dementia outcomes.

Core question
-------------
Does genetic susceptibility modify the association between gastrointestinal
(GI) disease and AD/dementia risk?

Primary screening model
-----------------------
Person-time-adjusted Poisson GLM:

    event ~ GI + genetic_modifier + GI:genetic_modifier + covariates
    offset(log(follow-up time))

Optional binary risk model
--------------------------
Logistic GLM with log follow-up time as an adjustment covariate.

Optional confirmatory model
---------------------------
Cox proportional hazards models for prioritized interaction signals.

Interpretation
--------------
Genetic variables are treated as susceptibility background and effect modifiers,
not as mediators.

Example
-------
python genetic_interaction_analysis.py \
    --master data/master_preprocessed_with_genetics.csv \
    --outdir results/genetic_interaction_analysis \
    --outcomes strict_ad,all_cause_dementia \
    --modifiers APOE4_carrier,AD_PRS_enhanced \
    --covariate-sets lean,expanded \
    --engines poisson \
    --cox-confirm
"""

from __future__ import annotations

import argparse
import json
import os
import time
import warnings
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

try:
    from lifelines import CoxPHFitter
    HAS_LIFELINES = True
except Exception:  # pragma: no cover - optional dependency
    CoxPHFitter = None
    HAS_LIFELINES = False


# -----------------------------------------------------------------------------
# Definitions
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

MAIN_EXPOSURES = [
    "exp_diverticular",
    "exp_other_functional_intestinal",
    "exp_ibd",
    "exp_other_chronic_intestinal",
    "exp_ibs",
]
SECONDARY_EXPOSURES = ["exp_appendiceal", "exp_malabsorption"]
ALL_EXPOSURES = MAIN_EXPOSURES + SECONDARY_EXPOSURES

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
EXPOSURE_LABELS = {
    "exp_diverticular": "Diverticular disease",
    "exp_other_functional_intestinal": "Other functional intestinal disorders",
    "exp_other_chronic_intestinal": "Other chronic intestinal disease",
    "exp_ibd": "Inflammatory bowel disease",
    "exp_ibs": "Irritable bowel syndrome",
    "exp_appendiceal": "Appendiceal disease",
    "exp_malabsorption": "Malabsorption",
}

OUTCOME_SPECS = {
    "strict_ad": {
        "time": "ad_time_years",
        "event": "ad_event",
        "label": "Incident AD",
        "primary": True,
    },
    "all_cause_dementia": {
        "time": "dementia_time_years",
        "event": "dementia_event",
        "label": "All-cause dementia",
        "primary": False,
    },
    "cognition_defined_mci": {
        "time": "mci_time_years",
        "event": "mci_event",
        "label": "Cognition-defined MCI-like status",
        "primary": False,
    },
    "mci_or_ad": {
        "time": "mci_ad_time_years",
        "event": "mci_ad_event",
        "label": "MCI-like status or AD",
        "primary": False,
    },
    "incident_mci_or_ad_legacy": {
        "time": "mci_or_ad_time_years",
        "event": "mci_or_ad_event",
        "label": "Legacy MCI-or-AD composite",
        "primary": False,
    },
}

APOE_CANDIDATES = [
    "apoe4_carrier", "apoe_e4_carrier", "apoe4", "apoe_e4", "apoe_e4_carrier_model"
]
APOE_COUNT_CANDIDATES = ["apoe_e4_count", "apoe4_count", "apoe_e4_count_z"]
APOE_GENOTYPE_CANDIDATES = ["apoe_genotype", "apoe_genotype.x", "apoe_genotype.y"]
AD_PRS_STD_CANDIDATES = ["ad_prs_std", "ad_prs", "prs_ad", "ad_prs_z", "p26206_i0", "p26206"]
AD_PRS_ENH_CANDIDATES = ["ad_prs_enhanced_std", "ad_prs_enhanced", "ad_prs_enh", "p26207_i0", "p26207"]

LEAN_COVAR_CANDIDATES = [
    "age_at_baseline", "sex", "ethnicity_5cat", "edu_level", "townsend_i0"
]
EXPANDED_COVAR_CANDIDATES = [
    "age_at_baseline", "sex", "ethnicity_5cat", "edu_level", "townsend_i0",
    "bmi_i0", "whr_i0", "smoke_status_3cat", "smoke_pack_i0", "alcohol_freq",
    "sleep_hours", "total_sedentary_hours", "activity_score_days", "diet_quality_proxy",
    "hx_diabetes", "hx_hypertension", "hx_hyperlipidemia", "hx_chd_cvd",
    "hx_stroke_tia", "hx_depression_anxiety", "hx_sleep_disorder", "hx_parkinson_other_nd",
    "med_statin", "med_antihypertensive", "med_antidiabetic",
    "med_ppi_gi_drug", "med_laxative", "med_antidepressant", "med_steroid_immunosuppressive",
    "bowel_cancer_screening", "parental_dementia_any", "parental_dementia_count",
]

GENETIC_PC_PREFIXES = ["p22009_", "p26201_"]
GENETIC_QC_CANDIDATES = [
    "p22000", "p22000_i0", "p22001", "p22001_i0",
    "p22006", "p22006_i0", "p22019", "p22019_i0",
    "p22027", "p22027_i0", "p26200", "p26200_i0",
    "prs_release_test_flag",
]

GI_PRS_CANDIDATES = {
    "exp_diverticular": ["prs_diverticular", "diverticular_prs", "prs_gi_diverticular"],
    "exp_ibd": ["prs_ibd", "ibd_prs", "prs_gi_ibd"],
    "exp_ibs": ["prs_ibs", "ibs_prs", "prs_gi_ibs"],
    "exp_other_functional_intestinal": ["prs_functional_intestinal", "prs_ofid", "ofid_prs"],
    "exp_malabsorption": ["prs_malabsorption", "malabsorption_prs"],
    "exp_appendiceal": ["prs_appendiceal", "appendiceal_prs"],
    "exp_other_chronic_intestinal": ["prs_chronic_intestinal", "prs_ocid", "ocid_prs"],
}


@dataclass
class RunConfig:
    master: str
    outdir: str
    outcomes: List[str]
    modifiers: List[str]
    covariate_sets: List[str]
    engines: List[str]
    no_stratified: bool
    no_shared_liability: bool
    resume: bool
    max_models: int
    max_iter: int
    min_n: int
    min_events: int
    min_exposed: int
    min_exposed_events: int
    min_n_strat: int
    min_events_strat: int
    max_pcs_lean: int
    max_pcs_expanded: int
    cox_confirm: bool
    max_cox_confirm: int
    cox_penalizer: float
    preflight: bool
    force_full_read: bool
    seed: int
    tabdir: str
    qcdir: str


# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

T0 = time.time()
MODEL_COUNTER = 0
COMPLETED_KEYS: set[str] = set()


def log(msg: str) -> None:
    elapsed = (time.time() - T0) / 60
    print(f"[{elapsed:7.2f} min] {msg}", flush=True)


@contextmanager
def timed(stage: str):
    t = time.time()
    log(f"START {stage}")
    yield
    log(f"END   {stage}: {time.time() - t:.2f} sec")


def safe_numeric(x) -> pd.Series:
    return pd.to_numeric(x, errors="coerce")


def ensure_binary(x) -> pd.Series:
    x = safe_numeric(x)
    return x.where(x.isin([0, 1]), np.nan)


def winsorize(x, q=(0.005, 0.995)) -> pd.Series:
    x = safe_numeric(x).copy()
    lo, hi = x.quantile(q[0]), x.quantile(q[1])
    if pd.notna(lo) and pd.notna(hi):
        x.loc[x < lo] = lo
        x.loc[x > hi] = hi
    return x


def zscore(x) -> pd.Series:
    x = safe_numeric(x)
    mu = x.mean(skipna=True)
    sd = x.std(skipna=True)
    if pd.isna(sd) or sd == 0:
        return pd.Series(np.nan, index=x.index)
    return (x - mu) / sd


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


def detect_column(candidates: Sequence[str], columns: Iterable[str]) -> Optional[str]:
    colset = set(columns)
    for c in candidates:
        if c in colset:
            return c
    return None


def clean_string_name(x) -> str:
    out = str(x)
    for a, b in [
        (" ", "_"), ("/", "_"), ("\\", "_"), ("-", "_"), (":", "_"), (";", "_"),
        ("(", "_"), (")", "_"), ("[", "_"), ("]", "_"), (",", "_"), (".", "_"),
    ]:
        out = out.replace(a, b)
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_")


def save_table(df: Optional[pd.DataFrame], name: str, tabdir: str) -> None:
    path = os.path.join(tabdir, name)
    if df is None:
        df = pd.DataFrame()
    df.to_csv(path, index=False)
    log(f"[SAVE] {path} rows={len(df)}")


def append_table(rows, name: str, tabdir: str) -> None:
    if rows is None:
        return
    if isinstance(rows, dict):
        rows = [rows]
    if len(rows) == 0:
        return

    if name == "fast_model_diagnostics_live.csv":
        path = os.path.join(tabdir, "fast_model_diagnostics_live.jsonl")
        with open(path, "a", encoding="utf-8") as f:
            for r in rows:
                f.write(json.dumps(r, ensure_ascii=False) + "\n")
        return

    path = os.path.join(tabdir, name)
    df = pd.DataFrame(rows)
    header = not os.path.exists(path)
    df.to_csv(path, mode="a", index=False, header=header)


def read_live_diagnostics_jsonl(path: str) -> pd.DataFrame:
    rows = []
    if not os.path.exists(path):
        return pd.DataFrame()
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(json.loads(line))
            except Exception:
                continue
    return pd.DataFrame(rows)


def fdr_add(df: pd.DataFrame, pcol: str, group_cols: Sequence[str], outcol: str = "fdr") -> pd.DataFrame:
    if df is None or len(df) == 0 or pcol not in df.columns:
        return df
    out = df.copy()
    out[outcol] = np.nan
    if group_cols:
        for _, idx in out.groupby(list(group_cols), dropna=False).groups.items():
            out.loc[idx, outcol] = fdr_bh(out.loc[idx, pcol].values)
    else:
        out[outcol] = fdr_bh(out[pcol].values)
    return out


# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------

def coalesce_pair(df: pd.DataFrame, base: str) -> pd.DataFrame:
    xcol, ycol = f"{base}.x", f"{base}.y"
    if xcol in df.columns and ycol in df.columns:
        df[base] = df[xcol].combine_first(df[ycol])
    elif xcol in df.columns and base not in df.columns:
        df[base] = df[xcol]
    elif ycol in df.columns and base not in df.columns:
        df[base] = df[ycol]
    return df


def normalize_apoe_columns(df: pd.DataFrame) -> pd.DataFrame:
    for base in ["apoe4_carrier", "apoe_genotype", "rs429358", "rs7412"]:
        df = coalesce_pair(df, base)

    if "apoe_genotype" not in df.columns:
        geno_like = [c for c in df.columns if "apoe_genotype" in c.lower()]
        vals = None
        for c in geno_like:
            cur = df[c].astype("string")
            vals = cur if vals is None else vals.combine_first(cur)
        if vals is not None:
            df["apoe_genotype"] = vals

    if "apoe4_carrier" not in df.columns or safe_numeric(df["apoe4_carrier"]).notna().sum() == 0:
        apoe_like = [c for c in df.columns if "apoe4_carrier" in c.lower()]
        vals = None
        for c in apoe_like:
            cur = safe_numeric(df[c])
            vals = cur if vals is None else vals.combine_first(cur)
        if vals is not None:
            df["apoe4_carrier"] = vals

    if "apoe4_carrier" not in df.columns or safe_numeric(df["apoe4_carrier"]).notna().sum() == 0:
        if "apoe_genotype" in df.columns:
            g = df["apoe_genotype"].astype(str).str.lower().str.strip()
            carrier = pd.Series(np.nan, index=df.index, dtype=float)
            notna = ~g.isin(["", "nan", "na", "none"])
            carrier.loc[notna] = 0
            carrier.loc[notna & g.str.contains("e4", na=False)] = 1
            df["apoe4_carrier"] = carrier

    if "apoe_e4_count" not in df.columns and "apoe_genotype" in df.columns:
        g = df["apoe_genotype"].astype(str).str.lower()
        cnt = pd.Series(np.nan, index=df.index, dtype=float)
        notna = ~g.isin(["", "nan", "na", "none"])
        cnt.loc[notna] = g.loc[notna].str.count("e4").astype(float)
        df["apoe_e4_count"] = cnt

    return df


def max_binary_across(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    cols = [c for c in cols if c in df.columns]
    if not cols:
        return pd.Series(np.nan, index=df.index)
    tmp = pd.concat([safe_numeric(df[c]) for c in cols], axis=1)
    out = tmp.max(axis=1, skipna=True)
    out[tmp.isna().all(axis=1)] = np.nan
    return out


def add_grouped_exposures(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    for grouped_exp, components in GROUPED_GI_MAP.items():
        if grouped_exp not in d.columns:
            d[grouped_exp] = max_binary_across(d, components)
        else:
            d[grouped_exp] = ensure_binary(d[grouped_exp])
    return d


def quartile(x) -> pd.Series:
    x = safe_numeric(x)
    try:
        q = pd.qcut(x, 4, labels=False, duplicates="drop")
        return q.astype(float) + 1
    except Exception:
        return pd.Series(np.nan, index=x.index)


def tertile_high_low(x) -> pd.Series:
    x = safe_numeric(x)
    lo = x.quantile(1 / 3)
    hi = x.quantile(2 / 3)
    out = pd.Series(np.nan, index=x.index, dtype=float)
    out.loc[x <= lo] = 0
    out.loc[x >= hi] = 1
    return out


def outcome_available(df: pd.DataFrame, spec: dict) -> bool:
    return spec["time"] in df.columns and spec["event"] in df.columns


def build_selected_columns(master_path: str, outcomes: Sequence[str]) -> Tuple[List[str], List[str], List[str]]:
    header = pd.read_csv(master_path, nrows=0).columns.tolist()
    hset = set(header)
    selected = {"eid"}

    for outcome in outcomes:
        spec = OUTCOME_SPECS.get(outcome)
        if spec:
            selected.update([spec["time"], spec["event"]])

    selected.update(ALL_EXPOSURES)
    for comps in GROUPED_GI_MAP.values():
        selected.update(comps)

    for lst in [APOE_CANDIDATES, APOE_COUNT_CANDIDATES, APOE_GENOTYPE_CANDIDATES, AD_PRS_STD_CANDIDATES, AD_PRS_ENH_CANDIDATES]:
        selected.update(lst)

    selected.update(LEAN_COVAR_CANDIDATES)
    selected.update(EXPANDED_COVAR_CANDIDATES)
    selected.update(GENETIC_QC_CANDIDATES)

    for v in GI_PRS_CANDIDATES.values():
        selected.update(v)

    pc_cols = [c for c in header if any(c.startswith(pref) for pref in GENETIC_PC_PREFIXES)]
    selected.update(pc_cols)
    selected.update([c for c in header if "apoe" in c.lower() or c.lower() in ["rs429358", "rs7412"]])

    usecols = [c for c in header if c in selected]
    missing_important = sorted([c for c in selected if c not in hset and c in (set(ALL_EXPOSURES) | set(LEAN_COVAR_CANDIDATES))])
    return header, usecols, missing_important


def build_covariate_matrix(
    df: pd.DataFrame,
    covar_candidates: Sequence[str],
    pc_limit: int = 10,
    prefix: str = "cov",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    covar_candidates = [c for c in covar_candidates if c in df.columns]
    pc_cols = sorted([c for c in covar_candidates if any(c.startswith(pref) for pref in GENETIC_PC_PREFIXES)])[:pc_limit]
    non_pc_cols = [c for c in covar_candidates if c not in pc_cols]
    covars = list(dict.fromkeys(non_pc_cols + pc_cols))

    mats = []
    actions = []
    for c in covars:
        s = df[c]
        raw_nonmiss = int(s.notna().sum())
        if raw_nonmiss == 0:
            actions.append({"original_column": c, "action": "drop_all_missing", "n_created": 0})
            continue

        num = safe_numeric(s)
        num_nonmiss = int(num.notna().sum())
        retention = num_nonmiss / raw_nonmiss if raw_nonmiss else 0

        if retention >= 0.80:
            vals = set(pd.unique(num.dropna()))
            if len(vals) <= 2 and vals.issubset({0, 1}):
                out = num.fillna(0).astype(float)
                if out.nunique(dropna=True) > 1:
                    mats.append(pd.DataFrame({f"{prefix}__{clean_string_name(c)}": out}, index=df.index))
                    actions.append({"original_column": c, "action": "binary_numeric_kept", "n_created": 1})
                else:
                    actions.append({"original_column": c, "action": "drop_constant_binary", "n_created": 0})
            else:
                out = zscore(winsorize(num))
                med = out.median(skipna=True)
                out = out.fillna(0 if pd.isna(med) else med)
                if out.nunique(dropna=True) > 1 and out.std(skipna=True) > 1e-8:
                    mats.append(pd.DataFrame({f"{prefix}__{clean_string_name(c)}": out}, index=df.index))
                    actions.append({"original_column": c, "action": "continuous_zscore", "n_created": 1})
                else:
                    actions.append({"original_column": c, "action": "drop_constant_continuous", "n_created": 0})
        else:
            ss = s.astype("string").fillna("Missing").str.strip()
            nunique = int(ss.nunique(dropna=True))
            if 2 <= nunique <= 30:
                dummies = pd.get_dummies(ss, prefix=f"{prefix}__{clean_string_name(c)}", drop_first=True, dtype=float)
                keep = [dc for dc in dummies.columns if dummies[dc].nunique(dropna=True) > 1]
                if keep:
                    mats.append(dummies[keep])
                    actions.append({"original_column": c, "action": "categorical_one_hot", "n_created": len(keep)})
                else:
                    actions.append({"original_column": c, "action": "drop_categorical_constant", "n_created": 0})
            else:
                actions.append({"original_column": c, "action": f"drop_high_cardinality_or_low_numeric_retention_{nunique}", "n_created": 0})

    mat = pd.concat(mats, axis=1) if mats else pd.DataFrame(index=df.index)

    drop_cols = []
    seen = {}
    for c in mat.columns:
        h = pd.util.hash_pandas_object(mat[c], index=False).sum()
        if h in seen and mat[c].equals(mat[seen[h]]):
            drop_cols.append(c)
        else:
            seen[h] = c
    if drop_cols:
        mat = mat.drop(columns=drop_cols)

    return mat.astype(float, copy=False), pd.DataFrame(actions)


def drop_problematic_columns(model_df: pd.DataFrame, protected: Sequence[str]) -> Tuple[pd.DataFrame, List[str], str]:
    d = model_df.copy().replace([np.inf, -np.inf], np.nan)
    protected = set(protected)
    dropped = []

    protected_present = [c for c in protected if c in d.columns]
    d = d.dropna(subset=protected_present)

    for c in list(d.columns):
        if c in protected:
            continue
        if d[c].isna().all():
            d = d.drop(columns=[c])
            dropped.append(f"{c}:all_na")
            continue
        if d[c].isna().any():
            med = d[c].median(skipna=True)
            d[c] = d[c].fillna(0.0 if pd.isna(med) else med)
        nunq = d[c].nunique(dropna=True)
        std = d[c].std(skipna=True)
        if nunq <= 1 or (pd.notna(std) and std < 1e-10):
            d = d.drop(columns=[c])
            dropped.append(f"{c}:constant")

    for c in protected:
        if c in d.columns and c not in ["time_years", "event", "offset_log_time"]:
            if d[c].nunique(dropna=True) <= 1:
                return d, dropped, f"protected_predictor_constant:{c}"

    return d, dropped, "ok"


# -----------------------------------------------------------------------------
# Model functions
# -----------------------------------------------------------------------------

def model_key(label_dict: dict, predictors: Sequence[str], engine: str) -> str:
    parts = [
        engine,
        label_dict.get("analysis", ""),
        label_dict.get("covariate_set", ""),
        label_dict.get("outcome", ""),
        label_dict.get("exposure", ""),
        label_dict.get("modifier", ""),
        label_dict.get("stratum", ""),
        "|".join(predictors),
    ]
    return "::".join(map(str, parts))


def fit_fast_event_model(
    df: pd.DataFrame,
    time_col: str,
    event_col: str,
    predictors: Sequence[str],
    covars: Sequence[str],
    label_dict: dict,
    config: RunConfig,
    engine: str = "poisson",
    min_n: Optional[int] = None,
    min_events: Optional[int] = None,
) -> Tuple[Optional[List[dict]], dict]:
    global MODEL_COUNTER, COMPLETED_KEYS

    min_n = config.min_n if min_n is None else min_n
    min_events = config.min_events if min_events is None else min_events
    key = model_key(label_dict, predictors, engine)
    if config.resume and key in COMPLETED_KEYS:
        log(f"[SKIP-DONE] {key}")
        return [], {"model_key": key, **label_dict, "engine": engine, "status": "skipped_completed"}

    if config.max_models and MODEL_COUNTER >= config.max_models:
        return None, {"model_key": key, **label_dict, "engine": engine, "status": "max_models_reached"}

    MODEL_COUNTER += 1
    t_model = time.time()
    log(f"[MODEL-START #{MODEL_COUNTER}] {key}")

    diagnostics = {
        "model_key": key,
        **label_dict,
        "engine": engine,
        "time_col": time_col,
        "event_col": event_col,
        "predictors": ";".join(predictors),
        "n_covars_input": len(covars),
        "status": "",
        "error": "",
    }

    missing = [c for c in [time_col, event_col] + list(predictors) if c not in df.columns]
    if missing:
        diagnostics.update({"status": "missing_required_columns", "missing": ";".join(missing)})
        append_table(diagnostics, "fast_model_diagnostics_live.csv", config.tabdir)
        return [], diagnostics

    cols = list(dict.fromkeys([time_col, event_col] + list(predictors) + [c for c in covars if c in df.columns]))
    sub = df[cols].copy()
    for c in cols:
        sub[c] = safe_numeric(sub[c])

    sub = sub.rename(columns={time_col: "time_years", event_col: "event"})
    sub = sub.loc[sub["time_years"].notna() & (sub["time_years"] > 0) & sub["event"].notna()].copy()

    protected = ["time_years", "event"] + list(predictors)
    sub, dropped, clean_status = drop_problematic_columns(sub, protected)
    diagnostics["dropped_covars"] = ";".join(dropped)
    diagnostics["clean_status"] = clean_status
    diagnostics["n_after_clean"] = int(len(sub))
    diagnostics["n_event_after_clean"] = int(sub["event"].sum()) if len(sub) else 0
    diagnostics["n_covars_after_clean"] = int(len([c for c in sub.columns if c not in protected]))

    if clean_status != "ok":
        diagnostics["status"] = clean_status
        diagnostics["elapsed_sec"] = round(time.time() - t_model, 3)
        append_table(diagnostics, "fast_model_diagnostics_live.csv", config.tabdir)
        return [], diagnostics
    if len(sub) < min_n:
        diagnostics["status"] = "n_below_threshold"
        diagnostics["elapsed_sec"] = round(time.time() - t_model, 3)
        append_table(diagnostics, "fast_model_diagnostics_live.csv", config.tabdir)
        return [], diagnostics
    if int(sub["event"].sum()) < min_events:
        diagnostics["status"] = "events_below_threshold"
        diagnostics["elapsed_sec"] = round(time.time() - t_model, 3)
        append_table(diagnostics, "fast_model_diagnostics_live.csv", config.tabdir)
        return [], diagnostics

    for p in predictors:
        if p in sub.columns:
            vals = set(pd.unique(sub[p].dropna()))
            if vals.issubset({0, 1}):
                diagnostics[f"n_{p}_1"] = int((sub[p] == 1).sum())
                diagnostics[f"events_{p}_1"] = int(sub.loc[sub[p] == 1, "event"].sum())

    y = sub["event"].astype(float)
    xcols = [c for c in sub.columns if c not in ["time_years", "event"]]
    X = sub[xcols].astype(float)
    X = sm.add_constant(X, has_constant="add")

    try:
        if engine == "poisson":
            offset = np.log(np.clip(sub["time_years"].astype(float).values, 1e-6, None))
            fit = sm.GLM(y, X, family=sm.families.Poisson(), offset=offset).fit(maxiter=config.max_iter, disp=0)
            effect_name = "irr"
        elif engine == "logit":
            X2 = X.copy()
            X2["log_followup_time"] = np.log(np.clip(sub["time_years"].astype(float).values, 1e-6, None))
            fit = sm.GLM(y, X2, family=sm.families.Binomial()).fit(maxiter=config.max_iter, disp=0)
            X = X2
            effect_name = "or"
        else:
            raise ValueError(f"Unknown engine: {engine}")

        rows = []
        for term in predictors:
            if term not in fit.params.index:
                continue
            coef = float(fit.params[term])
            se = float(fit.bse[term])
            p = float(fit.pvalues[term])
            rows.append({
                **label_dict,
                "engine": engine,
                "term": term,
                "coef": coef,
                "se": se,
                "effect": float(np.exp(coef)),
                "effect_name": effect_name,
                "ci_low": float(np.exp(coef - 1.96 * se)),
                "ci_high": float(np.exp(coef + 1.96 * se)),
                "p_value": p,
                "n": int(len(sub)),
                "n_event": int(sub["event"].sum()),
                "n_covariates": int(X.shape[1] - 1 - len(predictors)),
                "model_key": key,
            })

        diagnostics["status"] = "fit_ok"
        diagnostics["elapsed_sec"] = round(time.time() - t_model, 3)
        diagnostics["n_fit"] = int(len(sub))
        diagnostics["n_event_fit"] = int(sub["event"].sum())
        append_table(diagnostics, "fast_model_diagnostics_live.csv", config.tabdir)
        log(f"[MODEL-DONE  #{MODEL_COUNTER}] status=fit_ok elapsed={diagnostics['elapsed_sec']} sec")
        return rows, diagnostics

    except Exception as e:
        diagnostics["status"] = "fit_failed"
        diagnostics["error"] = str(e)[:1000]
        diagnostics["elapsed_sec"] = round(time.time() - t_model, 3)
        append_table(diagnostics, "fast_model_diagnostics_live.csv", config.tabdir)
        log(f"[MODEL-FAIL  #{MODEL_COUNTER}] elapsed={diagnostics['elapsed_sec']} sec error={diagnostics['error'][:120]}")
        return [], diagnostics


def fit_cox_confirm(
    df: pd.DataFrame,
    time_col: str,
    event_col: str,
    predictors: Sequence[str],
    covars: Sequence[str],
    label_dict: dict,
    config: RunConfig,
) -> List[dict]:
    if not HAS_LIFELINES:
        return []
    t_model = time.time()
    log(f"[COX-START] {label_dict.get('analysis')} {label_dict.get('outcome')} {label_dict.get('exposure')} {label_dict.get('modifier')}")
    needed = list(dict.fromkeys([time_col, event_col] + list(predictors) + [c for c in covars if c in df.columns]))
    sub = df[needed].copy()
    for c in needed:
        sub[c] = safe_numeric(sub[c])
    sub = sub.rename(columns={time_col: "time_years", event_col: "event"})
    sub = sub.loc[sub["time_years"].notna() & (sub["time_years"] > 0) & sub["event"].notna()].copy()
    sub, _, clean_status = drop_problematic_columns(sub, ["time_years", "event"] + list(predictors))
    if clean_status != "ok" or len(sub) < config.min_n or int(sub["event"].sum()) < config.min_events:
        return []
    try:
        cph = CoxPHFitter(penalizer=config.cox_penalizer)
        cph.fit(sub, duration_col="time_years", event_col="event", robust=False, show_progress=False)
        out = []
        for term in predictors:
            if term in cph.summary.index:
                s = cph.summary.loc[term]
                out.append({
                    **label_dict,
                    "term": term,
                    "coef": float(s["coef"]),
                    "se": float(s["se(coef)"]),
                    "hr": float(np.exp(s["coef"])),
                    "ci_low": float(np.exp(s["coef lower 95%"])),
                    "ci_high": float(np.exp(s["coef upper 95%"])),
                    "p_value": float(s["p"]),
                    "n": int(len(sub)),
                    "n_event": int(sub["event"].sum()),
                    "elapsed_sec": round(time.time() - t_model, 3),
                    "status": "fit_ok",
                })
        log(f"[COX-DONE] elapsed={time.time() - t_model:.2f} sec")
        return out
    except Exception as e:
        log(f"[COX-FAIL] {str(e)[:160]}")
        return []


# -----------------------------------------------------------------------------
# Main analysis
# -----------------------------------------------------------------------------

def parse_args() -> RunConfig:
    parser = argparse.ArgumentParser(description="Fast genetic susceptibility and GI-by-genetic interaction screening.")
    parser.add_argument("--master", required=True, help="Master preprocessed table with GI exposures, outcomes, covariates, and genetics.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--outcomes", default="strict_ad,all_cause_dementia")
    parser.add_argument("--modifiers", default="APOE4_carrier,AD_PRS_enhanced,APOE4_count,AD_PRS_standard")
    parser.add_argument("--covariate-sets", default="lean,expanded", help="Comma-separated: lean,expanded.")
    parser.add_argument("--engines", default="poisson", help="Comma-separated: poisson,logit.")
    parser.add_argument("--no-stratified", action="store_true")
    parser.add_argument("--no-shared-liability", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--max-models", type=int, default=0)
    parser.add_argument("--max-iter", type=int, default=80)
    parser.add_argument("--min-n", type=int, default=1000)
    parser.add_argument("--min-events", type=int, default=20)
    parser.add_argument("--min-exposed", type=int, default=20)
    parser.add_argument("--min-exposed-events", type=int, default=3)
    parser.add_argument("--min-n-strat", type=int, default=300)
    parser.add_argument("--min-events-strat", type=int, default=8)
    parser.add_argument("--max-pcs-lean", type=int, default=10)
    parser.add_argument("--max-pcs-expanded", type=int, default=20)
    parser.add_argument("--cox-confirm", action="store_true")
    parser.add_argument("--max-cox-confirm", type=int, default=30)
    parser.add_argument("--cox-penalizer", type=float, default=0.02)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--force-full-read", action="store_true")
    parser.add_argument("--seed", type=int, default=20260530)
    args = parser.parse_args()

    outdir = args.outdir
    tabdir = os.path.join(outdir, "tables")
    qcdir = os.path.join(outdir, "qc")
    os.makedirs(tabdir, exist_ok=True)
    os.makedirs(qcdir, exist_ok=True)

    return RunConfig(
        master=args.master,
        outdir=outdir,
        outcomes=[x.strip() for x in args.outcomes.split(",") if x.strip()],
        modifiers=[x.strip() for x in args.modifiers.split(",") if x.strip()],
        covariate_sets=[x.strip() for x in args.covariate_sets.split(",") if x.strip()],
        engines=[x.strip() for x in args.engines.split(",") if x.strip()],
        no_stratified=args.no_stratified,
        no_shared_liability=args.no_shared_liability,
        resume=args.resume,
        max_models=args.max_models,
        max_iter=args.max_iter,
        min_n=args.min_n,
        min_events=args.min_events,
        min_exposed=args.min_exposed,
        min_exposed_events=args.min_exposed_events,
        min_n_strat=args.min_n_strat,
        min_events_strat=args.min_events_strat,
        max_pcs_lean=args.max_pcs_lean,
        max_pcs_expanded=args.max_pcs_expanded,
        cox_confirm=args.cox_confirm,
        max_cox_confirm=args.max_cox_confirm,
        cox_penalizer=args.cox_penalizer,
        preflight=args.preflight,
        force_full_read=args.force_full_read,
        seed=args.seed,
        tabdir=tabdir,
        qcdir=qcdir,
    )


def run_analysis(config: RunConfig) -> None:
    global COMPLETED_KEYS
    np.random.seed(config.seed)

    log("=" * 90)
    log("[SCRIPT] Genetic susceptibility and GI-by-genetic interaction screen")
    log(f"[MASTER] {config.master}")
    log(f"[OUTDIR] {config.outdir}")
    log(f"[OUTCOMES] {config.outcomes}")
    log(f"[MODIFIERS] {config.modifiers}")
    log(f"[COVARIATE_SETS] {config.covariate_sets}")
    log(f"[ENGINES] {config.engines}")
    log("=" * 90)

    with timed("read_header"):
        header, usecols, missing_important = build_selected_columns(config.master, config.outcomes)
    if config.force_full_read:
        usecols = header

    pd.DataFrame({
        "item": ["n_header_cols", "n_selected_cols", "missing_important"],
        "value": [len(header), len(usecols), ";".join(missing_important)],
    }).to_csv(os.path.join(config.qcdir, "read_column_selection_summary.csv"), index=False)

    with timed("read_master_csv"):
        df = pd.read_csv(config.master, usecols=usecols, low_memory=False)
    log(f"[LOAD-DONE] shape={df.shape}")

    if "eid" in df.columns:
        df["eid"] = df["eid"].astype(str).str.strip()

    with timed("prepare_variables"):
        df = normalize_apoe_columns(df)
        df = add_grouped_exposures(df)

        apoe_col = detect_column(APOE_CANDIDATES, df.columns)
        apoe_count_col = detect_column(APOE_COUNT_CANDIDATES, df.columns)
        apoe_genotype_col = detect_column(APOE_GENOTYPE_CANDIDATES, df.columns)
        ad_prs_std_col = detect_column(AD_PRS_STD_CANDIDATES, df.columns)
        ad_prs_enh_col = detect_column(AD_PRS_ENH_CANDIDATES, df.columns)

        if apoe_col is not None:
            df[apoe_col] = ensure_binary(df[apoe_col])
        if apoe_count_col is not None:
            df[apoe_count_col] = safe_numeric(df[apoe_count_col])
        if ad_prs_std_col is not None:
            df[ad_prs_std_col] = zscore(winsorize(df[ad_prs_std_col]))
            df["ad_prs_std_quartile"] = quartile(df[ad_prs_std_col])
            df["ad_prs_std_high_vs_low"] = tertile_high_low(df[ad_prs_std_col])
        if ad_prs_enh_col is not None:
            df[ad_prs_enh_col] = zscore(winsorize(df[ad_prs_enh_col]))
            df["ad_prs_enh_quartile"] = quartile(df[ad_prs_enh_col])
            df["ad_prs_enh_high_vs_low"] = tertile_high_low(df[ad_prs_enh_col])

    available_outcomes = {k: v for k, v in OUTCOME_SPECS.items() if k in config.outcomes and outcome_available(df, v)}
    available_exposures = [e for e in ALL_EXPOSURES if e in df.columns]

    susceptibility_map = {}
    if apoe_col is not None:
        susceptibility_map["APOE4_carrier"] = {"col": apoe_col, "type": "binary", "label": "APOE e4 carrier"}
    if apoe_count_col is not None:
        susceptibility_map["APOE4_count"] = {"col": apoe_count_col, "type": "continuous", "label": "APOE e4 count"}
    if ad_prs_std_col is not None:
        susceptibility_map["AD_PRS_standard"] = {"col": ad_prs_std_col, "type": "continuous", "label": "AD PRS standard", "quartile": "ad_prs_std_quartile", "highlow": "ad_prs_std_high_vs_low"}
    if ad_prs_enh_col is not None:
        susceptibility_map["AD_PRS_enhanced"] = {"col": ad_prs_enh_col, "type": "continuous", "label": "AD PRS enhanced", "quartile": "ad_prs_enh_quartile", "highlow": "ad_prs_enh_high_vs_low"}

    modifiers_to_run = [m for m in config.modifiers if m in susceptibility_map]

    with timed("build_covariates"):
        genetic_pc_cols = [c for c in df.columns if any(c.startswith(pref) for pref in GENETIC_PC_PREFIXES)]
        genetic_qc_cols = [c for c in GENETIC_QC_CANDIDATES if c in df.columns]
        lean_candidates = [c for c in LEAN_COVAR_CANDIDATES if c in df.columns] + genetic_pc_cols
        expanded_candidates = [c for c in EXPANDED_COVAR_CANDIDATES if c in df.columns] + genetic_pc_cols + genetic_qc_cols

        lean_mat, lean_actions = build_covariate_matrix(df, lean_candidates, pc_limit=config.max_pcs_lean, prefix="lean")
        expanded_mat, expanded_actions = build_covariate_matrix(df, expanded_candidates, pc_limit=config.max_pcs_expanded, prefix="expd")

        lean_actions["covariate_set"] = "lean"
        expanded_actions["covariate_set"] = "expanded"
        covar_actions = pd.concat([lean_actions, expanded_actions], ignore_index=True)

        covars_to_concat = []
        if "lean" in config.covariate_sets:
            covars_to_concat.append(lean_mat)
        if "expanded" in config.covariate_sets:
            covars_to_concat.append(expanded_mat)
        df = pd.concat([df.reset_index(drop=True)] + [x.reset_index(drop=True) for x in covars_to_concat], axis=1).copy()

        covar_sets = {"lean": list(lean_mat.columns), "expanded": list(expanded_mat.columns)}

    save_table(covar_actions, "covariate_encoding_actions.csv", config.tabdir)
    save_table(pd.DataFrame([{"covariate_set": k, "n_covariates": len(v), "covariates": ";".join(v)} for k, v in covar_sets.items()]), "covariate_sets.csv", config.tabdir)

    for e in available_exposures:
        df[e] = ensure_binary(df[e])
    for spec in available_outcomes.values():
        df[spec["time"]] = safe_numeric(df[spec["time"]])
        df[spec["event"]] = ensure_binary(df[spec["event"]])
    for m in modifiers_to_run:
        df[susceptibility_map[m]["col"]] = safe_numeric(df[susceptibility_map[m]["col"]])

    log(f"[CHECK] n={len(df)}")
    log(f"[CHECK] available_outcomes={list(available_outcomes.keys())}")
    log(f"[CHECK] available_exposures={available_exposures}")
    log(f"[CHECK] modifiers_to_run={modifiers_to_run}")
    log(f"[CHECK] lean_covars={len(covar_sets['lean'])}, expanded_covars={len(covar_sets['expanded'])}")

    var_diag_rows = []
    diag_cols = list(set(
        available_exposures
        + [x["time"] for x in available_outcomes.values()]
        + [x["event"] for x in available_outcomes.values()]
        + [susceptibility_map[m]["col"] for m in modifiers_to_run]
    ))
    for c in diag_cols:
        if c in df.columns:
            sx = safe_numeric(df[c])
            var_diag_rows.append({
                "column": c,
                "nonmissing": int(sx.notna().sum()),
                "unique_n": int(sx.nunique(dropna=True)),
                "mean": float(sx.mean(skipna=True)) if sx.notna().sum() else np.nan,
                "sd": float(sx.std(skipna=True)) if sx.notna().sum() else np.nan,
            })
    save_table(pd.DataFrame(var_diag_rows), "variable_diagnostics.csv", config.tabdir)

    pd.DataFrame([{
        "n_total": len(df),
        "available_outcomes": ";".join(available_outcomes.keys()),
        "available_exposures": ";".join(available_exposures),
        "apoe_col": apoe_col,
        "apoe_count_col": apoe_count_col,
        "apoe_genotype_col": apoe_genotype_col,
        "ad_prs_std_col": ad_prs_std_col,
        "ad_prs_enh_col": ad_prs_enh_col,
        "modifiers_run": ";".join(modifiers_to_run),
        "engines": ";".join(config.engines),
        "covariate_sets": ";".join(config.covariate_sets),
        "lean_covariates_n": len(covar_sets["lean"]),
        "expanded_covariates_n": len(covar_sets["expanded"]),
    }]).to_csv(os.path.join(config.qcdir, "preflight_summary.csv"), index=False)

    if config.preflight:
        log("[DONE] Preflight requested. No models were fitted.")
        return

    diag_path = os.path.join(config.tabdir, "fast_model_diagnostics_live.csv")
    diag_jsonl_path = os.path.join(config.tabdir, "fast_model_diagnostics_live.jsonl")
    if config.resume and (os.path.exists(diag_jsonl_path) or os.path.exists(diag_path)):
        try:
            if os.path.exists(diag_jsonl_path):
                old_diag = read_live_diagnostics_jsonl(diag_jsonl_path)
            else:
                old_diag = pd.read_csv(diag_path, engine="python", on_bad_lines="skip")
            if "model_key" in old_diag.columns and "status" in old_diag.columns:
                COMPLETED_KEYS = set(old_diag.loc[old_diag["status"].isin(["fit_ok", "fit_failed", "n_below_threshold", "events_below_threshold"]), "model_key"].dropna().astype(str))
                log(f"[RESUME] completed model keys={len(COMPLETED_KEYS)}")
        except Exception as e:
            log(f"[RESUME-WARN] failed to read existing diagnostics: {e}")

    # B4a Genetic main effects
    log("[RUN] Genetic main effects")
    main_rows = []
    diagnostics = []
    for engine in config.engines:
        for covar_set_name in config.covariate_sets:
            covars = covar_sets.get(covar_set_name, [])
            for outcome_name, spec in available_outcomes.items():
                for modifier_name in modifiers_to_run:
                    gcol = susceptibility_map[modifier_name]["col"]
                    other_susc = [susceptibility_map[m]["col"] for m in modifiers_to_run if m != modifier_name and susceptibility_map[m]["col"] in df.columns]
                    if modifier_name.startswith("APOE"):
                        other_susc = [x for x in other_susc if "apoe" not in x.lower()]
                    if modifier_name.startswith("AD_PRS"):
                        other_susc = [x for x in other_susc if "prs" not in x.lower()]
                    covars_use = list(dict.fromkeys(covars + other_susc))
                    rows, diag = fit_fast_event_model(
                        df, spec["time"], spec["event"], [gcol], covars_use,
                        label_dict={
                            "analysis": "genetic_main_effect",
                            "covariate_set": covar_set_name,
                            "outcome": outcome_name,
                            "outcome_label": spec["label"],
                            "exposure": "",
                            "exposure_abbr": "",
                            "exposure_label": "",
                            "modifier": modifier_name,
                            "modifier_label": susceptibility_map[modifier_name]["label"],
                        },
                        config=config,
                        engine=engine,
                    )
                    diagnostics.append(diag)
                    if rows is None:
                        break
                    for r in rows:
                        r["predictor"] = modifier_name
                        r["predictor_col"] = gcol
                        main_rows.append(r)
                        append_table(r, "genetic_main_effects_live.csv", config.tabdir)
                    if config.max_models and MODEL_COUNTER >= config.max_models:
                        break
                if config.max_models and MODEL_COUNTER >= config.max_models:
                    break
            if config.max_models and MODEL_COUNTER >= config.max_models:
                break
        if config.max_models and MODEL_COUNTER >= config.max_models:
            break

    main_df = pd.DataFrame(main_rows)
    if len(main_df):
        main_df = fdr_add(main_df, "p_value", ["engine", "outcome", "covariate_set"], "fdr")
    save_table(main_df, "genetic_main_effects.csv", config.tabdir)

    # B4b GI x genetic interactions
    log("[RUN] GI-by-genetic interactions")
    interaction_rows = []
    if not (config.max_models and MODEL_COUNTER >= config.max_models):
        for engine in config.engines:
            for covar_set_name in config.covariate_sets:
                covars = covar_sets.get(covar_set_name, [])
                for outcome_name, spec in available_outcomes.items():
                    for exposure in available_exposures:
                        for modifier_name in modifiers_to_run:
                            gcol = susceptibility_map[modifier_name]["col"]
                            inter_col = f"inter__{clean_string_name(exposure)}__{clean_string_name(modifier_name)}"
                            df[inter_col] = safe_numeric(df[exposure]) * safe_numeric(df[gcol])
                            other_susc = [susceptibility_map[m]["col"] for m in modifiers_to_run if m != modifier_name and susceptibility_map[m]["col"] in df.columns]
                            if modifier_name.startswith("APOE"):
                                other_susc = [x for x in other_susc if "apoe" not in x.lower()]
                            if modifier_name.startswith("AD_PRS"):
                                other_susc = [x for x in other_susc if "prs" not in x.lower()]
                            covars_use = list(dict.fromkeys(covars + other_susc))
                            rows, diag = fit_fast_event_model(
                                df, spec["time"], spec["event"], [exposure, gcol, inter_col], covars_use,
                                label_dict={
                                    "analysis": "gi_genetic_interaction",
                                    "covariate_set": covar_set_name,
                                    "outcome": outcome_name,
                                    "outcome_label": spec["label"],
                                    "exposure": exposure,
                                    "exposure_abbr": EXPOSURE_ABBR.get(exposure, exposure),
                                    "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                                    "modifier": modifier_name,
                                    "modifier_label": susceptibility_map[modifier_name]["label"],
                                    "interaction_col": inter_col,
                                },
                                config=config,
                                engine=engine,
                            )
                            diagnostics.append(diag)
                            if rows is None:
                                break
                            for r in rows:
                                if r["term"] == inter_col:
                                    r["interaction_term"] = inter_col
                                    r["effect_interaction"] = r["effect"]
                                    r["ci_low_interaction"] = r["ci_low"]
                                    r["ci_high_interaction"] = r["ci_high"]
                                    r["p_interaction"] = r["p_value"]
                                    interaction_rows.append(r)
                                    append_table(r, "gi_genetic_interactions_live.csv", config.tabdir)
                            if config.max_models and MODEL_COUNTER >= config.max_models:
                                break
                        if config.max_models and MODEL_COUNTER >= config.max_models:
                            break
                    if config.max_models and MODEL_COUNTER >= config.max_models:
                        break
                if config.max_models and MODEL_COUNTER >= config.max_models:
                    break
            if config.max_models and MODEL_COUNTER >= config.max_models:
                break

    interaction_df = pd.DataFrame(interaction_rows)
    if len(interaction_df):
        interaction_df = fdr_add(interaction_df, "p_interaction", ["engine", "outcome", "covariate_set"], "fdr_interaction")
    save_table(interaction_df, "gi_genetic_interactions.csv", config.tabdir)

    # B4c Stratified GI-outcome associations
    strat_rows = []
    if not config.no_stratified and not (config.max_models and MODEL_COUNTER >= config.max_models):
        log("[RUN] Stratified GI-outcome associations")
        for engine in config.engines:
            for outcome_name, spec in available_outcomes.items():
                for exposure in available_exposures:
                    for modifier_name in modifiers_to_run:
                        info = susceptibility_map[modifier_name]
                        gcol = info["col"]
                        covars = covar_sets.get("lean", [])
                        strata = []
                        strat_col = None
                        if modifier_name == "APOE4_carrier" and info["type"] == "binary":
                            strata = [(0, "APOE e4 non-carrier"), (1, "APOE e4 carrier")]
                            strat_col = gcol
                        elif modifier_name in ["AD_PRS_standard", "AD_PRS_enhanced"]:
                            qcol = info.get("quartile")
                            if qcol and qcol in df.columns:
                                strata = [(1, "Q1"), (2, "Q2"), (3, "Q3"), (4, "Q4")]
                                strat_col = qcol
                        if not strata or strat_col is None:
                            continue
                        for lv, lv_label in strata:
                            sub = df.loc[safe_numeric(df[strat_col]) == lv].copy()
                            covars_use = [c for c in covars if c not in [exposure, gcol, strat_col]]
                            rows, diag = fit_fast_event_model(
                                sub, spec["time"], spec["event"], [exposure], covars_use,
                                label_dict={
                                    "analysis": "stratified_gi_outcome_by_genetics",
                                    "covariate_set": "lean",
                                    "outcome": outcome_name,
                                    "outcome_label": spec["label"],
                                    "exposure": exposure,
                                    "exposure_abbr": EXPOSURE_ABBR.get(exposure, exposure),
                                    "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                                    "modifier": modifier_name,
                                    "modifier_label": info["label"],
                                    "stratifier_col": strat_col,
                                    "stratum": lv_label,
                                    "stratum_value": lv,
                                },
                                config=config,
                                engine=engine,
                                min_n=config.min_n_strat,
                                min_events=config.min_events_strat,
                            )
                            diagnostics.append(diag)
                            if rows is None:
                                break
                            for r in rows:
                                r["stratum_effect"] = r["effect"]
                                strat_rows.append(r)
                                append_table(r, "stratified_gi_outcome_by_genetics_live.csv", config.tabdir)
                            if config.max_models and MODEL_COUNTER >= config.max_models:
                                break
                        if config.max_models and MODEL_COUNTER >= config.max_models:
                            break
                    if config.max_models and MODEL_COUNTER >= config.max_models:
                        break
                if config.max_models and MODEL_COUNTER >= config.max_models:
                    break
    strat_df = pd.DataFrame(strat_rows)
    if len(strat_df):
        strat_df = fdr_add(strat_df, "p_value", ["engine", "outcome", "modifier"], "fdr")
    save_table(strat_df, "stratified_gi_outcome_by_genetics.csv", config.tabdir)

    # B4d Joint GI-genetic risk gradient
    log("[RUN] Joint GI-genetic risk gradient")
    joint_rows = []

    def add_joint_rows(outcome_name: str, spec: dict, exposure: str, modifier_name: str) -> None:
        info = susceptibility_map[modifier_name]
        gcol = info["col"]
        if modifier_name == "APOE4_carrier":
            genetic_status = ensure_binary(df[gcol])
            status_label = {0: "Genetic low", 1: "Genetic high"}
        elif modifier_name in ["AD_PRS_standard", "AD_PRS_enhanced"]:
            hl = info.get("highlow")
            if not hl or hl not in df.columns:
                return
            genetic_status = safe_numeric(df[hl])
            status_label = {0: "Genetic low", 1: "Genetic high"}
        else:
            return

        e = ensure_binary(df[exposure])
        event = ensure_binary(df[spec["event"]])
        timev = safe_numeric(df[spec["time"]])
        base = pd.DataFrame({"exposure": e, "genetic_status": genetic_status, "event": event, "time": timev}).dropna()
        base = base.loc[base["time"] > 0].copy()
        if len(base) == 0:
            return
        for gi_val in [0, 1]:
            for g_val in [0, 1]:
                ss = base.loc[(base["exposure"] == gi_val) & (base["genetic_status"] == g_val)]
                if len(ss) == 0:
                    continue
                n = int(len(ss))
                ev = int(ss["event"].sum())
                person_years = float(ss["time"].sum())
                joint_rows.append({
                    "analysis": "joint_gi_genetic_risk_gradient",
                    "outcome": outcome_name,
                    "outcome_label": spec["label"],
                    "exposure": exposure,
                    "exposure_abbr": EXPOSURE_ABBR.get(exposure, exposure),
                    "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                    "modifier": modifier_name,
                    "modifier_label": info["label"],
                    "joint_group": f"GI{gi_val}_GEN{g_val}",
                    "gi_status": "GI+" if gi_val == 1 else "GI-",
                    "genetic_status": status_label.get(g_val, str(g_val)),
                    "n": n,
                    "n_event": ev,
                    "event_rate": ev / n if n else np.nan,
                    "person_years": person_years,
                    "incidence_rate_per_1000py": ev / person_years * 1000 if person_years > 0 else np.nan,
                })

    for outcome_name, spec in available_outcomes.items():
        for exposure in available_exposures:
            for modifier_name in modifiers_to_run:
                if modifier_name in ["APOE4_carrier", "AD_PRS_enhanced", "AD_PRS_standard"]:
                    add_joint_rows(outcome_name, spec, exposure, modifier_name)
    joint_df = pd.DataFrame(joint_rows)
    save_table(joint_df, "joint_gi_genetic_risk_gradient.csv", config.tabdir)

    # B4e shared-liability placeholder/diagnostics
    shared_rows = []
    shared_diag = []
    if not config.no_shared_liability:
        log("[RUN] GI PRS shared-liability support diagnostics")
        for exposure, prs_candidates in GI_PRS_CANDIDATES.items():
            if exposure not in df.columns:
                shared_diag.append({"exposure": exposure, "status": "missing_exposure"})
                continue
            prs_col = detect_column(prs_candidates, df.columns)
            if prs_col is None:
                shared_diag.append({"exposure": exposure, "status": "missing_gi_prs"})
                continue
            df[prs_col] = zscore(winsorize(df[prs_col]))
            # Intentionally diagnostic-only unless GI-specific PRS are present and a
            # separate case-control model is explicitly planned.
            shared_diag.append({"exposure": exposure, "prs_col": prs_col, "status": "prs_detected_not_modeled_by_default"})
    shared_df = pd.DataFrame(shared_rows)
    save_table(shared_df, "gi_prs_shared_liability_support.csv", config.tabdir)
    save_table(pd.DataFrame(shared_diag), "gi_prs_shared_liability_diagnostics.csv", config.tabdir)

    # Optional Cox confirmation
    cox_df = pd.DataFrame()
    if config.cox_confirm and len(interaction_df) > 0:
        log("[RUN] Optional Cox confirmation for prioritized interactions")
        candidates = interaction_df.copy()
        candidates = candidates.loc[candidates["covariate_set"].eq("lean")].copy()
        candidates = candidates.sort_values("p_interaction").head(config.max_cox_confirm)
        cox_rows = []
        for _, r in candidates.iterrows():
            outcome_name = r["outcome"]
            exposure = r["exposure"]
            modifier_name = r["modifier"]
            spec = OUTCOME_SPECS.get(outcome_name)
            if spec is None or modifier_name not in susceptibility_map:
                continue
            gcol = susceptibility_map[modifier_name]["col"]
            inter_col = r["interaction_term"]
            if inter_col not in df.columns:
                df[inter_col] = safe_numeric(df[exposure]) * safe_numeric(df[gcol])
            rows = fit_cox_confirm(
                df, spec["time"], spec["event"], [exposure, gcol, inter_col], covar_sets.get("lean", []),
                label_dict={
                    "analysis": "cox_confirmatory_top_interaction",
                    "outcome": outcome_name,
                    "outcome_label": spec["label"],
                    "exposure": exposure,
                    "exposure_abbr": EXPOSURE_ABBR.get(exposure, exposure),
                    "exposure_label": EXPOSURE_LABELS.get(exposure, exposure),
                    "modifier": modifier_name,
                    "modifier_label": susceptibility_map[modifier_name]["label"],
                    "interaction_col": inter_col,
                    "screen_engine": r.get("engine"),
                    "screen_effect_interaction": r.get("effect_interaction"),
                    "screen_p_interaction": r.get("p_interaction"),
                },
                config=config,
            )
            for rr in rows:
                if rr["term"] == inter_col:
                    rr["hr_interaction"] = rr["hr"]
                    rr["p_interaction"] = rr["p_value"]
                    cox_rows.append(rr)
                    append_table(rr, "cox_confirmatory_top_interactions_live.csv", config.tabdir)
        cox_df = pd.DataFrame(cox_rows)
    save_table(cox_df, "cox_confirmatory_top_interactions.csv", config.tabdir)

    # Final diagnostics and report
    diag_jsonl_path = os.path.join(config.tabdir, "fast_model_diagnostics_live.jsonl")
    if os.path.exists(diag_jsonl_path):
        diag_df = read_live_diagnostics_jsonl(diag_jsonl_path)
    elif len(diagnostics) > 0:
        diag_df = pd.DataFrame(diagnostics)
    else:
        diag_df = pd.DataFrame()
    save_table(diag_df, "fast_model_diagnostics_final.csv", config.tabdir)

    if len(diag_df):
        status_cols = [c for c in ["analysis", "engine", "covariate_set", "outcome", "status"] if c in diag_df.columns]
        status_summary = diag_df.groupby(status_cols, dropna=False).size().reset_index(name="n") if status_cols else pd.DataFrame()
        timing_summary = diag_df.copy()
        if "elapsed_sec" in timing_summary.columns:
            timing_summary["elapsed_sec"] = safe_numeric(timing_summary["elapsed_sec"])
            timing_summary = timing_summary.sort_values("elapsed_sec", ascending=False)
    else:
        status_summary = pd.DataFrame()
        timing_summary = pd.DataFrame()

    save_table(status_summary, "model_status_summary.csv", config.tabdir)
    save_table(timing_summary.head(100), "model_timing_top100.csv", config.tabdir)

    cohort_summary = pd.DataFrame([{
        "n_total": len(df),
        "available_outcomes": ";".join(available_outcomes.keys()),
        "available_exposures": ";".join(available_exposures),
        "modifiers_run": ";".join(modifiers_to_run),
        "engines": ";".join(config.engines),
        "covariate_sets": ";".join(config.covariate_sets),
        "main_effect_rows": len(main_df),
        "interaction_rows": len(interaction_df),
        "stratified_rows": len(strat_df),
        "joint_gradient_rows": len(joint_df),
        "cox_confirmatory_rows": len(cox_df),
        "elapsed_minutes": round((time.time() - T0) / 60, 2),
    }])
    save_table(cohort_summary, "cohort_and_result_summary.csv", config.tabdir)

    metadata = {
        "script": "genetic_interaction_analysis.py",
        "master_path": config.master,
        "outdir": config.outdir,
        "interpretation": "genetic susceptibility context and effect modification, not mediation",
        "design": {
            "fast_screen": "Poisson/logistic GLM for interaction screening",
            "poisson_model": "event ~ exposure + genetic + exposure:genetic + covariates + offset(log follow-up time)",
            "cox_confirmatory": "optional Cox models only for prioritized top interactions",
        },
        "detected_columns": {
            "apoe_col": apoe_col,
            "apoe_count_col": apoe_count_col,
            "apoe_genotype_col": apoe_genotype_col,
            "ad_prs_std_col": ad_prs_std_col,
            "ad_prs_enh_col": ad_prs_enh_col,
        },
        "run_settings": {
            "outcomes": config.outcomes,
            "modifiers": config.modifiers,
            "covariate_sets": config.covariate_sets,
            "engines": config.engines,
            "max_models": config.max_models,
            "cox_confirm": config.cox_confirm,
        },
    }
    with open(os.path.join(config.outdir, "run_metadata.json"), "w", encoding="utf-8") as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    report = [
        "Genetic susceptibility and GI-by-genetic interaction screen finished.",
        "=" * 80,
        f"N total = {len(df):,}",
        f"Outcomes run = {list(available_outcomes.keys())}",
        f"Exposures run = {available_exposures}",
        f"Modifiers run = {modifiers_to_run}",
        f"Engines = {config.engines}",
        f"Covariate sets = {config.covariate_sets}",
        "",
        "Main outputs:",
        " - tables/genetic_main_effects.csv",
        " - tables/gi_genetic_interactions.csv",
        " - tables/stratified_gi_outcome_by_genetics.csv",
        " - tables/joint_gi_genetic_risk_gradient.csv",
        " - tables/cox_confirmatory_top_interactions.csv",
        "",
        "Interpretation:",
        " - Genetic variables are used as susceptibility/effect-modification factors, not mediators.",
        " - Fast screening uses person-time-adjusted Poisson GLM unless another engine is selected.",
        " - The interaction term estimates departure from multiplicativity on the incidence-rate scale for Poisson models.",
        " - Optional Cox confirmation is reserved for prioritized signals.",
        "",
        "Rows:",
        f" - Genetic main effects: {len(main_df)}",
        f" - GI x genetic interactions: {len(interaction_df)}",
        f" - Stratified models: {len(strat_df)}",
        f" - Joint gradients: {len(joint_df)}",
        f" - Cox confirmations: {len(cox_df)}",
        f" - Elapsed minutes: {round((time.time() - T0) / 60, 2)}",
    ]
    if len(status_summary):
        report.extend(["", "Model status summary:", status_summary.to_string(index=False)])
    with open(os.path.join(config.outdir, "analysis_report.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(report))

    log("=" * 90)
    log("[DONE] Outputs saved to:")
    log(config.outdir)
    log("[READ FIRST]")
    log(os.path.join(config.outdir, "analysis_report.txt"))
    log(os.path.join(config.tabdir, "cohort_and_result_summary.csv"))
    log(os.path.join(config.tabdir, "gi_genetic_interactions.csv"))
    log(os.path.join(config.tabdir, "joint_gi_genetic_risk_gradient.csv"))
    log("=" * 90)


def main() -> None:
    config = parse_args()
    run_analysis(config)


if __name__ == "__main__":
    main()
