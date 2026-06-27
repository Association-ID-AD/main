Multimodal trajectories of intestinal antecedents and Alzheimer's disease risk
================================================================================

Project overview
----------------
This repository provides the analysis code for the study:

    Multimodal trajectories of intestinal antecedents and Alzheimer's disease
    risk across multinational longitudinal cohorts

The study evaluates whether clinically recorded intestinal diseases can serve as
routine-care antecedent signals for Alzheimer's disease (AD) risk enrichment.
The analysis integrates longitudinal disease records, incident AD and dementia
outcomes, blood biochemistry, NMR metabolomics, cognitive phenotypes, brain
imaging-derived phenotypes, and genetic susceptibility profiles.

The primary discovery cohort is UK Biobank. External cross-cohort support
analyses were performed using harmonized ageing cohorts including HRS, CHARLS,
LASI, and SHARE. The full multi-cohort framework contains 565,307 participants,
including 366,466 eligible UK Biobank participants and 198,841 participants from
external ageing cohorts.

The main scientific questions are:

1. Are intestinal diseases associated with future AD and dementia risk?
2. Are these associations disease-specific rather than a generic digestive
   disease signal?
3. Which metabolic, brain imaging, cognitive, and genetic profiles contextualize
   the intestinal disease-AD association?
4. Can these multimodal signals be translated into individualized 1-20 year AD
   risk trajectories and warning windows?

The code is organized as a sequence of analysis scripts. Each script can be run
independently after the required preprocessed input tables have been generated.


Repository structure
--------------------
A recommended repository layout is:

    .
    ├── README.txt
    ├── scripts/
    │   ├── gi_ad_analysis.py
    │   ├── metabolic_bridge_analysis.py
    │   ├── metabolic_neural_bridge_analysis.py
    │   ├── metabolic_cognition_bridge_analysis.py
    │   ├── genetic_interaction_analysis.py
    │   └── risk_warning_window_patch.py
    ├── data/
    │   └── README_data.txt
    └── results/

The scripts are designed to write analysis outputs to user-specified output
directories. Raw cohort data are not included in this repository because access
to UK Biobank, HRS, CHARLS, LASI, and SHARE is controlled by the respective
cohort data-access policies.


Required input files
--------------------
The public scripts assume that cohort-specific preprocessing has already produced
the following analysis-ready tables.

1. master_preprocessed_with_genetics.csv

   Participant-level master table containing:
   - participant ID: eid
   - baseline date: time0
   - intestinal disease exposure indicators
   - intestinal disease dates or durations where available
   - AD, dementia, and MCI-related outcomes
   - follow-up times
   - baseline demographics and socioeconomic covariates
   - lifestyle, clinical, medication, and family-history covariates
   - APOE and AD polygenic risk score variables where available
   - genetic principal components and genetic quality-control fields where used

2. metabolic_prepared.csv

   Participant-level metabolic table containing:
   - eid
   - blood biochemistry markers
   - NMR metabolomic markers
   - optional metabolic assessment date fields

3. brain_idp_prepared.csv

   Participant-level brain imaging table containing:
   - eid
   - imaging_date
   - structural MRI-derived phenotypes
   - diffusion MRI-derived phenotypes
   - other UK Biobank imaging-derived phenotypes after quality control

4. cognition_prepared.csv

   Participant-level cognitive phenotype table containing:
   - eid
   - cognition_date
   - cognitive task-derived phenotypes
   - optional cognition instance indicator

5. Optional label-map CSV

   A two-column CSV file with:
   - feature
   - label

   This can be used to replace raw UK Biobank feature IDs with readable marker
   labels in output tables.

The exact preprocessing scripts that create these files are project-specific and
are not included here unless they are explicitly added to the repository.


Intestinal disease definitions
------------------------------
The main analysis uses eight grouped intestinal antecedent categories.

    DIV   exp_diverticular                  K57
    IBD   exp_ibd                           K50, K51
    IBS   exp_ibs                           K58
    OCID  exp_other_chronic_intestinal      K52, K55, K56, K63
    APD   exp_appendiceal                   K35-K38
    MAL   exp_malabsorption                 K90
    OFID  exp_other_functional_intestinal   K59
    ARD   exp_anorectal                     K60-K64

The scripts reconstruct grouped exposures from ICD-level exposure columns when
the grouped variables are not already present in the master table.

For example:
    exp_ibd = max(exp_k50, exp_k51)
    exp_appendiceal = max(exp_k35, exp_k36, exp_k37, exp_k38)

For time-ordered bridge analyses, exposure dates are used when available. The
scripts search for date columns using common naming patterns such as:

    date_k57, dt_k57, date_diverticular, dt_diverticular

If exposure dates are available, GI disease is required to occur before imaging
or cognition assessment in the corresponding bridge analyses.


Software requirements
---------------------
The scripts are written for Python 3.9 or later. The main Python dependencies are:

    numpy
    pandas
    scipy
    statsmodels
    lifelines
    scikit-learn
    matplotlib
    torch
    torchtuples
    pycox

Not all packages are required for every script. For example, pycox is required
only when DeepHit-based risk modeling is enabled in gi_ad_analysis.py.

A minimal installation command is:

    pip install numpy pandas scipy statsmodels lifelines scikit-learn matplotlib

For DeepHit support:

    pip install torch torchtuples pycox

For reproducibility, users are encouraged to create a dedicated conda or virtual
environment and record package versions.


Script 1: gi_ad_analysis.py
---------------------------
Purpose:
    Main intestinal disease-AD/dementia longitudinal association analysis.

Core analyses:
    - construction of grouped intestinal antecedent exposures
    - construction or use of AD, dementia, MCI-like, and MCI-or-AD outcomes
    - staged Cox proportional hazards models
    - optional baseline logistic support analysis
    - optional DeepHit-based absolute risk prediction
    - reverse-causation sensitivity analyses
    - recent-exposure exclusion sensitivity analyses
    - detection-bias sensitivity using bowel cancer screening
    - subgroup and interaction analyses

Recommended command:

    python scripts/gi_ad_analysis.py \
        --input data/master_preprocessed_with_genetics.csv \
        --outdir results/gi_ad_analysis

To skip DeepHit:

    python scripts/gi_ad_analysis.py \
        --input data/master_preprocessed_with_genetics.csv \
        --outdir results/gi_ad_analysis \
        --skip-deephit

To run DeepHit:

    python scripts/gi_ad_analysis.py \
        --input data/master_preprocessed_with_genetics.csv \
        --outdir results/gi_ad_analysis \
        --run-deephit

Main output directories:
    results/gi_ad_analysis/main/
    results/gi_ad_analysis/sensitivity/
    results/gi_ad_analysis/subgroup/
    results/gi_ad_analysis/deephit/
    results/gi_ad_analysis/qc/

Representative outputs:
    main/main_strict_ad_cox_results.csv
    main/gi_ad_main_cox_all.csv
    sensitivity/reverse_causation_cox_all.csv
    sensitivity/exposure_recency_cox_all.csv
    subgroup/subgroup_strict_ad_grouped_cox.csv
    deephit/deephit_risk_contrast_all.csv
    qc/outcome_qc.csv
    run_metadata.json
    analysis_report.txt

Model interpretation:
    Cox models estimate hazard ratios for each intestinal antecedent relative
    to participants without the corresponding disease. DeepHit outputs are
    supplementary horizon-specific model-based absolute risk contrasts and are
    not used as primary etiological estimates.


Script 2: metabolic_bridge_analysis.py
--------------------------------------
Purpose:
    Systemic metabolic bridge analysis linking intestinal antecedents to AD risk.

Scientific path:
    GI disease -> blood/NMR metabolic marker -> AD risk

Core analyses:
    - Stage 1: GI exposure to individual metabolic marker using robust OLS
    - Stage 2: metabolic marker to incident AD using Cox models
    - horizon-specific logistic support analysis
    - all-cause dementia support analysis
    - candidate metabolic bridge marker selection
    - restricted g-computation at the 10-year AD horizon
    - PCA-based metabolic module construction for selected marker sets
    - volcano and bridge-effect figures

Recommended command:

    python scripts/metabolic_bridge_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --metabolic data/metabolic_prepared.csv \
        --outdir results/metabolic_bridge_analysis \
        --label-map data/metabolic_feature_labels.csv

The label map is optional.

Representative outputs:
    tables/screen1_all.csv
    tables/screen2_cox_all.csv
    tables/screen2_logit_all.csv
    tables/candidate_markers_all.csv
    tables/bridge_marker_all.csv
    tables/bridge_module_all.csv
    tables/metabolic_module_scores_for_downstream.csv
    tables/metabolic_bridge_summary_table.csv
    figures/volcano_stage1_*.png
    figures/bridge_marker_*.png
    qc/covariate_blocks.csv
    run_metadata.json
    analysis_report.txt

Bridge interpretation:
    A metabolic bridge marker is prioritized when it is associated with both
    intestinal antecedent status and AD-related risk after the specified
    screening thresholds. Restricted g-computation estimates total, direct, and
    indirect risk contrasts at a fixed prediction horizon. These estimates are
    interpreted as pathway-consistent observational signatures, not definitive
    causal mediation.


Script 3: metabolic_neural_bridge_analysis.py
---------------------------------------------
Purpose:
    Marker-level metabolic-neural bridge analysis linking intestinal antecedents
    to brain imaging-derived phenotypes through metabolic markers.

Scientific path:
    GI disease -> metabolic marker -> brain imaging-derived phenotype

Core analyses:
    - direct GI to brain IDP screening
    - import and standardization of metabolic bridge candidates from
      metabolic_bridge_analysis.py
    - marker-level GI -> marker -> IDP product-of-coefficients bridge
    - marker-to-IDP screen adjusted for GI exposure and covariates
    - generation of significant bridge paths
    - brain map source table generation
    - IDP -> AD and IDP -> dementia support Cox analyses
    - AD overlay onto metabolic-neural bridge paths

Recommended command:

    python scripts/metabolic_neural_bridge_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --brain data/brain_idp_prepared.csv \
        --metabolic data/metabolic_prepared.csv \
        --metabolic-bridge-dir results/metabolic_bridge_analysis \
        --outdir results/metabolic_neural_bridge_analysis

Optional IDP panel file:

    python scripts/metabolic_neural_bridge_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --brain data/brain_idp_prepared.csv \
        --metabolic data/metabolic_prepared.csv \
        --metabolic-bridge-dir results/metabolic_bridge_analysis \
        --outdir results/metabolic_neural_bridge_analysis \
        --panel-file data/selected_idps.txt

Representative outputs:
    tables/direct_gi_to_idp_screen_full.csv
    tables/idp_candidates_by_exposure.csv
    tables/marker_candidates_from_metabolic_bridge.csv
    tables/marker_to_idp_screen.csv
    tables/marker_idp_bridge_effects.csv
    tables/marker_idp_bridge_effects_significant.csv
    tables/figure_source_paths.csv
    tables/brain_map_source_data.csv
    tables/idp_to_ad_cox_all.csv
    tables/idp_to_dementia_cox_all.csv
    tables/marker_idp_bridge_with_ad_overlay.csv
    qc/covariate_blocks.csv
    run_metadata.json
    analysis_report.txt

Temporality rule:
    If exposure dates are available, exposed participants are retained only when
    the GI exposure date precedes the imaging date. Unexposed participants are
    retained as reference participants.

Bridge interpretation:
    The path effect is alpha * beta, where alpha is the coefficient for
    GI -> metabolic marker and beta is the coefficient for marker -> IDP
    adjusted for GI status and covariates. This is a marker-level bridge signal,
    not proof of causal mediation.


Script 4: metabolic_cognition_bridge_analysis.py
------------------------------------------------
Purpose:
    Marker-level metabolic-cognition bridge analysis linking intestinal
    antecedents to cognitive phenotypes through metabolic markers.

Scientific path:
    GI disease -> metabolic marker -> cognition phenotype

Core analyses:
    - Stage 1: GI to metabolic marker robust OLS
    - Stage 2: metabolic marker to cognition robust OLS adjusted for GI
    - alpha * beta path-effect table
    - Sobel statistic and false discovery rate correction
    - Sankey/alluvial source table
    - optional continuous-outcome g-computation mediation for selected paths
    - cognition -> AD support Cox analysis
    - cognition -> dementia support Cox analysis
    - cognition -> 10-year AD logistic support analysis
    - summary plots and QC tables

Recommended command:

    python scripts/metabolic_cognition_bridge_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --cognition data/cognition_prepared.csv \
        --metabolic data/metabolic_prepared.csv \
        --outdir results/metabolic_cognition_bridge_analysis

Resume from existing stage tables:

    python scripts/metabolic_cognition_bridge_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --cognition data/cognition_prepared.csv \
        --metabolic data/metabolic_prepared.csv \
        --outdir results/metabolic_cognition_bridge_analysis \
        --resume-from-stage-tables

Stop after path table generation:

    python scripts/metabolic_cognition_bridge_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --cognition data/cognition_prepared.csv \
        --metabolic data/metabolic_prepared.csv \
        --outdir results/metabolic_cognition_bridge_analysis \
        --resume-from-stage-tables \
        --stop-after-paths

Representative outputs:
    tables/stage1_gi_to_marker_all.csv
    tables/stage2_marker_to_cognition_all.csv
    tables/marker_cognition_bridge_paths_all.csv
    tables/marker_cognition_bridge_paths_significant.csv
    tables/sankey_source_paths.csv
    tables/gi_marker_cognition_mediation_all.csv
    tables/support_cognition_to_ad_cox.csv
    tables/support_cognition_to_dementia_cox.csv
    tables/support_cognition_to_ad_logit_10y.csv
    tables/marker_cognition_summary_table.csv
    figures/top_marker_cognition_paths_product_forest.png
    qc/covariate_blocks.csv
    run_metadata.json
    analysis_report.txt

Cognitive phenotype handling:
    The script uses raw cognitive phenotypes and also constructs selected PCA
    composites when the required source tests are available, including:
    - cognition_speed_pc1
    - cognition_accuracy_pc1
    - cognition_error_pc1

Temporality rule:
    If GI exposure dates are available, exposed participants are retained only
    when the GI exposure date occurs on or before the cognition assessment date.


Script 5: genetic_interaction_analysis.py
-----------------------------------------
Purpose:
    Genetic susceptibility and effect-modification analysis.

Scientific question:
    Does AD genetic susceptibility modify or stratify the association between
    intestinal antecedents and AD/dementia risk?

Core analyses:
    - AD genetic main effects
    - GI x genetic susceptibility interaction screening
    - stratified GI-outcome models by APOE or AD PRS strata
    - joint GI-genetic risk gradients
    - optional Cox confirmation for prioritized interaction signals
    - live model diagnostics and timing summaries

Primary fast model:
    The default fast screening model is a person-time-adjusted Poisson GLM:

        event ~ GI + genetic_modifier + GI:genetic_modifier + covariates
        offset(log(follow-up time))

    Exponentiated coefficients are reported as incidence rate ratios.

Recommended command:

    python scripts/genetic_interaction_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --outdir results/genetic_interaction_analysis

Run only a subset of outcomes or modifiers:

    python scripts/genetic_interaction_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --outdir results/genetic_interaction_analysis \
        --outcomes strict_ad \
        --modifiers APOE4_carrier,AD_PRS_enhanced

Run optional Cox confirmation:

    python scripts/genetic_interaction_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --outdir results/genetic_interaction_analysis \
        --cox-confirm \
        --max-cox-confirm 30

Resume a long run:

    python scripts/genetic_interaction_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --outdir results/genetic_interaction_analysis \
        --resume

Representative outputs:
    tables/genetic_main_effects.csv
    tables/gi_genetic_interactions.csv
    tables/stratified_gi_outcome_by_genetics.csv
    tables/joint_gi_genetic_risk_gradient.csv
    tables/cox_confirmatory_top_interactions.csv
    tables/fast_model_diagnostics_final.csv
    tables/model_status_summary.csv
    tables/cohort_and_result_summary.csv
    qc/preflight_summary.csv
    run_metadata.json
    analysis_report.txt

Interpretation:
    This module treats genetics as susceptibility background and effect
    modifier, not as a mediator. The multiplicative interaction term evaluates
    departure from multiplicativity on the incidence-rate scale. Joint gradient
    tables are used to summarize absolute incidence-rate amplification across
    combined GI and genetic strata.


Script 6: risk_warning_window_patch.py
--------------------------------------
Purpose:
    Patch an existing dynamic risk timing script to strengthen individualized
    AD warning-window analysis.

This script is a code patch utility rather than a primary analysis script.

It modifies a risk-timing script by:
    - fixing an off-by-one issue in risk-acceleration window derivation
    - adding GI burden and GI disease indicators to individual risk-curve output
    - creating any_gi and gi_burden_group variables
    - replacing a limited subgroup summary with GI-aware warning-window summaries
    - saving both warning_window_subgroup_summary.csv and the backward-compatible
      subgroup_high_risk_summary.csv alias

Recommended command:

    python scripts/risk_warning_window_patch.py \
        --input scripts/risk_timing_stratification.py \
        --output scripts/risk_timing_stratification_warning_windows.py

Dry-run check:

    python scripts/risk_warning_window_patch.py \
        --input scripts/risk_timing_stratification.py \
        --dry-run

Dry-run with diff preview:

    python scripts/risk_warning_window_patch.py \
        --input scripts/risk_timing_stratification.py \
        --dry-run \
        --show-diff

Outputs:
    A patched Python script. The patched script, when run, should produce
    GI-aware individual risk-curve and warning-window subgroup tables.


Recommended analysis order
--------------------------
The scripts are intended to be run in the following order:

1. Main longitudinal association analysis

    python scripts/gi_ad_analysis.py \
        --input data/master_preprocessed_with_genetics.csv \
        --outdir results/gi_ad_analysis

2. Metabolic bridge analysis

    python scripts/metabolic_bridge_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --metabolic data/metabolic_prepared.csv \
        --outdir results/metabolic_bridge_analysis

3. Metabolic-neural bridge analysis

    python scripts/metabolic_neural_bridge_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --brain data/brain_idp_prepared.csv \
        --metabolic data/metabolic_prepared.csv \
        --metabolic-bridge-dir results/metabolic_bridge_analysis \
        --outdir results/metabolic_neural_bridge_analysis

4. Metabolic-cognition bridge analysis

    python scripts/metabolic_cognition_bridge_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --cognition data/cognition_prepared.csv \
        --metabolic data/metabolic_prepared.csv \
        --outdir results/metabolic_cognition_bridge_analysis

5. Genetic susceptibility and interaction analysis

    python scripts/genetic_interaction_analysis.py \
        --master data/master_preprocessed_with_genetics.csv \
        --outdir results/genetic_interaction_analysis

6. Optional risk warning-window patch

    python scripts/risk_warning_window_patch.py \
        --input scripts/risk_timing_stratification.py \
        --output scripts/risk_timing_stratification_warning_windows.py

7. Run the patched risk-timing script if dynamic risk prediction is included in
   the public repository.


Statistical details
-------------------
Longitudinal association models:
    Cox proportional hazards models are used for disease-specific intestinal
    antecedent associations with incident AD and dementia outcomes. A small
    penalization term is used for numerical stability. Hazard ratios and 95%
    confidence intervals are reported.

Metabolic bridge models:
    Stage 1 uses robust OLS:
        M = alpha_0 + alpha * GI + covariates + error

    Stage 2 for AD uses Cox models:
        h(t) = h0(t) exp(theta * M + delta * GI + covariates)

    Candidate bridge markers are selected when evidence supports both
    GI-to-marker and marker-to-AD relationships under specified thresholds.

Metabolic-neural and metabolic-cognition bridge models:
    Stage 2 uses robust OLS:
        P = beta_0 + beta * M + delta * GI + covariates + error

    Product-of-coefficients bridge effect:
        bridge effect = alpha * beta

    Delta-method standard error:
        SE(alpha * beta) =
            sqrt(beta^2 * SE(alpha)^2 + alpha^2 * SE(beta)^2)

    Sobel Z statistic:
        Z = alpha * beta / SE(alpha * beta)

Genetic interaction models:
    The fast primary model is person-time-adjusted Poisson GLM:
        log E(Y) = log(T) + beta_G * GI + beta_Z * Z
                   + beta_GZ * GI * Z + covariates

    beta_GZ estimates multiplicative interaction on the incidence-rate scale.

Dynamic risk timing:
    Individual cumulative risk curves are estimated from 1 to 20 years.
    Warning onset is the first year at which cumulative risk enters the
    training-derived warning zone. The risk-acceleration window is the 3-year
    interval with the largest increase in annual predicted risk increments.


Quality-control outputs
-----------------------
Each script writes QC and metadata outputs. Users should inspect these files
before interpreting main results.

Common QC outputs include:
    - cohort summary
    - variable diagnostics
    - covariate block definitions
    - exposure and outcome availability
    - model diagnostics
    - run metadata
    - human-readable analysis report

Recommended first files to inspect:
    results/*/analysis_report.txt
    results/*/run_metadata.json
    results/*/qc/*.csv


Reproducibility notes
---------------------
1. Random seeds are set in scripts where random procedures are used.
2. Long-running scripts write live diagnostics or intermediate tables where
   appropriate.
3. Several scripts support resume or partial rerun options.
4. Model results depend on the exact preprocessing, feature availability,
   missingness handling, and cohort access version.
5. The code assumes that predictors are available at or before the relevant
   prediction or assessment date. Users should avoid adding future information
   to the master table.


Data availability
-----------------
The repository does not redistribute individual-level cohort data.

Data access should be requested from the original cohort platforms:

    UK Biobank: https://ukbiobank.dnanexus.com
    CHARLS:    https://charls.pku.edu.cn
    HRS:       https://hrsdata.isr.umich.edu/
    LASI:      https://lasi-india.org/
    SHARE:     https://releases.sharedataportal.eu/

Users must comply with each cohort's data-use agreement, privacy policy, and
publication requirements.


Code availability
-----------------
The code used for this study is intended to be shared at:

    https://github.com/Association-ID-AD/main

Before public release, users should confirm that:
    - all personal/local paths have been removed
    - no restricted data or participant-level output is included
    - no cohort-protected fields are redistributed
    - environment and dependency information is documented
    - scripts have been tested using example or synthetic inputs where possible


Citation
--------
If using this code, please cite the associated manuscript:

    Liang Z, Tan T, Yao Y, Lam CK, Dai L, Lian Z, Shen D.
    Multimodal trajectories of intestinal antecedents and Alzheimer's disease
    risk across multinational longitudinal cohorts.

License
-------
Please add the intended software license before public release. If no license is
provided, downstream users do not automatically receive permission to reuse,
modify, or redistribute the code.

Recommended options include:
    - MIT License for permissive academic/software reuse
    - Apache License 2.0 for permissive reuse with explicit patent language
    - GPLv3 if derivative code must remain open-source

Contact
-------
For questions about the study and code, please contact the lead contact listed
in the manuscript resource availability section.

