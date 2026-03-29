# Source Code for "Dynamic risk stratification using time-varying anion gap in sepsis-induced myocardial injury: insights from a joint modeling framework"

This repository contains the R scripts used for data extraction, preprocessing, statistical analysis, and visualization in the study.

## Repository Structure

The code is organized into two main directories corresponding to the two datasets used in the study:

### 1. `mimic_iv/`
Contains the scripts used for the primary analysis based on the MIMIC-IV database.
*   **Data Extraction & Preprocessing**: Scripts starting with `step0_` to `step3_` handle the extraction of baseline and longitudinal data, cohort definition (SIMI), and missing data imputation using Random Forest.
*   **Statistical Modeling**: Scripts starting with `step4_` perform the Bayesian joint modeling (JMbayes2) to evaluate the time-varying association between anion gap/lactate trajectories and 28-day mortality.
*   **Validation & Sensitivity**: Scripts starting with `step5_` and `step6_` conduct model diagnostics, sensitivity analyses (complete-case), and subgroup analyses.
*   **Causal Inference**: Scripts like `step7_causal_inference_ow_dr.r` implement Overlap Weighting (OW) and Doubly Robust (DR) estimation.

### 2. `zigong/`
Contains the scripts used for external validation based on the Zigong Fourth People's Hospital database.
*   The workflow mirrors the primary analysis, adapting the data extraction and formatting for the Zigong dataset while applying the identical Bayesian joint modeling framework to confirm the robustness of the findings.

## Requirements

The analysis was performed using R. Key packages required include:
*   `tidyverse` (Data manipulation and visualization)
*   `survival` (Survival analysis)
*   `JMbayes2` (Bayesian joint modeling)
*   `nlme` (Linear mixed-effects models)
*   `missForest` (Missing data imputation)
*   `WeightIt` (Propensity score weighting)
*   `rms` (Restricted cubic splines)

## Note on Data Privacy

To protect patient privacy and comply with data use agreements, this repository only contains the analysis code. The raw datasets (MIMIC-IV and Zigong) are not included.
*   **MIMIC-IV** is publicly available via PhysioNet (https://doi.org/10.13026/rnrx-jm88) to credentialed users.
*   **Zigong ICU Infection database** is available via PhysioNet (https://doi.org/10.13026/xpt9-z726).
