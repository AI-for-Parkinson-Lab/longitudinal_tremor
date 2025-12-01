# Longitudinal tremor
This repository contains code related to the scientific publication: <b>'Daily-life tremor measures are sensitive to progression in early Parkinson disease'</b>.

## Overview
The repository contains code to transform weekly tremor measures into results presented in the scientific publication. The Matlab scripts are enumerated to ensure proper serial execution
and improve readability:
* `step1_participant_selection.m`: Apply inclusion/exclusion criteria and create subgroups of participants;
* `step2_descriptives.m`: Extract demographic and clinical characteristics of the selected participants;
* `step3_l1trend_filtering.m`: Apply l1 trend filtering to the biweekly tremor measures, perform interpolation to extracted trends, and plot individual examples of trends;
* `step4_IPCW.m`: Extract weights for inverse probability-of-censoring weighting;
* `step5_group_level_figures.m`: Create plots of (changes in) sensor-derived tremor measures over time on a group level;
* `step6a_sensitivity_to_change.m`: Assess the sensitivity to change of sensor-derived tremor measures and clinical scores;
* `step6b_SRM_time_plots.m`: Show the standardized response mean over time for the sensor-derived tremor measures (Supplementary analysis);
* `step6c_sensitivity_supplement.m`: Assess the sensitivity to change in matched medicated and unweighted unmedicated groups (Supplementary analysis);
* `step7_correlation_patient_reported_tremor.m`: Assess the correlation between changes in patient-reported tremor score and sensor-derived tremor measures/clinical tremor scores;
* `step8_regression_analysis.m`: Assess the effect of medication, disease duration and watch side on changes in sensor-derived tremor measures;
* `step9_responsiveness_dopa_treatment_initiation.m`: Assess the responsiveness of sensor-derived tremor measures to dopaminergic treatment initiation;
  
Note that, while executing `step4_IPCW.m`, we performed a survival analysis in R (see `Survival_analysis_tremor.R`) to estimate the probability to remain unmedicated (uncensored) for each participant.
Based on these probabilites, inverse probability-of-censoring weights were computed in `step4_IPCW.m`.

## Contact
For questions, please reach out to nienke.timmermans@radboudumc.nl.
