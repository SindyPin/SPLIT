# SPLIT Supplementary Table S21: validation and sensitivity analyses

All analyses here run on the stored SPLIT outputs (no new PlaNet predictions).

Inputs expected (relative to the SPLIT project folder):
- SM/2_Drug_Targets.csv, SM/1_Causal_Genes.csv
- 6_CG_Drugs_PlaNet/ALL_Scripts/drugs_to_test.txt (the 1,776 drug panel)
- 8_FINAL_RESULTS/NCT04809974/Ranking/*, 8_FINAL_RESULTS/NCT04809974/Avoid/*
- 6_CG_Drugs_PlaNet/4_Score_Ranking/NCT04880161/* (respiratory ranking, avoid lists, enriched tables)
- bbw_union.csv: FDA boxed warning labels pooled from the enrichment xlsx files (1,202 drugs)
- data/lc_registry_drugs.csv: drug interventions in ClinicalTrials.gov Long COVID studies, accessed 2026-10-11

Scripts:
- analysis.py   panel reconciliation, S vs target count, baselines, boxed warning benchmark,
                weight / threshold / perturbation sensitivity, E diagnostics, cross context agreement,
                oncology enrichment  -> results/results.json, results/*.csv, figures/*.png
- analysis2.py  scale guarded composite (tau = 0.05), registry base rate, revised top 10 tables
                -> results/results_part2.json, results/panel_*_guarded.csv
- planet_experiments_todo.py  job lists and evaluators for experiments that need new PlaNet runs
                (negative gene set controls, leakage stress test, calibration, third context)

Key outputs:
- results/panel_cog_guarded.csv, results/panel_resp_guarded.csv: per drug S, E, AE, composite,
  guarded composite, q10 flags, perturbation membership frequency, boxed warning label, oncology flag
- results/weight_membership_*.csv: fraction of the 231 weight triples under which each drug is in q10
- results/stability_membership_*.csv: fraction of 1,000 perturbation replicates in q10

Python 3.11+, pandas, numpy, scipy, scikit-learn, matplotlib, openpyxl.
