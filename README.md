# ✂️ SPLIT: Safety Prioritization for Long-COVID Drug Repurposing via a Causal Integrated Targeting Framework

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

> **Integrating Causal Driver Genes with Clinical Knowledge Graphs for Safety-First Drug Repurposing in Long COVID**

A computational framework that combines causal gene discovery with clinical knowledge graph learning to identify drugs that should be **avoided** (deprioritized) for Long COVID clinical trials based on unfavorable predicted safety, toxicity burden, and efficacy profiles.

---

## 🎯 Overview

![SPLIT Framework Workflow](Workflow.png)

**SPLIT** addresses a critical gap in drug repurposing: rather than simply ranking candidates by predicted efficacy, we prioritize **safety-first deprioritization**—systematically identifying drugs with unacceptable risk-benefit profiles before they enter costly clinical trials.

### Key Features

- **Multi-method causal gene discovery**: Integrates Mendelian Randomization (MR), Control Theory (CT), and Differential Causal Effects (DCE) to identify 1,725 Long COVID driver genes
- **Knowledge graph-based prediction**: Uses PlaNet to generate cohort-anchored predictions of safety (S), adverse events (AE), and efficacy (E)
- **Counterfactual trial design**: Evaluates drugs within fixed cohort contexts by substituting interventions while preserving protocol structure
- **Multi-tier deprioritization**: Generates complementary avoid lists (composite, threshold-based, counterfactual-delta) with configurable stringency
- **Pharmacology enrichment**: Contextualizes predictions with wet-lab evidence from PubChem, ChEMBL, Open Targets, and DailyMed
- **Cross-validation with active trials**: Compares avoid list drugs against ongoing Long COVID clinical and preclinical programs

---

## 📊 Key Results

### Causal Gene Landscape

| Evidence Source | Genes |
|-----------------|-------|
| Control Theory (CT) | 1,641 |
| Mendelian Randomization (MR) | 50 |
| Differential Causal Effects (DCE) | 45 |
| **Multi-method convergent** | **11** |
| **Total unique genes** | **1,725** |

### Safety-Coverage Trade-off

A critical finding: **maximizing causal gene coverage does not improve efficacy but substantially compromises safety**.

| Drug | Causal Genes | Safety (S) | Efficacy (E) |
|------|--------------|------------|--------------|
| Niagen (Original) | — | 0.549 | 0.519 |
| Leronlimab | 1 | 0.441 | 0.521 |
| Vortioxetine | 11 | 0.332 | 0.521 |
| Ritonavir | 70 | 0.215 | 0.527 |

### Deprioritization Outputs

From 1,625 screened drugs, we generated multi-tier avoid lists for two clinical cohorts:

#### NCT04809974 (Cognitive Cohort)

| List Type | Threshold | Drugs (n) | Description |
|-----------|-----------|-----------|-------------|
| Primary q05 | Worst 5% | 82 | Core deprioritized set |
| Primary q10 | Worst 10% | 163 | Default threshold |
| Primary q20 | Worst 20% | 325 | Extended set |
| Threshold-based | Any extreme | 855 | ≥1 extreme flag |
| Multi-extreme | ≥2 flags | 121 | Highest concern |
| Delta (vs Niagen) | Any worse | 456 | Worse than original |

**Total unique drugs**: 883

#### NCT04880161 (Respiratory Cohort)

| List Type | Threshold | Drugs (n) | Description |
|-----------|-----------|-----------|-------------|
| Primary q05 | Worst 5% | 82 | Core deprioritized set |
| Primary q10 | Worst 10% | 163 | Default threshold |
| Primary q20 | Worst 20% | 325 | Extended set |
| Threshold-based | Any extreme | 803 | ≥1 extreme flag |
| Delta (vs Ampion) | Any worse | 412 | Worse than original |

**Total unique drugs**: 824

**113 drugs** appear in all three avoid list types for NCT04809974, representing the highest-confidence deprioritization targets.

---

## 🔬 Data Files

### Core Data

| File | Description | Rows |
|------|-------------|------|
| `1_Causal_Genes.csv` | Long COVID driver genes with evidence source labels | 1,725 |
| `2_Drug_Targets.csv` | Drug-to-causal-gene mappings from 5 databases | 19,172 |
| `3_Ground_Truth.csv` | PlaNet predictions for template + counterfactuals | 4 scenarios |

### Avoid Lists

#### NCT04809974 (Cognitive Cohort)

| File | Description | Drugs |
|------|-------------|-------|
| `NCT04809974_avoid_primary_composite_q10.csv` | Worst 10% by composite score | 163 |
| `NCT04809974_avoid_secondary_guard_q20.csv` | Any extreme metric flag | 855 |
| `NCT04809974_avoid_delta_q10.csv` | Worse than Niagen on any dimension | 456 |
| `NCT04809974_avoid_sensitivity_primary_q05.csv` | Worst 5% (most stringent) | 82 |
| `NCT04809974_avoid_sensitivity_primary_q10.csv` | Worst 10% primary | 163 |
| `NCT04809974_avoid_sensitivity_primary_q20.csv` | Worst 20% (extended) | 325 |
| `NCT04809974_avoid_sensitivity_secondary_q05.csv` | Worst 5% secondary | 232 |
| `NCT04809974_avoid_sensitivity_secondary_q10.csv` | Worst 10% secondary | 456 |
| `NCT04809974_avoid_sensitivity_secondary_q20.csv` | Worst 20% secondary | 850 |

#### NCT04880161 (Respiratory Cohort)

| File | Description | Drugs |
|------|-------------|-------|
| `NCT04880161_avoid_primary_composite_q10.csv` | Worst 10% by composite score | 163 |
| `NCT04880161_avoid_secondary_guard_q20.csv` | Any extreme metric flag | 803 |
| `NCT04880161_avoid_delta_q10.csv` | Worse than Ampion on any dimension | 412 |
| `NCT04880161_avoid_sensitivity_primary_q05.csv` | Worst 5% (most stringent) | 82 |
| `NCT04880161_avoid_sensitivity_primary_q10.csv` | Worst 10% primary | 163 |
| `NCT04880161_avoid_sensitivity_primary_q20.csv` | Worst 20% (extended) | 325 |
| `NCT04880161_avoid_sensitivity_secondary_q05.csv` | Worst 5% secondary | 232 |
| `NCT04880161_avoid_sensitivity_secondary_q10.csv` | Worst 10% secondary | 456 |
| `NCT04880161_avoid_sensitivity_secondary_q20.csv` | Worst 20% secondary | 803 |

### Enrichment Tables

| File | Description | Candidates |
|------|-------------|------------|
| `NCT04809974_enriched_pharmacology.xlsx` | Cognitive cohort enrichment | 42 |
| `NCT04880161_enriched_pharmacology.xlsx` | Respiratory cohort enrichment | 12 |

---

## 🧪 Cross-Reference with Active Long COVID Trials

A key validation of our framework: several drugs appearing on avoid lists are simultaneously under investigation in active Long COVID clinical trials, supported by preclinical evidence. This underscores that avoid lists function as **hypothesis-generating filters** for specific cohort-drug combinations, not universal exclusion criteria.

### Drugs on Avoid Lists with Active Long COVID Clinical Programs

#### NCT04809974 (Cognitive Cohort) — 4 Drugs

| Drug | Avoid Lists | Long COVID Trial | Phase | Preclinical Evidence |
|------|-------------|------------------|-------|---------------------|
| **Baricitinib** | All 9 files | REVERSE-LC (NCT05858515) | Phase 3 | In vitro: inhibits SARS-CoV-2 entry via NAK blockade in 3D liver organoids; Mouse: reverses neurocognitive deficits in SCID model |
| **Naltrexone (LDN)** | 6 files | RECOVER-TLC pivotal trial | Phase 2/3 | In vitro: restores TRPM3 ion channel function in NK cells from Long COVID patients (Sasso 2024) |
| **Ibudilast** | 4 files | Phase 2-3 RCT | Phase 2-3 | In vitro: PDE4/PDE10 inhibition reduces pro-inflammatory cytokines; BBB penetrant |
| **Ritonavir** | 2 files | RECOVER-VITAL (NCT05576662) | Phase 3 | In vitro: nirmatrelvir IC₅₀ ~19nM against 3CL protease; effective in Vero E6, Calu-3 cells |

#### NCT04880161 (Respiratory Cohort) — 9 Drugs

| Drug | Avoid Lists | Long COVID Trial | Phase | Preclinical Evidence |
|------|-------------|------------------|-------|---------------------|
| **Baricitinib** | 4 files | REVERSE-LC (NCT05858515) | Phase 3 | JAK1/2 inhibition blocks IL-6, IFN signaling |
| **Ibudilast** | 4 files | Phase 2-3 RCT | Phase 2-3 | Reduces TNF-α, IL-1β, IL-6 in vitro |
| **Rosuvastatin** | 7 files | Observational studies | Observational | Endothelial anti-inflammatory effects |
| **Metformin** | 2 files | COVID-OUT (Phase 3 completed) | Phase 3 | Reduces SARS-CoV-2 replication by 99% in Caco2/Calu3 cells; AMPK→mTOR inhibition |
| **Ivabradine** | 3 files | RECOVER-AUTONOMIC | Phase 2/3 | Selective If channel inhibition in sinoatrial node cells |
| **Modafinil** | 2 files | RECOVER-SLEEP | Phase 2/3 | Wake-promoting via dopaminergic/histaminergic pathways |
| **Guanfacine** | 3 files | Open-label studies | Phase 2 | α2A-agonism strengthens prefrontal cortex connectivity in primate models |
| **Atorvastatin** | 3 files | Observational studies | Observational | Reduces endothelial inflammation markers |
| **Ritonavir** | 2 files | RECOVER-VITAL | Phase 3 | Hamster/mouse models confirm viral load reduction |

### Key Insight

The presence of drugs like **baricitinib** (flagged across all avoid lists for NCT04809974) as the leading immunomodulator candidate in the $12M REVERSE-LC trial illustrates that:

1. **Context matters**: A drug flagged for avoidance in combination with nicotinamide riboside may be beneficial as a standalone intervention
2. **Avoid lists are cohort-specific**: Safety-efficacy profiles depend on the specific trial population and co-interventions
3. **Preclinical validation exists**: These drugs have robust wet-lab evidence supporting their investigation in Long COVID

---

## 🏥 Clinical Cohorts

### Primary Cohort: NCT04809974 (Cognitive)

| Characteristic | Value |
|----------------|-------|
| Drug | Niagen (Nicotinamide Riboside) |
| Sponsor | University of Pennsylvania |
| Target Domain | Cognitive (brain fog) |
| Sample Size | 72 |
| Age (mean ± SD) | 45.9 ± 13.4 |
| Female | 72.4% |
| Race (White) | 91% |
| **Unique avoid list drugs** | **883** |

### Sensitivity Cohort: NCT04880161 (Respiratory)

| Characteristic | Value |
|----------------|-------|
| Drug | Ampion (Inhaled Biologic) |
| Sponsor | Ampio Pharmaceuticals |
| Target Domain | Respiratory |
| Sample Size | 32 |
| Age (mean ± SD) | 52.1 ± 13.0 |
| Female | 56.3% |
| Race (White/Black) | 75% / 25% |
| **Unique avoid list drugs** | **824** |

---

## 🎯 Top Deprioritized Drugs

### NCT04809974 (Cognitive Cohort) - Top 10 to Avoid

| Rank | Drug | Composite | Safety | Efficacy | AE Score | Rationale |
|------|------|-----------|--------|----------|----------|-----------|
| 1 | Upadacitinib | 0.159 | 0.241 | 0.515 | 2.27 | JAK inhibitor; immunosuppressive |
| 2 | Secnidazole | 0.164 | 0.124 | 0.516 | 1.14 | Lowest safety score |
| 3 | Fluticasone propionate | 0.185 | 0.127 | 0.516 | 0.82 | Very low safety |
| 4 | Pivampicillin | 0.188 | 0.320 | 0.514 | 2.18 | Low efficacy; high AE |
| 5 | Empagliflozin | 0.191 | 0.183 | 0.517 | 1.66 | SGLT2 inhibitor; metabolic risks |
| 6 | Fidaxomicin | 0.191 | 0.200 | 0.517 | 2.03 | High AE burden |
| 7 | Ethionamide | 0.201 | 0.196 | 0.517 | 1.87 | Hepatotoxicity risk |
| 8 | Bithionol | 0.204 | 0.299 | 0.515 | 2.55 | Highest AE burden |
| 9 | Thiethylperazine | 0.205 | 0.219 | 0.517 | 1.97 | Extrapyramidal risk |
| 10 | Tavaborole | 0.210 | 0.168 | 0.519 | 2.25 | Low safety; high AE |

### NCT04880161 (Respiratory Cohort) - Top 10 to Avoid

| Rank | Drug | Composite | Safety | Efficacy | AE Score | Rationale |
|------|------|-----------|--------|----------|----------|-----------|
| 1 | Pazopanib | 0.142 | 0.183 | 0.512 | 2.89 | TKI; highest AE burden |
| 2 | Azacitidine | 0.156 | 0.201 | 0.513 | 2.67 | Hypomethylating agent |
| 3 | Carboplatin | 0.161 | 0.178 | 0.514 | 2.45 | Chemotherapy toxicity |
| 4 | Rosuvastatin | 0.168 | 0.110 | 0.516 | 1.12 | Very low safety score |
| 5 | Paclitaxel | 0.172 | 0.192 | 0.514 | 2.38 | Chemotherapy; neuropathy |
| 6 | Decitabine | 0.175 | 0.205 | 0.513 | 2.51 | Myelosuppression risk |
| 7 | Temozolomide | 0.179 | 0.211 | 0.514 | 2.33 | CNS penetrant; bone marrow |
| 8 | Upadacitinib | 0.182 | 0.238 | 0.515 | 2.19 | JAK inhibitor |
| 9 | Venetoclax | 0.185 | 0.245 | 0.515 | 2.28 | BCL-2 inhibitor; TLS risk |
| 10 | Ibrutinib | 0.188 | 0.251 | 0.515 | 2.21 | BTK inhibitor; bleeding |

### Drug Classes Enriched in Avoid Lists

| Class | NCT04809974 | NCT04880161 | Common Concerns |
|-------|-------------|-------------|-----------------|
| **Antimicrobials** | 15 (9.2%) | 12 (7.4%) | Nitroimidazoles, antitubercular agents |
| **CNS/Psychiatric** | 11 (6.7%) | 8 (4.9%) | Phenothiazines, anticonvulsants |
| **JAK inhibitors** | 7 (4.3%) | 5 (3.1%) | Immunosuppression, infections |
| **Antidiabetics** | 6 (3.7%) | 5 (3.1%) | SGLT2/DPP-4 inhibitors |
| **Targeted oncology** | 4 (2.5%) | 18 (11.0%) | TKIs, PARP inhibitors, immunotherapy |
| **Corticosteroids** | 3 (1.8%) | 4 (2.5%) | Fluticasone formulations |

---

## 📖 Methods Summary

### Stage A: Causal Gene Discovery

1. **Mendelian Randomization (MR)**: Identifies genes with genetically-anchored causal effects on Long COVID risk using GWAS + GTEx eQTL data
2. **Control Theory (CT)**: Identifies network-indispensable nodes whose removal disrupts system controllability
3. **Differential Causal Effects (DCE)**: Quantifies pathway-level gene-gene interaction differences between Long COVID cases and controls

### Stage B: Drug-Gene Mapping

- Maps 1,725 causal genes to pharmacological modulators using DrugBank, ChEMBL, TTD, DrugCentral, and GtoPdb
- Screening panel: 1,776 drugs with ≥11 causal gene targets

### Stage C: PlaNet Prediction

- Generates cohort-anchored predictions of Safety (S), Adverse Events (AE), and Efficacy (E)
- Uses counterfactual trial design to evaluate drug substitutions within fixed cohort contexts

### Stage D: Safety-First Deprioritization

- **Primary composite**: Ranks by weighted combination of S, E, and inverse AE
- **Threshold-based**: Flags drugs with extreme values on individual metrics
- **Counterfactual-delta**: Identifies drugs worse than the original trial intervention

### Stage E: Cross-Reference Validation

- Compares avoid list drugs against active Long COVID clinical trials (RECOVER, REVERSE-LC, etc.)
- Annotates with preclinical/wet lab evidence (in vitro, animal models)
- Confirms hypothesis-generating (not exclusionary) nature of avoid lists

---

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- [PlaNet](https://github.com/mims-harvard/PlaNet) for clinical knowledge graph learning
- [GTEx Consortium](https://gtexportal.org/) for eQTL data
- [ClinicalTrials.gov](https://clinicaltrials.gov/) for trial data
- [DrugBank](https://go.drugbank.com/), [ChEMBL](https://www.ebi.ac.uk/chembl/), [Open Targets](https://www.opentargets.org/) for drug-target data
- [RECOVER Initiative](https://recovercovid.org/) for Long COVID clinical trial context

---

## 📧 Contact

For questions or collaborations, please open an issue or contact:

- **Lead Author**: Sindy Pinero
- **Project Link**: [https://github.com/yourusername/SPLIT](https://github.com/yourusername/SPLIT)

---

<p align="center">
  <b>SPLIT: Because knowing what to avoid is as important as knowing what to pursue.</b>
</p>
