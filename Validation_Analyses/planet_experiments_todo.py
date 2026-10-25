"""Experiments that need new PlaNet predictions (not run in the ESWA revision).

Each function writes the list of (trial_context, drug) jobs it needs, in the same
format consumed by 4_Experiment/predict_all_trials.py, and a second function
evaluates the outputs once the predictions exist. Run order:

    python planet_experiments_todo.py plan   -> writes jobs/*.txt
    <run PlaNet on the job lists with the existing PBS pipeline>
    python planet_experiments_todo.py eval   -> reads predictions/, writes results/

E1  Negative control for the causal layer: rebuild the screening panel from
    (a) 20 random gene sets of 1,725 genes drawn from the 20,000 gene background and
    (b) an association only gene set (differentially expressed genes, no causal filter),
    score the drugs that are new to each panel, and compare q10 lists (Jaccard) with
    the causal panel list. Expected cost: about 1,500 new drugs per random set.
E2  Leakage stress test: rescore the 1,624 panel drugs with (a) three alternative
    placeholder wordings for the stripped outcome fields and (b) the drug field
    shuffled across records. Predictions must be invariant to (a) and must change
    under (b).
E3  Calibration: score completed Long COVID and post viral trials with known
    adverse event rates (ClinicalTrials.gov results section) and report calibration
    curves for S and the AE probabilities.
E4  Third context: a non respiratory, non oncology Long COVID trial (for example a
    fatigue or autonomic trial) and a deliberately mismatched oncology trial, to test
    whether the respiratory oncology pattern is a training distribution effect.
"""
import sys, os, re, json
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(2026)
norm = lambda s: re.sub(r'[^A-Z0-9]', '', str(s).upper())


def plan():
    os.makedirs(os.path.join(HERE, 'jobs'), exist_ok=True)
    dt = pd.read_csv(os.path.join(HERE, 'data', '2_Drug_Targets.csv'))
    # E1a: random gene sets. 2_Drug_Targets.csv only lists causal targets, so the
    # full drug to gene map (LG_Causal_Genes/gene_lists/*_results) is needed here.
    with open(os.path.join(HERE, 'jobs', 'E1_README.txt'), 'w') as f:
        f.write('Rebuild drug panels from LG_Causal_Genes/gene_lists/* using 20 random gene sets '
                'of size 1725 (seed 2026) and the association only set; keep drugs with >= 11 targets; '
                'list drugs not already in drugs_to_test.txt.\n')
    # E2: placeholder variants and shuffled drug field
    panel = [l.strip() for l in open(os.path.join(HERE, 'data', 'drugs_to_test.txt')) if l.strip()]
    variants = {'placeholder_A': 'Outcome data not available.',
                'placeholder_B': 'Results withheld for this analysis.',
                'placeholder_C': 'No outcome information is provided.'}
    json.dump(variants, open(os.path.join(HERE, 'jobs', 'E2_placeholders.json'), 'w'), indent=1)
    shuffled = list(panel); rng.shuffle(shuffled)
    pd.DataFrame({'record_drug': panel, 'scored_as': shuffled}).to_csv(
        os.path.join(HERE, 'jobs', 'E2_shuffled_drug_field.csv'), index=False)
    # E3 and E4: trial ids to parse
    open(os.path.join(HERE, 'jobs', 'E3_E4_trials.txt'), 'w').write(
        '# E3 calibration: completed trials with posted results\nNCT05047952\nNCT04604704\nNCT05430152\n'
        '# E4 third context (fatigue or autonomic) and mismatched oncology context\nNCT06404086\nNCT05481177\n')
    print('job lists written to', os.path.join(HERE, 'jobs'))


def eval():
    from scipy import stats
    pred = os.path.join(HERE, 'predictions')
    if not os.path.isdir(pred):
        sys.exit('predictions/ not found: run PlaNet on the job lists first')
    # E2a invariance: Spearman and max absolute difference of S across placeholder variants
    frames = {k: pd.read_csv(os.path.join(pred, f'{k}_ranked_drugs_by_composite_SAE.tsv'), sep='\t')
              for k in ['placeholder_A', 'placeholder_B', 'placeholder_C'] if
              os.path.exists(os.path.join(pred, f'{k}_ranked_drugs_by_composite_SAE.tsv'))}
    out = {}
    keys = list(frames)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = frames[keys[i]].set_index('drug'), frames[keys[j]].set_index('drug')
            c = a.index.intersection(b.index)
            out[f'{keys[i]}_vs_{keys[j]}'] = dict(
                spearman_S=float(stats.spearmanr(a.loc[c, 'safety_score'], b.loc[c, 'safety_score'])[0]),
                max_abs_dS=float((a.loc[c, 'safety_score'] - b.loc[c, 'safety_score']).abs().max()))
    json.dump(out, open(os.path.join(HERE, 'results', 'E2_invariance.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    {'plan': plan, 'eval': eval}[sys.argv[1]]()
