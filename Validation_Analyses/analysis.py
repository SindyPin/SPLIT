"""SPLIT additional analyses for the ESWA revision.
Runs entirely on existing SPLIT outputs (no new PlaNet predictions).
Outputs: results/results.json, results/*.csv, figures/*.png
"""
import re, json, itertools, warnings
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, average_precision_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
rng = np.random.default_rng(2026)
import os
os.makedirs('results', exist_ok=True); os.makedirs('figures', exist_ok=True)
R = {}

norm = lambda s: re.sub(r'[^A-Z0-9]', '', str(s).upper())
W = dict(S=0.4, E=0.4, AE=0.2)

# ---------------------------------------------------------------- data
dt = pd.read_csv('data/SM/2_Drug_Targets.csv'); dt['key'] = dt.Drug.map(norm)
panel = [norm(l.strip()) for l in open('data/6_CG_Drugs_PlaNet/ALL_Scripts/drugs_to_test.txt') if l.strip()]
pt = dt[dt.key.isin(panel) & (dt.Gene_Count >= 11)].drop_duplicates('key')
gc_map = dict(zip(pt.key, pt.Gene_Count))
genes_map = dict(zip(pt.key, pt.Genes))
panel_keys = list(gc_map)

def gc_lookup(k):
    if k in gc_map: return gc_map[k]
    c = [x for x in panel_keys if x.startswith(k)]
    if len(c) == 1: return gc_map[c[0]]
    if len(c) > 1:  # truncated name, ambiguous: take the longest common
        return gc_map[sorted(c, key=len)[0]]
    return np.nan

paths = {
 'cog': 'data/8_FINAL_RESULTS/NCT04809974/Ranking/NCT04809974_ranked_drugs_by_composite_SAE.tsv',
 'resp': 'data/6_CG_Drugs_PlaNet/4_Score_Ranking/NCT04880161/NCT04880161_ranked_drugs_by_composite_SAE.tsv'}
avoid_paths = {
 'cog': 'data/8_FINAL_RESULTS/NCT04809974/Avoid/NCT04809974_avoid_primary_composite_q10.csv',
 'resp': 'data/6_CG_Drugs_PlaNet/4_Score_Ranking/NCT04880161/NCT04880161_avoid_primary_composite_q10.csv'}
base_paths = {
 'cog': None,
 'resp': 'data/6_CG_Drugs_PlaNet/4_Score_Ranking/NCT04880161/Master/NCT04880161_real_baseline_summary.csv'}

D = {}
for c, p in paths.items():
    d = pd.read_csv(p, sep='\t')
    d['key'] = d.drug.map(norm)
    d['gc'] = d.key.map(gc_lookup)
    d['log_gc'] = np.log10(d.gc)
    av = pd.read_csv(avoid_paths[c]); avk = set(av.drug.map(norm))
    d['avoid_q10'] = d.key.isin(avk)
    d = d.set_index('key'); D[c] = d[~d.index.duplicated()]
common = D['cog'].index.intersection(D['resp'].index)
R['panel'] = dict(mapped_drugs=19172, passed_k11=len(panel), scored_per_cohort={c: int(len(D[c])) for c in D},
                  common_scored=int(len(common)), gc_matched={c: int(D[c].gc.notna().sum()) for c in D},
                  avoid_q10_size={c: int(D[c].avoid_q10.sum()) for c in D})

# boxed warning labels
bbw = pd.read_csv('bbw_union.csv'); bbw['key'] = bbw.drug.map(norm)
bbw['bbw'] = (bbw.black_box_warnings.str.startswith('Has')).astype(int)
bbw_map = dict(zip(bbw.key, bbw.bbw))
for c in D: D[c]['bbw'] = D[c].index.map(bbw_map)

def boot_ci(fn, n, B=1000):
    vals = []
    for _ in range(B):
        idx = rng.integers(0, n, n)
        try: vals.append(fn(idx))
        except Exception: pass
    return [float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))]

# ---------------------------------------------------------------- 1 S vs gene count, full panel
R['s_vs_gc'] = {}
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, c in zip(axes, D):
    d = D[c].dropna(subset=['gc'])
    rs, ps = stats.spearmanr(d.gc, d.safety_score)
    rp, pp = stats.pearsonr(d.log_gc, d.safety_score)
    ci = boot_ci(lambda i: stats.spearmanr(d.gc.values[i], d.safety_score.values[i])[0], len(d))
    out = dict(n=int(len(d)), spearman_S=float(rs), spearman_S_ci=ci, p_S=float(ps),
               pearson_logS=float(rp),
               spearman_AE=float(stats.spearmanr(d.gc, d.ae_score)[0]),
               spearman_E=float(stats.spearmanr(d.gc, d.efficacy_score)[0]),
               spearman_C=float(stats.spearmanr(d.gc, d.composite_score)[0]))
    # within strata
    bins = pd.cut(d.gc, [10, 15, 20, 30, 50, 100, 300])
    out['S_by_gc_bin'] = d.groupby(bins).safety_score.agg(['count', 'median', 'std']).reset_index().astype(str).to_dict('records')
    R['s_vs_gc'][c] = out
    ax.scatter(d.gc, d.safety_score, s=6, alpha=.35, color='#1f4e79')
    ax.set_xscale('log'); ax.set_xlabel('Causal gene targets per drug (log scale)'); ax.set_ylabel('Predicted S')
    ax.set_title(f"{'Cognitive' if c=='cog' else 'Respiratory'} cohort: Spearman rho = {rs:.2f} (n = {len(d)})", fontsize=10)
# four drug reference
four = {'Leronlimab': 1, 'Vortioxetine': 11, 'Ritonavir': 70}
d = D['cog']; xs, ys = [], []
for k, g in four.items():
    kk = norm(k)
    if kk in d.index: xs.append(g); ys.append(d.loc[kk, 'safety_score'])
R['s_vs_gc']['four_drug_pearson_cog'] = float(stats.pearsonr(xs, ys)[0]) if len(xs) >= 3 else None
R['s_vs_gc']['four_drug_points'] = dict(zip([k for k in four if norm(k) in d.index], [float(y) for y in ys]))
plt.tight_layout(); plt.savefig('figures/fig_S_vs_targets_fullpanel.png', dpi=200); plt.close()

# ---------------------------------------------------------------- 2 baselines for the avoid list
def auc_safe(y, s):
    y = np.asarray(y); s = np.asarray(s); m = ~np.isnan(s) & ~np.isnan(y)
    return float(roc_auc_score(y[m], s[m])) if len(set(y[m])) == 2 else np.nan
R['avoid_baselines'] = {}
for c in D:
    d = D[c]; y = d.avoid_q10.astype(int).values
    out = {}
    for name, s in [('neg_gene_count', -d.gc), ('neg_S', -d.safety_score), ('neg_E', -d.efficacy_score), ('AE', d.ae_score)]:
        out['auroc_' + name] = auc_safe(y, s.values)
    # Jaccard of bottom 10% by gene count (most targets) with avoid list
    dd = d.dropna(subset=['gc']); n10 = int(round(0.1 * len(dd)))
    top_gc = set(dd.sort_values('gc', ascending=False).index[:n10]); av = set(dd.index[dd.avoid_q10])
    out['jaccard_topGC_vs_avoid'] = len(top_gc & av) / len(top_gc | av)
    out['overlap_topGC_vs_avoid'] = int(len(top_gc & av)); out['n10'] = n10
    out['expected_overlap_random'] = n10 * len(av) / len(dd)
    R['avoid_baselines'][c] = out

# ---------------------------------------------------------------- 3 boxed warning benchmark
R['bbw'] = {}
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, c in zip(axes, D):
    d = D[c].dropna(subset=['bbw']); y = d.bbw.astype(int).values
    out = dict(n_labelled=int(len(d)), n_bbw=int(y.sum()), prevalence=float(y.mean()))
    preds = {'S (lower = risk)': -d.safety_score.values, 'AE burden': d.ae_score.values,
             'Composite (lower = risk)': -d.composite_score.values, 'E (lower = risk)': -d.efficacy_score.values,
             'Gene count': d.gc.values}
    for name, s in preds.items():
        m = ~np.isnan(s)
        a = roc_auc_score(y[m], s[m]); ap = average_precision_score(y[m], s[m])
        ci = boot_ci(lambda i: roc_auc_score(y[m][i], s[m][i]), m.sum(), 500)
        out[name] = dict(auroc=float(a), auroc_ci=ci, auprc=float(ap))
    # residual S beyond gene count
    dd = d.dropna(subset=['gc']); slope, inter = np.polyfit(dd.log_gc, dd.safety_score, 1)
    resid = dd.safety_score - (slope * dd.log_gc + inter)
    yy = dd.bbw.astype(int).values
    out['S_residual_given_gc'] = dict(auroc=float(roc_auc_score(yy, -resid.values)),
                                      auroc_ci=boot_ci(lambda i: roc_auc_score(yy[i], -resid.values[i]), len(dd), 500))
    # enrichment of BBW in the q10 avoid list
    tab = pd.crosstab(d.avoid_q10, d.bbw); orr, pf = stats.fisher_exact(tab.values[::-1, ::-1]) if tab.shape == (2, 2) else (np.nan, np.nan)
    out['bbw_rate_in_avoid'] = float(d[d.avoid_q10].bbw.mean()); out['bbw_rate_not_avoid'] = float(d[~d.avoid_q10].bbw.mean())
    out['fisher_OR'] = float(orr); out['fisher_p'] = float(pf)
    R['bbw'][c] = out
    names = ['S (lower = risk)', 'AE burden', 'Composite (lower = risk)', 'Gene count', 'E (lower = risk)']
    vals = [out[n]['auroc'] for n in names]; lo = [out[n]['auroc_ci'][0] for n in names]; hi = [out[n]['auroc_ci'][1] for n in names]
    ax.barh(names, vals, xerr=[np.array(vals) - lo, np.array(hi) - vals], color='#1f4e79', alpha=.85)
    ax.axvline(0.5, ls='--', c='k', lw=.8); ax.set_xlim(0.3, 0.85); ax.set_xlabel('AUROC for FDA boxed warning')
    ax.set_title(f"{'Cognitive' if c=='cog' else 'Respiratory'} (n = {len(d)}, {int(y.sum())} boxed)", fontsize=10)
plt.tight_layout(); plt.savefig('figures/fig_boxed_warning_benchmark.png', dpi=200); plt.close()

# ---------------------------------------------------------------- composite recomputation helpers
def minmax(x): return (x - x.min()) / (x.max() - x.min())
def composite(d, w, noise=None):
    S, E, A = d.safety_score.values.copy(), d.efficacy_score.values.copy(), d.ae_score.values.copy()
    if noise is not None:
        S = S + rng.normal(0, noise * (S.max() - S.min()), len(S))
        E = E + rng.normal(0, noise * (E.max() - E.min()), len(E))
        A = A + rng.normal(0, noise * (A.max() - A.min()), len(A))
    return w['S'] * minmax(S) + w['E'] * minmax(E) + w['AE'] * (1 - minmax(A))
def q10set(d, comp):
    n = int(round(0.1 * len(d))); return set(d.index[np.argsort(comp)[:n]])
def jacc(a, b): return len(a & b) / len(a | b)

# check recomputed composite reproduces the stored one
R['composite_check'] = {c: float(stats.spearmanr(composite(D[c], W), D[c].composite_score)[0]) for c in D}

# ---------------------------------------------------------------- 4 weight sensitivity
R['weights'] = {}
grid = [(a / 20, b / 20, 1 - a / 20 - b / 20) for a in range(21) for b in range(21 - a)]
for c in D:
    d = D[c]; base = q10set(d, composite(d, W)); js = []; cnt = pd.Series(0, index=d.index)
    for ws, we, wa in grid:
        s = q10set(d, composite(d, dict(S=ws, E=we, AE=wa))); js.append((ws, we, wa, jacc(base, s))); cnt[list(s)] += 1
    jsdf = pd.DataFrame(js, columns=['wS', 'wE', 'wAE', 'jaccard']); jsdf.to_csv(f'results/weight_sensitivity_{c}.csv', index=False)
    frac = cnt / len(grid)
    R['weights'][c] = dict(n_weightings=len(grid), jaccard_median=float(jsdf.jaccard.median()), jaccard_min=float(jsdf.jaccard.min()),
                           jaccard_q25=float(jsdf.jaccard.quantile(.25)),
                           jaccard_median_nontrivial=float(jsdf[(jsdf.wS >= .2) & (jsdf.wE >= .2) & (jsdf.wAE >= .1)].jaccard.median()),
                           core_ge95=int((frac[list(base)] >= .95).sum()), core_ge50=int((frac[list(base)] >= .5).sum()),
                           base_size=len(base))
    frac.rename('fraction_weightings_in_q10').to_csv(f'results/weight_membership_{c}.csv')

# ---------------------------------------------------------------- 5 k threshold sensitivity
R['k_threshold'] = {}
for c in D:
    d0 = D[c].dropna(subset=['gc']); base = q10set(d0, composite(d0, W)); rows = []
    for k in [11, 15, 20, 30, 50]:
        d = d0[d0.gc >= k]; s = q10set(d, composite(d, W)); b = base & set(d.index)
        rows.append(dict(k=k, panel=len(d), q10=len(s), jaccard_vs_default=jacc(b, s) if b | s else np.nan,
                         retained_fraction=len(b & s) / len(b) if b else np.nan))
    R['k_threshold'][c] = rows

# ---------------------------------------------------------------- 6 perturbation stability
R['stability'] = {}
for c in D:
    d = D[c]; base = q10set(d, composite(d, W)); cnt = pd.Series(0, index=d.index); B = 1000
    for _ in range(B):
        cnt[list(q10set(d, composite(d, W, noise=0.05)))] += 1
    frac = cnt / B
    R['stability'][c] = dict(B=B, noise='5% of score range', core_ge95=int((frac[list(base)] >= .95).sum()),
                             core_ge80=int((frac[list(base)] >= .8).sum()), fragile_lt50=int((frac[list(base)] < .5).sum()),
                             base_size=len(base), mean_membership_of_flagged=float(frac[list(base)].mean()))
    frac.rename('membership_frequency').to_csv(f'results/stability_membership_{c}.csv')

# ---------------------------------------------------------------- 7 E compression and composite variance decomposition
R['E_diagnostics'] = {}
for c in D:
    d = D[c]; E = d.efficacy_score
    comp = dict(S=W['S'] * minmax(d.safety_score), E=W['E'] * minmax(d.efficacy_score), AE=W['AE'] * (1 - minmax(d.ae_score)))
    C = comp['S'] + comp['E'] + comp['AE']
    # share of composite variance attributable to each term (covariance decomposition)
    share = {k: float(np.cov(v, C)[0, 1] / C.var()) for k, v in comp.items()}
    R['E_diagnostics'][c] = dict(E_min=float(E.min()), E_max=float(E.max()), E_range=float(E.max() - E.min()), E_sd=float(E.std()),
                                 E_cv=float(E.std() / E.mean()), S_sd=float(d.safety_score.std()), AE_sd=float(d.ae_score.std()),
                                 spearman_E_S=float(stats.spearmanr(E, d.safety_score)[0]), spearman_E_AE=float(stats.spearmanr(E, d.ae_score)[0]),
                                 spearman_E_gc=float(stats.spearmanr(E, d.gc, nan_policy='omit')[0]),
                                 composite_variance_share=share,
                                 spearman_C_S=float(stats.spearmanr(C, d.safety_score)[0]), spearman_C_E=float(stats.spearmanr(C, E)[0]),
                                 spearman_C_AE=float(stats.spearmanr(C, d.ae_score)[0]))
# E compression: raw PlaNet efficacy for the two trial contexts
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.hist(D['cog'].efficacy_score, bins=60, alpha=.7, label='Cognitive (NCT04809974)', color='#c0392b')
ax.hist(D['resp'].efficacy_score, bins=60, alpha=.7, label='Respiratory (NCT04880161)', color='#1f4e79')
ax.set_xlabel('Predicted E = P(drug > placebo)'); ax.set_ylabel('Drugs'); ax.legend(fontsize=8); plt.tight_layout()
plt.savefig('figures/fig_E_distributions.png', dpi=200); plt.close()

# ---------------------------------------------------------------- 8 cross cohort agreement
a, b = D['cog'].loc[common], D['resp'].loc[common]
ov = int((a.avoid_q10 & b.avoid_q10).sum()); n = len(common); k1, k2 = int(a.avoid_q10.sum()), int(b.avoid_q10.sum())
exp = k1 * k2 / n
p_hyper = float(stats.hypergeom.sf(ov - 1, n, k1, k2))
R['cross_cohort'] = dict(n=n, shared_q10=ov, expected_random=exp, hypergeom_p_ge=p_hyper,
                         spearman_S=float(stats.spearmanr(a.safety_score, b.safety_score)[0]),
                         spearman_E=float(stats.spearmanr(a.efficacy_score, b.efficacy_score)[0]),
                         spearman_AE=float(stats.spearmanr(a.ae_score, b.ae_score)[0]),
                         spearman_C=float(stats.spearmanr(a.composite_score, b.composite_score)[0]),
                         shared_drugs=sorted(a.index[a.avoid_q10 & b.avoid_q10].tolist()))

# ---------------------------------------------------------------- 9 oncology fraction (name suffix + indications)
onc_re = re.compile(r'neoplasm|carcinoma|leukemia|leukaemia|lymphoma|cancer|melanoma|myeloma|sarcoma|tumor|tumour|glioma', re.I)
ind = {}
for p in ['data/6_CG_Drugs_PlaNet/4_Score_Ranking/NCT04880161/NCT04880161_avoid_master_table_enriched_pharmacology.xlsx',
          'data/6_CG_Drugs_PlaNet/4_Score_Ranking/NCT04880161/NCT04880161_avoid_secondary_guard_q20_enriched_pharmacology.xlsx',
          'data/8_FINAL_RESULTS/NCT04809974/Enrichment/NCT04809974_avoid_secondary_guard_q20_enriched_pharmacology.xlsx']:
    x = pd.read_excel(p)
    for k, v in zip(x.drug.map(norm), x.approved_indications):
        if isinstance(v, str) and k not in ind: ind[k] = v
def is_onc(key):
    v = ind.get(key)
    if isinstance(v, str) and onc_re.search(v): return 1
    if re.search(r'(TINIB|CICLIB|RAFENIB|PLATIN|TAXEL|RUBICIN|MUSTINE|CITABINE|CITIDINE|TECAN|POSIDE|ZOMIB|TRABINE|ZOLOMIDE|BLASTINE|CRISTINE)$', key): return 1
    return 0
R['oncology'] = {}
for c in D:
    d = D[c]; d['onc'] = [is_onc(k) for k in d.index]
    tab = pd.crosstab(d.avoid_q10, d.onc)
    orr, pf = stats.fisher_exact(tab.values[::-1, ::-1])
    R['oncology'][c] = dict(onc_in_panel=float(d.onc.mean()), onc_in_avoid=float(d[d.avoid_q10].onc.mean()),
                            onc_in_rest=float(d[~d.avoid_q10].onc.mean()), fisher_OR=float(orr), fisher_p=float(pf),
                            n_onc_panel=int(d.onc.sum()))

# ---------------------------------------------------------------- save
for c in D:
    D[c].reset_index()[['key', 'drug', 'gc', 'safety_score', 'efficacy_score', 'ae_score', 'composite_score', 'avoid_q10', 'bbw', 'onc']].to_csv(f'results/panel_{c}.csv', index=False)
json.dump(R, open('results/results.json', 'w'), indent=1, default=str)
print(json.dumps(R, indent=1, default=str)[:12000])
