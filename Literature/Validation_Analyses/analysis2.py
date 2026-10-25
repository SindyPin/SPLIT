"""Part 2: scale aware (guarded) composite, trial drug base rate, revised tables."""
import re, json, numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
rng = np.random.default_rng(7)
norm = lambda s: re.sub(r'[^A-Z0-9]', '', str(s).upper())
R = json.load(open('results/results.json'))
P = {c: pd.read_csv(f'results/panel_{c}.csv').set_index('key') for c in ['cog', 'resp']}
TAU = 0.05  # minimum range for a dimension to be treated as informative

def minmax(x): return (x - x.min()) / (x.max() - x.min())
def guarded_composite(d, w=dict(S=.4, E=.4, AE=.2), tau=TAU):
    comps = {'S': d.safety_score, 'E': d.efficacy_score, 'AE': d.ae_score}
    rng_ = {k: float(v.max() - v.min()) for k, v in comps.items()}
    active = {k: (rng_[k] >= tau) for k in comps}
    wsum = sum(w[k] for k in comps if active[k])
    C = np.zeros(len(d)); used = {}
    for k, v in comps.items():
        if not active[k]: continue
        nv = minmax(v) if k != 'AE' else 1 - minmax(v)
        C += (w[k] / wsum) * nv.values; used[k] = w[k] / wsum
    return pd.Series(C, index=d.index), used, rng_
def q10(d, comp): n = int(round(0.1 * len(d))); return set(comp.sort_values().index[:n])

out = {}
G = {}
for c, d in P.items():
    C, used, rg = guarded_composite(d)
    d['composite_guarded'] = C; d['avoid_q10_guarded'] = d.index.isin(q10(d, C))
    orig = set(d.index[d.avoid_q10]); new = set(d.index[d.avoid_q10_guarded])
    info = dict(effective_weights=used, ranges=rg, overlap_with_original=len(orig & new), jaccard=len(orig & new) / len(orig | new),
                spearman_vs_original=float(stats.spearmanr(C, d.composite_score)[0]))
    # composition of the guarded avoid list
    a = d[d.avoid_q10_guarded]
    info['pct_lowS'] = float((a.safety_score < .3).mean()); info['pct_lowE'] = float((a.efficacy_score < .5).mean()); info['pct_highAE'] = float((a.ae_score > 2).mean())
    info['median_S_avoid'] = float(a.safety_score.median()); info['median_AE_avoid'] = float(a.ae_score.median())
    info['median_S_panel'] = float(d.safety_score.median()); info['median_AE_panel'] = float(d.ae_score.median())
    # boxed warning check
    dd = d.dropna(subset=['bbw']); y = dd.bbw.astype(int).values
    info['bbw_auroc_guarded'] = float(roc_auc_score(y, -dd.composite_guarded.values))
    info['bbw_rate_in_guarded_avoid'] = float(dd[dd.avoid_q10_guarded].bbw.mean()); info['bbw_rate_rest'] = float(dd[~dd.avoid_q10_guarded].bbw.mean())
    tab = pd.crosstab(dd.avoid_q10_guarded, dd.bbw); orr, pf = stats.fisher_exact(tab.values[::-1, ::-1]); info['bbw_fisher_OR'] = float(orr); info['bbw_fisher_p'] = float(pf)
    # oncology
    info['onc_in_guarded_avoid'] = float(a.onc.mean())
    # stability under perturbation
    base = new; cnt = pd.Series(0, index=d.index)
    for _ in range(1000):
        dn = d.copy()
        for col in ['safety_score', 'efficacy_score', 'ae_score']:
            dn[col] = dn[col] + rng.normal(0, 0.05 * (d[col].max() - d[col].min()), len(d))
        Cn, _, _ = guarded_composite(dn); cnt[list(q10(dn, Cn))] += 1
    frac = cnt / 1000; info['stable_core_ge95'] = int((frac[list(base)] >= .95).sum()); info['fragile_lt50'] = int((frac[list(base)] < .5).sum())
    info['mean_membership'] = float(frac[list(base)].mean())
    d['membership_guarded'] = frac
    a = d[d.avoid_q10_guarded]
    top = a.sort_values('composite_guarded').head(10)[['drug', 'composite_guarded', 'safety_score', 'efficacy_score', 'ae_score', 'membership_guarded']]
    info['top10'] = top.round(3).to_dict('records')
    G[c] = info; P[c] = d
# cross cohort under guard
a, b = P['cog'], P['resp'].loc[P['cog'].index]
sh = a.index[a.avoid_q10_guarded & b.avoid_q10_guarded].tolist()
n = len(a); k1 = int(a.avoid_q10_guarded.sum()); k2 = int(b.avoid_q10_guarded.sum())
G['cross'] = dict(shared=len(sh), expected=k1 * k2 / n, p=float(stats.hypergeom.sf(len(sh) - 1, n, k1, k2)),
                  spearman_C=float(stats.spearmanr(a.composite_guarded, b.composite_guarded)[0]), shared_drugs=sorted(sh))
out['guarded'] = G

# ---------------------------------------------------------------- trial drug base rate
reg = pd.read_csv('data/lc_registry_drugs.csv', comment='#'); reg['key'] = reg.drug.map(norm)
panel_keys = list(P['cog'].index)
def match(k):
    if k in panel_keys: return k
    c = [x for x in panel_keys if x.startswith(k) or k.startswith(x) and len(x) >= 12]
    return c[0] if len(c) == 1 else (None if not c else c[0])
reg['panel_key'] = reg.key.map(match)
inpanel = reg.dropna(subset=['panel_key']).drop_duplicates('panel_key')
T = {}
for lab, sub in [('all_registry', inpanel), ('active_only', inpanel[inpanel.active == 1])]:
    keys = sub.panel_key.tolist(); r = {}
    for c in ['cog', 'resp']:
        d = P[c].loc[keys]
        for col in ['avoid_q10', 'avoid_q10_guarded']:
            k = int(d[col].sum()); N = len(keys); p0 = P[c][col].mean()
            r[f'{c}_{col}'] = dict(flagged=k, n=N, expected=N * p0, binom_p_ge=float(stats.binom.sf(k - 1, N, p0)))
    either = P['cog'].loc[keys].avoid_q10 | P['resp'].loc[keys].avoid_q10
    p_either = (P['cog'].avoid_q10 | P['resp'].loc[P['cog'].index].avoid_q10).mean()
    r['either_cohort_q10'] = dict(flagged=int(either.sum()), n=len(keys), expected=len(keys) * p_either, binom_p_ge=float(stats.binom.sf(int(either.sum()) - 1, len(keys), p_either)))
    r['drugs_in_panel'] = sorted(keys); r['flagged_drugs_either'] = sorted(either.index[either].tolist())
    T[lab] = r
T['registry_drugs_total'] = int(len(reg)); T['registry_drugs_in_panel'] = int(len(inpanel)); T['not_in_panel'] = sorted(reg[reg.panel_key.isna()].drug.tolist())
out['trial_base_rate'] = T

# figure: variance share and guarded vs original composition
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
for ax, c in zip(axes, ['cog', 'resp']):
    sh = R['E_diagnostics'][c]['composite_variance_share']
    ax.bar(['S', 'E', 'AE'], [sh['S'], sh['E'], sh['AE']], color=['#1f4e79', '#c0392b', '#7f8c8d'])
    ax.axhline(0, c='k', lw=.6); ax.set_ylabel('Share of composite variance'); ax.set_ylim(-0.1, 0.85)
    ax.set_title(('Cognitive' if c == 'cog' else 'Respiratory') + f" cohort (E range = {R['E_diagnostics'][c]['E_range']:.3f})", fontsize=10)
plt.tight_layout(); plt.savefig('figures/fig_composite_variance_share.png', dpi=200); plt.close()

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
for ax, c in zip(axes, ['cog', 'resp']):
    ws = pd.read_csv(f'results/weight_sensitivity_{c}.csv')
    sc = ax.scatter(ws.wS, ws.wE, c=ws.jaccard, cmap='viridis', s=28, vmin=0, vmax=1)
    ax.scatter([.4], [.4], marker='*', s=160, c='red', edgecolor='k', label='default (0.4, 0.4, 0.2)'); ax.legend(fontsize=8, loc='upper right')
    ax.set_xlabel('weight on S'); ax.set_ylabel('weight on E'); ax.set_title(('Cognitive' if c == 'cog' else 'Respiratory') + ': Jaccard of q10 list vs default', fontsize=10)
    plt.colorbar(sc, ax=ax, label='Jaccard')
plt.tight_layout(); plt.savefig('figures/fig_weight_sensitivity.png', dpi=200); plt.close()

for c in P: P[c].reset_index().to_csv(f'results/panel_{c}_guarded.csv', index=False)
json.dump(out, open('results/results_part2.json', 'w'), indent=1, default=str)
print(json.dumps(out, indent=1, default=str)[:9000])
