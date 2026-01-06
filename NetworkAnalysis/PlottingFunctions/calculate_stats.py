import numpy as np
import pandas as pd
from scipy import stats
from typing import List, Dict
import itertools

from statsmodels.stats.multicomp import pairwise_tukeyhsd

def welch_degrees_of_freedom(data1: np.ndarray, data2: np.ndarray) -> float:
    s1_sq, s2_sq = np.var(data1, ddof=1), np.var(data2, ddof=1)
    n1, n2      = len(data1), len(data2)

    num = (s1_sq/n1 + s2_sq/n2) ** 2
    den = ((s1_sq/n1)**2) / (n1 - 1) + ((s2_sq/n2)**2) / (n2 - 1)
    return num / den


def perform_statistical_analysis(
    df: pd.DataFrame,
    data_column: str,
    group_column: str,
    groups: List[str] = None,
    alpha: float = 0.05
) -> Dict:

    if groups is None:
        groups = sorted(df[group_column].dropna().unique())

    results = {'anova': {}, 'pairwise_tests': {}, 'group_stats': {}, 'tukey': None}

    # ── Descriptive stats ────────────────────────────────────────────────────
    group_data = {g: df[df[group_column] == g][data_column].dropna().values
                  for g in groups}

    for g, vals in group_data.items():
        results['group_stats'][g] = {
            'n':     len(vals),
            'mean':  np.mean(vals),
            'std':   np.std(vals, ddof=1),
            'sem':   stats.sem(vals) if len(vals) > 1 else np.nan,
            'median': np.median(vals),
            '95ci':  stats.t.interval(
                        0.95, len(vals)-1,
                        loc=np.mean(vals), scale=stats.sem(vals)
                     ) if len(vals) > 1 else (np.nan, np.nan)
        }

    # ── One-way ANOVA ────────────────────────────────────────────────────────
    all_vals = [group_data[g] for g in groups]
    f_stat, p_anova = stats.f_oneway(*all_vals)
    results['anova'] = {
        'f_statistic': f_stat,
        'p_value':     p_anova,
        'df_between':  len(groups) - 1,
        'df_within':   sum(len(v) for v in all_vals) - len(groups)
    }

    # ── Tukey HSD post-hoc (family-wise α) ───────────────────────────────────
    tukey = pairwise_tukeyhsd(
        endog=df[data_column].dropna().values,
        groups=df[group_column].dropna().values,
        alpha=alpha
    )
    results['tukey'] = tukey        # store full object for printing later

    # ── Pairwise Welch t-tests (optional, kept for effect sizes) ─────────────
    for g1, g2 in itertools.combinations(groups, 2):
        v1, v2 = group_data[g1], group_data[g2]
        t_stat, p_val = stats.ttest_ind(v1, v2, equal_var=False)
        df_welch      = welch_degrees_of_freedom(v1, v2)

        # Cohen’s d (unequal-n pooled)
        n1, n2 = len(v1), len(v2)
        pooled_sd = np.sqrt(((n1-1)*np.var(v1, ddof=1) +
                             (n2-1)*np.var(v2, ddof=1)) / (n1+n2-2))
        cohens_d  = (np.mean(v1) - np.mean(v2)) / pooled_sd

        results['pairwise_tests'][f'{g1}_vs_{g2}'] = {
            't_statistic':      t_stat,
            'p_value':          p_val,
            'degrees_of_freedom': df_welch,
            'cohens_d':         cohens_d,
            'mean_difference':  np.mean(v1) - np.mean(v2)
        }

    return results


def print_statistical_summary(results: Dict, alpha: float = 0.05) -> None:
    # ANOVA
    dfb = results['anova']['df_between']
    dfw = results['anova']['df_within']
    print("\n=== ANOVA ===")
    print(f"F({dfb}, {dfw}) = {results['anova']['f_statistic']:.4f}  "
          f"p = {results['anova']['p_value']:.4g}")

    # Descriptives
    print("\n=== Group Descriptives ===")
    for g, s in results['group_stats'].items():
        print(f"{g:>10s}:  n={s['n']:3d}  "
              f"mean±SEM={s['mean']:.3g}±{s['sem']:.3g}  "
              f"95%CI=[{s['95ci'][0]:.3g}, {s['95ci'][1]:.3g}]")

    # Tukey HSD
    print("\n=== Tukey HSD (α = {:.2f}) ===".format(alpha))
    print(results['tukey'].summary())        # pretty ASCII table

    # Welch t (effect sizes)
    print("\n=== Pairwise Welch t-tests (effect sizes) ===")
    for cmp, s in results['pairwise_tests'].items():
        mark = '*' if s['p_value'] < alpha else 'ns'
        df_s = f"{s['degrees_of_freedom']:.1f}"
        print(f"{cmp:>15s}:  t({df_s})={s['t_statistic']:.3g}  "
              f"p={s['p_value']:.3g} {mark}  "
              f"d={s['cohens_d']:.3g}  "
              f"Δ={s['mean_difference']:.3g}")