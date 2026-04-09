import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

from scipy import stats
from itertools import combinations
from matplotlib.backends.backend_pdf import PdfPages

import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd


class MEAPlotter:
    def __init__(self, group_order=None, palette=None, div_order=None):
        """
        Publication-ish defaults + centralized config.
        """
        # --- FONT CONFIGURATION ---
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans', 'sans-serif']
        plt.rcParams['svg.fonttype'] = 'none'
        plt.rcParams['pdf.fonttype'] = 42

        plt.rcParams['axes.titlesize'] = 9
        plt.rcParams['axes.labelsize'] = 8
        plt.rcParams['xtick.labelsize'] = 6
        plt.rcParams['ytick.labelsize'] = 6
        plt.rcParams['legend.fontsize'] = 7
        plt.rcParams['legend.title_fontsize'] = 8

        # --- TICK APPEARANCE ---
        plt.rcParams['xtick.major.size'] = 3
        plt.rcParams['ytick.major.size'] = 3
        plt.rcParams['xtick.major.width'] = 0.8
        plt.rcParams['ytick.major.width'] = 0.8

        # --- GLOBAL FONT WEIGHT ---
        plt.rcParams['axes.titleweight'] = 'bold'
        plt.rcParams['axes.labelweight'] = 'bold'
        plt.rcParams['font.weight'] = 'bold'

        # --- DEFAULTS ---
        self.group_order_default = list(group_order) if group_order is not None else None
        self.div_order_default = list(div_order) if div_order is not None else None
        self.palette_default = palette  # "Set2" OR {"WT":"#..",...}

    # --------------------------
    # helpers to resolve defaults
    # --------------------------
    def _resolve_group_order(self, df, col, group_order):
        if group_order is not None:
            return list(group_order)
        if self.group_order_default is not None:
            return list(self.group_order_default)
        return sorted(df[col].dropna().unique())

    def _resolve_div_order(self, df, col, div_order):
        if div_order is not None:
            return list(div_order)
        if self.div_order_default is not None:
            return list(self.div_order_default)
        vals = pd.to_numeric(df[col], errors="coerce").dropna().astype(int).unique()
        return sorted(vals)

    def _resolve_palette(self, palette):
        if palette is not None:
            return palette
        return self.palette_default

    def _prepare_data(self, df, cols, div_col=None):
        """
        Common cleaning helper.
        """
        d = df.dropna(subset=cols).copy()

        if div_col is not None:
            d[div_col] = pd.to_numeric(d[div_col], errors="coerce")
            d = d.dropna(subset=[div_col])
            d[div_col] = d[div_col].astype(int)

        return d

    # --------------------------
    # stats helpers
    # --------------------------
    @staticmethod
    def _stars_from_p(p):
        if not np.isfinite(p):
            return "ns"
        if p < 0.001:
            return "***"
        if p < 0.01:
            return "**"
        if p < 0.05:
            return "*"
        return "ns"

    @staticmethod
    def _safe_sem(x):
        x = pd.Series(x).dropna()
        return x.sem() if len(x) > 1 else np.nan

    @staticmethod
    def _safe_cohens_d(d1, d2):
        d1 = pd.Series(d1).dropna()
        d2 = pd.Series(d2).dropna()

        if len(d1) < 2 or len(d2) < 2:
            return np.nan

        sd1 = d1.std(ddof=1)
        sd2 = d2.std(ddof=1)
        pooled_std = np.sqrt((sd1**2 + sd2**2) / 2)

        if not np.isfinite(pooled_std) or pooled_std <= 0:
            return np.nan

        return (d1.mean() - d2.mean()) / pooled_std

    def get_style(
        self,
        figsize=(6, 4),
        title="Title",
        xlabel="X label",
        ylabel="Y label",
        xticks=None,
        xticklabels=None
    ):
        """
        Returns a styled (fig, ax) with the MEAPlotter defaults applied.
        No data is plotted. Intended for quick notebook experimentation.
        """
        fig, ax = plt.subplots(figsize=figsize)

        ax.set_title(title, fontsize=9, fontweight="bold", pad=6)
        ax.set_xlabel(xlabel, fontsize=8, fontweight="bold")
        ax.set_ylabel(ylabel, fontsize=8, fontweight="bold")

        ax.tick_params(axis="x", labelsize=6, width=0.8, length=3)
        ax.tick_params(axis="y", labelsize=6, width=0.8, length=3)

        if xticks is not None:
            ax.set_xticks(xticks)
        if xticklabels is not None:
            ax.set_xticklabels(xticklabels)

        sns.despine(ax=ax)
        return fig, ax

    # --------------------------
    # centralized stats dispatcher
    # --------------------------
    def _run_group_stats(self, df, group_col, y, group_order=None, alpha=0.05):
        """
        Central inferential logic.

        Rules:
        - <2 groups with data -> no test
        - exactly 2 groups -> independent t-test
        - >=3 groups -> one-way ANOVA, then Tukey HSD if omnibus is significant

        Returns dict with:
            test_type
            omnibus_stat
            omnibus_p
            pairwise_df
        """
        d = self._prepare_data(df, [group_col, y])

        if d.empty:
            return {
                "test_type": None,
                "omnibus_stat": np.nan,
                "omnibus_p": np.nan,
                "pairwise_df": pd.DataFrame()
            }

        group_order = self._resolve_group_order(d, group_col, group_order)
        present_groups = []
        for g in group_order:
            n = d.loc[d[group_col] == g, y].dropna().shape[0]
            if n > 0:
                present_groups.append(g)

        if len(present_groups) < 2:
            return {
                "test_type": None,
                "omnibus_stat": np.nan,
                "omnibus_p": np.nan,
                "pairwise_df": pd.DataFrame()
            }

        # ---------------- 2 groups -> t-test ----------------
        if len(present_groups) == 2:
            g1, g2 = present_groups
            d1 = d.loc[d[group_col] == g1, y].dropna()
            d2 = d.loc[d[group_col] == g2, y].dropna()

            if len(d1) <= 1 or len(d2) <= 1:
                pairwise_df = pd.DataFrame([{
                    "group1": g1,
                    "group2": g2,
                    "comparison": f"{g1} vs {g2}",
                    "p_adj": np.nan,
                    "raw_p": np.nan,
                    "reject": False,
                    "stat": np.nan,
                    "Grp1_Stats": f"{d1.mean():.2f} ± {self._safe_sem(d1):.2f} (n={len(d1)})" if len(d1) > 0 else f"NA (n={len(d1)})",
                    "Grp2_Stats": f"{d2.mean():.2f} ± {self._safe_sem(d2):.2f} (n={len(d2)})" if len(d2) > 0 else f"NA (n={len(d2)})",
                    "Cohen's d": np.nan
                }])
                return {
                    "test_type": "t-test",
                    "omnibus_stat": np.nan,
                    "omnibus_p": np.nan,
                    "pairwise_df": pairwise_df
                }

            # Standard 2-group independent t-test
            t_stat, p_val = stats.ttest_ind(d1, d2, equal_var=True, nan_policy="omit")
            cohens_d = self._safe_cohens_d(d1, d2)

            pairwise_df = pd.DataFrame([{
                "group1": g1,
                "group2": g2,
                "comparison": f"{g1} vs {g2}",
                "p_adj": p_val,
                "raw_p": p_val,
                "reject": bool(p_val < alpha) if np.isfinite(p_val) else False,
                "stat": t_stat,
                "Grp1_Stats": f"{d1.mean():.2f} ± {self._safe_sem(d1):.2f} (n={len(d1)})",
                "Grp2_Stats": f"{d2.mean():.2f} ± {self._safe_sem(d2):.2f} (n={len(d2)})",
                "Cohen's d": cohens_d
            }])

            return {
                "test_type": "t-test",
                "omnibus_stat": t_stat,
                "omnibus_p": p_val,
                "pairwise_df": pairwise_df
            }

        # ---------------- 3+ groups -> one-way ANOVA ----------------
        # safer formula quoting
        yq = f'Q("{y}")'
        gq = f'Q("{group_col}")'
        formula = f"{yq} ~ C({gq})"

        try:
            model = ols(formula, data=d).fit()
            anova_table = sm.stats.anova_lm(model, typ=2)
            row_name = f"C({gq})"
            omnibus_stat = anova_table.loc[row_name, "F"]
            omnibus_p = anova_table.loc[row_name, "PR(>F)"]
        except Exception:
            omnibus_stat = np.nan
            omnibus_p = np.nan

        pairwise_rows = []

        if np.isfinite(omnibus_p) and omnibus_p < alpha:
            tukey = pairwise_tukeyhsd(endog=d[y], groups=d[group_col], alpha=alpha)
            tukey_df = pd.DataFrame(
                tukey._results_table.data[1:],
                columns=tukey._results_table.data[0]
            )

            tukey_df = tukey_df.rename(columns={
                "p-adj": "p_adj",
                "meandiff": "stat"
            })

            for _, row in tukey_df.iterrows():
                g1 = row["group1"]
                g2 = row["group2"]

                d1 = d.loc[d[group_col] == g1, y].dropna()
                d2 = d.loc[d[group_col] == g2, y].dropna()

                pairwise_rows.append({
                    "group1": g1,
                    "group2": g2,
                    "comparison": f"{g1} vs {g2}",
                    "p_adj": float(row["p_adj"]) if pd.notna(row["p_adj"]) else np.nan,
                    "raw_p": np.nan,
                    "reject": bool(row["reject"]),
                    "stat": row["stat"],
                    "Grp1_Stats": f"{d1.mean():.2f} ± {self._safe_sem(d1):.2f} (n={len(d1)})" if len(d1) > 0 else f"NA (n={len(d1)})",
                    "Grp2_Stats": f"{d2.mean():.2f} ± {self._safe_sem(d2):.2f} (n={len(d2)})" if len(d2) > 0 else f"NA (n={len(d2)})",
                    "Cohen's d": self._safe_cohens_d(d1, d2)
                })
        else:
            # still return rows for all pairwise combinations so tables stay informative
            pairs = list(combinations(present_groups, 2))
            for g1, g2 in pairs:
                d1 = d.loc[d[group_col] == g1, y].dropna()
                d2 = d.loc[d[group_col] == g2, y].dropna()

                pairwise_rows.append({
                    "group1": g1,
                    "group2": g2,
                    "comparison": f"{g1} vs {g2}",
                    "p_adj": np.nan,
                    "raw_p": np.nan,
                    "reject": False,
                    "stat": np.nan,
                    "Grp1_Stats": f"{d1.mean():.2f} ± {self._safe_sem(d1):.2f} (n={len(d1)})" if len(d1) > 0 else f"NA (n={len(d1)})",
                    "Grp2_Stats": f"{d2.mean():.2f} ± {self._safe_sem(d2):.2f} (n={len(d2)})" if len(d2) > 0 else f"NA (n={len(d2)})",
                    "Cohen's d": self._safe_cohens_d(d1, d2)
                })

        pairwise_df = pd.DataFrame(pairwise_rows)

        return {
            "test_type": "one-way ANOVA",
            "omnibus_stat": omnibus_stat,
            "omnibus_p": omnibus_p,
            "pairwise_df": pairwise_df
        }

    def calculate_stats(self, df, x, y, order=None, alpha=0.05):
        """
        Single-panel public stats API.
        Centralized version:
        - 2 groups -> t-test
        - 3+ groups -> one-way ANOVA + Tukey
        """
        d = self._prepare_data(df, [x, y])
        order = self._resolve_group_order(d, x, order)

        result = self._run_group_stats(
            d,
            group_col=x,
            y=y,
            group_order=order,
            alpha=alpha
        )

        pairwise_df = result["pairwise_df"].copy()
        if pairwise_df.empty:
            return pd.DataFrame(columns=[
                "Comparison", "Grp1_Stats", "Grp2_Stats",
                "Stat", "p-val", "Sig", "Cohen's d",
                "Test", "Omnibus Stat", "Omnibus p-val"
            ])

        out_rows = []
        for _, row in pairwise_df.iterrows():
            p = row.get("p_adj", np.nan)
            out_rows.append({
                "Comparison": row["comparison"],
                "Grp1_Stats": row["Grp1_Stats"],
                "Grp2_Stats": row["Grp2_Stats"],
                "Stat": row["stat"],
                "p-val": p,
                "Sig": self._stars_from_p(p),
                "Cohen's d": row["Cohen's d"],
                "Test": result["test_type"],
                "Omnibus Stat": result["omnibus_stat"],
                "Omnibus p-val": result["omnibus_p"]
            })

        return pd.DataFrame(out_rows)

    def calculate_stats_by_div(self, df, div_col, group_col, y, div_order=None, group_order=None, alpha=0.05):
        """
        Per-DIV public stats API.
        Centralized version:
        - 2 groups -> t-test
        - 3+ groups -> one-way ANOVA + Tukey
        """
        d = self._prepare_data(df, [div_col, group_col, y], div_col=div_col)
        div_order = self._resolve_div_order(d, div_col, div_order)
        group_order = self._resolve_group_order(d, group_col, group_order)

        rows = []

        print(f"\n{'='*20} CENTRALIZED STATS: {y} by {div_col} {'='*20}")

        for div_val in div_order:
            div_data = d[d[div_col] == div_val].copy()
            if div_data.empty:
                continue

            result = self._run_group_stats(
                div_data,
                group_col=group_col,
                y=y,
                group_order=group_order,
                alpha=alpha
            )

            print(f"\n--- DIV {div_val} ---")
            print(f"Test: {result['test_type']}, omnibus p = {result['omnibus_p']}")

            pairwise_df = result["pairwise_df"]

            if pairwise_df.empty:
                rows.append({
                    "DIV": div_val,
                    "Comparison": "None",
                    "Grp1_Stats": "",
                    "Grp2_Stats": "",
                    "Stat": np.nan,
                    "p-val": np.nan,
                    "Sig": "ns",
                    "Cohen's d": np.nan,
                    "Test": result["test_type"],
                    "Omnibus Stat": result["omnibus_stat"],
                    "Omnibus p-val": result["omnibus_p"],
                    "group1": np.nan,
                    "group2": np.nan,
                    "reject": False
                })
                continue

            for _, row in pairwise_df.iterrows():
                p = row.get("p_adj", np.nan)
                stars = self._stars_from_p(p)

                print(f"  {row['comparison']}: p={p}, {stars}")

                rows.append({
                    "DIV": div_val,
                    "Comparison": row["comparison"],
                    "Grp1_Stats": row["Grp1_Stats"],
                    "Grp2_Stats": row["Grp2_Stats"],
                    "Stat": row["stat"],
                    "p-val": p,
                    "Sig": stars,
                    "Cohen's d": row["Cohen's d"],
                    "Test": result["test_type"],
                    "Omnibus Stat": result["omnibus_stat"],
                    "Omnibus p-val": result["omnibus_p"],
                    "group1": row["group1"],
                    "group2": row["group2"],
                    "reject": bool(row["reject"])
                })

        return pd.DataFrame(rows)

    # --------------------------
    # backward-compatible Welch APIs
    # --------------------------
    def calculate_stats_welch(self, df, group_col, y, group_order=None):
        """
        Backward-compatible Welch t-test (uncorrected), pairwise across groups.
        Preserved so existing notebook code does not break.
        """
        data = df.dropna(subset=[group_col, y]).copy()
        group_order = self._resolve_group_order(data, group_col, group_order)
        pairs = list(combinations(group_order, 2))
        rows = []

        print(f"\n{'='*20} STATS (Welch, uncorrected): {y} {'='*20}")

        for g1, g2 in pairs:
            d1 = data.loc[data[group_col] == g1, y].dropna()
            d2 = data.loc[data[group_col] == g2, y].dropna()

            n1, n2 = len(d1), len(d2)
            if n1 <= 1 or n2 <= 1:
                print(f"{g1} vs {g2}: not enough data (n1={n1}, n2={n2})")
                rows.append({
                    "Comparison": f"{g1} vs {g2}",
                    "Grp1_Stats": f"NA (n={n1})" if n1 == 0 else f"{d1.mean():.2f} ± {d1.sem():.2f} (n={n1})",
                    "Grp2_Stats": f"NA (n={n2})" if n2 == 0 else f"{d2.mean():.2f} ± {d2.sem():.2f} (n={n2})",
                    "t-stat": np.nan,
                    "p-val": np.nan,
                    "Sig": "ns",
                    "Cohen's d": np.nan
                })
                continue

            t_stat, p_val = stats.ttest_ind(d1, d2, equal_var=False, nan_policy="omit")
            mean1, mean2 = d1.mean(), d2.mean()
            sem1, sem2 = d1.sem(), d2.sem()
            cohens_d = self._safe_cohens_d(d1, d2)
            stars = self._stars_from_p(p_val)

            print(f"{g1} (n={n1}) vs {g2} (n={n2}): p={p_val:.4e} ({stars}), d={cohens_d:.3f}")

            rows.append({
                "Comparison": f"{g1} vs {g2}",
                "Grp1_Stats": f"{mean1:.2f} ± {sem1:.2f} (n={n1})",
                "Grp2_Stats": f"{mean2:.2f} ± {sem2:.2f} (n={n2})",
                "t-stat": t_stat,
                "p-val": p_val,
                "Sig": stars,
                "Cohen's d": cohens_d
            })

        return pd.DataFrame(rows)

    def calculate_stats_by_div_welch(self, df, div_col, group_col, y, div_order=None, group_order=None):
        """
        Backward-compatible Welch t-tests (uncorrected) per DIV.
        Preserved so existing notebook code does not break.
        """
        data = df.dropna(subset=[div_col, group_col, y]).copy()
        data[div_col] = pd.to_numeric(data[div_col], errors="coerce")
        data = data.dropna(subset=[div_col])
        data[div_col] = data[div_col].astype(int)

        div_order = self._resolve_div_order(data, div_col, div_order)
        group_order = self._resolve_group_order(data, group_col, group_order)
        pairs = list(combinations(group_order, 2))
        rows = []

        print(f"\n{'='*20} STATS (Welch, uncorrected): {y} by {div_col} {'='*20}")

        for div_val in div_order:
            div_data = data[data[div_col] == div_val]
            if div_data.empty:
                continue

            print(f"\n--- DIV {div_val} ---")
            for g1, g2 in pairs:
                d1 = div_data.loc[div_data[group_col] == g1, y].dropna()
                d2 = div_data.loc[div_data[group_col] == g2, y].dropna()

                n1, n2 = len(d1), len(d2)
                if n1 <= 1 or n2 <= 1:
                    print(f"  {g1} vs {g2}: not enough data (n1={n1}, n2={n2})")
                    rows.append({
                        "DIV": div_val,
                        "Comparison": f"{g1} vs {g2}",
                        "Grp1_Stats": f"NA (n={n1})" if n1 == 0 else f"{d1.mean():.2f} ± {d1.sem():.2f} (n={n1})",
                        "Grp2_Stats": f"NA (n={n2})" if n2 == 0 else f"{d2.mean():.2f} ± {d2.sem():.2f} (n={n2})",
                        "t-stat": np.nan,
                        "p-val": np.nan,
                        "Sig": "ns",
                        "Cohen's d": np.nan
                    })
                    continue

                t_stat, p_val = stats.ttest_ind(d1, d2, equal_var=False, nan_policy="omit")
                mean1, mean2 = d1.mean(), d2.mean()
                sem1, sem2 = d1.sem(), d2.sem()
                cohens_d = self._safe_cohens_d(d1, d2)
                stars = self._stars_from_p(p_val)

                print(f"  {g1} (n={n1}) vs {g2} (n={n2}): p={p_val:.4e} ({stars}), d={cohens_d:.3f}")

                rows.append({
                    "DIV": div_val,
                    "Comparison": f"{g1} vs {g2}",
                    "Grp1_Stats": f"{mean1:.2f} ± {sem1:.2f} (n={n1})",
                    "Grp2_Stats": f"{mean2:.2f} ± {sem2:.2f} (n={n2})",
                    "t-stat": t_stat,
                    "p-val": p_val,
                    "Sig": stars,
                    "Cohen's d": cohens_d
                })

        return pd.DataFrame(rows)

    # --------------------------
    # annotation helper
    # --------------------------
    def _annotate_stats_brackets(self, ax, stats_df, div_order, div_to_x, offsets, y_top_reserved_frac=0.20):
        """
        Draw significance brackets from centralized stats output.
        Expects stats_df with columns:
            DIV, group1, group2, reject, p-val, Sig
        """
        if stats_df is None or len(stats_df) == 0:
            return ax

        ymin_ax, ymax_ax = ax.get_ylim()
        yr = ymax_ax - ymin_ax
        base_top = ymax_ax - y_top_reserved_frac * yr
        bracket_step = 0.06 * yr
        max_annotation_y = base_top

        for div_val in div_order:
            div_tests = stats_df[
                (stats_df["DIV"] == div_val) &
                (stats_df["reject"] == True)
            ].copy()

            if div_tests.empty:
                continue

            level = 0
            for _, t in div_tests.iterrows():
                g1 = t["group1"]
                g2 = t["group2"]
                stars = t["Sig"]

                if stars == "ns":
                    continue
                if pd.isna(g1) or pd.isna(g2):
                    continue
                if g1 not in offsets or g2 not in offsets:
                    continue

                x1 = div_to_x[div_val] + offsets[g1]
                x2 = div_to_x[div_val] + offsets[g2]
                bracket_y = base_top + level * bracket_step

                ax.plot(
                    [x1, x1, x2, x2],
                    [
                        bracket_y,
                        bracket_y + 0.25 * bracket_step,
                        bracket_y + 0.25 * bracket_step,
                        bracket_y
                    ],
                    c="black",
                    lw=1.0,
                    zorder=7
                )

                ax.text(
                    (x1 + x2) / 2,
                    bracket_y + 0.30 * bracket_step,
                    stars,
                    ha="center",
                    va="bottom",
                    fontsize=6,
                    zorder=8
                )

                max_annotation_y = max(max_annotation_y, bracket_y + 0.55 * bracket_step)
                level += 1

        if max_annotation_y > ax.get_ylim()[1]:
            ax.set_ylim(ax.get_ylim()[0], max_annotation_y + 0.04 * yr)

        return ax

    # --------------------------
    # plots
    # --------------------------
    def plot_bars(self, df, x, y, order=None, palette=None, title=None, ax=None, annotate=False):
        """
        Simple grouped bar + scatter (single panel).
        annotate=False by default.
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6))

        data = df.dropna(subset=[x, y]).copy()
        order = self._resolve_group_order(data, x, order)
        palette = self._resolve_palette(palette)

        sns.barplot(
            data=data, x=x, y=y, order=order, palette=palette,
            errorbar='se', capsize=0.12, alpha=0.75,
            edgecolor='black', linewidth=0.8, ax=ax
        )

        sns.stripplot(
            data=data, x=x, y=y, order=order,
            color='black', alpha=0.7, jitter=0.18, size=5,
            edgecolor='black', linewidth=0.3, ax=ax
        )

        ax.set_title(title if title else y, fontweight='bold', pad=10)
        ax.set_xlabel("")
        ax.set_ylabel(y.replace('_', ' '), fontweight='bold')
        sns.despine(ax=ax)

        # Optional single-panel annotation using centralized stats
        if annotate:
            stats_df = self.calculate_stats(df, x=x, y=y, order=order, alpha=0.05)

            if len(stats_df) > 0:
                ymin_ax, ymax_ax = ax.get_ylim()
                yr = ymax_ax - ymin_ax
                base_y = ymax_ax + 0.02 * yr
                step = 0.08 * yr

                pos_map = {g: i for i, g in enumerate(order)}
                level = 0

                for _, row in stats_df.iterrows():
                    p = row["p-val"]
                    stars = row["Sig"]

                    if stars == "ns" or not np.isfinite(p):
                        continue

                    g1, g2 = row["Comparison"].split(" vs ")
                    x1, x2 = pos_map[g1], pos_map[g2]
                    yb = base_y + level * step

                    ax.plot([x1, x1, x2, x2], [yb, yb + 0.2 * step, yb + 0.2 * step, yb], c="black", lw=1.1)
                    ax.text((x1 + x2) / 2, yb + 0.25 * step, stars, ha="center", va="bottom", fontsize=9)
                    level += 1

                ax.set_ylim(ymin_ax, base_y + (level + 1) * step)

        return ax

    def plot_line_sem_by_div(
        self,
        df,
        div_col,
        group_col,
        y,
        div_order=None,
        group_order=None,
        palette=None,
        title=None,
        ylabel=None,
        xlabel=None,
        ax=None,
        fig_width=6,
        aspect=1.0,
        scatter=True,
        sem_fill=True,
        annotate=True,
        # x-position control
        group_offset_width=0.25,
        scatter_jitter=0.04,
        # clipping / outlier display
        clip_upper=True,
        clip_method="quantile",
        clip_q=0.98,
        clip_k=1.5,
        show_upper_outliers=True,
        # legend
        show_legend=True,
        legend_outside=True,
        legend_loc="upper right"
    ):
        """
        Developmental trajectory plot:
        - mean line
        - SEM fill
        - individual well scatter
        - genotype-specific x-offset within each DIV
        - upper clipped outliers shown as lollipops/triangles
        - centralized stats:
            * 2 groups -> t-test
            * 3+ groups -> one-way ANOVA + Tukey
        """
        d = self._prepare_data(df, [div_col, group_col, y], div_col=div_col)

        div_order = self._resolve_div_order(d, div_col, div_order)
        group_order = self._resolve_group_order(d, group_col, group_order)
        palette = self._resolve_palette(palette)
        div_to_x = {div: i for i, div in enumerate(div_order)}

        if ax is None:
            fig_height = fig_width * aspect
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        # ---------------- Colors ----------------
        if isinstance(palette, dict):
            color_map = {g: palette[g] for g in group_order}
        else:
            colors = sns.color_palette(palette, len(group_order))
            color_map = {g: colors[i] for i, g in enumerate(group_order)}

        # ---------------- X offsets ----------------
        n_groups = len(group_order)
        if n_groups == 1:
            offsets = {group_order[0]: 0.0}
        else:
            raw_offsets = np.linspace(-group_offset_width, group_offset_width, n_groups)
            offsets = {g: raw_offsets[i] for i, g in enumerate(group_order)}

        # ---------------- Robust upper clip ----------------
        cap = None
        if clip_upper:
            yvals = d[y].dropna()
            if len(yvals) > 0:
                if clip_method == "iqr":
                    q1, q3 = yvals.quantile([0.25, 0.75])
                    iqr = float(q3 - q1)
                    if not np.isfinite(iqr) or iqr == 0:
                        cap = float(yvals.quantile(clip_q))
                    else:
                        cap = float(q3 + clip_k * iqr)
                else:
                    cap = float(yvals.quantile(clip_q))

        # ---------------- Plot each group ----------------
        for g in group_order:
            gdata = d[d[group_col] == g].copy()
            color = color_map[g]
            xoff = offsets[g]

            means = []
            sems = []
            xline = []

            for div in div_order:
                vals = gdata.loc[gdata[div_col] == div, y].dropna()
                means.append(vals.mean() if len(vals) > 0 else np.nan)
                sems.append(vals.sem() if len(vals) > 1 else 0.0)
                xline.append(div_to_x[div] + xoff)

            means = np.array(means, dtype=float)
            sems = np.array(sems, dtype=float)
            xline = np.array(xline, dtype=float)

            if sem_fill:
                ax.fill_between(
                    xline,
                    means - sems,
                    means + sems,
                    color=color,
                    alpha=0.25,
                    linewidth=0,
                    zorder=1
                )

            ax.plot(
                xline,
                means,
                color=color,
                linewidth=1.0,
                marker="o",
                markersize=2,
                label=g,
                zorder=3
            )

            if scatter and len(gdata) > 0:
                xvals = gdata[div_col].map(div_to_x).astype(float).values + xoff
                xvals = xvals + np.random.uniform(-scatter_jitter, scatter_jitter, size=len(xvals))

                if cap is not None:
                    main_mask = gdata[y].values <= cap
                else:
                    main_mask = np.ones(len(gdata), dtype=bool)

                ax.scatter(
                    xvals[main_mask],
                    gdata[y].values[main_mask],
                    color=color,
                    alpha=0.6,
                    s=10,
                    edgecolor="black",
                    linewidth=0.25,
                    zorder=2
                )

        # ---------------- Apply clipped y-limit ----------------
        if cap is not None:
            ymin = float(d[y].min())
            pad = 0.08 * max(cap - ymin, 1e-9)
            ax.set_ylim(ymin, cap + pad)

        # ---------------- Upper outliers ----------------
        if show_upper_outliers and cap is not None:
            upper_out = d[d[y] > cap].copy()
            if len(upper_out) > 0:
                ymin_ax, ymax_ax = ax.get_ylim()
                yr = ymax_ax - ymin_ax

                y_tip = ymax_ax - 0.01 * yr
                y_stem0 = y_tip - 0.04 * yr

                for _, r in upper_out.iterrows():
                    div_val = r[div_col]
                    g = r[group_col]
                    if div_val not in div_order or g not in group_order:
                        continue

                    x = float(div_to_x[div_val]) + offsets[g]

                    ax.plot([x, x], [y_stem0, y_tip], color=color_map[g], lw=0.9, zorder=5)
                    ax.scatter(
                        [x], [y_tip],
                        marker="^", s=26,
                        color=color_map[g],
                        edgecolor="black",
                        linewidth=0.25,
                        zorder=6
                    )

        # ---------------- Centralized stats + annotation ----------------
        if annotate:
            stats_df = self.calculate_stats_by_div(
                df=d,
                div_col=div_col,
                group_col=group_col,
                y=y,
                div_order=div_order,
                group_order=group_order,
                alpha=0.05
            )
            self._annotate_stats_brackets(ax, stats_df, div_order, div_to_x, offsets)

        # ---------------- Axes / ticks / labels ----------------
        ax.set_xticks(range(len(div_order)))
        ax.set_xticklabels([str(v) for v in div_order])
        ax.set_xlim(-0.5, len(div_order) - 0.5)

        ax.set_xlabel(xlabel if xlabel else div_col, fontweight="bold", fontsize=8)
        ax.set_ylabel(ylabel if ylabel else y.replace("_", " "), fontweight="bold", fontsize=8)
        ax.set_title(title if title else f"{y} trajectory", fontweight="bold", pad=6, fontsize=9)

        ax.tick_params(axis="x", labelsize=6, width=0.8, length=3)
        ax.tick_params(axis="y", labelsize=6, width=0.8, length=3)

        if show_legend:
            if legend_outside:
                ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, title_fontsize=8)
            else:
                ax.legend(frameon=False, loc=legend_loc, fontsize=7, title_fontsize=8)

        sns.despine(ax=ax)
        return ax

    def plot_bars_by_div(
        self,
        df,
        div_col,
        group_col,
        y,
        div_order=None,
        group_order=None,
        palette=None,
        title=None,
        ax=None,
        fig_width=6,
        aspect=1.0,
        annotate=True,
        # clipping / outliers
        clip_upper=True,
        clip_method="quantile",
        clip_k=1.5,
        clip_q=0.98,
        show_upper_outliers=True,
        # legend control
        legend_loc="upper right",
        legend_outside=False
    ):
        """
        Grouped bars by DIV with scatter.
        Centralized stats:
            * 2 groups -> t-test
            * 3+ groups -> one-way ANOVA + Tukey
        """
        data = self._prepare_data(df, [div_col, group_col, y], div_col=div_col)

        div_order = self._resolve_div_order(data, div_col, div_order)
        group_order = self._resolve_group_order(data, group_col, group_order)
        palette = self._resolve_palette(palette)

        if ax is None:
            fig_height = fig_width * aspect
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        sns.barplot(
            data=data,
            x=div_col, y=y,
            hue=group_col,
            order=div_order,
            hue_order=group_order,
            palette=palette,
            errorbar="se",
            capsize=0.10,
            alpha=0.75,
            edgecolor='black',
            linewidth=0.8,
            ax=ax
        )

        sns.stripplot(
            data=data,
            x=div_col, y=y,
            hue=group_col,
            order=div_order,
            hue_order=group_order,
            dodge=True,
            jitter=0.15,
            alpha=0.75,
            size=5,
            edgecolor='black',
            linewidth=0.3,
            ax=ax
        )

        # ---- legend cleanup ----
        handles, labels = ax.get_legend_handles_labels()
        wanted = [str(g) for g in group_order]
        kept_handles, kept_labels, seen = [], [], set()
        for h, lab in zip(handles, labels):
            if lab in wanted and lab not in seen:
                kept_handles.append(h)
                kept_labels.append(lab)
                seen.add(lab)

        if legend_outside:
            ax.legend(
                kept_handles, kept_labels, title=group_col,
                frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left"
            )
        else:
            ax.legend(kept_handles, kept_labels, title=group_col, frameon=False, loc=legend_loc)

        # ---- clip cap ----
        cap = None
        if clip_upper:
            yvals = data[y].dropna()
            if len(yvals) > 0:
                if clip_method == "iqr":
                    q1, q3 = yvals.quantile([0.25, 0.75])
                    iqr = float(q3 - q1)
                    if not np.isfinite(iqr) or iqr == 0:
                        cap = float(yvals.quantile(clip_q))
                    else:
                        cap = float(q3 + clip_k * iqr)
                else:
                    cap = float(yvals.quantile(clip_q))

                ymin = float(yvals.min())
                pad = 0.06 * max(cap - ymin, 1e-9)
                ax.set_ylim(ymin, cap + pad)

        # ---- upper outliers ----
        if show_upper_outliers and cap is not None:
            upper_out = data[data[y] > cap]
            if len(upper_out) > 0:
                n_groups = len(group_order)
                bar_width = 0.8 / n_groups

                ymin_ax, ymax_ax = ax.get_ylim()
                yr = ymax_ax - ymin_ax

                y_tip = ymax_ax - 0.01 * yr
                y_stem0 = y_tip - 0.06 * yr

                for _, r in upper_out.iterrows():
                    div_val = r[div_col]
                    g = r[group_col]
                    if div_val not in div_order or g not in group_order:
                        continue

                    div_idx = div_order.index(div_val)
                    g_idx = group_order.index(g)
                    x = div_idx + (g_idx - n_groups / 2 + 0.5) * bar_width

                    ax.plot([x, x], [y_stem0, y_tip], color="black", lw=0.9, zorder=10)
                    ax.scatter([x], [y_tip], marker="^", s=42, color="black", zorder=11)

        # ---- centralized stats annotations ----
        if annotate:
            stats_df = self.calculate_stats_by_div(
                df=data,
                div_col=div_col,
                group_col=group_col,
                y=y,
                div_order=div_order,
                group_order=group_order,
                alpha=0.05
            )

            # build offsets to match grouped bars
            n_groups = len(group_order)
            bar_width = 0.8 / n_groups
            offsets = {
                g: (i - n_groups / 2 + 0.5) * bar_width
                for i, g in enumerate(group_order)
            }
            div_to_x = {div: i for i, div in enumerate(div_order)}

            self._annotate_stats_brackets(ax, stats_df, div_order, div_to_x, offsets)

        ax.set_title(title if title else f"{y} by {div_col}", fontweight="bold", pad=10)
        ax.set_xlabel("")
        ax.set_ylabel(y.replace("_", " "), fontweight="bold")
        ax.tick_params(axis="both", labelsize=11)
        sns.despine(ax=ax)

        return ax

    # --------------------------
    # wrapper: save pdf
    # --------------------------
    def save_pdf(
        self,
        df,
        y,
        mode="by_div",               # "by_div" or "by_group"
        div_col="DIV",
        group_col="NeuronType",
        x=None,                      # only for mode="by_group"
        div_order=None,
        group_order=None,
        palette=None,
        title=None,
        filename="Report.pdf",
        pdf_obj=None,
        # plot options forwarded
        annotate=True,
        clip_upper=True,
        clip_method="quantile",
        clip_q=0.98,
        clip_k=1.5,
        show_upper_outliers=True,
        legend_outside=False
    ):
        """
        Orchestrator:
        - runs centralized stats
        - makes plot
        - adds stats table
        - saves a single PDF page

        Stats rule:
        - 2 groups -> t-test
        - 3+ groups -> one-way ANOVA + Tukey
        """
        if not filename.endswith(".pdf"):
            filename += ".pdf"

        palette = self._resolve_palette(palette)

        # ---- compute stats ----
        if mode == "by_div":
            stats_df = self.calculate_stats_by_div(
                df=df,
                div_col=div_col,
                group_col=group_col,
                y=y,
                div_order=div_order,
                group_order=group_order,
                alpha=0.05
            )
        elif mode == "by_group":
            if x is None:
                raise ValueError("For mode='by_group', you must pass x=<group column> (e.g., 'NeuronType').")
            stats_df = self.calculate_stats(
                df=df,
                x=x,
                y=y,
                order=group_order,
                alpha=0.05
            )
        else:
            raise ValueError("mode must be 'by_div' or 'by_group'")

        # ---- build fig ----
        fig = plt.figure(figsize=(12, 11))
        gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.0])
        ax_plot = fig.add_subplot(gs[0])
        ax_table = fig.add_subplot(gs[1])
        ax_table.axis("off")

        # ---- plot ----
        if mode == "by_div":
            self.plot_bars_by_div(
                df=df,
                div_col=div_col,
                group_col=group_col,
                y=y,
                div_order=div_order,
                group_order=group_order,
                palette=palette,
                title=title,
                ax=ax_plot,
                annotate=annotate,
                clip_upper=clip_upper,
                clip_method=clip_method,
                clip_q=clip_q,
                clip_k=clip_k,
                show_upper_outliers=show_upper_outliers,
                legend_outside=legend_outside
            )
        else:
            self.plot_bars(
                df=df,
                x=x,
                y=y,
                order=group_order,
                palette=palette,
                title=title,
                ax=ax_plot,
                annotate=annotate
            )

        # ---- table ----
        if stats_df is None or len(stats_df) == 0:
            ax_table.text(0.5, 0.5, "No stats to display", ha="center", va="center")
        else:
            # choose columns
            if mode == "by_div" and "DIV" in stats_df.columns:
                cols = [
                    "DIV", "Test", "Omnibus Stat", "Omnibus p-val",
                    "Comparison", "Grp1_Stats", "Grp2_Stats",
                    "Stat", "p-val", "Sig", "Cohen's d"
                ]
            else:
                cols = [
                    "Test", "Omnibus Stat", "Omnibus p-val",
                    "Comparison", "Grp1_Stats", "Grp2_Stats",
                    "Stat", "p-val", "Sig", "Cohen's d"
                ]

            # format
            cell_text = []
            for _, row in stats_df.iterrows():
                row_vals = []
                for c in cols:
                    v = row.get(c, "")
                    if isinstance(v, float) and np.isnan(v):
                        row_vals.append("NA")
                    elif c in ["Stat", "Omnibus Stat"] and isinstance(v, (float, np.floating)):
                        row_vals.append(f"{v:.2f}")
                    elif c in ["p-val", "Omnibus p-val"] and isinstance(v, (float, np.floating)):
                        row_vals.append(f"{v:.3e}")
                    elif c == "Cohen's d" and isinstance(v, (float, np.floating)):
                        row_vals.append(f"{v:.2f}" if np.isfinite(v) else "NA")
                    else:
                        row_vals.append(str(v))
                cell_text.append(row_vals)

            table = ax_table.table(
                cellText=cell_text,
                colLabels=cols,
                loc="center",
                cellLoc="center"
            )
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1.05, 1.55)

            for (r, c), cell in table.get_celld().items():
                if r == 0:
                    cell.set_text_props(weight="bold")

        # ---- layout control ----
        right = 0.82 if legend_outside else 0.95
        fig.subplots_adjust(top=0.92, bottom=0.05, left=0.07, right=right, hspace=0.15)

        # ---- save ----
        if pdf_obj is not None:
            pdf_obj.savefig(fig)
        else:
            with PdfPages(filename) as pdf:
                pdf.savefig(fig)
            print(f"✅ PDF saved to: {filename}")

        plt.close(fig)