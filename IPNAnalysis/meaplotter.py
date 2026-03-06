import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations
from matplotlib.backends.backend_pdf import PdfPages


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

    # --------------------------
    # stats
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

    def calculate_stats_welch(self, df, group_col, y, group_order=None):
        """
        Welch t-test (uncorrected), pairwise across groups (single x-axis panel).
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
            pooled_std = np.sqrt((d1.std(ddof=1)**2 + d2.std(ddof=1)**2) / 2)
            cohens_d = (mean1 - mean2) / pooled_std if np.isfinite(pooled_std) and pooled_std > 0 else np.nan
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
        Welch t-tests (uncorrected) per DIV, all pairwise group comparisons.
        Returns DataFrame with a DIV column for table output.
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
                pooled_std = np.sqrt((d1.std(ddof=1)**2 + d2.std(ddof=1)**2) / 2)
                cohens_d = (mean1 - mean2) / pooled_std if np.isfinite(pooled_std) and pooled_std > 0 else np.nan
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
    # plots
    # --------------------------
    def plot_bars(self, df, x, y, order=None, palette=None, title=None, ax=None, annotate=False):
        """
        Simple grouped bar + scatter (single panel).
        annotate=False by default (since you're not focusing on this right now).
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
        scatter=True,
        sem_fill=True,
        annotate=True,
        # x-position control
        group_offset_width=0.25,   # spread groups slightly within each DIV
        scatter_jitter=0.04,       # small extra jitter around each group offset
        # clipping / outlier display
        clip_upper=True,
        clip_method="quantile",    # "quantile" or "iqr"
        clip_q=0.98,
        clip_k=1.5,
        show_upper_outliers=True,
        # legend
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
        - uncorrected Welch t-test annotations per DIV
        """

        # ---------------- Data prep ----------------
        d = df.dropna(subset=[div_col, group_col, y]).copy()
        d[div_col] = pd.to_numeric(d[div_col], errors="coerce")
        d = d.dropna(subset=[div_col])
        d[div_col] = d[div_col].astype(int)

        div_order = self._resolve_div_order(d, div_col, div_order)
        group_order = self._resolve_group_order(d, group_col, group_order)
        palette = self._resolve_palette(palette)
        div_to_x = {div: i for i, div in enumerate(div_order)}

        if ax is None:
            fig, ax = plt.subplots(figsize=(max(8, len(div_order) * 1.0 + 2), 5))

        # ---------------- Colors ----------------
        if isinstance(palette, dict):
            color_map = {g: palette[g] for g in group_order}
        else:
            colors = sns.color_palette(palette, len(group_order))
            color_map = {g: colors[i] for i, g in enumerate(group_order)}

        # ---------------- X offsets for groups within each DIV ----------------
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

            # SEM fill
            if sem_fill:
                ax.fill_between(
                    xline,
                    means - sems,
                    means + sems,
                    color=color,
                    alpha=0.18,
                    linewidth=0,
                    zorder=1
                )

            # Mean line
            ax.plot(
                xline,
                means,
                color=color,
                linewidth=1.5,
                marker="o",
                markersize=3.5,
                label=g,
                zorder=3
            )

            # Scatter wells
            if scatter and len(gdata) > 0:
                xvals = gdata[div_col].map(div_to_x).astype(float).values + xoff
                xvals = xvals + np.random.uniform(-scatter_jitter, scatter_jitter, size=len(xvals))

                if cap is not None:
                    # only non-upper-outliers in main scatter
                    main_mask = gdata[y].values <= cap
                else:
                    main_mask = np.ones(len(gdata), dtype=bool)

                ax.scatter(
                    xvals[main_mask],
                    gdata[y].values[main_mask],
                    color=color,
                    alpha=0.70,
                    s=16,
                    edgecolor="black",
                    linewidth=0.25,
                    zorder=2
                )

        # ---------------- Apply clipped y-limit ----------------
        if cap is not None:
            ymin = float(d[y].min())
            pad = 0.08 * max(cap - ymin, 1e-9)
            ax.set_ylim(ymin, cap + pad)

        # ---------------- Upper outliers as lollipops ----------------
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

        # ---------------- Uncorrected Welch stats per DIV ----------------
        all_tests = []
        pairs = list(combinations(range(len(group_order)), 2))

        for div_val in div_order:
            div_data = d[d[div_col] == div_val]
            for idx1, idx2 in pairs:
                g1, g2 = group_order[idx1], group_order[idx2]
                d1 = div_data.loc[div_data[group_col] == g1, y].dropna()
                d2 = div_data.loc[div_data[group_col] == g2, y].dropna()

                if len(d1) > 1 and len(d2) > 1:
                    t_stat, p_val = stats.ttest_ind(d1, d2, equal_var=False, nan_policy="omit")
                    all_tests.append({
                        "div_val": div_val,
                        "idx1": idx1,
                        "idx2": idx2,
                        "g1": g1,
                        "g2": g2,
                        "p": p_val
                    })

        # ---------------- Statistical brackets ----------------
        if annotate and len(all_tests) > 0:
            ymin_ax, ymax_ax = ax.get_ylim()
            yr = ymax_ax - ymin_ax

            # keep bracket region below outlier triangles
            base_top = ymax_ax - 0.20 * yr
            bracket_step = 0.06 * yr
            max_annotation_y = base_top

            for div_val in div_order:
                div_tests = [t for t in all_tests if t["div_val"] == div_val]
                if not div_tests:
                    continue

                level = 0
                for t in div_tests:
                    p = t["p"]
                    if not np.isfinite(p):
                        continue
                    if p < 0.001:
                        stars = "***"
                    elif p < 0.01:
                        stars = "**"
                    elif p < 0.05:
                        stars = "*"
                    else:
                        continue

                    g1 = t["g1"]
                    g2 = t["g2"]

                    x1 = div_to_x[div_val] + offsets[g1]
                    x2 = div_to_x[div_val] + offsets[g2]

                    bracket_y = base_top + level * bracket_step

                    ax.plot(
                        [x1, x1, x2, x2],
                        [bracket_y,
                        bracket_y + 0.25 * bracket_step,
                        bracket_y + 0.25 * bracket_step,
                        bracket_y],
                        c="black",
                        lw=1.2,
                        zorder=7
                    )

                    ax.text(
                        (x1 + x2) / 2,
                        bracket_y + 0.30 * bracket_step,
                        stars,
                        ha="center",
                        va="bottom",
                        fontsize=9,
                        zorder=8
                    )

                    max_annotation_y = max(max_annotation_y, bracket_y + 0.55 * bracket_step)
                    level += 1

            if max_annotation_y > ax.get_ylim()[1]:
                ax.set_ylim(ax.get_ylim()[0], max_annotation_y + 0.04 * yr)

        # ---------------- Axes / ticks / labels ----------------
        #ax.set_xticks(div_order)
        ax.set_xticks(range(len(div_order)))
        ax.set_xticklabels([str(v) for v in div_order])
        ax.set_xlim(-0.5,len(div_order)-0.5)

        ax.set_xlabel(xlabel if xlabel else div_col, fontweight="bold", fontsize=8)
        ax.set_ylabel(ylabel if ylabel else y.replace("_", " "), fontweight="bold", fontsize=8)
        ax.set_title(title if title else f"{y} trajectory", fontweight="bold", pad=6, fontsize=9)

        ax.tick_params(axis="x", labelsize=6, width=0.8, length=3)
        ax.tick_params(axis="y", labelsize=6, width=0.8, length=3)

        if legend_outside:
            ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, title_fontsize=8)
        else:
            ax.legend(frameon=False, loc=legend_loc, fontsize=7, title_fontsize=8)

        sns.despine(ax=ax)
        return ax


    def plot_bars_by_div(
        self, df, div_col, group_col, y,
        div_order=None, group_order=None,
        palette=None, title=None, ax=None, annotate=True,
        # clipping / outliers
        clip_upper=True,
        clip_method="quantile",  # "iqr" or "quantile"
        clip_k=1.5,
        clip_q=0.98,
        show_upper_outliers=True,
        # legend control
        legend_loc="upper right",      # inside by default
        legend_outside=False
    ):
        """
        Grouped bars by DIV with scatter.
        Upper clipping prevents one insane well from flattening the plot.
        Upper outliers shown as clipped lollipop markers (triangles) at top.
        Welch stats are handled by the wrapper; here we only draw brackets from p-values passed in.
        """
        data = df.dropna(subset=[div_col, group_col, y]).copy()
        data[div_col] = pd.to_numeric(data[div_col], errors="coerce")
        data = data.dropna(subset=[div_col])
        data[div_col] = data[div_col].astype(int)

        div_order = self._resolve_div_order(data, div_col, div_order)
        group_order = self._resolve_group_order(data, group_col, group_order)
        palette = self._resolve_palette(palette)

        if ax is None:
            fig, ax = plt.subplots(figsize=(max(10, len(div_order) * 1.1 + 2), 6))

        # ---- plot bars + points ----
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
            edgecolor="black",
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
            edgecolor="black",
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

        # ---- compute clip cap ----
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
                pad = 0.06 * max(cap - ymin, 1e-9)   # SMALL pad to avoid huge whitespace
                ax.set_ylim(ymin, cap + pad)

        # ---- show upper outliers as clipped triangles at top ----
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
                    x = div_idx + (g_idx - n_groups/2 + 0.5) * bar_width

                    ax.plot([x, x], [y_stem0, y_tip], color="black", lw=0.9, zorder=10)
                    ax.scatter([x], [y_tip], marker="^", s=42, color="black", zorder=11)

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
        - runs Welch stats (uncorrected)
        - makes plot
        - adds stats table
        - saves a single PDF page

        NOTE: No multiple-comparison correction (as requested).
        """
        if not filename.endswith(".pdf"):
            filename += ".pdf"

        palette = self._resolve_palette(palette)

        # ---- compute stats ----
        if mode == "by_div":
            stats_df = self.calculate_stats_by_div_welch(
                df, div_col=div_col, group_col=group_col, y=y,
                div_order=div_order, group_order=group_order
            )
        elif mode == "by_group":
            if x is None:
                raise ValueError("For mode='by_group', you must pass x=<group column> (e.g., 'NeuronType').")
            stats_df = self.calculate_stats_welch(df, group_col=x, y=y, group_order=group_order)
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
                df=df, x=x, y=y,
                order=group_order,
                palette=palette,
                title=title,
                ax=ax_plot
            )

        # ---- table ----
        if stats_df is None or len(stats_df) == 0:
            ax_table.text(0.5, 0.5, "No stats to display", ha="center", va="center")
        else:
            # choose columns
            if "DIV" in stats_df.columns:
                cols = ["DIV", "Comparison", "Grp1_Stats", "Grp2_Stats", "t-stat", "p-val", "Sig", "Cohen's d"]
            else:
                cols = ["Comparison", "Grp1_Stats", "Grp2_Stats", "t-stat", "p-val", "Sig", "Cohen's d"]

            # format
            cell_text = []
            for _, row in stats_df.iterrows():
                row_vals = []
                for c in cols:
                    v = row.get(c, "")
                    if isinstance(v, float) and np.isnan(v):
                        row_vals.append("NA")
                    elif c == "t-stat" and isinstance(v, (float, np.floating)):
                        row_vals.append(f"{v:.2f}")
                    elif c == "p-val" and isinstance(v, (float, np.floating)):
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

        # ---- layout control (NO tight_layout; it fights legends) ----
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