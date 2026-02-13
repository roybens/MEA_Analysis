import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations
from statannotations.Annotator import Annotator
from matplotlib.backends.backend_pdf import PdfPages

class MEAPlotter:
    def __init__(self):
        """
        Initializes plotting with Publication-Quality settings.
        """
        # --- FONT CONFIGURATION ---
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans', 'sans-serif']
        plt.rcParams['svg.fonttype'] = 'none' 
        plt.rcParams['pdf.fonttype'] = 42

    def calculate_stats(self, df, x, y, order=None):
        """
        STEP 1: Pure Analysis.
        Returns a DataFrame containing T-stats, P-values, and Cohen's d.
        """
        data = df.copy()
        if order is None: order = sorted(data[x].unique())
        
        pairs = list(combinations(order, 2))
        stats_list = []
        
        print(f"\n{'='*20} STATISTICAL ANALYSIS: {y} {'='*20}")
        
        for group1, group2 in pairs:
            d1 = data[data[x] == group1][y]
            d2 = data[data[x] == group2][y]
            
            # Welch's t-test
            t_stat, p_val = stats.ttest_ind(d1, d2, equal_var=False)
            
            # Cohen's d
            pooled_std = np.sqrt((d1.std()**2 + d2.std()**2) / 2)
            cohens_d = (d1.mean() - d2.mean()) / pooled_std
            
            # Stars
            stars = "ns"
            if p_val < 0.001: stars = "***"
            elif p_val < 0.01: stars = "**"
            elif p_val < 0.05: stars = "*"

            stats_list.append({
                "Comparison": f"{group1} vs {group2}",
                "t-stat": t_stat,
                "p-val": p_val,
                "Sig": stars,
                "Cohen's d": cohens_d
            })
            
            print(f"{group1} vs {group2}: p={p_val:.4e} ({stars}), d={cohens_d:.3f}")
            
        return pd.DataFrame(stats_list)

    def plot_bars(self, df, x, y, order=None, palette=None, title=None, ax=None):
        """
        STEP 2: Pure Visualization.
        Draws the Bar+Strip plot and significance brackets on a given axis.
        Returns the axis.
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6))

        # 1. Bar Plot
        sns.barplot(data=df, x=x, y=y, order=order, palette=palette,
                    errorbar='se', capsize=0.15, alpha=0.8, 
                    edgecolor='black', linewidth=1.5, ax=ax)

        # 2. Strip Plot
        sns.stripplot(data=df, x=x, y=y, order=order, palette=palette,
                      color='black', alpha=0.6, jitter=0.2, size=5, 
                      edgecolor='black', linewidth=0.5, ax=ax)

        # 3. Add Significance Brackets (Calculates its own stats for placement)
        if order and len(order) > 1:
            pairs = list(combinations(order, 2))
            try:
                annotator = Annotator(ax, pairs, data=df, x=x, y=y, order=order)
                annotator.configure(test='t-test_welch', text_format='star', loc='outside', verbose=False)
                annotator.apply_and_annotate()
            except Exception as e:
                print(f"(!) Could not draw brackets: {e}")

        # Formatting
        ax.set_title(title if title else y, fontweight='bold', pad=15)
        ax.set_xlabel("")
        ax.set_ylabel(y.replace('_', ' '), fontweight='bold')
        sns.despine(ax=ax)
        
        return ax

    def save_pdf(self, df, x, y, order, stats_df, palette=None, title=None, filename="Report.pdf"):
        """
        STEP 3: Report Generation.
        Combines Plot (Top) and Stats Table (Bottom) into one PDF.
        """
        # Create a figure with 2 rows (Plot is 3x taller than Table)
        fig = plt.figure(figsize=(8, 11)) 
        gs = fig.add_gridspec(2, 1, height_ratios=[3, 1])

        # Draw Plot in Top Section
        ax_plot = fig.add_subplot(gs[0])
        self.plot_bars(df, x, y, order, palette, title, ax=ax_plot)

        # Draw Table in Bottom Section
        ax_table = fig.add_subplot(gs[1])
        ax_table.axis('off')
        
        # Prepare table text
        cell_text = []
        for _, row in stats_df.iterrows():
            # FIX: Use .format() instead of f-string for key with quote
            cohen_val = "{:.3f}".format(row["Cohen's d"])
            
            cell_text.append([
                row['Comparison'],
                f"{row['t-stat']:.3f}",
                f"{row['p-val']:.4e}",
                row['Sig'],
                cohen_val
            ])
            
        col_labels = ["Comparison", "t-stat", "p-value", "Sig", "Cohen's d"]
        
        # Add table
        table = ax_table.table(cellText=cell_text, colLabels=col_labels, 
                               loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.2, 1.5)

        plt.tight_layout()
        
        # Save
        # Auto-append .pdf if missing
        if not filename.endswith('.pdf'): filename += ".pdf"

        with PdfPages(filename) as pdf:
            pdf.savefig(fig)
        print(f"✅ PDF saved to: {filename}")
        plt.close(fig)