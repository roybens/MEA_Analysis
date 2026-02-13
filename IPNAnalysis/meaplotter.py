import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations
from statannotations.Annotator import Annotator

class MEAPlotter:
    def __init__(self):
        """
        Initializes the plotting engine with Publication-Quality settings.
        No data processing happens here.
        """
        
        # 1. TEXT AS TEXT (Crucial for Inkscape/Illustrator)
        plt.rcParams['svg.fonttype'] = 'none' 
        
        # 2. PDF COMPATIBILITY (Type 42 = TrueType)
        # This embeds the font so it looks correct on any computer
        plt.rcParams['pdf.fonttype'] = 42
        
        # 3. FONT FALLBACK (Solves your Linux error)
        # It tries Arial first (Publication Standard). 
        # If missing, it uses DejaVu Sans (Linux Standard), which Inkscape handles well.
        plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans', 'sans-serif']
        plt.rcParams['font.family'] = 'sans-serif'


    def analyze(self, df, x, y, order=None, palette=None, title=None, save_path=None):
        """
        Exact replication of your analysis snippet:
        1. Pairwise Welch's T-tests + Cohen's d
        2. Bar Plot + Strip Plot (with stars)
        3. Descriptive Summary Table
        """
        # Create a working copy
        data = df.copy()
        
        # 1. Setup Order
        if order is None:
            order = sorted(data[x].unique())
        
        # Filter data to strictly match the order list (removes unused categories)
        data = data[data[x].isin(order)]

        # ==========================================
        # PART A: STATISTICS (Your Exact Logic)
        # ==========================================
        print(f"\n{'='*20} STATISTICAL ANALYSIS: {y} {'='*20}")
        print("Pairwise Welch's t-tests (Equal Variance = False):")
        print("-" * 60)

        pairs = list(combinations(order, 2))
        stats_results = []
        
        for group1, group2 in pairs:
            d1 = data[data[x] == group1][y]
            d2 = data[data[x] == group2][y]
            
            # Welch's t-test
            t_stat, p_val = stats.ttest_ind(d1, d2, equal_var=False)
            
            # Cohen's d (Using the formula from your snippet)
            mean1, mean2 = d1.mean(), d2.mean()
            std1, std2 = d1.std(), d2.std()
            n1, n2 = len(d1), len(d2)
            
            # Pooled STD (Unweighted average variance, as per your snippet)
            pooled_std = np.sqrt((std1**2 + std2**2) / 2)
            cohens_d = (mean1 - mean2) / pooled_std
            
            stats_results.append({
                'Comparison': f'{group1} vs {group2}',
                't-statistic': t_stat,
                'p-value': p_val,
                'Cohen\'s d': cohens_d,
                'Mean1': mean1,
                'Mean2': mean2
            })
            
            # Print Log
            sig = "ns"
            if p_val < 0.001: sig = "***"
            elif p_val < 0.01: sig = "**"
            elif p_val < 0.05: sig = "*"
            
            print(f"{group1} vs {group2}:")
            print(f"  t({n1+n2-2}) = {t_stat:.3f}, p = {p_val:.4e} ({sig})")
            print(f"  Diff: {mean1:.1f} - {mean2:.1f} = {mean1-mean2:.1f}")
            print(f"  Cohen's d = {cohens_d:.3f}\n")

        # ==========================================
        # PART B: VISUALIZATION
        # ==========================================
        fig, ax = plt.subplots(figsize=(6, 6))

        # 1. Bar Plot (Mean +/- SE)
        sns.barplot(
            data=data, x=x, y=y, order=order, palette=palette,
            errorbar='se', capsize=0.15, alpha=0.8, 
            edgecolor='black', linewidth=1.5, ax=ax
        )

        # 2. Strip Plot (Individual Points)
        sns.stripplot(
            data=data, x=x, y=y, order=order, palette=palette,
            color='black', alpha=0.6, jitter=0.2, size=5,
            edgecolor='black', linewidth=0.5, ax=ax
        )

        # 3. Add Significance Brackets (Automated)
        try:
            annotator = Annotator(ax, pairs, data=data, x=x, y=y, order=order)
            annotator.configure(test='t-test_welch', text_format='star', loc='outside', verbose=False)
            annotator.apply_and_annotate()
        except Exception:
            print("(!) Could not draw stats brackets (possibly insufficient data).")

        # Formatting
        ax.set_title(title if title else y, fontweight='bold', pad=15)
        ax.set_xlabel("")
        ax.set_ylabel(y.replace('_', ' '), fontweight='bold')
        sns.despine()
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to: {save_path}")
            
        plt.show()

        # ==========================================
        # PART C: SUMMARY TABLE
        # ==========================================
        print(f"\n{'='*20} DESCRIPTIVE STATISTICS {'='*20}")
        summary = data.groupby(x)[y].agg(['count', 'mean', 'std', 'sem', 'median'])
        print(summary)
        
        return pd.DataFrame(stats_results)